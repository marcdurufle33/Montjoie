#ifndef MONTJOIE_FILE_WEDGE_CLASSICAL_CXX

namespace Montjoie
{
  
  //! default constructor
  WedgeClassical::WedgeClassical() : WedgeReference<1>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  size_t WedgeClassical::GetMemorySize() const
  {
    size_t taille = WedgeReference<1>::GetMemorySize();
    taille += CoordinateDofs2D_quad.GetMemorySize() + CoordinateDofs2D_tri.GetMemorySize()
      + CoordinateDofs3D.GetMemorySize() + NumDofs2D_quad.GetMemorySize()
      + NumDofs2D_tri.GetMemorySize();
    taille += NumDofs3D.GetMemorySize();
    taille += function_basis_tri.GetMemorySize();
    taille += lob_basis.GetMemorySize() + lob.GetMemorySize();
    taille += CoordinateDof.GetMemorySize() + NumDofTri.GetMemorySize();
    taille += G_Gquad.GetMemorySize() + GL_Gquad.GetMemorySize() + dG_Gquad.GetMemorySize()
      + GL_G.GetMemorySize() + G_G.GetMemorySize() + Gquad_GL.GetMemorySize();
    taille += G_Geom.GetMemorySize() + dG_Geom.GetMemorySize() + ChGeom.GetMemorySize()
      + GL_Geom.GetMemorySize() + RhGeom.GetMemorySize() + dGL_Geom.GetMemorySize();
    taille += ChTri.GetMemorySize() + RhTri.GetMemorySize() + ChConstTri.GetMemorySize()
      + RhConstTri.GetMemorySize() + RhNodalTri.GetMemorySize() + 
      RhTriX.GetMemorySize() + RhTriY.GetMemorySize();
    
    taille += RhConstTriX.GetMemorySize() + RhConstTriY.GetMemorySize();
    taille += ValGauss0.GetMemorySize() + ValGauss1.GetMemorySize()
      + DerivGauss0.GetMemorySize() + DerivGauss1.GetMemorySize();
    
    taille += GradPhiWeight.GetMemorySize() + stiffness_matrix1D.GetMemorySize();
    return taille;
  }
  
  
  //! constructing finite element
  void WedgeClassical
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                           int rtri, int r_quad, int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    if (gauss_z == -1)
      gauss_z = Globatto<Real_wp>::QUADRATURE_LOBATTO;
    
    WedgeReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                              rtri, r_quad, type_surf_tri, type_surf_quad, gauss_z);

    ConstructFunctions(gauss_z, type_surf_quad);
    
    this->Fb_geom.ComputeCoefficientTransformation();    
  }
  

  //! construction of basis functions
  void WedgeClassical::ConstructFunctions(int gauss_z, int type_surf_quad)
  {
    this->nb_dof_quad = (order+1)*(order+1);
    this->nb_dof_tri = (order+1)*(order+2)/2;

    // along z, we are using either Gauss or Gauss-Lobatto points (H1 or DG)
    if (gauss_z == -1)
      gauss_z = lob_basis.QUADRATURE_LOBATTO;
    
    lob_basis.ConstructQuadrature(order, gauss_z);

    int r = order;
    int order_geom = this->GetGeometryOrder();
    
    VectReal_wp points_dof1d;
    VectR2 points_dof2d_tri, points_dof2d_quad;
    VectR3 points_dof3d;    

    if (order == order_geom)
      {
	//this->elt_geom.dof_equal_nodal = true;
        nb_dof_loc = this->GetNbPointsNodalElt();
        
	points_dof2d_tri = this->PointsNodal2D_tri();
        this->SetPointsDof2D_tri(this->PointsNodal2D_tri());
                
        FacesDof = this->Fb_geom.GetNodalNumber();
        NumDofs3D = this->Fb_geom.GetNumNodes3D();
	CoordinateDofs3D = this->Fb_geom.GetCoordinateNodes3D();
        NumDofs2D_quad = this->GetNumNodes2D_quad();
	CoordinateDofs2D_quad = this->GetCoordinateNodes2D_quad();
        NumDofs2D_tri = this->GetNumNodes2D_tri();
	CoordinateDofs2D_tri = this->GetCoordinateNodes2D_tri();
	
        NumDofTri = Fb_geom.GetNumNodesTri();
        CoordinateDof = this->Fb_geom.GetCoordinateNodes3D();
        
        function_basis_tri.ConstructFiniteElement(r, r, r, -1, r, type_surf_quad);
      }
    else
      {
	//this->elt_geom.dof_equal_nodal = false;
        MeshNumbering<Dimension3>::ConstructPrismaticNumbering(r, NumDofs3D, CoordinateDofs3D);
        MeshNumbering<Dimension2>::
          ConstructQuadrilateralNumbering(r, NumDofs2D_quad, CoordinateDofs2D_quad);
        MeshNumbering<Dimension2>::
          ConstructTriangularNumbering(r, NumDofs2D_tri, CoordinateDofs2D_tri);
    
        FacesDof.Reallocate((r+1)*(r+1), 5);
        FacesDof.Fill(-1);
        // quadrangular faces
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            {
              int node = NumDofs2D_quad(i,j);
              FacesDof(node, 1) = NumDofs3D(i, 0, j);
              FacesDof(node, 2) = NumDofs3D(r-i, i, j);
              FacesDof(node, 3) = NumDofs3D(0, i, j);
            }
        
        // triangular faces
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            {
              int node = NumDofs2D_tri(i,j);
              FacesDof(node,0) = NumDofs3D(i, j, 0);
              FacesDof(node,4) = NumDofs3D(i, j, r);
            }
        
	WedgeGeomReference::
	  ConstructLobattoPoints(r, points_dof1d, points_dof2d_tri, NumDofs2D_tri,
				 points_dof2d_quad, NumDofs2D_quad, points_dof3d, NumDofs3D);

	this->SetPointsDof1D(points_dof1d);
	this->SetPointsDofND(points_dof3d);
        
        function_basis_tri.ConstructFiniteElement(r, r, r, -1, r, type_surf_quad);
        
        nb_dof_loc = points_dof3d.GetM();
        CoordinateDof.Reallocate(nb_dof_loc, 2);
        NumDofTri.Reallocate(points_dof2d_tri.GetM(), order+1);
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r-i; j++)
            for (int k = 0; k <= r; k++)
              {
                int node = NumDofs3D(i, j, k);
                CoordinateDof(node, 0) = NumDofs2D_tri(i, j);
                CoordinateDof(node, 1) = k;
                NumDofTri(NumDofs2D_tri(i, j), k) = node;
              }
        
      }
    
    nb_dof_boundaries = (order+1)*(order+2)+3*order*(order-1);
    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    points_dof1d = function_basis_tri.PointsNodal1D();
    this->SetPointsDof1D(function_basis_tri.PointsNodal1D());
    points_dof3d.Reallocate((order+1)*(order+2)*(order+1)/2);
    points_dof2d_quad.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  points_dof2d_quad(NumDofs2D_quad(i, j)).Init(points_dof1d(i), lob_basis.Points(j));
	  for (int k = 0; k <= order-i; k++)
	    {
	      int ndof = NumDofs2D_tri(i, k);
	      Real_wp x = points_dof2d_tri(ndof)(0);
	      Real_wp y = points_dof2d_tri(ndof)(1);
	      points_dof3d(NumDofTri(ndof, j)).Init(x, y, lob_basis.Points(j));
	    }
	}
    
    this->SetPointsDof2D_tri(points_dof2d_tri);
    this->SetPointsDof2D_quad(points_dof2d_quad);    
    this->SetPointsDofND(points_dof3d);

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
    
    lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);
    ValGauss0.Reallocate(order+1);
    ValGauss1.Reallocate(order+1);
    DerivGauss0.Reallocate(order+1);
    DerivGauss1.Reallocate(order+1);
    GradPhiWeight.Reallocate(order+1, order+1);
    lob_basis.ComputeGradPhi(1e3*epsilon_machine);
    for (int i = 0; i <= order; i++)
      {
	for (int j = 0; j <= order; j++)
          GradPhiWeight(i, j) = lob_basis.GradPhi(i, j)*lob_basis.Weights(j);
		
	ValGauss0(i) = lob_basis.EvaluatePhi(i, Real_wp(0));
	ValGauss1(i) = lob_basis.EvaluatePhi(i, Real_wp(1));

        DerivGauss0(i) = lob_basis.EvaluatePhiGrad(i, Real_wp(0));
	DerivGauss1(i) = lob_basis.EvaluatePhiGrad(i, Real_wp(1));
      }

    G_Geom.Reallocate(order+1, order_geom+1);
    GL_Geom.Reallocate(order+1, order_geom+1);
    dGL_Geom.Reallocate(order+1, order_geom+1);
    dG_Geom.Reallocate(order+1, order_geom+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_geom; j++)
        {
          G_Geom(i, j) = lob_basis.EvaluatePhi(i, this->PointsNodal1D(j));
          GL_Geom(i, j) = lob.EvaluatePhi(i, this->PointsNodal1D(j));
          dGL_Geom(i, j) = lob.EvaluatePhiGrad(i, this->PointsNodal1D(j));
          dG_Geom(i, j) = lob_basis.EvaluatePhiGrad(i, this->PointsNodal1D(j));
        }

    GL_Gquad.Reallocate(order+1, order_quad+1);
    G_Gquad.Reallocate(order+1, order_quad+1);
    dG_Gquad.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        {
          GL_Gquad(i, j) = lob.EvaluatePhi(i, this->Points1D(j));
          G_Gquad(i, j) = lob_basis.EvaluatePhi(i, this->Points1D(j));
          dG_Gquad(i, j) = lob_basis.EvaluatePhiGrad(i, this->Points1D(j));
        }
    
    Globatto<Real_wp> gauss_basis;
    gauss_basis.ConstructQuadrature(order, gauss_basis.QUADRATURE_GAUSS);
    GL_G.Reallocate(order+1, order_quad+1);
    G_G.Reallocate(order+1, order_quad+1);
    Gquad_GL.Reallocate(order+1, order_quad+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order_quad; j++)
        {
          GL_G(i,j) = lob.EvaluatePhi(i, gauss_basis.Points(j));
          G_G(i,j) = lob_basis.EvaluatePhi(i, gauss_basis.Points(j));
          Gquad_GL(i, j) = gauss_basis.EvaluatePhi(i, lob.Points(j));
        }

    int Ntri = points_dof2d_tri.GetM();
    VectReal_wp phi;     VectR2 grad_phi;
    const VectR2& points_nodal2d_tri = this->PointsNodal2D_tri();
    const VectR2& points2d_tri = this->Points2D_tri();
    const VectReal_wp& weights2d_tri = this->Weights2D_tri();
    ChGeom.Reallocate(Ntri, points_nodal2d_tri.GetM());
    RhGeom.Reallocate(Ntri, 2*points_nodal2d_tri.GetM());
    for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
      {
        function_basis_tri.ComputeValuesPhiRef(points_nodal2d_tri(j), phi);
        function_basis_tri.ComputeGradientPhiRef(points_nodal2d_tri(j), grad_phi);
        for (int i = 0; i < Ntri; i++)
          {
            ChGeom(i, j) = phi(i);
            RhGeom(i, 2*j) = grad_phi(i)(0);
            RhGeom(i, 2*j+1) = grad_phi(i)(1);
          }
            
      }
    
    int Nquad = points2d_tri.GetM();
    ChTri.Reallocate(Ntri, Nquad);
    RhTri.Reallocate(Ntri, 2*Nquad);
    RhTriX.Reallocate(Ntri, 2*Nquad);
    RhTriY.Reallocate(Ntri, 2*Nquad);
    RhNodalTri.Reallocate(Ntri, 2*Ntri);
    ChConstTri.Reallocate(Ntri, Ntri);
    RhConstTri.Reallocate(Ntri, 2*Ntri);
    RhConstTriX.Reallocate(Ntri, 2*Ntri);
    RhConstTriY.Reallocate(Ntri, 2*Ntri);
    ChTri.Fill(0); RhTri.Fill(0);
    ChConstTri.Fill(0); RhConstTri.Fill(0);
    RhConstTriX.Fill(0); RhConstTriY.Fill(0);
    mass_matrix.Reallocate(Ntri, Ntri);
    mass_matrix.Fill(0);

    for (int k = 0; k < Ntri; k++)
      {
        function_basis_tri.ComputeGradientPhiNodalRef(points_dof2d_tri(k), grad_phi);
	for (int i = 0; i < Ntri; i++)
          {
            RhNodalTri(i, 2*k) = grad_phi(i)(0);
            RhNodalTri(i, 2*k+1) = grad_phi(i)(1);
          }
      }
    
    for (int k = 0; k < Nquad; k++)
      {
	function_basis_tri.ComputeValuesPhiNodalRef(points2d_tri(k), phi);
	function_basis_tri.ComputeGradientPhiNodalRef(points2d_tri(k), grad_phi);
	for (int i = 0; i < Ntri; i++)
	  {
	    ChTri(i, k) = phi(i);
	    RhTri(i, 2*k) = grad_phi(i)(0);
	    RhTri(i, 2*k+1) = grad_phi(i)(1);
	    RhTriX(i, 2*k) = grad_phi(i)(0);
	    RhTriX(i, 2*k+1) = 0;
	    RhTriY(i, 2*k) = 0;
	    RhTriY(i, 2*k+1) = grad_phi(i)(1);
	    for (int j = 0; j < Ntri; j++)
	      {
		ChConstTri(i, j) += weights2d_tri(k)*phi(i)*phi(j);
		RhConstTri(i, 2*j) += weights2d_tri(k)*phi(j)*grad_phi(i)(0);
		RhConstTri(i, 2*j+1) += weights2d_tri(k)*phi(j)*grad_phi(i)(1);
                RhConstTriX(i, 2*j) += weights2d_tri(k)*phi(j)*grad_phi(i)(0);
                RhConstTriY(i, 2*j+1) += weights2d_tri(k)*phi(j)*grad_phi(i)(1);
                if (j >= i)
                  mass_matrix(i, j) += weights2d_tri(k)*phi(i)*phi(j);
	      }
	  }
      }
    
    mass_matrix_chol = mass_matrix;
    GetCholesky(mass_matrix_chol);
    
    stiffness_matrix1D.Reallocate(order+1, order+1);
    stiffness_matrix1D.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        for (int k = 0; k <= order; k++)
          stiffness_matrix1D(i, j) += weights1d_z(k)*lob_basis.GradPhi(i, k)*lob_basis.GradPhi(j, k);
    
    const_stiff_matrix.Reallocate(2*Ntri, 2*Ntri);
    const_grad_matrix.Reallocate(Ntri, 2*Ntri);
    
    Real_wp vloc;
    for (int i = 0; i < Ntri; i++)
      for (int j = 0; j < Ntri; j++)
	{
          for (int m = 0; m < 2; m++)
            {
              vloc = 0;
              for (int k = 0; k < Nquad; k++)
                vloc += weights2d_tri(k)*RhTri(i, 2*k+m)*ChTri(j, k);
              
              const_grad_matrix(i, 2*j+m) = vloc;
            }
          
          for (int m = 0; m < 2; m++)
            for (int n = 0; n < 2; n++)
              {
                vloc = 0;
                for (int k = 0; k < Nquad; k++)
                  vloc += weights2d_tri(k)*RhTri(i, 2*k+m)*RhTri(j, 2*k+n);
                
                const_stiff_matrix(2*i + m, 2*j + n) = vloc;
              }
        }
    
    TriangleClassical* tri = new TriangleClassical();
    QuadrangleGauss* quad = new QuadrangleGauss();
    tri->ConstructFiniteElement(order, order, order_quad);    
    quad->ConstructFiniteElement(order, order_geom, order_quad, type_surf_quad);
    element_tri_surf = tri;
    element_quad_surf = quad;
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
  void WedgeClassical::
  ComputeGaussIntegralSurfaceGen(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    if (num_loc%4 != 0)
      {
        res.Fill(0);
	const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
	int Nquad = (order_quad+1)*(order_quad+1);
        Vector2 w1(Nquad); FillZero(w1);
        Vector2 w2(Nquad); FillZero(w2);
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order_quad; j++)
            for (int k = 0; k <= order_quad; k++)
              w1(i*(order_quad+1) + j) += GL_G(i, k)*feval(NumQuad2D(k, j));
        
        for (int i = 0; i <= order; i++)
          for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order_quad; k++)
              w2(NumDofs2D_quad(i, j)) += G_G(j, k)*w1(i*(order_quad+1) + k);
        
        for (int j = 0; j < (order+1)*(order+1); j++)
          {
            int num_dof = FacesDof(j, num_loc);
            res(num_dof) = w2(j);
          }
      }
    else
      ComputeIntegralSurfaceRef(feval, res, num_loc);
  }


  //! computation of U on nodal points
  /*!
    \param[in] Un components of U on dofs
    \param[out] Unode values of U on nodal points
   */
  template<class Vector1, class Vector2>
  void WedgeClassical::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    int order_geom = this->GetGeometryOrder();
    if ( order == order_geom)
      {
        if (abs(1.0-ValGauss0(0)) < 1e4*epsilon_machine)
          {
            // Lobatto case, nodal and dofs points are the same
            Unode = Un;
          }
        else
          {
            // Gauss case
            Unode.Fill(0);
	    int Ntri = this->PointsNodal2D_tri().GetM();
            for (int j = 0; j < Ntri; j++)
              for (int i = 0; i <= order; i++)
                {
                  int node = NumDofTri(j, i);
                  for (int k = 0; k <= order; k++)
                    Unode(node) += G_Geom(k, i)*Un(NumDofTri(j, k));
                }
          }
      }
    else
      {
        Unode.Fill(0);
	int Ntri = this->PointsNodal2D_tri().GetM();
        Vector1 Utri_node(Ntri);
        Vector1 Utri_dof(this->nb_dof_tri);
        Utri_node.Fill(0);
	const Matrix<int>& NumNodesTri = Fb_geom.GetNumNodesTri();
        for (int i = 0; i <= order; i++)
          {
            for (int j = 0; j < this->nb_dof_tri; j++)
              Utri_dof(j) = Un(NumDofTri(j, i));
            
            Mlt(SeldonTrans, ChGeom, Utri_dof, Utri_node);
            
            for (int j = 0; j < Ntri; j++)
              for (int k = 0; k <= order_geom; k++)
                Unode(NumNodesTri(j, k)) += G_Geom(i, k)*Utri_node(j);
            }
      }
  }
  
  
  //! computation of u on nodal points of a face
  /*!
    \param[in] Un dof components of u
    \param[out] Unode values of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1, class Vector2>
  void WedgeClassical
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Fill(0);
    typename Vector2::value_type vloc;
    int Ntri = this->PointsNodal2D_tri().GetM();
    int order_geom = this->GetGeometryOrder();
    const Matrix<int>& NumNodes2D_quad = this->GetNumNodes2D_quad();
    switch (num_loc)
      {
      case 0 :
        {
          Vector2 Udof(this->nb_dof_tri);
          Udof.Fill(0);
          for (int i = 0; i <= order; i++)
            for (int j = 0; j < Udof.GetM(); j++)
              Udof(j) += ValGauss0(i)*Un(NumDofTri(j, i));
          
          for (int i = 0; i < Ntri; i++)
            for (int j = 0; j < this->nb_dof_tri; j++)
              Unode(i) += ChGeom(j, i)*Udof(j);
          
        }
        break;
      case 1 :
        {
          Vector2 Utri(order_geom + 1);
          for (int i = 0; i <= order; i++)
            {
              Utri.Fill(0);
              for (int j = 0; j <= order; j++)
                {
                  int j2 = function_basis_tri.GetNodalNumber(0, j);
                  vloc = Un(NumDofTri(j2, i));
                  
                  for (int k = 0; k <= order_geom; k++)
                    Utri(k) += GL_Geom(j, k)*vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                for (int k = 0; k <= order_geom; k++)
                  Unode(NumNodes2D_quad(j, k)) += G_Geom(i, k)*Utri(j);
            }                            
        }
        break;
      case 2 :
        {
          Vector2 Utri(order_geom+1);
          for (int i = 0; i <= order; i++)
            {
              Utri.Fill(0);
              for (int j = 0; j <= order; j++)
                {
                  int j2 = function_basis_tri.GetNodalNumber(1, j);
                  vloc = Un(NumDofTri(j2, i));
                  
                  for (int k = 0; k <= order_geom; k++)
                    Utri(k) += GL_Geom(j, k)*vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                for (int k = 0; k <= order_geom; k++)
                  Unode(NumNodes2D_quad(j, k)) += G_Geom(i, k)*Utri(j);
            }                            
        }
        break;
      case 3 :
        {
          Vector2 Utri(order_geom+1);
          for (int i = 0; i <= order; i++)
            {
              Utri.Fill(0);
              for (int j = 0; j <= order; j++)
                {
                  int j2 = function_basis_tri.GetNodalNumber(2, order-j);
                  vloc = Un(NumDofTri(j2, i));
                  
                  for (int k = 0; k <= order_geom; k++)
                    Utri(k) += GL_Geom(j, k)*vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                for (int k = 0; k <= order_geom; k++)
                  Unode(NumNodes2D_quad(j, k)) += G_Geom(i, k)*Utri(j);
            }                            
        }
        break;
      case 4 :
        {
          Vector2 Udof(this->nb_dof_tri);
          Udof.Fill(0);
          for (int i = 0; i <= order; i++)
            for (int j = 0; j < Udof.GetM(); j++)
              Udof(j) += G_Geom(i, order_geom)*Un(NumDofTri(j, i));
          
          for (int i = 0; i < Ntri; i++)
            for (int j = 0; j < this->nb_dof_tri; j++)
              Unode(i) += ChGeom(j, i)*Udof(j);
          
        }
        break;
      }
        
  }
  
  
  //! computation of grad u on nodal points of a face
  /*!
    \param[in] Un dof components of u
    \param[out] Unode gradient of u on nodal points of the face
    \param[in] num_loc local position of the face in the element
  */  
  template<class Vector1, class Vector2>
  void WedgeClassical
  ::ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Fill(0);
    typename Vector2::value_type vloc;
    int Ntri = this->PointsNodal2D_tri().GetM();
    int order_geom = this->GetGeometryOrder();
    const Matrix<int>& NumNodes2D_quad = this->GetNumNodes2D_quad();
    const TriangleGeomReference& function_geom_tri = Fb_geom.GetTriangularSurfaceFiniteElement();
    switch (num_loc)
      {
      case 0 :
        {
          Vector2 Ux(this->nb_dof_tri);
          Ux.Fill(0);
          for (int i = 0; i <= order; i++)
            for (int j = 0; j < Ux.GetM(); j++)
              Ux(j) += DerivGauss0(i)*Un(NumDofTri(j, i));
          
          for (int i = 0; i < Ntri; i++)
            for (int j = 0; j < this->nb_dof_tri; j++)
              Unode(3*i+2) += ChGeom(j, i)*Ux(j);
          
          Vector2 dUx(2*Ntri);
          dUx.Fill(0);
          for (int i = 0; i <= order; i++)
            {
              for (int j = 0; j < Ux.GetM(); j++)
                Ux(j) = Un(NumDofTri(j, i));
              
              Mlt(SeldonTrans, RhGeom, Ux, dUx);
                            
              for (int j = 0; j < Ntri; j++)
                {
                  Unode(3*j) += ValGauss0(i)*dUx(2*j);
                  Unode(3*j+1) += ValGauss0(i)*dUx(2*j+1);
                }
            }
        }
        break;
      case 1 :
        {
          Vector2 Utri(order_geom + 1), dUtri(order_geom+1), dUy(order_geom + 1);
          for (int i = 0; i <= order; i++)
            {
              Utri.Fill(0); dUtri.Fill(0);
              for (int j = 0; j <= order; j++)
                {
                  int j2 = function_basis_tri.GetNodalNumber(0, j);
                  vloc = Un(NumDofTri(j2, i));
                  
                  for (int k = 0; k <= order_geom; k++)
                    {
                      Utri(k) += GL_Geom(j, k)*vloc;
                      dUtri(k) += dGL_Geom(j, k)*vloc;
                    }
                }
              
              for (int j = 0; j <= order_geom; j++)
                {
                  int j2 = function_geom_tri.GetNodalNumber(0, j);
                  vloc = 0;
                  for (int k = 0; k < Ntri; k++)
                    vloc += RhGeom(k, 2*j2+ 1)*Un(NumDofTri(k, i));
                  
                  dUy(j) = vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                for (int k = 0; k <= order_geom; k++)
                  {
                    int p = 3*NumNodes2D_quad(j, k);
                    Unode(p) += G_Geom(i, k)*dUtri(j);
                    Unode(p+1) += G_Geom(i, k)*dUy(j);
                    Unode(p+2) += dG_Geom(i, k)*Utri(j);
                  }
            }                            
        }
        break;
      case 2 :
        {
          Vector2 Utri(order_geom + 1), dUx(order_geom+1), dUy(order_geom + 1);
          for (int i = 0; i <= order; i++)
            {
              Utri.Fill(0);
              for (int j = 0; j <= order; j++)
                {
                  int j2 = function_basis_tri.GetNodalNumber(1, j);
                  vloc = Un(NumDofTri(j2, i));
                  
                  for (int k = 0; k <= order_geom; k++)
                    Utri(k) += GL_Geom(j, k)*vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                {
                  int j2 = function_geom_tri.GetNodalNumber(1, j);
                  vloc = 0;
                  for (int k = 0; k < Ntri; k++)
                    vloc += RhGeom(k, 2*j2)*Un(NumDofTri(k, i));
                  
                  dUx(j) = vloc;

                  vloc = 0;
                  for (int k = 0; k < Ntri; k++)
                    vloc += RhGeom(k, 2*j2+ 1)*Un(NumDofTri(k, i));
                  
                  dUy(j) = vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                for (int k = 0; k <= order_geom; k++)
                  {
                    int p = 3*NumNodes2D_quad(j, k);
                    Unode(p) += G_Geom(i, k)*dUx(j);
                    Unode(p+1) += G_Geom(i, k)*dUy(j);
                    Unode(p+2) += dG_Geom(i, k)*Utri(j);
                  }
            }                            
        }
        break;
      case 3 :
        {
          Vector2 Utri(order_geom + 1), dUtri(order_geom+1), dUx(order_geom + 1);
          for (int i = 0; i <= order; i++)
            {
              Utri.Fill(0); dUtri.Fill(0);
              for (int j = 0; j <= order; j++)
                {
                  int j2 = function_basis_tri.GetNodalNumber(2, order-j);
                  vloc = Un(NumDofTri(j2, i));
                  
                  for (int k = 0; k <= order_geom; k++)
                    {
                      Utri(k) += GL_Geom(j, k)*vloc;
                      dUtri(k) += dGL_Geom(j, k)*vloc;
                    }
                }
              
              for (int j = 0; j <= order_geom; j++)
                {
                  int j2 = function_geom_tri.GetNodalNumber(2, order_geom-j);
                  vloc = 0;
                  for (int k = 0; k < Ntri; k++)
                    vloc += RhGeom(k, 2*j2)*Un(NumDofTri(k, i));
                  
                  dUx(j) = vloc;
                }
              
              for (int j = 0; j <= order_geom; j++)
                for (int k = 0; k <= order_geom; k++)
                  {
                    int p = 3*NumNodes2D_quad(j, k);
                    Unode(p) += G_Geom(i, k)*dUx(j);
                    Unode(p+1) += G_Geom(i, k)*dUtri(j);
                    Unode(p+2) += dG_Geom(i, k)*Utri(j);
                  }
            }                            
        }
        break;
      case 4 :
        {
          Vector2 Ux(this->nb_dof_tri);
          Ux.Fill(0);
          for (int i = 0; i <= order; i++)
            for (int j = 0; j < Ux.GetM(); j++)
              Ux(j) += DerivGauss1(i)*Un(NumDofTri(j, i));
          
          for (int i = 0; i < Ntri; i++)
            for (int j = 0; j < this->nb_dof_tri; j++)
              Unode(3*i+2) += ChGeom(j, i)*Ux(j);
          
          Vector2 dUx(2*Ntri);
          dUx.Fill(0);
          for (int i = 0; i <= order; i++)
            {
              for (int j = 0; j < Ux.GetM(); j++)
                Ux(j) = Un(NumDofTri(j, i));
              
              Mlt(SeldonTrans, RhGeom, Ux, dUx);
                            
              for (int j = 0; j < Ntri; j++)
                {
                  Unode(3*j) += ValGauss1(i)*dUx(2*j);
                  Unode(3*j+1) += ValGauss1(i)*dUx(2*j+1);
                }
            }
        }
        break;
      }
  }

  
  //! Vh is overwritten by L^-1 Vh or L^-T Vh where M = L L^t is the mass matrix
  template<class Vector1>
  void WedgeClassical::SolveCholeskyGen(const SeldonTranspose& TransA, Vector1& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    Vector1 Ux(Ntri);
    if (TransA.Trans())
      {
        for (int k = 0; k <= order; k++)
          {	
            Real_wp coef = 1.0/sqrt(lob_basis.Weights(k));
            for (int i = 0; i < Ntri; i++)
              Ux(i) = coef*Vh(NumDofTri(i, k));
            
            Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, Ux);
            
            for (int i = 0; i < Ntri; i++)
              Vh(NumDofTri(i, k)) = Ux(i);
          }
      }
    else
      {
        for (int k = 0; k <= order; k++)
          {	
            for (int i = 0; i < Ntri; i++)
              Ux(i) = Vh(NumDofTri(i, k));
            
            Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, Ux);
            Real_wp coef = 1.0/sqrt(lob_basis.Weights(k));
            
            for (int i = 0; i < Ntri; i++)
              Vh(NumDofTri(i, k)) = coef*Ux(i);
          }
      }

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
  void WedgeClassical::ApplyRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    int Nquad = this->Points2D_tri().GetM();
    Vector2 Ux(Ntri), dUx(Ntri), U_quad(Nquad), dU_quad(2*Nquad);
    Vh.Fill(0);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Nquad; i++)
	  {
	    U_quad(i) = 0;
	    for (int m = 0; m <= order; m++)
	      U_quad(i) += lob_basis.GradPhi(k, m)*Uh(3*NumQuadTri(i, m)+2);
	    
	    dU_quad(2*i) = Uh(3*NumQuadTri(i, k));
	    dU_quad(2*i+1) = Uh(3*NumQuadTri(i, k)+1);
	  }
	
	Mlt(ChTri, U_quad, Ux);
	Mlt(RhTri, dU_quad, dUx);
	
	for (int i = 0; i < Ntri; i++)
	  Vh(NumDofTri(i, k)) = Ux(i)+dUx(i);
      } 
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
  void WedgeClassical::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    int Nquad = this->Points2D_tri().GetM();
    Vector2 Ux(Ntri), U_quad(Nquad), dU_quad(2*Nquad);
    Vh.Fill(0);
    for (int k = 0; k <= order; k++)
      {
	for (int i = 0; i < Ntri; i++)
	  Ux(i) = Uh(NumDofTri(i, k));
	
	Mlt(SeldonTrans, ChTri, Ux, U_quad);
	Mlt(SeldonTrans, RhTri, Ux, dU_quad);
	
	for (int i = 0; i < Nquad; i++)
	  {
	    for (int m = 0; m <= order; m++)
	      Vh(3*NumQuadTri(i, m)+2) += lob_basis.GradPhi(k, m)*U_quad(i);
	    
	    Vh(3*NumQuadTri(i, k)) = dU_quad(2*i);
	    Vh(3*NumQuadTri(i, k)+1) = dU_quad(2*i+1);
	  }
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
    Vx = Rh^1 Uh, Vy = Rh^2 Uh, Vz = Rh^3 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
    where (Rh)^3_{i,j} = d/dz(\varphi_i(\xi_j) )
  */  
  template<class Vector1, class Vector2>
  void WedgeClassical
  ::ApplyRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const
  {
    int Ntri = this->nb_dof_tri;
    int Nquad = this->Points2D_tri().GetM();
    Vector2 Uz(Ntri), dUx(Ntri), dUy(Ntri), U_quad(Nquad), dU_quad(2*Nquad);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Nquad; i++)
	  {
	    U_quad(i) = 0;
	    for (int m = 0; m <= order; m++)
	      U_quad(i) += lob_basis.GradPhi(k, m)*Uh(3*NumQuadTri(i, m)+2);
	    
	    dU_quad(2*i) = Uh(3*NumQuadTri(i, k));
	    dU_quad(2*i+1) = Uh(3*NumQuadTri(i, k)+1);
	  }
	
	Mlt(ChTri, U_quad, Uz);
	Mlt(RhTriX, dU_quad, dUx);
	Mlt(RhTriY, dU_quad, dUy);
	
	for (int i = 0; i < Ntri; i++)
	  {
            Vx(NumDofTri(i, k)) = dUx(i);
            Vy(NumDofTri(i, k)) = dUy(i);
            Vz(NumDofTri(i, k)) = Uz(i);
          }
      } 
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
  void WedgeClassical::ApplyChGen(const Vector1& Uh, Vector2& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    int Nquad = this->Points2D_tri().GetM();
    Vector2 Ux(Ntri), U_quad(Nquad);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Nquad; i++)
	  U_quad(i) = Uh(NumQuadTri(i, k));
		
	Mlt(ChTri, U_quad, Ux);
		
	for (int i = 0; i < Ntri; i++)
	  Vh(NumDofTri(i, k)) = Ux(i);
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
  void WedgeClassical::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    int Nquad = NumQuadTri.GetM();
    Vector2 Ux(Ntri), U_quad(Nquad);
    for (int k = 0; k <= order; k++)
      {
	for (int i = 0; i < Ntri; i++)
	  Ux(i) = Uh(NumDofTri(i, k));
	
	Mlt(SeldonTrans, ChTri, Ux, U_quad);
	
	for (int i = 0; i < Nquad; i++)
	  Vh(NumQuadTri(i, k)) = U_quad(i);
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
  void WedgeClassical::ApplyShTransposeGen(int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    bool lobatto = (abs(1.0-ValGauss0(0)) < 1e4*epsilon_machine);
    typename Vector1::value_type vloc;
    switch (num_loc)
      {
      case 0 :
        {
          Vector1 U0(this->nb_dof_tri);          
    
          if (lobatto)
            {
              for (int i = 0; i < U0.GetM(); i++)
                U0(i) = Uh(NumDofTri(i, 0));
            }
          else
            {
              for (int i = 0; i < U0.GetM(); i++)
                {
                  U0(i) = 0;
                  for (int k = 0; k <= order; k++)
                    U0(i) += ValGauss0(k)*Uh(NumDofTri(i, k));
                }
            }
          
          if ((r == 0) || (r == order))
	    Mlt(SeldonTrans, ChTri, U0, Vh);
	  else
	    Mlt(SeldonTrans, ProjOperatorTriOrder(r), U0, Vh);
        }
        break;
      case 1 :
	if ((r != 0) && (r != order))
	  {
            Vh.Fill(0);
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  vloc = Uh(NumDofTri(node, j));
		  int ndof = NumDofs2D_quad(i, j);
                  for (int k = 0; k < Ah.GetN(); k++)
		    Vh(k) += Ah(ndof, k)*vloc;
                }
	  }
        else if (lobatto)
          {            
            Vector1 Ux((order+1)*(order+1));
            Ux.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  vloc = Uh(NumDofTri(node, j));
                  for (int k = 0; k <= order; k++)
                    Ux(NumDofs2D_quad(i, k)) += G_Gquad(j, k)*vloc;
                }
            
            Vh.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int p = NumDofs2D_quad(i, j);
                  for (int k = 0; k <= order; k++)
                    Vh(NumDofs2D_quad(k, j)) += GL_Gquad(i, k)*Ux(p);
                }
          }
        else
          {
            Vh.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  vloc = Uh(NumDofTri(node, j));
                  for (int k = 0; k <= order; k++)
                    Vh(NumDofs2D_quad(k, j)) += GL_Gquad(i, k)*vloc;
                }
          }
        break;
      case 2 :
	if ((r != 0) && (r != order))
	  {
            Vh.Fill(0);
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  vloc = Uh(NumDofTri(node, j));
		  int ndof = NumDofs2D_quad(i, j);
                  for (int k = 0; k < Ah.GetN(); k++)
		    Vh(k) += Ah(ndof, k)*vloc;
                }
	  }
        else if (lobatto)
          {            
            Vector1 Ux((order+1)*(order+1));
            Ux.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  vloc = Uh(NumDofTri(node, j));
                  for (int k = 0; k <= order; k++)
                    Ux(NumDofs2D_quad(i, k)) += G_Gquad(j, k)*vloc;
                }
            
            Vh.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int p = NumDofs2D_quad(i, j);
                  for (int k = 0; k <= order; k++)
                    Vh(NumDofs2D_quad(k, j)) += GL_Gquad(i, k)*Ux(p);
                }
          }
        else
          {
            Vh.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  vloc = Uh(NumDofTri(node, j));
                  for (int k = 0; k <= order; k++)
                    Vh(NumDofs2D_quad(k, j)) += GL_Gquad(i, k)*vloc;
                }
          }
        break;
      case 3 :
	if ((r != 0) && (r != order))
	  {
            Vh.Fill(0);
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  vloc = Uh(NumDofTri(node, j));
		  int ndof = NumDofs2D_quad(i, j);
                  for (int k = 0; k < Ah.GetN(); k++)
		    Vh(k) += Ah(ndof, k)*vloc;
                }
	  }
        else if (lobatto)
          {            
            Vector1 Ux((order+1)*(order+1));
            Ux.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  vloc = Uh(NumDofTri(node, j));
                  for (int k = 0; k <= order; k++)
                    Ux(NumDofs2D_quad(i, k)) += G_Gquad(j, k)*vloc;
                }
            
            Vh.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int p = NumDofs2D_quad(i, j);
                  for (int k = 0; k <= order; k++)
                    Vh(NumDofs2D_quad(k, j)) += GL_Gquad(i, k)*Ux(p);
                }
          }
        else
          {
            Vh.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  vloc = Uh(NumDofTri(node, j));
                  for (int k = 0; k <= order; k++)
                    Vh(NumDofs2D_quad(k, j)) += GL_Gquad(i, k)*vloc;
                }
          }
        break;
      case 4 :
        {
          Vector1 U0(this->nb_dof_tri);          
          
          if (lobatto)
            {
              for (int i = 0; i < U0.GetM(); i++)
                U0(i) = Uh(NumDofTri(i, order));
            }
          else
            {
              for (int i = 0; i < U0.GetM(); i++)
                {
                  U0(i) = 0;
                  for (int k = 0; k <= order; k++)
                    U0(i) += ValGauss1(k)*Uh(NumDofTri(i, k));
                }
            }
          
          if ((r == 0) || (r == order))
	    Mlt(SeldonTrans, ChTri, U0, Vh);
	  else
	    Mlt(SeldonTrans, ProjOperatorTriOrder(r), U0, Vh);
        }
        break;
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
  void WedgeClassical
  ::ApplyShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& Vh, int r) const
  {
    bool lobatto = (abs(1.0-ValGauss0(0)) < 1e4*epsilon_machine);
    typename Vector1::value_type vloc;
    
    switch (num_loc)
      {
      case 0 :
        {
          Vector1 U0(this->nb_dof_tri);          
          if ((r == 0) || (r == order))
	    Mlt(ChTri, Uh, U0);
	  else
	    Mlt(ProjOperatorTriOrder(r), Uh, U0);
	  
          if (lobatto)
            {
              for (int i = 0; i < U0.GetM(); i++)
                Vh(NumDofTri(i, 0)) += alpha*U0(i);
            }
          else
            {
              for (int i = 0; i < U0.GetM(); i++)
                for (int k = 0; k <= order; k++)
                  Vh(NumDofTri(i, k)) += alpha*ValGauss0(k)*U0(i);
            }
        }
        break;
      case 1 :
	if ((r != 0) && (r != order))
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
		  int ndof = NumDofs2D_quad(i, j);
                  for (int k = 0; k < Ah.GetN(); k++)
                    vloc += Ah(ndof, k)*Uh(k);
                  
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
	  }
	else if (lobatto)
          {
            Vector1 Ux((order+1)*(order+1));
            Ux.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int p = NumDofs2D_quad(i, j);
                  for (int k = 0; k <= order; k++)
                    Ux(p) += GL_Gquad(i, k)*Uh(NumDofs2D_quad(k, j));
                }
            
            vloc = 0;
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
                  for (int k = 0; k <= order; k++)
                    vloc += G_Gquad(j, k)*Ux(NumDofs2D_quad(i, k));
                  
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
          }
        else
          {
            vloc = 0;
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
                  for (int k = 0; k <= order; k++)
                    vloc += GL_Gquad(i, k)*Uh(NumDofs2D_quad(k, j));
                  
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
          }
        break;
      case 2 :
	if ((r != 0) && (r != order))
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	    for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
		  int ndof = NumDofs2D_quad(i, j);
                  for (int k = 0; k < Ah.GetN(); k++)
                    vloc += Ah(ndof, k)*Uh(k);
                  
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
	  }
	else if (lobatto)
          {
            Vector1 Ux((order+1)*(order+1));
            Ux.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int p = NumDofs2D_quad(i, j);
                  for (int k = 0; k <= order; k++)
                    Ux(p) += GL_Gquad(i, k)*Uh(NumDofs2D_quad(k, j));
                }
            
            vloc = 0;
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
                  for (int k = 0; k <= order; k++)
                    vloc += G_Gquad(j, k)*Ux(NumDofs2D_quad(i, k));
                  
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
          }
        else
          {
            vloc = 0;
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
                  for (int k = 0; k <= order; k++)
                    vloc += GL_Gquad(i, k)*Uh(NumDofs2D_quad(k, j));
                  
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
          }
        break;
      case 3 :
	if ((r != 0) && (r != order))
	  {
	    const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
	    for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
		  int ndof = NumDofs2D_quad(i, j);
                  for (int k = 0; k < Ah.GetN(); k++)
                    vloc += Ah(ndof, k)*Uh(k);
                  
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
	  }
        else if (lobatto)
          {
            Vector1 Ux((order+1)*(order+1));
            Ux.Fill(0);
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  int p = NumDofs2D_quad(i, j);
                  for (int k = 0; k <= order; k++)
                    Ux(p) += GL_Gquad(i, k)*Uh(NumDofs2D_quad(k, j));
                }
            
            vloc = 0;
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
                  for (int k = 0; k <= order; k++)
                    vloc += G_Gquad(j, k)*Ux(NumDofs2D_quad(i, k));
                  
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
          }
        else
          {
            vloc = 0;
            for (int i = 0; i <= order; i++)
              for (int j = 0; j <= order; j++)
                {
                  vloc = 0;
                  for (int k = 0; k <= order; k++)
                    vloc += GL_Gquad(i, k)*Uh(NumDofs2D_quad(k, j));
                  
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  Vh(NumDofTri(node, j)) += alpha*vloc;
                }
          }
        break;
      case 4 :
        {
          Vector1 U0(this->nb_dof_tri);          
          if ((r == 0) || (r == order))
	    Mlt(ChTri, Uh, U0);
	  else
	    Mlt(ProjOperatorTriOrder(r), Uh, U0);
	  
          if (lobatto)
            {
              for (int i = 0; i < U0.GetM(); i++)
                Vh(NumDofTri(i, order)) += alpha*U0(i);
            }
          else
            {
              for (int i = 0; i < U0.GetM(); i++)
                for (int k = 0; k <= order; k++)
                  Vh(NumDofTri(i, k)) += alpha*ValGauss1(k)*U0(i);
            }
        }
        break;
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
  void WedgeClassical
  ::ApplyNablaShTransposeGen(int num_loc, const Vector1& res, Vector2& Vh, int r) const
  {
    typename Vector1::value_type vloc, vx, vy, vz;
    bool variable_order = false;
    if ((r != 0) && (r != order_quad))
      variable_order = true;

    Vector2 feval;
    if (variable_order)
      feval.Reallocate(3*GetNbDofBoundary(num_loc));
    else
      feval.Reallocate(3*GetNbQuadBoundary(num_loc));
    
    feval.Fill(0);
    
    switch (num_loc)
      {
      case 0 :
        {
          Vector1 U0(feval.GetM()/3), dU0(2*feval.GetM()/3);
          Vector1 Uquad(this->nb_dof_tri), dUquad(this->nb_dof_tri);
          Uquad.Fill(0); dUquad.Fill(0);
          for (int i = 0; i < Uquad.GetM(); i++)
            {
              for (int k = 0; k <= order; k++)
                {
                  dUquad(i) += ValGauss0(k)*res(NumDofTri(i, k));
                  Uquad(i) += DerivGauss0(k)*res(NumDofTri(i, k));
                }
            }
          
          if (variable_order)
            {
              Copy(Uquad, U0);
              Mlt(SeldonTrans, RhNodalTri, dUquad, dU0);
            }
          else
            {
              Mlt(SeldonTrans, ChTri, Uquad, U0);
              Mlt(SeldonTrans, RhTri, dUquad, dU0);
            }
          
          for (int i = 0; i < U0.GetM(); i++)
            {
              feval(3*i) = dU0(2*i);
              feval(3*i+1) = dU0(2*i+1);
              feval(3*i+2) = U0(i);
            }
    
        }
        break;
      case 1 :
        {
          if (variable_order)
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  for (int j = 0; j <= order; j++)
                    {
                      vx = 0, vy = 0, vz = 0;
                      for (int k = 0; k < RhNodalTri.GetM(); k++)
                        {
                          vx += RhNodalTri(k, 2*node)*res(NumDofTri(k, j));
                          vy += RhNodalTri(k, 2*node+1)*res(NumDofTri(k, j));
                        }
                      
                      for (int k = 0; k <= order; k++)
                        vz += lob_basis.GradPhi(k, j)*res(NumDofTri(node, k));
                      
                      feval(3*NumDofs2D_quad(i, j)) = vx;
                      feval(3*NumDofs2D_quad(i, j)+1) = vy;
                      feval(3*NumDofs2D_quad(i, j)+2) = vz;
                    }
                }     
            }
          else
            {
              Vector1 Ux(2*(order+1)), Udof(this->nb_dof_tri);
              for (int i = 0; i <= order; i++)
                {
                  Udof.Fill(0);
                  for (int j = 0; j < Udof.GetM(); j++)
                    for (int k = 0; k <= order; k++)
                      Udof(j) += res(NumDofTri(j, k))*G_Gquad(k, i);
                  
                  function_basis_tri.ApplyNablaShTranspose(0, Udof, Ux);
                  
                  for (int j = 0; j <= order; j++)
                    {
                      feval(3*NumDofs2D_quad(j, i)) = Ux(2*j);
                      feval(3*NumDofs2D_quad(j, i)+1) = Ux(2*j+1);
                    }
                }
              
              Udof.Reallocate((order+1)*(order+1));
              Udof.Fill(0);
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int node = function_basis_tri.GetNodalNumber(0, i);
                    vloc = res(NumDofTri(node, j));
                    
                    for (int k = 0; k <= order; k++)
                      Udof(NumDofs2D_quad(i, k)) += dG_Gquad(j, k)*vloc;
                  }
              
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int p = NumDofs2D_quad(i, j);
                    for (int k = 0; k <= order; k++)
                      feval(3*NumDofs2D_quad(k, j)+2) += GL_Gquad(i, k)*Udof(p);
                  }          
            }
        }
        break;
      case 2 :
        {
          if (variable_order)
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  for (int j = 0; j <= order; j++)
                    {
                      vx = 0, vy = 0, vz = 0;
                      for (int k = 0; k < RhNodalTri.GetM(); k++)
                        {
                          vx += RhNodalTri(k, 2*node)*res(NumDofTri(k, j));
                          vy += RhNodalTri(k, 2*node+1)*res(NumDofTri(k, j));
                        }
                      
                      for (int k = 0; k <= order; k++)
                        vz += lob_basis.GradPhi(k, j)*res(NumDofTri(node, k));
                      
                      feval(3*NumDofs2D_quad(i, j)) = vx;
                      feval(3*NumDofs2D_quad(i, j)+1) = vy;
                      feval(3*NumDofs2D_quad(i, j)+2) = vz;
                    }
                }     
            }
          else
            {
              Vector1 Ux(2*(order+1)), Udof(this->nb_dof_tri);
              for (int i = 0; i <= order; i++)
                {
                  Udof.Fill(0);
                  for (int j = 0; j < Udof.GetM(); j++)
                    for (int k = 0; k <= order; k++)
                      Udof(j) += res(NumDofTri(j, k))*G_Gquad(k, i);
                  
                  function_basis_tri.ApplyNablaShTranspose(1, Udof, Ux);
                  
                  for (int j = 0; j <= order; j++)
                    {
                      feval(3*NumDofs2D_quad(j, i)) = Ux(2*j);
                      feval(3*NumDofs2D_quad(j, i)+1) = Ux(2*j+1);
                    }                            
                }
              
              Udof.Reallocate((order+1)*(order+1));
              Udof.Fill(0);
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int node = function_basis_tri.GetNodalNumber(1, i);
                    vloc = res(NumDofTri(node, j));
                    
                    for (int k = 0; k <= order; k++)
                      Udof(NumDofs2D_quad(i, k)) += dG_Gquad(j, k)*vloc;
                  }
              
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int p = NumDofs2D_quad(i, j);
                    for (int k = 0; k <= order; k++)
                      feval(3*NumDofs2D_quad(k, j)+2) += GL_Gquad(i, k)*Udof(p);
                  }
              
            }
        }
        break;
      case 3 :
        {
          if (variable_order)
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  for (int j = 0; j <= order; j++)
                    {
                      vx = 0, vy = 0, vz = 0;
                      for (int k = 0; k < RhNodalTri.GetM(); k++)
                        {
                          vx += RhNodalTri(k, 2*node)*res(NumDofTri(k, j));
                          vy += RhNodalTri(k, 2*node+1)*res(NumDofTri(k, j));
                        }
                      
                      for (int k = 0; k <= order; k++)
                        vz += lob_basis.GradPhi(k, j)*res(NumDofTri(node, k));
                      
                      feval(3*NumDofs2D_quad(i, j)) = vx;
                      feval(3*NumDofs2D_quad(i, j)+1) = vy;
                      feval(3*NumDofs2D_quad(i, j)+2) = vz;
                    }
                }     
            }
          else
            {
              Vector1 Ux(2*(order+1)), Udof(this->nb_dof_tri);
              for (int i = 0; i <= order; i++)
                {
                  Udof.Fill(0);
                  for (int j = 0; j < Udof.GetM(); j++)
                    for (int k = 0; k <= order; k++)
                      Udof(j) += res(NumDofTri(j, k))*G_Gquad(k, i);
                  
                  function_basis_tri.ApplyNablaShTranspose(2, Udof, Ux);
                  
                  for (int j = 0; j <= order; j++)
                    {
                      feval(3*NumDofs2D_quad(j, i)) = Ux(2*(order-j));
                      feval(3*NumDofs2D_quad(j, i)+1) = Ux(2*(order-j)+1);
                    }
                }
              
              Udof.Reallocate((order+1)*(order+1));
              Udof.Fill(0);
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int node = function_basis_tri.GetNodalNumber(2, order-i);
                    vloc = res(NumDofTri(node, j));
                    
                    for (int k = 0; k <= order; k++)
                      Udof(NumDofs2D_quad(i, k)) += dG_Gquad(j, k)*vloc;
                  }
              
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int p = NumDofs2D_quad(i, j);
                    for (int k = 0; k <= order; k++)
                      feval(3*NumDofs2D_quad(k, j)+2) += GL_Gquad(i, k)*Udof(p);
                  }
              
            }
        }
        break;
      case 4 :
        {
          Vector1 U0(feval.GetM()/3), dU0(2*feval.GetM()/3);
          Vector1 Uquad(this->nb_dof_tri), dUquad(this->nb_dof_tri);
          Uquad.Fill(0); dUquad.Fill(0);
          for (int i = 0; i < Uquad.GetM(); i++)
            {
              for (int k = 0; k <= order; k++)
                {
                  dUquad(i) += ValGauss1(k)*res(NumDofTri(i, k));
                  Uquad(i) += DerivGauss1(k)*res(NumDofTri(i, k));
                }
            }
          
          if (variable_order)
            {
              Copy(Uquad, U0);
              Mlt(SeldonTrans, RhNodalTri, dUquad, dU0);
            }
          else
            {
              Mlt(SeldonTrans, ChTri, Uquad, U0);
              Mlt(SeldonTrans, RhTri, dUquad, dU0);
            }
          
          for (int i = 0; i < U0.GetM(); i++)
            {
              feval(3*i) = dU0(2*i);
              feval(3*i+1) = dU0(2*i+1);
              feval(3*i+2) = U0(i);
            }
        }
      }
    
    if (variable_order)
      {
        typename Vector1::value_type val;
        Vh.Fill(0);
        if (num_loc%4 == 0)
          {
            const Matrix<Real_wp>& Ah = ProjOperatorTriOrder(r);
            for (int i = 0; i < Ah.GetM(); i++)
              {
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    val = Ah(i, j);
                    Vh(3*j) += val*feval(3*i);
                    Vh(3*j+1) += val*feval(3*i+1);
                    Vh(3*j+2) += val*feval(3*i+2);
                  }
              }
          }    
        else
          {
            const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            for (int i = 0; i < Ah.GetM(); i++)
              {
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    val = Ah(i, j);
                    Vh(3*j) += val*feval(3*i);
                    Vh(3*j+1) += val*feval(3*i+1);
                    Vh(3*j+2) += val*feval(3*i+2);
                  }
              }
          }
      }
    else
      Copy(feval, Vh);
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
  void WedgeClassical
  ::ApplyNablaShGen(const T0& alpha, int num_loc, const Vector1& Uh, Vector2& res, int r) const
  {
    typename Vector1::value_type vloc, val, vx, vy, vz;
    Vector1 feval;
    
    bool variable_order = false;
    if ((r != 0) && (r != order_quad))
      variable_order = true;
    
    if (variable_order)
      {
        feval.Reallocate(3*GetNbDofBoundary(num_loc));
        feval.Fill(0);
        if (num_loc%4 == 0)
          {
            const Matrix<Real_wp>& Ah = ProjOperatorTriOrder(r);
            for (int i = 0; i < Ah.GetM(); i++)
              {
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    val = Ah(i, j);
                    feval(3*i) += val*Uh(3*j);
                    feval(3*i+1) += val*Uh(3*j+1);
                    feval(3*i+2) += val*Uh(3*j+2);
                  }
              }
          }    
        else
          {
            const Matrix<Real_wp>& Ah = ProjOperatorQuadOrder(r);
            for (int i = 0; i < Ah.GetM(); i++)
              {
                for (int j = 0; j < Ah.GetN(); j++)
                  {
                    val = Ah(i, j);
                    feval(3*i) += val*Uh(3*j);
                    feval(3*i+1) += val*Uh(3*j+1);
                    feval(3*i+2) += val*Uh(3*j+2);
                  }
              }
          }
      }
    else
      {
        feval.Reallocate(3*GetNbQuadBoundary(num_loc));
        Copy(Uh, feval);
      }
    
    switch (num_loc)
      {
      case 0 :
        {
          Vector1 U0(feval.GetM()/3), dU0(2*feval.GetM()/3);
          for (int i = 0; i < U0.GetM(); i++)
            {
              dU0(2*i) = alpha*feval(3*i);
              dU0(2*i+1) = alpha*feval(3*i+1);
              U0(i) = alpha*feval(3*i+2);
            }
          
          Vector1 Uquad(this->nb_dof_tri), dUquad(this->nb_dof_tri);
          if (variable_order)
            {
              Copy(U0, Uquad);
              Mlt(RhNodalTri, dU0, dUquad);
            }
          else
            {
              Mlt(ChTri, U0, Uquad);
              Mlt(RhTri, dU0, dUquad);
            }
          
          for (int i = 0; i < Uquad.GetM(); i++)
            {
              for (int k = 0; k <= order; k++)
                res(NumDofTri(i, k)) += ValGauss0(k)*dUquad(i) + DerivGauss0(k)*Uquad(i);
            }
        }
        break;
      case 1 :
        {
          if (variable_order)
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  for (int j = 0; j <= order; j++)
                    {
                      vx = alpha*feval(3*NumDofs2D_quad(i, j));
                      vy = alpha*feval(3*NumDofs2D_quad(i, j)+1);
                      vz = alpha*feval(3*NumDofs2D_quad(i, j)+2);

                      for (int k = 0; k < RhNodalTri.GetM(); k++)
                        {
                          res(NumDofTri(k, j)) += RhNodalTri(k, 2*node)*vx;
                          res(NumDofTri(k, j)) += RhNodalTri(k, 2*node+1)*vy;
                        }
                      
                      for (int k = 0; k <= order; k++)
                        res(NumDofTri(node, k)) += lob_basis.GradPhi(k, j)*vz;
                      
                    }
                } 
            }
          else
            {
              Vector1 Ux(2*(order+1)), Udof(this->nb_dof_tri);
              for (int i = 0; i <= order; i++)
                {
                  for (int j = 0; j <= order; j++)
                    {
                      Ux(2*j) = alpha*feval(3*NumDofs2D_quad(j, i));
                      Ux(2*j+1) = alpha*feval(3*NumDofs2D_quad(j, i)+1);
                    }
                  
                  function_basis_tri.ComputeIntegralSurfaceGradientRef(Ux, Udof, 0);
                  for (int j = 0; j < Udof.GetM(); j++)
                    for (int k = 0; k <= order; k++)
                      res(NumDofTri(j, k)) += Udof(j)*G_Gquad(k, i);
                }
              
              Udof.Reallocate((order+1)*(order+1));
              Udof.Fill(0);
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int p = NumDofs2D_quad(i, j);
                    for (int k = 0; k <= order; k++)
                      Udof(p) += alpha*GL_Gquad(i, k)*feval(3*NumDofs2D_quad(k, j)+2);
                  }
              
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    vloc = 0;
                    for (int k = 0; k <= order; k++)
                      vloc += dG_Gquad(j, k)*Udof(NumDofs2D_quad(i, k));
                    
                    int node = function_basis_tri.GetNodalNumber(0, i);
                    res(NumDofTri(node, j)) += vloc;
                  }
            }
        }
        break;
      case 2 :
        {
          if (variable_order)
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  for (int j = 0; j <= order; j++)
                    {
                      vx = alpha*feval(3*NumDofs2D_quad(i, j));
                      vy = alpha*feval(3*NumDofs2D_quad(i, j)+1);
                      vz = alpha*feval(3*NumDofs2D_quad(i, j)+2);

                      for (int k = 0; k < RhNodalTri.GetM(); k++)
                        {
                          res(NumDofTri(k, j)) += RhNodalTri(k, 2*node)*vx;
                          res(NumDofTri(k, j)) += RhNodalTri(k, 2*node+1)*vy;
                        }
                      
                      for (int k = 0; k <= order; k++)
                        res(NumDofTri(node, k)) += lob_basis.GradPhi(k, j)*vz;
                      
                    }
                } 
            }
          else
            {
              Vector1 Ux(2*(order+1)), Udof(this->nb_dof_tri);
              for (int i = 0; i <= order; i++)
                {
                  for (int j = 0; j <= order; j++)
                    {
                      Ux(2*j) = alpha*feval(3*NumDofs2D_quad(j, i));
                      Ux(2*j+1) = alpha*feval(3*NumDofs2D_quad(j, i)+1);
                    }
                  
                  function_basis_tri.ComputeIntegralSurfaceGradientRef(Ux, Udof, 1);
                  for (int j = 0; j < Udof.GetM(); j++)
                    for (int k = 0; k <= order; k++)
                      res(NumDofTri(j, k)) += Udof(j)*G_Gquad(k, i);
                }
              
              Udof.Reallocate((order+1)*(order+1));
              Udof.Fill(0);
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int p = NumDofs2D_quad(i, j);
                    for (int k = 0; k <= order; k++)
                      Udof(p) += alpha*GL_Gquad(i, k)*feval(3*NumDofs2D_quad(k, j)+2);
                  }
              
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    vloc = 0;
                    for (int k = 0; k <= order; k++)
                      vloc += dG_Gquad(j, k)*Udof(NumDofs2D_quad(i, k));
                    
                    int node = function_basis_tri.GetNodalNumber(1, i);
                    res(NumDofTri(node, j)) += vloc;
                  }
            }
        }
        break;
      case 3 :
        {
          if (variable_order)
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  for (int j = 0; j <= order; j++)
                    {
                      vx = alpha*feval(3*NumDofs2D_quad(i, j));
                      vy = alpha*feval(3*NumDofs2D_quad(i, j)+1);
                      vz = alpha*feval(3*NumDofs2D_quad(i, j)+2);
                      
                      for (int k = 0; k < RhNodalTri.GetM(); k++)
                        {
                          res(NumDofTri(k, j)) += RhNodalTri(k, 2*node)*vx;
                          res(NumDofTri(k, j)) += RhNodalTri(k, 2*node+1)*vy;
                        }
                      
                      for (int k = 0; k <= order; k++)
                        res(NumDofTri(node, k)) += lob_basis.GradPhi(k, j)*vz;
                      
                    }
                } 
            }
          else
            {
              Vector1 Ux(2*(order+1)), Udof(this->nb_dof_tri);
              for (int i = 0; i <= order; i++)
                {
                  for (int j = 0; j <= order; j++)
                    {
                      Ux(2*(order-j)) = alpha*feval(3*NumDofs2D_quad(j, i));
                      Ux(2*(order-j)+1) = alpha*feval(3*NumDofs2D_quad(j, i)+1);
                    }
                  
                  function_basis_tri.ComputeIntegralSurfaceGradientRef(Ux, Udof, 2);
                  for (int j = 0; j < Udof.GetM(); j++)
                    for (int k = 0; k <= order; k++)
                      res(NumDofTri(j, k)) += Udof(j)*G_Gquad(k, i);
                }
              
              Udof.Reallocate((order+1)*(order+1));
              Udof.Fill(0);
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    int p = NumDofs2D_quad(i, j);
                    for (int k = 0; k <= order; k++)
                      Udof(p) += alpha*GL_Gquad(i, k)*feval(3*NumDofs2D_quad(k, j)+2);
                  }
              
              for (int i = 0; i <= order; i++)
                for (int j = 0; j <= order; j++)
                  {
                    vloc = 0;
                    for (int k = 0; k <= order; k++)
                      vloc += dG_Gquad(j, k)*Udof(NumDofs2D_quad(i, k));
                    
                    int node = function_basis_tri.GetNodalNumber(2, order-i);
                    res(NumDofTri(node, j)) += vloc;
                  }
            }
        }
        break;
      case 4 :
        {
          Vector1 U0(feval.GetM()/3), dU0(2*feval.GetM()/3);
          for (int i = 0; i < U0.GetM(); i++)
            {
              dU0(2*i) = alpha*feval(3*i);
              dU0(2*i+1) = alpha*feval(3*i+1);
              U0(i) = alpha*feval(3*i+2);
            }
          
          Vector1 Uquad(this->nb_dof_tri), dUquad(this->nb_dof_tri);
          if (variable_order)
            {
              Copy(U0, Uquad);
              Mlt(RhNodalTri, dU0, dUquad);
            }
          else
            {
              Mlt(ChTri, U0, Uquad);
              Mlt(RhTri, dU0, dUquad);
            }
          
          for (int i = 0; i < Uquad.GetM(); i++)
            {
              for (int k = 0; k <= order; k++)
                res(NumDofTri(i, k)) += ValGauss1(k)*dUquad(i) + DerivGauss1(k)*Uquad(i);
            }
        }
      }
  }


  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void WedgeClassical::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    Vector2 Ux(Ntri), dUx(Ntri), U_quad(Ntri), dU_quad(2*Ntri);
    Vh.Fill(0);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Ntri; i++)
	  {
	    U_quad(i) = 0;
	    for (int m = 0; m <= order; m++)
	      U_quad(i) += GradPhiWeight(k, m)*Uh(3*NumDofTri(i, m)+2);
	    
	    dU_quad(2*i) = lob_basis.Weights(k)*Uh(3*NumDofTri(i, k));
	    dU_quad(2*i+1) = lob_basis.Weights(k)*Uh(3*NumDofTri(i, k)+1);
	  }
	
	Mlt(ChConstTri, U_quad, Ux);
	Mlt(RhConstTri, dU_quad, dUx);
	
	for (int i = 0; i < Ntri; i++)
	  Vh(NumDofTri(i, k)) = Ux(i)+dUx(i);
      } 
  }


  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    (Rh^3)_{i, j} = \int_K \varphi_j d/dz(\varphi_i) dx
    then Vx = Rh^1 Uh, Vy = Rh^2 Uh, Vz = Rh^3 Uh
  */
  template<class Vector1, class Vector2>
  void WedgeClassical
  ::ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy, Vector2& Vz) const
  {
    int Ntri = this->nb_dof_tri;
    Vector2 Uz(Ntri), dUx(Ntri), dUy(Ntri), U_quad(Ntri), dU_quad(2*Ntri);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Ntri; i++)
	  {
	    U_quad(i) = 0;
	    for (int m = 0; m <= order; m++)
	      U_quad(i) += GradPhiWeight(k, m)*Uh(3*NumDofTri(i, m)+2);
	    
	    dU_quad(2*i) = lob_basis.Weights(k)*Uh(3*NumDofTri(i, k));
	    dU_quad(2*i+1) = lob_basis.Weights(k)*Uh(3*NumDofTri(i, k)+1);
	  }
	
	Mlt(ChConstTri, U_quad, Uz);
	Mlt(RhConstTriX, dU_quad, dUx);
	Mlt(RhConstTriY, dU_quad, dUy);
	
	for (int i = 0; i < Ntri; i++)
	  {
            Vx(NumDofTri(i, k)) = dUx(i);
            Vy(NumDofTri(i, k)) = dUy(i);
            Vz(NumDofTri(i, k)) = Uz(i);
          }
      } 
  }
  
  
  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void WedgeClassical::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    Vector2 Ux(Ntri), U_quad(Ntri), dU_quad(2*Ntri);
    Vh.Fill(0);
    for (int k = 0; k <= order; k++)
      {
	for (int i = 0; i < Ntri; i++)
	  Ux(i) = Uh(NumDofTri(i, k));
	
	Mlt(SeldonTrans, ChConstTri, Ux, U_quad);
	Mlt(SeldonTrans, RhConstTri, Ux, dU_quad);
	
	for (int i = 0; i < Ntri; i++)
	  {
	    for (int m = 0; m <= order; m++)
	      Vh(3*NumDofTri(i, m)+2) += GradPhiWeight(k, m)*U_quad(i);
	    
	    Vh(3*NumDofTri(i, k)) = lob_basis.Weights(k)*dU_quad(2*i);
	    Vh(3*NumDofTri(i, k)+1) = lob_basis.Weights(k)*dU_quad(2*i+1);
	  }
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
  void WedgeClassical::AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    for (int i2 = 0; i2 < this->nb_dof_tri; i2++)
      for (int j2 = 0; j2 < this->nb_dof_tri; j2++)
        {
          T val = mass*mass_matrix(i2, j2);
          for (int k = 0; k <= order; k++)
            {
              int i = NumDofTri(i2, k);
              int j = NumDofTri(j2, k);
	      A.AddInteraction(m+i, n+j, val*weights1d_z(k));
            }
        }
  }
    
  
  //! we add elementary matrix for a constant jacobian
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass, C, D, E coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p + \sum_{p,q}  C(p, q) S^{p,q}
    where S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
    where R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void WedgeClassical
  ::AddConstantElemMatrixGen(int m, int n, const T& mass, const TinyMatrix<T, Prop, 3, 3>& C,
			     const TinyVector<T, 3>& D, const TinyVector<T, 3>& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    T vloc, vloc2, vloc3, vloc4;
    for (int i2 = 0; i2 < this->nb_dof_tri; i2++)
      for (int j2 = 0; j2 < this->nb_dof_tri; j2++)
        {
          if (!null_term(0))
	    {
	      vloc = mass*mass_matrix(i2, j2);
	      for (int k = 0; k <= order; k++)
		{
		  int i = NumDofTri(i2, k);
		  int j = NumDofTri(j2, k);
		  A.AddInteraction(m+i, n+j, vloc*weights1d_z(k));
		}
	    }
	  
          if (!null_term(2))
	    {
	      vloc = D(0)*const_grad_matrix(i2, 2*j2) + D(1)*const_grad_matrix(i2, 2*j2+1);
	      for (int k = 0; k <= order; k++)
		{
		  int i = NumDofTri(i2, k);
		  int j = NumDofTri(j2, k);
		  A.AddInteraction(i+m, n+j, vloc*weights1d_z(k));
		  
		  vloc2 = D(2)*mass_matrix(i2, j2);
		  
		  for (int l = 0; l <= order; l++)
		    {
		      j = NumDofTri(j2, l);
		      A.AddInteraction(m+i, n+j, vloc2*weights1d_z(l)*lob_basis.GradPhi(k, l));
		    }
		}
	    }

          if (!null_term(3))
	    {
	      vloc = E(0)*const_grad_matrix(j2, 2*i2) + E(1)*const_grad_matrix(j2, 2*i2+1);
	      for (int k = 0; k <= order; k++)
		{
		  int i = NumDofTri(i2, k);
		  int j = NumDofTri(j2, k);
		  A.AddInteraction(i+m, n+j, vloc*weights1d_z(k));
		  
		  vloc2 = E(2)*weights1d_z(k)*mass_matrix(i2, j2);
		  
		  for (int l = 0; l <= order; l++)
		    {
		      j = NumDofTri(j2, l);
		      A.AddInteraction(m+i, n+j, vloc2*lob_basis.GradPhi(l, k));
		    }
		}
	    }
	  
	  if (null_term(1))
	    continue;

          vloc = C(0, 0)*const_stiff_matrix(2*j2, 2*i2) + C(0, 1)*const_stiff_matrix(2*j2+1, 2*i2)
            + C(1, 0)*const_stiff_matrix(2*j2, 2*i2+1)
            + C(1, 1)*const_stiff_matrix(2*j2+1, 2*i2+1);
          
          vloc2 = C(2, 0)*const_grad_matrix(j2, 2*i2) + C(2, 1)*const_grad_matrix(j2, 2*i2+1);
          vloc3 = C(0, 2)*const_grad_matrix(i2, 2*j2) + C(1, 2)*const_grad_matrix(i2, 2*j2+1);
          
          for (int k = 0; k <= order; k++)
            {
              int i = NumDofTri(i2, k);
              int j = NumDofTri(j2, k);
	      A.AddInteraction(i+m, n+j, vloc*weights1d_z(k));
              
              vloc4 = C(2, 2)*mass_matrix(i2, j2);
              
              for (int l = 0; l <= order; l++)
                {
                  j = NumDofTri(j2, l);
		  A.AddInteraction(m+i, n+j, vloc2*weights1d_z(l)*lob_basis.GradPhi(k, l)
				   + vloc3*weights1d_z(k)*lob_basis.GradPhi(l, k)
				   + vloc4*stiffness_matrix1D(l, k));
		}
	    }
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
  void WedgeClassical
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
  void WedgeClassical::AddVariableMassMatrixGen(int off_row, int off_col,
						const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    ElementReference<Dimension3, 1>::AddVariableMassMatrixOpt(off_row, off_col, A, mat);
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point local point where functions are evaluated
    \param[out] res values of basis functions on point
  */
  void WedgeClassical::ComputeValuesPhiRef(const R3& point, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    VectReal_wp Pz(r+1), phi2d;
    for (int i = 0; i <= r; i++)
      Pz(i) = lob_basis.EvaluatePhi(i, point(2));
    
    R2 pt_tri(point(0), point(1));
    function_basis_tri.ComputeValuesPhiNodalRef(pt_tri, phi2d);   
    
    // final result
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node2d = NumDofs2D_tri(i, j);
	  for (int k = 0; k <= r; k++)
	    phi(NumDofs3D(i, j, k)) = phi2d(node2d)*Pz(k);
	}
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on point
  */
  void WedgeClassical::ComputeGradientPhiRef(const R3& point, VectR3& dphi) const
  {
    dphi.Reallocate(nb_dof_loc);
    int r = order;
    VectReal_wp Pz(r+1), dPz(r+1), phi2d; VectR2 grad_phi2d;
    for (int i = 0; i <= r; i++)
      {
	Pz(i) = lob_basis.EvaluatePhi(i, point(2));
	dPz(i) = lob_basis.EvaluatePhiGrad(i, point(2));
      }
    
    R2 pt_tri(point(0), point(1));
    function_basis_tri.ComputeValuesPhiNodalRef(pt_tri, phi2d);   
    function_basis_tri.ComputeGradientPhiNodalRef(pt_tri, grad_phi2d);   
    
    // final result
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  int node2d = NumDofs2D_tri(i, j);
	  for (int k = 0; k <= r; k++)
	    dphi(NumDofs3D(i, j, k)).Init(grad_phi2d(node2d)(0)*Pz(k),
                                          grad_phi2d(node2d)(1)*Pz(k), phi2d(node2d)*dPz(k));
	}    
  }
  
  
  //! retrieves values of a single basis functions on all quadrature points
  void WedgeClassical::GetValueSinglePhiQuadrature(int k, VectReal_wp& phi) const
  {
    int k1 = CoordinateDof(k, 0);
    int k2 = CoordinateDof(k, 1);
    phi.Zero();
    for (int i = 0; i < NumQuadTri.GetM(); i++)
      phi(NumQuadTri(i, k2)) = ChTri(k1, i); 
  }
  

  //! retrieves values of a single basis functions on all quadrature points
  void WedgeClassical
  ::GetGradientSinglePhiQuadrature(int k, VectReal_wp& phi, VectR3& grad_phi) const
  {
    int k1 = CoordinateDof(k, 0);
    int k2 = CoordinateDof(k, 1);
    phi.Zero(); grad_phi.Zero();
    for (int i = 0; i < NumQuadTri.GetM(); i++)
      {
	int n = NumQuadTri(i, k2);
	phi(n) = ChTri(k1, i);
	grad_phi(n)(0) = RhTri(k1, 2*i);
	grad_phi(n)(1) = RhTri(k1, 2*i+1);
	for (int j = 0; j <= order; j++)
	  {
	    n = NumQuadTri(i, j);
	    grad_phi(n)(2) = lob_basis.GradPhi(k2, j)*ChTri(k1, i);
	  }
      }
  }


  //! Retrieving values of basis functions on a quadrature point
  /*!
    \param[in] k quadrature point number
    \param[out] phi values of basis functions on quadrature point k
   */  
  void WedgeClassical::GetValuePhiOnQuadraturePoint(int k, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc); phi.Fill(0);
    int k1 = CoordinateQuad(k, 0);
    int k2 = CoordinateQuad(k, 1);
    int num_loc = CoordinateQuad(k, 2);
    if (num_loc < 0)
      {
        for (int i = 0; i < nb_dof_tri; i++)
          phi(NumDofTri(i, k2)) = ChTri(i, k1);
      }
    else
      {
        switch (num_loc)
          {
          case 0:
            {
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  phi(NumDofTri(i, j)) = ChTri(i, k1)*ValGauss0(j);
              
            }
            break;
          case 1 :
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(0, i);
                  for (int j = 0; j <= order; j++)
                    phi(NumDofTri(node, j)) = GL_Gquad(i, k1)*G_Gquad(j, k2);
                }              
            }
            break;
          case 2 :
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(1, i);
                  for (int j = 0; j <= order; j++)
                    phi(NumDofTri(node, j)) = GL_Gquad(i, k1)*G_Gquad(j, k2);
                }
            }
            break;
          case 3 :
            {
              for (int i = 0; i <= order; i++)
                {
                  int node = function_basis_tri.GetNodalNumber(2, order-i);
                  for (int j = 0; j <= order; j++)
                    phi(NumDofTri(node, j)) = GL_Gquad(i, k1)*G_Gquad(j, k2);
                }
            }
            break;
          case 4 :
            {
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  phi(NumDofTri(i, j)) = ChTri(i, k1)*ValGauss1(j);
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
  void WedgeClassical::GetGradientPhiOnQuadraturePoint(int k, VectR3& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc); FillZero(grad_phi);
    int k1 = CoordinateQuad(k, 0);
    int k2 = CoordinateQuad(k, 1);
    int num_loc = CoordinateQuad(k, 2);
    if (num_loc < 0)
      {
        for (int i = 0; i < this->nb_dof_tri; i++)
          {
            grad_phi(NumDofTri(i, k2))(0) = RhTri(i, 2*k1);
            grad_phi(NumDofTri(i, k2))(1) = RhTri(i, 2*k1+1);
            for (int j = 0; j <= order; j++)
              grad_phi(NumDofTri(i, j))(2) = lob_basis.GradPhi(j, k2)*ChTri(i, k1);
          }
      }
    else
      {
        switch (num_loc)
          {
          case 0:
            {
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  grad_phi(NumDofTri(i, j)).Init(RhTri(i, 2*k1)*ValGauss0(j),
                                                 RhTri(i, 2*k1+1)*ValGauss0(j),
                                                 ChTri(i, k1)*DerivGauss0(j));
              
            }
            break;
          case 1 :
            {
              VectReal_wp phi(function_basis_tri.GetNbDof());
              VectR2 dphi(function_basis_tri.GetNbDof());
              int n = this->Points2D_tri().GetM() + k1;
              function_basis_tri.GetValuePhiOnQuadraturePoint(n, phi);
              function_basis_tri.GetGradientPhiOnQuadraturePoint(n, dphi);
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  grad_phi(NumDofTri(i, j)).Init(dphi(i)(0)*G_Gquad(j, k2),
                                                 dphi(i)(1)*G_Gquad(j, k2),
                                                 phi(i)*dG_Gquad(j, k2));
            }
            break;
          case 2 :
            {
              VectReal_wp phi(function_basis_tri.GetNbDof());
              VectR2 dphi(function_basis_tri.GetNbDof());
              int n = this->Points2D_tri().GetM() + (order+1) + k1;
              function_basis_tri.GetValuePhiOnQuadraturePoint(n, phi);
              function_basis_tri.GetGradientPhiOnQuadraturePoint(n, dphi);
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  grad_phi(NumDofTri(i, j)).Init(dphi(i)(0)*G_Gquad(j, k2),
                                                 dphi(i)(1)*G_Gquad(j, k2),
                                                 phi(i)*dG_Gquad(j, k2));
            }
            break;
          case 3 :
            {
              VectReal_wp phi(function_basis_tri.GetNbDof());
              VectR2 dphi(function_basis_tri.GetNbDof());
              int n = this->Points2D_tri().GetM() + 2*(order+1) + order-k1;
              function_basis_tri.GetValuePhiOnQuadraturePoint(n, phi);
              function_basis_tri.GetGradientPhiOnQuadraturePoint(n, dphi);
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  grad_phi(NumDofTri(i, j)).Init(dphi(i)(0)*G_Gquad(j, k2),
                                                 dphi(i)(1)*G_Gquad(j, k2),
                                                 phi(i)*dG_Gquad(j, k2));
            }
            break;
          case 4 :
            {
              for (int i = 0; i < this->nb_dof_tri; i++)
                for (int j = 0; j <= order; j++)
                  grad_phi(NumDofTri(i, j)).Init(RhTri(i, 2*k1)*ValGauss1(j),
                                                 RhTri(i, 2*k1+1)*ValGauss1(j),
                                                 ChTri(i, k1)*DerivGauss1(j));
            }
            break;
          }
      }
  }
  
  
  //! computation of projection between finite element of different orders  
  void WedgeClassical::
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
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorQuadOrder(r).Reallocate((order+1)*(order+1), Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    for (int i1 = 0; i1 <= order; i1++)
	      for (int i2 = 0; i2 <= order; i2++)
		{
		  int k = NumDofs2D_quad(i1, i2);
		  Real_wp yt = Pts(i)(j)(1), xt = Pts(i)(j)(0);
		  ProjOperatorQuadOrder(r)(k, j) = lob.EvaluatePhi(i1, xt)*lob_basis.EvaluatePhi(i2, yt);
		}	    
	  }
      }
  }
    

  //! Vh is overwritten by M Vh where M is the mass matrix
  void WedgeClassical::MltMassMatrix(VectReal_wp& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    VectReal_wp Ux(Ntri);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Ntri; i++)
	  Ux(i) = Vh(NumDofTri(i, k));
        
        Seldon::MltCholesky(SeldonTrans, mass_matrix_chol, Ux);
        Seldon::MltCholesky(SeldonNoTrans, mass_matrix_chol, Ux);
        
        for (int i = 0; i < Ntri; i++)
          Vh(NumDofTri(i, k)) = lob_basis.Weights(k)*Ux(i);
      } 
  }
  
  
  //! Vh is overwritten by M^-1 Vh where M is the mass matrix
  void WedgeClassical::SolveMassMatrix(VectReal_wp& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    VectReal_wp Ux(Ntri);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Ntri; i++)
	  Ux(i) = Vh(NumDofTri(i, k));
        
        Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, Ux);
        Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, Ux);
        Real_wp coef = 1.0/lob_basis.Weights(k);
        
        for (int i = 0; i < Ntri; i++)
          Vh(NumDofTri(i, k)) = coef*Ux(i);
      } 
  }
  

  //! Vh is overwritten by M^-1 Vh where M is the mass matrix
  void WedgeClassical::SolveMassMatrix(VectComplex_wp& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    VectReal_wp Ux(Ntri), Uy;
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Ntri; i++)
	  {
            Ux(i) = real(Vh(NumDofTri(i, k)));
            Uy(i) = imag(Vh(NumDofTri(i, k)));
          }
        
        Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, Ux);
        Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, Ux);

        Seldon::SolveCholesky(SeldonNoTrans, mass_matrix_chol, Uy);
        Seldon::SolveCholesky(SeldonTrans, mass_matrix_chol, Uy);
        
        Real_wp coef = 1.0/lob_basis.Weights(k);
        
        for (int i = 0; i < Ntri; i++)
          Vh(NumDofTri(i, k)) = coef*Complex_wp(Ux(i), Uy(i));
      } 
  }
  
  
  //! Vh is overwritten by M Vh where M is the mass matrix
  void WedgeClassical::MltMassMatrix(VectComplex_wp& Vh) const
  {
    int Ntri = this->nb_dof_tri;
    VectReal_wp Ux(Ntri), Uy(Ntri);
    for (int k = 0; k <= order; k++)
      {	
	for (int i = 0; i < Ntri; i++)
	  {
            Ux(i) = real(Vh(NumDofTri(i, k)));
            Uy(i) = imag(Vh(NumDofTri(i, k)));
          }
        
        Seldon::MltCholesky(SeldonTrans, mass_matrix_chol, Ux);
        Seldon::MltCholesky(SeldonNoTrans, mass_matrix_chol, Ux);

        Seldon::MltCholesky(SeldonTrans, mass_matrix_chol, Uy);
        Seldon::MltCholesky(SeldonNoTrans, mass_matrix_chol, Uy);
        
        for (int i = 0; i < Ntri; i++)
          Vh(NumDofTri(i, k)) = lob_basis.Weights(k)*Complex_wp(Ux(i), Uy(i));
      } 
  }


  //! displaying informations about WedgeClassical
  ostream& operator <<(ostream& out, const WedgeClassical& e)
  {
    out<<static_cast<const WedgeReference<1>&>(e);
    out<<"Number dof over the Wedge "<<e.nb_dof_loc<<endl;
    return out;
  }

  
  //! default constructor
  WedgeLobatto::WedgeLobatto() : WedgeClassical()
  {
    type_integration_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO;
  }
  
  
  //! construction of finite element
  void WedgeLobatto::
  ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			 int rtri, int r_quad, int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    WedgeClassical::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                           r, r, TriangleQuadrature::QUADRATURE_GAUSS,
                                           Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
  }

} // end namespace
  
#define MONTJOIE_FILE_WEDGE_CLASSICAL_CXX
#endif
