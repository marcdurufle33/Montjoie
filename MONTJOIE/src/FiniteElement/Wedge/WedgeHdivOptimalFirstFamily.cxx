#ifndef MONTJOIE_FILE_WEDGE_HDIV_OPTIMAL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  WedgeHdivOptimalFirstFamily::WedgeHdivOptimalFirstFamily() : WedgeReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void WedgeHdivOptimalFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return WedgeReference<3>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, (order+1)*(order+1));
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
    
    // dofs inside the wedge
    nmap.SetNbDofWedge(order, nb_dof_loc-nb_dof_boundaries);
    
    // rotation of dofs on faces
    Matrix<int> FacesDof_Rotation_Tri;
    MeshNumbering<Dimension3>::
      GetRotationTriangularFace(this->PointsDof2D_tri(), FacesDof_Rotation_Tri);
    
    nmap.SetFacesDofRotationTri(this->order, FacesDof_Rotation_Tri);
    
    Matrix<bool> OppositeSigneDof(6, order*(order+1)/2);
    OppositeSigneDof.Fill(false);
    for (int i = 3; i < 6; i++)
      for (int j = 0; j < order*(order+1)/2; j++)
	OppositeSigneDof(i, j) = true;
    
    nmap.SetSignDofRotationTri(order, OppositeSigneDof);
    
    Matrix<int> FacesDof_Rotation_Quad;
    MeshNumbering<Dimension3>::
      GetRotationQuadrilateralFace(this->PointsDof2D_quad(), FacesDof_Rotation_Quad);

    nmap.SetFacesDofRotationQuad(this->order, FacesDof_Rotation_Quad);
    
    OppositeSigneDof.Reallocate(8, (order+1)*(order+1));
    OppositeSigneDof.Fill(false);
    for (int i = 4; i < 8; i++)
      for (int j = 0; j < (order+1)*(order+1); j++)
	OppositeSigneDof(i, j) = true;
    
    nmap.SetSignDofRotationQuad(order, OppositeSigneDof);
  }


  size_t WedgeHdivOptimalFirstFamily::GetMemorySize() const
  {
    size_t taille = WedgeReference<3>::GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    taille += Fb_tri.GetMemorySize();
    taille += InverseBasisVDM.GetMemorySize() + LegendrePolynom.GetMemorySize();
    taille += Fb_triH1.GetMemorySize();
    taille += NumDofsXY.GetMemorySize() + NumDofsZ.GetMemorySize();
    taille += lob_lobatto.GetMemorySize() + lob_gauss.GetMemorySize();
    return taille;
  }

  
  //! construction of nearly orthogonal functions
  void WedgeHdivOptimalFirstFamily
  ::ConstructOrthogonalBasis(int r, int rgeom, int rquad, int type_quad)
  {
    if (rquad <= r+1)
      rquad = r+2;
    
    WedgeReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad,
					      r, r,
					      TriangleQuadrature::QUADRATURE_GAUSS,
					      Globatto<Real_wp>::QUADRATURE_GAUSS,
					      Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    VectReal_wp points1d, weights1d;
    ComputeGaussLegendre(points1d, weights1d, order);
    
    this->SetPoints1D(points1d);
    this->SetWeights1D(weights1d);

    nb_dof_tri = order*(order+1)/2;
    nb_dof_quad = (order+1)*(order+1);
    nb_dof_boundaries = 2*nb_dof_tri + 3*nb_dof_quad;
    nb_dof_loc = 3*order*(order+1)*(order+1)/2 + order*(order+1) + 3*(order+1)*(order+1);

    GetJacobiPolynomial(LegendrePolynom, r+2, Real_wp(0), Real_wp(0));
    
    // weight for orthonormalization
    InvWeightBasisRr.Reallocate(nb_dof_loc);
    InvWeightBasisRr.Fill(1);

    InvWeightBasisRr.Fill(1.0);
    VectReal_wp Scale(nb_dof_loc); Scale.Fill(0);
    VectR3 phi(nb_dof_loc);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
        for (int i = 0; i < nb_dof_loc; i++)
          Scale(i) += DotProd(phi(i), phi(i))*this->WeightsND(k);
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      InvWeightBasisRr(i) = 1.0/sqrt(Scale(i));
    
  }
  
  
  //! construction of finite element
  void WedgeHdivOptimalFirstFamily
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                         int rsurf_tri, int rsurf_quad,
                         int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    ConstructOrthogonalBasis(r, rgeom, rquad, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    this->nb_points_dof_inside = this->nb_dof_loc;
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    ConstructHdivElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace();

    TriangleClassical* Fb_tri = new TriangleClassical();
    Fb_tri->ConstructFiniteElement(order);
    element_tri_surf = Fb_tri;
    
    QuadrangleGauss* Fb_quad = new QuadrangleGauss();
    Fb_quad->ConstructFiniteElement(order);
    element_quad_surf = Fb_quad;
    
    this->num_dof_points_surf.Reallocate(5);
    for (int n = 0; n < 5; n++)
      {
	this->num_dof_points_surf(n).Reallocate(this->GetNbDofBoundary(n));
	for (int i = 0; i < this->GetNbDofBoundary(n); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    
  }
  
  
  //! construction of basis functions
  void WedgeHdivOptimalFirstFamily::ConstructFunctions()
  {
    Fb_tri.ConstructFiniteElement(order+1);
    if (order > 1)
      Fb_triH1.ConstructFiniteElement(order-1);
    
    VectReal_wp points_edge; VectR2 points_tri;
    TriangleGeomReference::ConstructLobattoPoints(order+2, 0,
						  points_edge, points_tri);

    VectR2 points_dof2d_tri, points_dof2d_quad;
    VectR3 points_dof3d;    
    int Ntri = order*(order+1)/2;
    points_dof2d_tri.Reallocate(Ntri);
    for (int i = 0; i < Ntri; i++)
      points_dof2d_tri(i) = points_tri(3*(order+2) + i);
    
    if (order > 1)
      {
	InverseBasisVDM.Reallocate(Ntri, Ntri);
	InverseBasisVDM.Fill(Real_wp(0));
	
	// evaluation of orthogonal polynomials	
	VectReal_wp phi;
	for (int k = 0; k < Ntri; k++)
	  {
	    Fb_triH1.ComputeValuesPhiOrthoRef(order-1, Fb_triH1.GetNumOrtho2D(),
					      Fb_triH1.GetInverseWeightPolynomial(),
					      points_dof2d_tri(k), phi);
	    
	    for (int j = 0; j < Ntri; j++)
	      InverseBasisVDM(j, k) = phi(j);
	  }
	
	GetInverse(InverseBasisVDM);
      }
    
    lob_gauss.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    lob_lobatto.ConstructQuadrature(order+2, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    points_dof2d_quad.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	points_dof2d_quad(i*(order+1) + j).Init(lob_gauss.Points(i), lob_gauss.Points(j));
      
    NumDofsXY.Reallocate(Fb_tri.GetNbDof(), order+1);
    NumDofsZ.Reallocate(Ntri, order+3);
    NumDofsXY.Fill(-1);
    NumDofsZ.Fill(-1);
    
    // dofs of the first triangular face
    for (int i = 0; i < Ntri; i++)
      NumDofsZ(i, 0) = i;
    
    // dofs of quadrilateral faces
    int Nquad = (order+1)*(order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  NumDofsXY(i, j) = Ntri + i*(order+1) + j;
	  NumDofsXY(order+1+i, j) = Ntri + Nquad + i*(order+1) + j;
	  NumDofsXY(2*(order+1)+i, j) = Ntri + 2*Nquad + (order-i)*(order+1) + j;
	}
    
    // dofs of the last triangular face
    for (int i = 0; i < Ntri; i++)
      NumDofsZ(i, order+2) = Ntri + 3*Nquad + i;
    
    // dofs inside
    int nb = 2*Ntri + 3*Nquad;
    for (int i = 3*(order+1); i < Fb_tri.GetNbDof(); i++) 
      for (int j = 0; j <= order; j++)
	NumDofsXY(i, j) = nb++;
    
    for (int i = 0; i < Ntri; i++)
      for (int j = 1; j <= order+1; j++)
	NumDofsZ(i, j) = nb++;
    
    points_dof3d.Reallocate(nb_dof_loc);
    for (int i = 0; i < Ntri; i++)
      for (int j = 0; j <= order+2; j++)
	points_dof3d(NumDofsZ(i, j)).Init(points_dof2d_tri(i)(0), points_dof2d_tri(i)(1),
					  lob_lobatto.Points(j));
    
    for (int i = 0; i < Fb_tri.GetNbDof(); i++) 
      for (int j = 0; j <= order; j++)
	points_dof3d(NumDofsXY(i, j)).Init(Fb_tri.PointsDofND(i)(0), Fb_tri.PointsDofND(i)(1),
					   lob_gauss.Points(j));

    this->SetPointsDof2D_tri(points_dof2d_tri);
    this->SetPointsDof2D_quad(points_dof2d_quad);
    this->SetPointsDofND(points_dof3d);            
  }
   
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void WedgeHdivOptimalFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    // values of triangular basis functions are computed
    R2 pt2D(point_loc(0), point_loc(1));
    VectR2 psi; VectReal_wp phi, phi_ortho;
    if (order > 1)
      {
	Fb_triH1.ComputeValuesPhiOrthoRef(order-1, Fb_triH1.GetNumOrtho2D(),
					  Fb_triH1.GetInverseWeightPolynomial(),
					  pt2D, phi_ortho);
	
	phi.Reallocate(order*(order+1)/2);
	phi.Fill(0);
	Mlt(InverseBasisVDM, phi_ortho, phi);
      }
    else
      {
	phi.Reallocate(1);
	phi.Fill(1.0);
      }

    Fb_tri.ComputeValuesPhiRef(pt2D, psi);
    
    for (int j = 0; j <= order; j++)
      {
	Real_wp val = lob_gauss.EvaluatePhi(j, point_loc(2)); 
	for (int i = 0; i < psi.GetM(); i++)
	  {
	    if ((i >= 2*(order+1)) && (i < 3*(order+1)))
	      res(NumDofsXY(i, j)).Init(-val*psi(i)(0), -val*psi(i)(1), 0.0);	
	    else
	      res(NumDofsXY(i, j)).Init(val*psi(i)(0), val*psi(i)(1), 0.0);
	  }
      }
    
    for (int j = 0; j <= order+2; j++)
      {
	Real_wp val = lob_lobatto.EvaluatePhi(j, point_loc(2)); 
	for (int i = 0; i < phi.GetM(); i++)
	  res(NumDofsZ(i, j)).Init(0.0, 0.0, val*phi(i));
      }

  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void WedgeHdivOptimalFirstFamily::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    res.Fill(0);
    
    // values of triangular basis functions are computed
    R2 pt2D(point_loc(0), point_loc(1));
    VectReal_wp div_psi; VectReal_wp phi, phi_ortho;
    if (order > 1)
      {
	Fb_triH1.ComputeValuesPhiOrthoRef(order-1, Fb_triH1.GetNumOrtho2D(),
					  Fb_triH1.GetInverseWeightPolynomial(),
					  pt2D, phi_ortho);
	
	phi.Reallocate(order*(order+1)/2);
	phi.Fill(0);
	Mlt(InverseBasisVDM, phi_ortho, phi);
      }
    else
      {
	phi.Reallocate(1);
	phi.Fill(1.0);
      }

    Fb_tri.ComputeDivPhiRef(pt2D, div_psi);
    
    for (int j = 0; j <= order; j++)
      {
	Real_wp val = lob_gauss.EvaluatePhi(j, point_loc(2)); 
	for (int i = 0; i < div_psi.GetM(); i++)
	  {
	    if ((i >= 2*(order+1)) && (i < 3*(order+1)))
	      res(NumDofsXY(i, j)) = -val*div_psi(i);
	    else
	      res(NumDofsXY(i, j)) = val*div_psi(i);
	  }
      }
    
    for (int j = 0; j <= order+2; j++)
      {
	Real_wp val = lob_lobatto.EvaluatePhiGrad(j, point_loc(2)); 
	for (int i = 0; i < phi.GetM(); i++)
	  res(NumDofsZ(i, j)) = val*phi(i);
      }
  }


  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void WedgeHdivOptimalFirstFamily
  ::ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    // part due to C_{r-1}^3
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a = 2.0*x - 1.0, b = 2.0*y - 1.0, c = 2.0*z - 1.0;
        
    VectReal_wp Px, Py, Pz;

    EvaluateJacobiPolynomial(LegendrePolynom, r+2, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, r+2, b, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, r+2, c, Pz);
    
    Real_wp vloc; 
    int num = 0;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	for (int k = 0; k <= r; k++)
	  {          
	    vloc = Px(i)*Py(j)*Pz(k);
	    
	    phi(num).Init(vloc, 0, 0);
	    phi(num+1).Init(0, vloc, 0);
	    num += 2;
	    
	    if (i+j==r)
	      {
		phi(num).Init(vloc*x, vloc*y, 0);
		num += 1;
	      }
	  }
    
    for (int i = 0; i <= r-1; i++)
      for (int j = 0; j <= r-1-i; j++)
	for (int k = 0; k <= r+2; k++)
	  {
	    vloc = Px(i)*Py(j)*Pz(k);
	    
	    phi(num).Init(0, 0, vloc);
	    num += 1;
	  }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < phi.GetM(); i++)
      phi(i) *= InvWeightBasisRr(i);
    
  }

  
  //! Evaluating divergence nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] div_phi divergence of nearly orthogonal functions on pointloc
  */
  void WedgeHdivOptimalFirstFamily::
  ComputeDivPhiOrthoRef(const R3& pointloc, VectReal_wp& div_phi) const
  {
    div_phi.Reallocate(nb_dof_loc);
    div_phi.Fill(0); abort();
  }
  
  
  //! Integration against basis functions on a face
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] num_loc face number
    (\omega, \xi_k) are here Gauss points and not Gauss-Lobatto points
  */  
  template<class Vector1, class Vector2>
  void WedgeHdivOptimalFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }


  void WedgeHdivOptimalFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    if ((num_loc == 0) || (num_loc == 3))
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  void WedgeHdivOptimalFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    if ((num_loc == 0) || (num_loc == 3))
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  //! displays details of class WedgeHdivOptimalFirstFamily
  ostream& operator <<(ostream& out, const WedgeHdivOptimalFirstFamily& e)
  {
    out << static_cast< const WedgeReference<3>& >(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_WEDGE_HDIV_OPTIMAL_FIRST_FAMILY_CXX
#endif
