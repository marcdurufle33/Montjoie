#ifndef MONTJOIE_FILE_HEXAHEDRON_HDIV_HP_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronHdivHpFirstFamily::HexahedronHdivHpFirstFamily() : HexahedronReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  

  //! returns size of memory used by the object
  size_t HexahedronHdivHpFirstFamily::GetMemorySize() const
  {
    size_t taille = HexahedronReference<3>::GetMemorySize();
    taille += NumDofsX.GetMemorySize() + NumDofsY.GetMemorySize() + NumDofsZ.GetMemorySize();
    
    taille += LegendrePolynom.GetMemorySize() + JacobiPolynom.GetMemorySize();    
    taille += ValLeg.GetMemorySize() + ValBubbleJac.GetMemorySize() 
      + DerBubbleJac.GetMemorySize();

    taille += ch1_node.GetMemorySize() + ch2_node.GetMemorySize() + ch3_node.GetMemorySize()
      + ch1_loc.GetMemorySize() + ch2_loc.GetMemorySize() + ch3_loc.GetMemorySize()
      + rh_loc.GetMemorySize();
    
    for (int k = 0; k < 6; k++)
      {
        taille += sh_loc(k).GetMemorySize() + sh_locX(k).GetMemorySize() + sh_locY(k).GetMemorySize()
          + sh_curl_loc(k).GetMemorySize() + sh_curl_locX(k).GetMemorySize() + sh_curl_locY(k).GetMemorySize()
          + sh_nodeX(k).GetMemorySize() + sh_nodeY(k).GetMemorySize() +
          sh_curl_nodeX(k).GetMemorySize() + sh_curl_nodeY(k).GetMemorySize();
      }
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        {
          taille += sparse_mass_matrix(i, j).GetMemorySize() +
            sparse_stiff_matrix(i, j).GetMemorySize();
        }
    
    taille += sparse_matrix_chol.GetMemorySize() + ValuePhi2D.GetMemorySize();
    
    return taille;
  }

  
  //! how to number mesh
  void HexahedronHdivHpFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return HexahedronReference<3>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, order*order);
    
    // dofs inside the tetrahedron
    nmap.SetNbDofHexahedron(order, nb_dof_loc-nb_dof_boundaries);
    
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    ElementReference<Dimension2, 3>::FindHdivRotationQuad(order, NumQuad2D,
							  ValuePhi2D, nmap);
  }

  
  //! construction of finite element
  void HexahedronHdivHpFirstFamily::ConstructFiniteElement(int r, int rgeom,
							   int rquad, int type_quad,
							   int rsurf_tri, int rsurf_quad,
							   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    // type_quad = HexahedronQuadrature::QUADRATURE_JACOBI1;    
    HexahedronReference<3>::ConstructFiniteElement(r, rgeom, rquad,
						   Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    
    ConstructHdivElementaryMatrix();
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    FindDofsOnFace();

    QuadrangleDgOrtho* Fb_quad = new QuadrangleDgOrtho();
    element_quad_surf = Fb_quad;
    Fb_quad->ConstructFiniteElement(order-1, order, order, Globatto<Real_wp>::QUADRATURE_GAUSS);
  }
  
  
  //! construction of basis functions
  void HexahedronHdivHpFirstFamily::ConstructFunctions()
  {
    nb_dof_quad = order*order;
    nb_dof_tri = 0;
    nb_dof_loc = 3*order*order*(order+1);
    nb_dof_boundaries = 6*order*order;
    
    VectR3 points_dof3d;
    points_dof3d = this->PointsND();

    this->nb_points_dof_inside = this->nb_points_quadrature_inside;
    this->num_dof_points_surf = this->num_quad_points_surf;
    
    this->SetPointsDof2D_quad(this->Points2D_quad());
    this->SetPointsDofND(points_dof3d);
    
    GetJacobiPolynomial(LegendrePolynom, order, Real_wp(0), Real_wp(0));
    GetJacobiPolynomial(JacobiPolynom, order, Real_wp(1), Real_wp(1));
    
    NumDofsX.Reallocate(order+1, order, order);
    NumDofsY.Reallocate(order, order+1, order);
    NumDofsZ.Reallocate(order, order, order+1);
    NumDofsX.Fill(-1); NumDofsY.Fill(-1); NumDofsZ.Fill(-1);
    
    // dofs on faces
    int N = order*order;
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
        {
          int node = i*order + j;
          NumDofsX(0, i, j) = node;
          NumDofsY(i, 0, j) = N + node;
          NumDofsZ(i, j, 0) = 2*N + node;
	  NumDofsZ(i, j, order) = 3*N + node;
	  NumDofsY(i, order, j) = 4*N + node;
	  NumDofsX(order, i, j) = 5*N + node;
        }
    
    // dofs inside the hexahedron
    int offset = nb_dof_boundaries;
    for (int diag = 0; diag <= order-1; diag++)
      for (int i = 0; i <= diag; i++)
	for (int j = 0; j <= diag; j++)
	  for (int k = 0; k < diag; k++)
	    if ((i == diag) || (j == diag) || (k == diag-1))
	      {
		NumDofsX(k+1, i, j) = offset++;
		NumDofsY(i, k+1, j) = offset++;
		NumDofsZ(i, j, k+1) = offset++;
	      }
    
    // storing P_j^(0,0) ( \xi_i)  and (1-x)*x*P_j^{1,1}(\xi_i)
    int r = order_quad;
    ValLeg.Reallocate(order, r+3);
    ValBubbleJac.Reallocate(order+1, r+3);
    DerBubbleJac.Reallocate(order+1, r+3);
    ValBubbleJac.Fill(0);
    DerBubbleJac.Fill(0);
    VectReal_wp Px, Jx, dJx;
    for (int i = 0; i <= r; i++)
      {
        EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2.0*this->Points1D(i)-1.0, Px);
        EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2.0*this->Points1D(i)-1.0, Jx, dJx);
        for (int j = 0; j < order; j++)
          ValLeg(j, i) = Px(j);
        
        ValBubbleJac(0, i) = 1.0-this->Points1D(i); DerBubbleJac(0, i) = -1.0;
        ValBubbleJac(order, i) = this->Points1D(i); DerBubbleJac(order, i) = 1.0;
        for (int j = 0; j < order-1; j++)
          {
            ValBubbleJac(j+1, i) = Jx(j)*this->Points1D(i)*(1.0-this->Points1D(i));
            DerBubbleJac(j+1, i) = (1.0 - 2.0*this->Points1D(i))*Jx(j)
              + 2.0*dJx(j)*this->Points1D(i)*(1.0-this->Points1D(i));
          }        
      }
    
    // values on extremities are stored at the end
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, Real_wp(-1), Px);
    for (int j = 0; j < order; j++)
      ValLeg(j, r+1) = Px(j);
    
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, Real_wp(1), Px);
    for (int j = 0; j < order; j++)
      ValLeg(j, r+2) = Px(j);
    
    ValBubbleJac(0, r+1) = 1.0; ValBubbleJac(order, r+1) = 0.0;
    DerBubbleJac(0, r+1) = -1.0; DerBubbleJac(order, r+1) = 1.0;
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, Real_wp(-1), Jx);
    for (int j = 0; j < order-1; j++)
      DerBubbleJac(j+1, r+1) = Jx(j);
    
    ValBubbleJac(0, r+2) = 0.0; ValBubbleJac(order, r+2) = 1.0;
    DerBubbleJac(0, r+2) = -1.0; DerBubbleJac(order, r+2) = 1.0;
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, Real_wp(1), Jx);
    for (int j = 0; j < order-1; j++)
      DerBubbleJac(j+1, r+2) = -Jx(j);
    
    // adding quadrature points of the face in CoordinateQuad3D
    CoordinateQuad3D.Resize(this->PointsND().GetM(), 3);
    N = nb_points_quadrature_inside;
    const Matrix<int>& NumQuad2D = this->GetNumQuad2D();
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          int p = N + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = r+1;
          CoordinateQuad3D(p, 1) = i;
          CoordinateQuad3D(p, 2) = j;

          p = N + (r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = r+1;
          CoordinateQuad3D(p, 2) = j;

          p = N + 2*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = j;
          CoordinateQuad3D(p, 2) = r+1;


          p = N + 3*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = j;
          CoordinateQuad3D(p, 2) = r+2;


          p = N + 4*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = i;
          CoordinateQuad3D(p, 1) = r+2;
          CoordinateQuad3D(p, 2) = j;

          p = N + 5*(r+1)*(r+1) + NumQuad2D(i, j);
          CoordinateQuad3D(p, 0) = r+2;
          CoordinateQuad3D(p, 1) = i;
          CoordinateQuad3D(p, 2) = j;
        }

    ValuePhi2D.Reallocate(order*order, this->Points2D_quad().GetM());
    for (int i = 0; i < order; i++)
      for (int j = 0; j < order; j++)
	for (int k1 = 0; k1 <= order; k1++)
	  for (int k2 = 0; k2 <= order; k2++)
	    ValuePhi2D(i*order + j, NumQuad2D(k1, k2)) = ValLeg(i, k1)*ValLeg(j, k2);    
  }
  
  
  //! constructing mass matrix
  void HexahedronHdivHpFirstFamily::ConstructMassMatrix()
  {
  }
  
  
  //! constructing stiffness matrix
  void HexahedronHdivHpFirstFamily::ConstructStiffnessMatrix()
  {
    lob_quad.ComputeGradPhi(1e3*epsilon_machine);
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point
  */
  void HexahedronHdivHpFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;

    VectReal_wp Px, Py, Pz, Jx, Jy, Jz, Jx2, Jy2, Jz2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y-1, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*z-1, Pz);

    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*x-1, Jx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*y-1, Jy);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*z-1, Jz);
    
    Jx2.Reallocate(order+1);
    Jy2.Reallocate(order+1);
    Jz2.Reallocate(order+1);
    
    Jx2(0) = 1.0-x; Jx2(order) = x;
    Jy2(0) = 1.0-y; Jy2(order) = y;
    Jz2(0) = 1.0-z; Jz2(order) = z;
    for (int i = 0; i < order-1; i++)
      {
        Jx2(i+1) = x*(1.0-x)*Jx(i);
        Jy2(i+1) = y*(1.0-y)*Jy(i);
        Jz2(i+1) = z*(1.0-z)*Jz(i);
      }
    
    for (int k = 0; k <= order; k++)
      for (int i = 0; i <= order-1; i++)
        for (int j = 0; j <= order-1; j++)
          {  
            node = NumDofsX(k, i, j);
	    phi(node)(0) = Jx2(k)*Py(i)*Pz(j);
            phi(node)(1) = 0;
            phi(node)(2) = 0;
            
	    node = NumDofsY(i, k, j);
            phi(node)(0) = 0;
            phi(node)(1) = -Jy2(k)*Px(i)*Pz(j);
            phi(node)(2) = 0;
	    
            node = NumDofsZ(i, j, k);
            phi(node)(0) = 0;
            phi(node)(1) = 0;
            phi(node)(2) = Px(i)*Py(j)*Jz2(k);
	  }	
    
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point
  */
  void HexahedronHdivHpFirstFamily::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    
    Real_wp x, y, z;
    x = point_loc(0);
    y = point_loc(1);
    z = point_loc(2);
    
    int node = 0;
    
    VectReal_wp Px, Py, Pz, Jx, Jy, Jz, dJx, dJy, dJz;
    VectReal_wp Jx2, Jy2, Jz2, dJx2, dJy2, dJz2;
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*x-1, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*y-1, Py);
    EvaluateJacobiPolynomial(LegendrePolynom, order-1, 2*z-1, Pz);
    
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*x-1, Jx, dJx);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*y-1, Jy, dJy);
    EvaluateJacobiPolynomial(JacobiPolynom, order-1, 2*z-1, Jz, dJz);
    
    Jx2.Reallocate(order+1); dJx2.Reallocate(order+1);
    Jy2.Reallocate(order+1); dJy2.Reallocate(order+1);
    Jz2.Reallocate(order+1); dJz2.Reallocate(order+1);
    
    Jx2(0) = 1.0-x; Jx2(order) = x; dJx2(0) = -1.0; dJx2(order) = 1.0;
    Jy2(0) = 1.0-y; Jy2(order) = y; dJy2(0) = -1.0; dJy2(order) = 1.0;
    Jz2(0) = 1.0-z; Jz2(order) = z; dJz2(0) = -1.0; dJz2(order) = 1.0;
    for (int i = 0; i < order-1; i++)
      {
        Jx2(i+1) = x*(1.0-x)*Jx(i);
        Jy2(i+1) = y*(1.0-y)*Jy(i);
        Jz2(i+1) = z*(1.0-z)*Jz(i);

        dJx2(i+1) = (1.0 - 2.0*x)*Jx(i) + 2.0*x*(1.0-x)*dJx(i);
        dJy2(i+1) = (1.0 - 2.0*y)*Jy(i) + 2.0*y*(1.0-y)*dJy(i);
        dJz2(i+1) = (1.0 - 2.0*z)*Jz(i) + 2.0*z*(1.0-z)*dJz(i);
      }
    
    for (int k = 0; k <= order; k++)
      for (int i = 0; i < order; i++)
        for (int j = 0; j < order; j++)
          { 
	    node = NumDofsX(k, i, j); 
            res(node) = dJx2(k)*Py(i)*Pz(j);
	    
	    node = NumDofsY(i, k, j);
            res(node) = -dJy2(k)*Px(i)*Pz(j);
	    
	    node = NumDofsZ(i, j, k);
            res(node) = dJz2(k)*Px(i)*Py(j);
	  }      
  }
  

  void HexahedronHdivHpFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    if (num_loc%2 == 0)
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  void HexahedronHdivHpFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    if (num_loc %= 2)
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  
  
  //! displays details of class HexahedronHcurlHpFirstFamily
  ostream& operator <<(ostream& out, const HexahedronHdivHpFirstFamily& e)
  {
    out<<static_cast<const HexahedronReference<3>&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_HEXAHEDRON_HDIV_HP_FIRST_FAMILY_CXX
#endif
