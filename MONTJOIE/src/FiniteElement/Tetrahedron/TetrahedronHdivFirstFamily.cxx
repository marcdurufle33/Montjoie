#ifndef MONTJOIE_FILE_TETRAHEDRON_HDIV_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronHdivFirstFamily::TetrahedronHdivFirstFamily() : TetrahedronReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
    
    type_basis = NODAL_GAUSS;
    //type_basis = NODAL_REGULAR;
  }
  

  size_t TetrahedronHdivFirstFamily::GetMemorySize() const
  {
    size_t taille = TetrahedronReference<3>::GetMemorySize();
    taille += InverseVDM_Edge.GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    taille += CoordinateDofs.GetMemorySize();
    return taille;
  }

  
  //! how to number mesh
  void TetrahedronHdivFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<3>::ConstructNumberMap(nmap, dg);
    
    // if r is order of the approximation
    // we have r dofs on each edge 
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    // nothing on quadrangular face !
    nmap.SetNbDofQuadrangle(order, 0);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
    
    // dofs inside the tetrahedron
    nmap.SetNbDofTetrahedron(order, nb_dof_loc-nb_dof_boundaries);
    nmap.SetNbDofHexahedron(order, 0);
    
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
    
    if (order > 1)
      {
	/* TriangleClassical Fb_tri;
	Fb_tri.ConstructFiniteElement(order-1, order-1, order);
	
	Fb_tri.SetInterpolationPoints(this->points_dof2d_tri);
	
	Matrix<Real_wp> ValuePhi2D(Fb_tri.GetNbDof(), this->points2d_tri.GetM());
	ValuePhi2D.Fill(0); VectReal_wp phi;
	for (int k = 0; k < this->points2d_tri.GetM(); k++)
	  {
	    Fb_tri.ComputeValuesPhiRef(this->points2d_tri(k), phi);
	    for (int i = 0; i < Fb_tri.GetNbDof(); i++)
	      ValuePhi2D(i, k) = phi(i);
	  }
	
	FaceReference<3>::FindHdivRotationTri(order, points2d_tri, weights2d_tri,
					      ValuePhi2D, nmap);
	*/
      }
  }


  //! construction of nearly orthogonal functions generating the finite element space
  void TetrahedronHdivFirstFamily
  ::ConstructOrthogonalBasis(int r, int rgeom, int rquad, int type_quad)
  {
    TetrahedronReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    int nb_dof_inside_tri = order*(order+1)/2;
    nb_dof_boundaries = 4*nb_dof_inside_tri;
    nb_dof_loc = order*(order+1)*(order+3)/2;
    nb_dof_tri = nb_dof_inside_tri;
    nb_dof_quad = 0;
    
    const Matrix<Real_wp>& CoefOddJacobi = Fb_geom.GetCoefOddJacobi();
    const Matrix<Real_wp>& CoefEvenJacobi = Fb_geom.GetCoefEvenJacobi();
    const VectReal_wp& CoefLegendre = Fb_geom.GetCoefLegendre();
    
    // weight for orthonormalization
    InvWeightBasisRr.Reallocate(nb_dof_loc);
    int num = 0;
    for (int i = 0; i < order; i++)
      {
        for (int j = 0; j <= order-1-i; j++)
          {
            for (int k = 0; k < order-i-j; k++)
              {
                Real_wp coef = CoefLegendre(i)*CoefOddJacobi(i, j)*CoefEvenJacobi(i+j, k);
                InvWeightBasisRr(num) = coef;
                InvWeightBasisRr(num+1) = coef;
                InvWeightBasisRr(num+2) = coef;
                num += 3;
              }
            
            InvWeightBasisRr(num) = CoefLegendre(i)*CoefOddJacobi(i, j)
              *CoefEvenJacobi(i+j, order-1-i-j);
            num ++;
          }
      }
    
    InvWeightBasisRr.Fill(1.0);
    VectReal_wp Scale(nb_dof_loc);
    VectR3 phi(nb_dof_loc); Scale.Fill(0);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
        for (int i = 0; i < nb_dof_loc; i++)
          Scale(i) += DotProd(phi(i), phi(i))*this->WeightsND(k);
      }
    
    for (int i = 0; i < nb_dof_loc; i++)
      InvWeightBasisRr(i) = 1.0/sqrt(Scale(i));

    // DISP(InvWeightBasisRr);

  }
  
  
  //! construction of finite element
  void TetrahedronHdivFirstFamily
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
    
    this->num_dof_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      {
	this->num_dof_points_surf(n).Reallocate(this->GetNbDofBoundary(n));
	for (int i = 0; i < this->GetNbDofBoundary(n); i++)
	  this->num_dof_points_surf(n)(i) = FacesDof(i, n);
      }    
  }
  
  
  //! construction of basis functions
  void TetrahedronHdivFirstFamily::ConstructFunctions()
  {
    // points where dofs are defined, and tangent on these points
    VectR3 phi, tangente_dof(nb_dof_loc), points_inside;
    VectR2 points_tri, points_quad;
    VectReal_wp points_lob, weights_lob;
    Matrix<int> NumTri, NumQuad; Array3D<int> NumNodes;        
    ComputeGaussLobatto(points_lob, weights_lob, order+1);

    VectR2 points_dof2d_tri;
    VectR3 points_dof3d; VectReal_wp points_dof1d;

    if (type_basis == NODAL_GAUSS)
      ComputeGaussLegendre(points_dof1d, weights_lob, order-1);
    else
      {
	points_dof1d.Reallocate(order);
	for (int i = 0; i < order; i++)
	  points_dof1d(i) = points_lob(i+1);
      }
    
    Fb_geom.ConstructLobattoPoints(order+2, points_lob, points_tri, points_inside);
    
    points_dof3d.Reallocate(nb_dof_loc);
    
    // dofs on faces
    int offset = 0; int Ntri = order*(order+1)/2;
    points_dof2d_tri.Reallocate(Ntri);
    for (int i = 0; i < order*(order+1)/2; i++)
      {
	Real_wp x = points_tri(3*(order+2) + i)(0);
	Real_wp y = points_tri(3*(order+2) + i)(1);
	points_dof2d_tri(i).Init(x, y);
        
	points_dof3d(offset).Init(x, y, 0);
	tangente_dof(offset).Init(0, 0, 1.0);
	
	points_dof3d(Ntri+offset).Init(x, 0, y);
	tangente_dof(Ntri+offset).Init(0, -1.0, 0);
	
	points_dof3d(2*Ntri+offset).Init(0, x, y);
	tangente_dof(2*Ntri+offset).Init(1.0, 0, 0);
        
	points_dof3d(3*Ntri+offset).Init(1.0-x-y, x, y);
	tangente_dof(3*Ntri+offset).Init(1.0, 1.0, 1.0);
	//Mlt(1.0/sqrt(3.0), tangente_dof(3*Ntri+offset));
	
	offset++;
      }
    
    // dofs inside
    offset += 3*Ntri;
    for (int i = 0; i < (order-1)*(order+1)*order/6; i++)
      {
	R3 pt = points_inside(2*(order+2)*(order+2)+2+i);
	points_dof3d(offset) = pt;
	points_dof3d(offset+1) = pt;
	points_dof3d(offset+2) = pt;
        
	tangente_dof(offset).Init(1, 0, 0);
	tangente_dof(offset+1).Init(0, 1, 0);
	tangente_dof(offset+2).Init(0, 0, 1);
	offset += 3;
      }

    this->SetPointsDof1D(points_dof1d);        
    this->SetPointsDof2D_tri(points_dof2d_tri);
    this->SetPointsDofND(points_dof3d);        

    Matrix<Real_wp> Mh(nb_dof_loc, nb_dof_loc);
    Mh.Fill(0);
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
	ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
	for (int i = 0; i < nb_dof_loc; i++)
	  for (int j = 0; j < nb_dof_loc; j++)
	    Mh(i, j) += this->WeightsND(k)*DotProd(phi(i), phi(j));
      }
    
    //Mh.WriteText("MassTetOpt"+to_str(order)+".dat");

    // computation of VDM = psi_j(xi_k) \cdot t_k
    InverseVDM_Edge.Reallocate(nb_dof_loc, nb_dof_loc);
    for (int k = 0; k < nb_dof_loc; k++)
      {
	ComputeValuesPhiOrthoRef(points_dof3d(k), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  InverseVDM_Edge(j, k) = DotProd(phi(j), tangente_dof(k));
      }
    
    //InverseVDM_Edge.WriteText("VdmTet"+to_str(order)+".dat");
    GetInverse(InverseVDM_Edge);
    
    // then checking phi_i(xi_k) \cdot t_k = delta_{i, k}
    for (int k = 0; k < nb_dof_loc; k++)
      {
	ComputeValuesPhiRef(points_dof3d(k), phi);            
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    Real_wp val = DotProd(phi(j), tangente_dof(k));
	    Real_wp val_exact = 0;
	    if (j == k)
	      val_exact = 1.0;
	    
	    if (abs(val - val_exact) > 1e4*epsilon_machine)
	      {
		DISP(val); DISP(val_exact);
		abort();
	      }
	  }
      }
  }
   
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TetrahedronHdivFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectR3 psi;
    ComputeValuesPhiOrthoRef(point_loc, psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	Add(InverseVDM_Edge(i, j), psi(j), res(i));
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void TetrahedronHdivFirstFamily
  ::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectReal_wp div_psi;
    ComputeDivPhiOrthoRef(point_loc, div_psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	res(i) += InverseVDM_Edge(i, j)*div_psi(j);
  }


  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void TetrahedronHdivFirstFamily::
  ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    // coordinates on the cube [-1, 1]
    Real_wp a(-1), b(-1), c(-1);
    if (abs(1.0-y-z) > epsilon_machine)
      a = 2.0*x/(1.0-y-z) - 1.0;
    
    if (abs(1.0-z) > epsilon_machine)
      b = 2.0*y/(1.0-z) - 1.0;
    
    c = 2.0*z - 1.0;

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    VectReal_wp Px;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px);
    Vector<VectReal_wp> Py(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, b, Py(i));
    
    Vector<VectReal_wp> Pz(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i));
    
    Real_wp valx(0), valy(0), vloc;
    Real_wp pow_oneMinusYZ = 1.0, pow_oneMinusZ = 1.0;
    int num = 0;
    for (int i = 0; i < r; i++)
      {
        // computation of valx = L_i(a) (1-y-z)^i
        if (i == 0)
          valx = 1.0;
        else if (i == 1)
          valx = 2.0*x + y + z - 1.0;
        else
          valx = pow_oneMinusYZ*Px(i);
        
        pow_oneMinusZ = 1.0;
        for (int j = 0; j < r-i; j++)
          {
            // computation of valy = P_j^(2i+1)(a) (1-z)^j
            if (j == 0)
              valy = 1.0;
            else if (j == 1)
              valy = 2.0*y + z - 1.0 - OddJacobiPolynom(i)(0,0)*(1.0-z);
            else
              valy = pow_oneMinusZ*Py(i)(j);
            
            for (int k = 0; k < r-i-j; k++)
              {
                vloc = valx*valy*Pz(i+j)(k);
                // part due to P_r-1^3
                phi(num).Init(vloc, 0, 0);
                phi(num+1).Init(0, vloc, 0);
                phi(num+2).Init(0, 0, vloc);
                num += 3;
              }
            
            // when vloc = x^i y^j z^r-1-i-j
            // part due to (vloc x, vloc y, vloc z)
            vloc = valx*valy*Pz(i+j)(r-1-i-j);
            phi(num).Init(vloc*x, vloc*y, vloc*z);
            num ++;
            
            pow_oneMinusZ *= 1.0-z;
          }
	
        pow_oneMinusYZ *= 1.0-y-z;
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
  void TetrahedronHdivFirstFamily::
  ComputeDivPhiOrthoRef(const R3& pointloc, VectReal_wp& div_phi) const
  {
    div_phi.Reallocate(nb_dof_loc);
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    // coordinates on the cube [-1, 1]
    Real_wp a(-1), b(-1), c(-1);
    Real_wp da_dx(0), da_dy(0), db_dy(0), db_dz(0);
    if (abs(1.0-y-z) > epsilon_machine)
      {
	a = 2.0*x/(1.0-y-z) - 1.0;
	da_dx = 2.0/(1.0-y-z);
	da_dy = 2.0*x/square(1.0-y-z);
      }
    
    if (abs(1.0-z) > epsilon_machine)
      {
	b = 2.0*y/(1.0-z) - 1.0;
	db_dy = 2.0/(1.0-z);
	db_dz = 2.0*y/square(1.0-z);
      }
    
    c = 2.0*z - 1.0;

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();
    const Vector<Matrix<Real_wp> >& OddJacobiPolynom = Fb_geom.GetOddJacobiPolynomial();
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    VectReal_wp Px, dPx;
    int r = order;
    EvaluateJacobiPolynomial(LegendrePolynom, r, a, Px, dPx);
    Vector<VectReal_wp> Py(r+1), dPy(r+1);
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(OddJacobiPolynom(i), r-i, b, Py(i), dPy(i));
    
    Vector<VectReal_wp> Pz(r+1), dPz(r+1);
    for (int i = 0; i <= order; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i), dPz(i));
    
    Real_wp valx(0), valy(0), valz(0), phi(0); R3 grad_phi;
    Real_wp dvalx_dx(0), dvalx_dy(0), dvaly_dy(0), dvaly_dz(0), dvalz_dz(0);
    Real_wp pow_oneMinusYZ = 1.0, pow_oneMinusZ = 1.0;    
    Real_wp powm1_oneMinusYZ = 0.0, powm1_oneMinusZ = 0.0;    
    int num = 0;
    for (int i = 0; i < r; i++)
      {
        if (i == 0)
          {
            valx = 1.0;
            dvalx_dx = 0;
            dvalx_dy = 0;
          }
        else if (i == 1)
          {
            valx = 2.0*x + y + z - 1.0;
            dvalx_dx = 2.0;
            dvalx_dy = 1.0;
          }
        else
          {
            valx = pow_oneMinusYZ*Px(i);
            dvalx_dx = pow_oneMinusYZ*dPx(i)*da_dx;
            dvalx_dy = pow_oneMinusYZ*dPx(i)*da_dy - i*powm1_oneMinusYZ*Px(i);
          }
        
        powm1_oneMinusZ = 0; pow_oneMinusZ = 1.0;
        for (int j = 0; j < r-i; j++)
          {
            if (j == 0)
              {
                valy = 1.0;
                dvaly_dy = 0;
                dvaly_dz = 0;
              }
            else if (j == 1)
              {
                valy = 2.0*y + z - 1.0 - OddJacobiPolynom(i)(0,0)*(1.0-z);
                dvaly_dy = 2.0;
                dvaly_dz = 1.0 + OddJacobiPolynom(i)(0,0);
              }
            else
              {
                valy = pow_oneMinusZ*Py(i)(j);
                dvaly_dy = pow_oneMinusZ*dPy(i)(j)*db_dy;
                dvaly_dz = pow_oneMinusZ*dPy(i)(j)*db_dz - j*powm1_oneMinusZ*Py(i)(j);
              }
                   
            for (int k = 0; k < r-i-j; k++)
              {
                valz = Pz(i+j)(k);
                dvalz_dz = 2.0*dPz(i+j)(k);
                phi = valx*valy*valz;
                grad_phi.Init(dvalx_dx*valy*valz, (dvalx_dy*valy + valx*dvaly_dy)*valz, 
                              dvalx_dy*valy*valz + valx*dvaly_dz*valz + valx*valy*dvalz_dz);
                
                // part due to P_r-1^3
                div_phi(num) = grad_phi(0);
                div_phi(num+1) = grad_phi(1);
                div_phi(num+2) = grad_phi(2);
                num += 3;
              }
            
            // when vloc = x^i y^j z^r-1-i-j
            // part due to (vloc x, vloc y, vloc z)
            div_phi(num) = 3.0*phi + grad_phi(0)*x + grad_phi(1)*y + grad_phi(2)*z;
            num ++;
            
            powm1_oneMinusZ = pow_oneMinusZ;
            pow_oneMinusZ *= 1.0-z;
          }
        
        powm1_oneMinusYZ = pow_oneMinusYZ;
	pow_oneMinusYZ *= 1.0-y-z;
      }
    
    // multiplying by inverse of weights in order to have orthonormality
    for (int i = 0; i < nb_dof_loc; i++)
      div_phi(i) *= InvWeightBasisRr(i);

  }


  void TetrahedronHdivFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    if ((num_loc == 0) || (num_loc == 2))
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  

  void TetrahedronHdivFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {
    if ((num_loc == 0) || (num_loc == 2))
      for (int i = 0; i < contrib.GetM(); i++)
	contrib(i) = -contrib(i);
  }
  
  
  //! displays details of class TetrahedronHdivFirstFamily
  ostream& operator <<(ostream& out, const TetrahedronHdivFirstFamily& e)
  {
    out << static_cast<const TetrahedronReference<3>& >(e);
    return  out;
  }

  
  //! default constructor
  TetrahedronHdivOptimalFirstFamily::TetrahedronHdivOptimalFirstFamily()
    : TetrahedronHdivFirstFamily()
  {
    type_basis = NODAL_LOBATTO;
  }


  //! displays details of class TetrahedronHcurlOptimalFirstFamily
  ostream& operator <<(ostream& out, const TetrahedronHdivOptimalFirstFamily& e)
  {
    out<<static_cast<const TetrahedronHdivFirstFamily&>(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_TETRAHEDRON_HDIV_FIRST_FAMILY_CXX
#endif
