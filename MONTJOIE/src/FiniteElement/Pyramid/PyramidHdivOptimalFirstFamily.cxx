#ifndef MONTJOIE_FILE_PYRAMID_HDIV_OPTIMAL_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidHdivOptimalFirstFamily::PyramidHdivOptimalFirstFamily() : PyramidReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void PyramidHdivOptimalFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<3>::ConstructNumberMap(nmap, dg);

    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, (order+1)*(order+1));
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
    
    // dofs inside the pyramid
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);
    
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
  

  size_t PyramidHdivOptimalFirstFamily::GetMemorySize() const
  {
    size_t taille = PyramidReference<3>::GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    taille += InverseBasisVDM.GetMemorySize();
    return taille;
  }

      
  //! construction nearly-orthogonal functions psi_i used to generate the correct space
  void PyramidHdivOptimalFirstFamily
  ::ConstructOrthogonalBasis(int r, int rgeom, int rquad, int type_quad)
  {
    if (rquad <= r)
      rquad = r+1;
    
    PyramidReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad, r, r,
                                                TriangleQuadrature::QUADRATURE_GAUSS,
						Globatto<Real_wp>::QUADRATURE_GAUSS);

    VectReal_wp points1d, weights1d;    
    ComputeGaussLegendre(points1d, weights1d, order);

    this->SetPoints1D(points1d);
    this->SetWeights1D(weights1d);

    nb_dof_tri = order*(order+1)/2;
    nb_dof_quad = (order+1)*(order+1);
    nb_dof_boundaries = 4*nb_dof_tri + nb_dof_quad;
    nb_dof_loc = order*(order+1)*(2*order+1)/2 + 2*order*(order+1) + (order+1)*(order+1);
    
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
  void PyramidHdivOptimalFirstFamily
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
  void PyramidHdivOptimalFirstFamily::ConstructFunctions()
  {
    VectReal_wp points_edge;
    VectR2 points_tri;
    TriangleGeomReference::ConstructLobattoPoints(order+2, 0,
						  points_edge, points_tri);
    
    VectR2 points_dof2d_tri, points_dof2d_quad;
    VectR3 points_dof3d;
    int Ntri = order*(order+1)/2;
    points_dof2d_tri.Reallocate(Ntri);
    for (int i = 0; i < Ntri; i++)
      points_dof2d_tri(i) = points_tri(3*(order+2) + i);
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
    
    points_dof2d_quad.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	points_dof2d_quad(i*(order+1) + j).Init(gauss.Points(i), gauss.Points(j));
    
    points_dof3d.Reallocate(nb_dof_loc);
    VectR3 tangente_dof(nb_dof_loc);
    int Nquad = (order+1)*(order+1);
    // dofs on quadrilateral base
    for (int i = 0; i < Nquad; i++)
      {
	Real_wp x = points_dof2d_quad(i)(0);
	Real_wp y = points_dof2d_quad(i)(1);
	points_dof3d(i).Init(2.0*x-1.0, 2.0*y-1.0, 0.0);
	tangente_dof(i).Init(0, 0, 1.0);
      }
    
    int nb = Nquad;
    // dofs on triangular faces
    for (int i = 0; i < Ntri; i++)
      {
	Real_wp x = points_dof2d_tri(i)(0);
	Real_wp y = points_dof2d_tri(i)(1);
	points_dof3d(nb).Init(2.0*x + y - 1.0, y - 1.0, y);
	tangente_dof(nb).Init(0, -1.0, 1.0);
	points_dof3d(Ntri+nb).Init(1.0-y, 2.0*x + y - 1.0, y);
	tangente_dof(Ntri+nb).Init(1.0, 0.0, 1.0);
	points_dof3d(2*Ntri+nb).Init(2.0*x + y - 1.0, 1.0 - y, y);
	tangente_dof(2*Ntri+nb).Init(0, -1.0, -1.0);
	points_dof3d(3*Ntri+nb).Init(y - 1.0, 2.0*x + y - 1.0, y);
	tangente_dof(3*Ntri+nb).Init(1.0, 0.0, -1.0);
	nb++;
      }
    
    nb += 3*Ntri;
    
    VectReal_wp points_lob; VectR2 points_quad;
    Matrix<int> NumNodesQuad(order+3, order+3);
    VectR3 points_inside; NumNodesQuad.Fill();
    Fb_geom.ConstructLobattoPoints(order+2, points_lob, points_tri,
				   points_quad, NumNodesQuad, points_inside);
    
    // dofs inside the pyramid
    for (int i = 3*(order+2)*(order+2)+2; i < points_inside.GetM(); i++)
      {
	points_dof3d(nb) = points_inside(i);
	tangente_dof(nb).Init(1, 0, 0);

	points_dof3d(nb+1) = points_inside(i);
	tangente_dof(nb+1).Init(0, 1, 0);

	points_dof3d(nb+2) = points_inside(i);
	tangente_dof(nb+2).Init(0, 0, 1);
	nb += 3;
      }

    this->SetPointsDof2D_tri(points_dof2d_tri);
    this->SetPointsDof2D_quad(points_dof2d_quad);
    this->SetPointsDofND(points_dof3d);        

    Matrix<Real_wp> Mh(nb_dof_loc, nb_dof_loc);
    Mh.Fill(0);     VectR3 phi;
    for (int k = 0; k < nb_points_quadrature_inside; k++)
      {
        ComputeValuesPhiOrthoRef(this->PointsND(k), phi);
        for (int i = 0; i < nb_dof_loc; i++)
          for (int j = 0; j < nb_dof_loc; j++)
            Mh(i, j) += this->WeightsND(k)*DotProd(phi(i), phi(j));
      }
    
    //Mh.WriteText("Mass.dat");
    
    InverseBasisVDM.Reallocate(nb_dof_loc, nb_dof_loc);
    for (int k = 0; k < nb_dof_loc; k++)
      {
        ComputeValuesPhiOrthoRef(points_dof3d(k), phi);
        for (int j = 0; j < nb_dof_loc; j++)
          InverseBasisVDM(j, k) = DotProd(phi(j), tangente_dof(k));
      }
    
    //InverseBasisVDM.WriteText("VDM.dat");
    GetInverse(InverseBasisVDM);
    
  }
   
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void PyramidHdivOptimalFirstFamily
  ::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectR3 psi;
    ComputeValuesPhiOrthoRef(point_loc, psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	Add(InverseBasisVDM(i, j), psi(j), res(i));
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void PyramidHdivOptimalFirstFamily
  ::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    VectReal_wp div_psi;
    ComputeDivPhiOrthoRef(point_loc, div_psi);
    FillZero(res);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	res(i) += InverseBasisVDM(i, j)*div_psi(j);
  }


  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void PyramidHdivOptimalFirstFamily
  ::ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    int r = order;
    // part due to C_{r-1}^3
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0), c = 2.0*z - 1.0;
    if (abs(1.0 - z) > epsilon_machine)
      {
        a = x/(1.0-z);
        b = y/(1.0-z);
      }

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();    
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    VectReal_wp Px, Py, pow_oneMinusZ(r+2);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, a, Px);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, b, Py);
    Vector<VectReal_wp> Pz(r+1);
    pow_oneMinusZ(0) = 1.0;
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i));
    
    for (int i = 0; i <= r; i++)
      pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
    
    Real_wp val, vloc; int num = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
          int m = max(i, j);
          if (m == 0)
            val = 1.0;
          else if (m == 1)
            {
              if (i == 0)
                val = y;
              else if (j == 0)
                val = x;
              else
                val = a*y;
            }
          else
            val = Px(i)*Py(j)*pow_oneMinusZ(m);
          
          for (int k = 0; k < r-m; k++)
            {
              vloc = val*Pz(m)(k);
              phi(num).Init(vloc, 0, 0);
              phi(num+1).Init(0, vloc, 0);
              phi(num+2).Init(0, 0, vloc);
              num += 3;
            }
        }
    
    for (int k = 0; k < r; k++)
      for (int m = 0; m <= k; m++)
	{
	  val = Px(k+1)*Py(m)*pow_oneMinusZ(k);
	  phi(num++).Init(val, 0, 0);
	  
	  val = Px(m)*Py(k+1)*pow_oneMinusZ(k);
	  phi(num++).Init(0, val, 0);
	  
	  val = Px(m)*Py(k+1)*pow_oneMinusZ(k);
	  phi(num++).Init(a*val, 0, -val);
	  
	  val = Px(k+1)*Py(m)*pow_oneMinusZ(k);
	  phi(num++).Init(0, b*val, -val);	  
	}
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
	  val = Px(i)*Py(j)*pow_oneMinusZ(r);
	  phi(num++).Init(a*val, b*val, -val);
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
  void PyramidHdivOptimalFirstFamily::
  ComputeDivPhiOrthoRef(const R3& pointloc, VectReal_wp& div_phi) const
  {
    div_phi.Reallocate(nb_dof_loc);

    int r = order;
    // part due to C_{r-1}^3
    Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
    Real_wp a(0), b(0), c = 2.0*z - 1.0, da_dx(0), da_dz(0), db_dy(0), db_dz(0);
    if (abs(1.0 - z) > epsilon_machine)
      {
        a = x/(1.0-z);
        b = y/(1.0-z);
	da_dx = 1.0/(1.0-z);
	da_dz = x/square(1.0-z);
	db_dy = da_dx;
	db_dz = y/square(1.0-z);
      }

    const Matrix<Real_wp>& LegendrePolynom = Fb_geom.GetLegendrePolynomial();    
    const Vector<Matrix<Real_wp> >& EvenJacobiPolynom = Fb_geom.GetEvenJacobiPolynomial();
    
    VectReal_wp Px, dPx, Py, dPy, pow_oneMinusZ(r+2), dpow_oneMinusZ(r+1);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, a, Px, dPx);
    EvaluateJacobiPolynomial(LegendrePolynom, r+1, b, Py, dPy);
    Vector<VectReal_wp> Pz(r+1), dPz(r+1);
    pow_oneMinusZ(0) = 1.0;
    dpow_oneMinusZ(0) = 0.0;
    for (int i = 0; i <= r; i++)
      EvaluateJacobiPolynomial(EvenJacobiPolynom(i), r-i, c, Pz(i), dPz(i));
    
    for (int i = 0; i <= r; i++)
      {
	pow_oneMinusZ(i+1) = pow_oneMinusZ(i)*(1.0-z);
	if (i > 0)
	  dpow_oneMinusZ(i) = -Real_wp(i)*pow_oneMinusZ(i-1);
      }
    
    Real_wp val, dv_dz; int num = 0; R3 grad;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
          int m = max(i, j);
          if (m == 0)
            {
	      val = 1.0;
	      grad.Fill(0.0);
	    }
          else if (m == 1)
            {
              if (i == 0)
                {
		  val = y;
		  grad.Init(0, 1, 0);
		}
              else if (j == 0)
                {
		  val = x;
		  grad.Init(1, 0, 0);
		}
              else
                {
		  val = a*y;
		  grad.Init(da_dx*y, a, da_dz*y);
		}
            }
          else
            {
	      val = Px(i)*Py(j)*pow_oneMinusZ(m);
	      grad(0) = da_dx*dPx(i)*Py(j)*pow_oneMinusZ(m);
	      grad(1) = db_dy*Px(i)*dPy(j)*pow_oneMinusZ(m);
	      grad(2) = pow_oneMinusZ(m)*(da_dz*dPx(i)*Py(j) + db_dz*Px(i)*dPy(j))
		+ dpow_oneMinusZ(m)*Px(i)*Py(j);
	    }
          
          for (int k = 0; k < r-m; k++)
            {
	      dv_dz = Pz(m)(k)*grad(2) + 2.0*val*dPz(m)(k);
	      
	      div_phi(num) = Pz(m)(k)*grad(0);
	      div_phi(num+1) = Pz(m)(k)*grad(1);
	      div_phi(num+2) = dv_dz;
              num += 3;
            }
        }
    
    for (int k = 0; k < r; k++)
      for (int m = 0; m <= k; m++)
	{
	  // val = Px(k+1)*Py(m)*pow_oneMinusZ(k);
	  div_phi(num++) = da_dx*dPx(k+1)*Py(m)*pow_oneMinusZ(k);
	  
	  // val = Px(m)*Py(k+1)*pow_oneMinusZ(k);
	  div_phi(num++) = db_dy*Px(m)*dPy(k+1)*pow_oneMinusZ(k);
	  
	  val = Px(m)*Py(k+1)*pow_oneMinusZ(k);
	  div_phi(num++) = da_dx*val + a*da_dx*dPx(m)*Py(k+1)*pow_oneMinusZ(k)
	    - pow_oneMinusZ(k)*(da_dz*dPx(m)*Py(k+1) + db_dz*Px(m)*dPy(k+1))
	    - dpow_oneMinusZ(k)*Px(m)*Py(k+1);
	  
	  val = Px(k+1)*Py(m)*pow_oneMinusZ(k);
	  div_phi(num++) = db_dy*val + b*db_dy*Px(k+1)*dPy(m)*pow_oneMinusZ(k)
	    - pow_oneMinusZ(k)*(da_dz*dPx(k+1)*Py(m) + db_dz*Px(k+1)*dPy(m))
	    - dpow_oneMinusZ(k)*Px(k+1)*Py(m);
	}
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
	  val = Px(i)*Py(j)*pow_oneMinusZ(r);
	  grad(0) = da_dx*dPx(i)*Py(j)*pow_oneMinusZ(r);
	  grad(1) = db_dy*Px(i)*dPy(j)*pow_oneMinusZ(r);
	  grad(2) = pow_oneMinusZ(r)*(da_dz*dPx(i)*Py(j) + db_dz*Px(i)*dPy(j))
	    + dpow_oneMinusZ(r)*Px(i)*Py(j);
	  
	  div_phi(num++) = val*(da_dx + db_dy) + a*grad(0) + b*grad(1) - grad(2);
	}
    
    // multiplying by inverse of weights
    for (int i = 0; i < div_phi.GetM(); i++)
      div_phi(i) *= InvWeightBasisRr(i);

  }
  

  void PyramidHdivOptimalFirstFamily
  ::ModifySignProjectionSurface(VectReal_wp& contrib, int num_loc) const
  {    
    Real_wp one(1);
    switch (num_loc)
      {
      case 0:
	Mlt(-one/4, contrib);
	break;
      case 1:
      case 2:
	Mlt(one/2, contrib);
	break;
      case 3:
      case 4:
	Mlt(-one/2, contrib);
	break;
      }
  }


  void PyramidHdivOptimalFirstFamily
  ::ModifySignProjectionSurface(VectComplex_wp& contrib, int num_loc) const
  {    
    Real_wp one(1);
    switch (num_loc)
      {
      case 0:
	Mlt(-one/4, contrib);
	break;
      case 1:
      case 2:
	Mlt(one/2, contrib);
	break;
      case 3:
      case 4:
	Mlt(-one/2, contrib);
	break;
      }
  }

  
  //! displays details of class PyramidHcurlFirstFamily
  ostream& operator <<(ostream& out, const PyramidHdivOptimalFirstFamily& e)
  {
    out << static_cast< const PyramidReference<3>& >(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HDIV_OPTIMAL_FIRST_FAMILY_CXX
#endif
