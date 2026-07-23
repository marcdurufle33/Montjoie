#ifndef MONTJOIE_FILE_PYRAMID_HDIV_FIRST_FAMILY_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidHdivFirstFamily::PyramidHdivFirstFamily() : PyramidReference<3>()
  {
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number mesh
  void PyramidHdivFirstFamily::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg == ElementReference_Base::DISCONTINUOUS)
      return PyramidReference<3>::ConstructNumberMap(nmap, dg);
    
    nmap.SetNbDofVertex(order, 0);
    nmap.SetNbDofEdge(order, 0);
    nmap.SetNbDofQuadrangle(order, order*order);
    // triangular face
    nmap.SetNbDofTriangle(order, order*(order+1)/2);
    
    // dofs inside the pyramid
    nmap.SetNbDofPyramid(order, nb_dof_loc-nb_dof_boundaries);
    
  }
  

  size_t PyramidHdivFirstFamily::GetMemorySize() const
  {
    size_t taille = PyramidReference<3>::GetMemorySize();
    taille += InvWeightBasisRr.GetMemorySize();
    return taille;
  }

    
  //! construction nearly-orthogonal functions psi_i used to generate the correct space
  void PyramidHdivFirstFamily::
  ConstructOrthogonalBasis(int r, int rgeom, int rquad, int type_quad)
  {
    if (rquad == 0)
      rquad = r;
    
    PyramidReference<3>::ConstructFiniteElement(r, rgeom, rquad, type_quad, r, r,
                                                TriangleQuadrature::QUADRATURE_GAUSS,
						Globatto<Real_wp>::QUADRATURE_LOBATTO);

    VectReal_wp points1d, weights1d;    
    ComputeGaussLegendre(points1d, weights1d, order);

    this->SetPoints1D(points1d);
    this->SetWeights1D(weights1d);

    nb_dof_tri = order*(order+1)/2;
    nb_dof_quad = order*order;
    nb_dof_boundaries = 4*nb_dof_tri + nb_dof_quad;
    // nb_dof_loc = order*(2*order*order+9*order+5)/2; 
    
    nb_dof_loc = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	for (int k = 0; k < r-max(i, j); k++)
	  {
	    nb_dof_loc += 3;
	    //if (i+j+k==r-1)
	    //nb_dof_loc++;
	  }
    
    for (int k = 0; k < r; k++)
      for (int m = 0; m <= k; m++)
	//if (k+m<=r-1)
	nb_dof_loc+=2;

    for (int k = 0; k < r; k++)
      for (int m = 0; m <= k; m++)
	if (k < r-1)
	  nb_dof_loc+=2;

    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	nb_dof_loc++;
    
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
  
  
  //! constructing finite element
  void PyramidHdivFirstFamily::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						      int rsurf_tri, int rsurf_quad,
						      int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    ConstructOrthogonalBasis(r, rgeom, rquad, type_quad);
    
    // computation of basis functions
    ConstructFunctions();
    
    // coefficients for fast computation of Fi on curved tets.
    this->Fb_geom.ComputeCoefficientTransformation();
    
    // construction of Value_PhiVec, Curl_Phi and elementary matrices
    //ConstructHdivElementaryMatrix(*this);
    
    // construction of FacesDof (we also check that dofs are correctly numbered)
    // FindDofsOnFace(*this);
  }
  
  
  //! construction of basis functions
  void PyramidHdivFirstFamily::ConstructFunctions()
  {
  }
   
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void PyramidHdivFirstFamily::ComputeValuesPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    ComputeValuesPhiOrthoRef(point_loc, res);
  }
  
  
  //! Evaluating divergence of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res divergence of basis functions on point_loc
  */
  void PyramidHdivFirstFamily::ComputeDivPhiRef(const R3& point_loc, VectReal_wp& res) const
  {
    res.Reallocate(nb_dof_loc);
    ComputeDivPhiOrthoRef(point_loc, res);
  }

  
  //! Evaluating nearly orthogonal functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of nearly orthogonal functions on pointloc
  */
  void PyramidHdivFirstFamily::ComputeValuesPhiOrthoRef(const R3& pointloc, VectR3& phi) const
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
    
    // part B_{r-1}^3
    Real_wp val, vloc; int num = 0;
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
        {
          int m = max(i, j);
	  // int m = i+j;
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
	      
	      /*if (i+j+k == r-1)
		{
		  phi(num).Init(x*vloc, y*vloc, z*vloc);
		  num++;
		  } */
            }
        }
    
    for (int k = 0; k < r; k++)
      for (int m = 0; m <= k; m++)
	//if (k+m<=r-1)
	  {
	    val = Px(k+1)*Py(m)*pow_oneMinusZ(k);
	    phi(num++).Init(val, 0, 0);
	    
	    val = Px(m)*Py(k+1)*pow_oneMinusZ(k);
	    phi(num++).Init(0, val, 0);
	  }
    
    for (int k = 0; k < r; k++)
      for (int m = 0; m <= k; m++)
	{
	  if (k < r-1)
	    {
	      val = Px(m)*Py(k+1)*pow_oneMinusZ(k);
	      phi(num++).Init(a*val, 0, -val);
	      
	      val = Px(k+1)*Py(m)*pow_oneMinusZ(k);
	      phi(num++).Init(0, b*val, -val);	  
	    }
	}
    
    
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
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
  void PyramidHdivFirstFamily::
  ComputeDivPhiOrthoRef(const R3& pointloc, VectReal_wp& div_phi) const
  {
    div_phi.Reallocate(nb_dof_loc);
    div_phi.Fill(0);
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
  void PyramidHdivFirstFamily::
  ComputeGaussIntegralSurfaceGen(const Vector1& feval, Vector2& res, int num_loc) const
  {
    ComputeIntegralSurfaceRef(feval, res, num_loc);
  }


  //! displays details of class PyramidHcurlFirstFamily
  ostream& operator <<(ostream& out, const PyramidHdivFirstFamily& e)
  {
    out << static_cast< const PyramidReference<3>& >(e);
    return  out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_HDIV_FIRST_FAMILY_CXX
#endif
