#ifndef MONTJOIE_FILE_TETRAHEDRON_MASSLUMPED_CXX

namespace Montjoie
{
  
  //! default constructor
  TetrahedronMassLumped::TetrahedronMassLumped() : TetrahedronReference<1>()
  {
    quadrature_equal_nodal = false;
    dof_equal_nodal = false;
    dof_equal_quadrature = true;
    
    mass_lumping = true;
    diagonal_mass = true;
  }

  
  //! how to number mesh
  void TetrahedronMassLumped::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    if (dg != ElementReference_Base::CONTINUOUS))
      return TetrahedronReference<1>::ConstructNumber(nmap, dg);
  
    nmap.SetNbDofVertex(order, 1);
    nmap.SetNbDofEdge(order, order-1);
    nmap.SetNbDofQuadrangle(order, 0);
    nmap.SetNbDofTriangle(order, PointsDof2D_Inside.GetM());
    nmap.SetNbDofTetrahedron(order, nb_dof_loc - nb_dof_boundaries);
    
    nmap.SetOppositeEdgesDofSymmetry(order, order-1);
    
    // rotation of dofs on faces
    Matrix<int> FacesDof_Rotation_Tri;
    MeshNumbering<Dimension3>::
      GetRotationTriangularFace(PointsDof2D_Inside, FacesDof_Rotation_Tri);
    
    nmap.SetFacesDofRotationTri(order, FacesDof_Rotation_Tri);
    
  }
  

  //! constructing finite element
  void TetrahedronMassLumped::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						     int rsurf_tri, int rsurf_quad,
						     int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronReference<1>::ConstructFiniteElement(r, rgeom, rquad, type_quad);
    
    ConstructFunctions();
    
    ConstructElementaryMatrix(*this);

    ComputeCoefficientTransformation();
  }


  //! returns finite element that is the restriction 
  //! of the current finite element on a triangular face
  const TriangleClassical& TetrahedronMassLumped::GetTriangularSurfaceFiniteElement() const
  {
    return function_basis_tri;
  }
  

  //! returns finite element that is the restriction 
  //! of the current finite element on a quadrilateral face
  const TriangleClassical& TetrahedronMassLumped::GetQuadrangularSurfaceFiniteElement() const
  {
    return function_basis_tri;
  }

  
  //! construction of basis functions
  void TetrahedronMassLumped::ConstructFunctions()
  {
    points3d = points_nodal3d;
    if (order == 1)
      {
	nb_dof_loc = 4;
	nb_dof_boundaries = 4;
	// we keep same points
	weights3d.Reallocate(4);
	Real_wp poids = Real_wp(1)/24;
	weights3d.Fill(poids);
      }
    else if (order == 2)
      {
	nb_dof_loc = 23;
	nb_dof_boundaries = 22;

	points3d.Resize(nb_dof_loc);
	weights3d.Reallocate(nb_dof_loc);
	FacesDof.Resize(9, 4);
	
	PointsDof2D_Inside.Reallocate(3);
	Real_wp alpha = (Real_wp(7) - sqrt(Real_wp(13)))/18;
	Real_wp weight_vertex = (Real_wp(13) - 3*sqrt(Real_wp(13)))/10080;
	Real_wp weight_edge = (Real_wp(4) - sqrt(Real_wp(13)))/315;
	Real_wp weight_alpha = (Real_wp(29) + 17*sqrt(Real_wp(13)))/10080;
	Real_wp weight_center = Real_wp(16)/315;
	Real_wp beta = Real_wp(1) - 2*alpha;
	PointsDof2D_Inside(0).Init(alpha, alpha);
	PointsDof2D_Inside(1).Init(beta, alpha);
	PointsDof2D_Inside(2).Init(alpha, beta);
	
	// points on first face z = 0
	points3d(10).Init(alpha, alpha, 0);
	points3d(11).Init(beta, alpha, 0);
	points3d(12).Init(alpha, beta, 0);
	FacesDof(6, 0) = 10; FacesDof(7, 0) = 11; FacesDof(8, 0) = 12;
	
	// points on second face y = 0
	points3d(13).Init(alpha, 0, alpha);
	points3d(14).Init(beta, 0, alpha);
	points3d(15).Init(alpha, 0, beta);
	FacesDof(6, 1) = 13; FacesDof(7, 1) = 14; FacesDof(8, 1) = 15;
	
	// points on third face x = 0
	points3d(16).Init(0, alpha, alpha);
	points3d(17).Init(0, beta, alpha);
	points3d(18).Init(0, alpha, beta);
	FacesDof(6, 2) = 16; FacesDof(7, 2) = 17; FacesDof(8, 2) = 18;
	
	// points on fourth face x + y + z = 1
	points3d(19).Init(beta, alpha, alpha);
	points3d(20).Init(alpha, beta, alpha);
	points3d(21).Init(alpha, alpha, beta);
	FacesDof(6, 3) = 19; FacesDof(7, 3) = 20; FacesDof(8, 3) = 21;
	
	// center of tet
	points3d(22).Init(0.25, 0.25, 0.25);
	
	for (int i = 0; i < 4; i++)
	  weights3d(i) = weight_vertex;
	
	for (int i = 0; i < 6; i++)
	  weights3d(4+i) = weight_edge;
	
	for (int i = 0; i < 12; i++)
	  weights3d(10+i) = weight_alpha;
	
	weights3d(22) = weight_center;
      }
    else 
      {
	cout<<"Not implemented "<<endl;
	abort();
      }
    
    points_dof3d = points3d;
    FillPositionDofBoundaries(FacesDof, this->power_two_face, this->PosDofOnFace);
    
    // we add quadrature points on the faces
    nb_points_quadrature_inside = nb_dof_loc;
    points3d.Resize(nb_dof_loc + 4*points2d_tri.GetM());
    int nb = nb_points_quadrature_inside;
    this->num_quad_points_surf.Reallocate(4);
    for (int n = 0; n < 4; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_tri.GetM());
    
    // Points of integration for the first face z = 0
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(0)(i) = nb;
	points3d(nb++).Init(points2d_tri(i)(0), points2d_tri(i)(1), Real_wp(0));
      }
    
    // Points of integration for the second face y = 0
    for (int i = 0;i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(1)(i) = nb;
	points3d(nb++).Init(points2d_tri(i)(0), Real_wp(0), points2d_tri(i)(1));
      }
    
    // Points of integration for the third face x = 0
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(2)(i) = nb;
	points3d(nb++).Init(Real_wp(0), points2d_tri(i)(0), points2d_tri(i)(1));
      }
    
    // Points of integration for the fourth face x+y+z = 1
    for (int i = 0; i < points2d_tri.GetM(); i++)
      {
	this->num_quad_points_surf(3)(i) = nb;
	points3d(nb++).Init(Real_wp(1)-points2d_tri(i)(0)-points2d_tri(i)(1),
			    points2d_tri(i)(0), points2d_tri(i)(1));
      }
    
    ComputeLagrangianFunctions();

    // on teste la formule de quadrature
    for (int p = 0; p <= 2*order; p++)
      {
	bool integration_ok = true;
	for (int i = 0; i <= p; i++)
	  for (int j = 0; j <= p-i; i++)
	    {
	      int k = p-i-j;
	      // on integre x^i y^j z^k
	      Real_wp val_exact = tgamma(Real_wp(i+1))*tgamma(Real_wp(j+1))
                *tgamma(Real_wp(k+1))/tgamma(Real_wp(p+4));
	      Real_wp val_app = 0;
	      for (int m = 0; m < nb_points_quadrature_inside; m++)
		val_app += weights3d(m)*pow(points3d(m)(0), i)*pow(points3d(m)(1), j)
                  *pow(points3d(m)(2), k);
	    
	      if (abs(val_app-val_exact) > 1e5*epsilon_machine)
		{
		  DISP(val_app); DISP(val_exact); DISP(i); DISP(j); DISP(k);
		  integration_ok = false;
		}
	    }
	
	if (!integration_ok)
	  {
	    cout<<"L'ordre "<<p<<" n'est pas integre exactement"<<endl;
	    break;
	  }
      }
    
  }
  
  
  //! computing basis functions
  void TetrahedronMassLumped::ComputeLagrangianFunctions()
  {
    Vector<MultivariatePolynomial<Real_wp> > Monomial;
    // specification of basis functions
    if (order == 1)
      {
	MultivariatePolynomial<Real_wp> x(3,1), y(3,1), z(3,1), L0, one(3,0);
	one(0,0,0) = 1.0; x(1,0,0) = 1.0; y(0,1,0) = 1.0; z(0,0,1) = 1.0;
	Monomial.Reallocate(4);
	Monomial(0) = one; Monomial(1) = x; Monomial(2) = y; Monomial(3) = z;
      }
    else if (order == 2)
      {
	// volume functions
	MultivariatePolynomial<Real_wp> x(3,1), y(3,1), z(3,1), L0, one(3,0);
	one(0,0,0) = 1.0; x(1,0,0) = 1.0; y(0,1,0) = 1.0; z(0,0,1) = 1.0;
	L0 = one - x - y - z;
	Monomial.Reallocate(23);
	// first P2
	Monomial(0) = one; Monomial(1) = x; Monomial(2) = y;
	Monomial(3) = z; Monomial(4) = x*x; Monomial(5) = y*y;
	Monomial(6) = z*z; Monomial(7) = x*y; Monomial(8) = x*z;
	Monomial(9) = y*z;
	// then x y (1-x-y-z) P1(x, y)
	Monomial(10) = x*y*L0; Monomial(11) = x*y*L0*x; Monomial(12) = x*y*L0*y;
	// then x z (1-x-y-z) P1(x, z)
	Monomial(13) = x*z*L0; Monomial(14) = x*z*L0*x; Monomial(15) = x*z*L0*z;
	// then y z (1-x-y-z) P1(y, z)
	Monomial(16) = y*z*L0; Monomial(17) = y*z*L0*y; Monomial(18) = y*z*L0*z;
	// then x y z P1(x, y)
	Monomial(19) = x*y*z; Monomial(20) = x*y*z*x; Monomial(21) = x*y*z*y;
	// then bubble function x y z (1-x-y-z)
	Monomial(22) = x*y*z*L0;
      }
    else
      {
	cout<<"This order has not been implemented "<<endl;
	abort();
      }
    
    // computation of VDM
    Matrix<Real_wp> VDM(nb_dof_loc, nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	VDM(i, j) = Monomial(i).Evaluate(points_dof3d(j));
    
    // inversion of VDM
    GetInverse(VDM);
    
    // basis functions
    BasisFunction.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      {
	BasisFunction(i).Zero();
	for (int j = 0; j < nb_dof_loc; j++)
	  Add(VDM(i, j), Monomial(j), BasisFunction(i));
      }

    // checking delta_ij
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  Real_wp phi = BasisFunction(i).Evaluate(points_dof3d(j));
	  Real_wp phi_ref = 0;
	  if (i == j)
	    phi_ref = 1.0;
	  
	  if (abs(phi-phi_ref) > 1e4*epsilon_machine)
	    {
	      DISP(i); DISP(j); DISP(phi);
	      abort();
	    }
	}
    
    // derivate
    BasisDx_Function.Reallocate(nb_dof_loc);
    BasisDy_Function.Reallocate(nb_dof_loc);
    BasisDz_Function.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      {
	DerivatePolynomial(BasisFunction(i), BasisDx_Function(i), 0);
	DerivatePolynomial(BasisFunction(i), BasisDy_Function(i), 1);
	DerivatePolynomial(BasisFunction(i), BasisDz_Function(i), 2);
      }
  }
  
  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res values of basis functions on point_loc
  */
  void TetrahedronMassLumped::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    phi.Reallocate(nb_dof_loc);
    for (int num_dof = 0; num_dof < phi.GetM(); num_dof++)
      phi(num_dof) = BasisFunction(num_dof).Evaluate(point_loc);
  }

  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] point_loc local point where functions are evaluated
    \param[out] res gradient of basis functions on point_loc
  */
  void TetrahedronMassLumped::ComputeGradientPhiRef(const R3& point_loc, VectR3& res) const
  {
    res.Reallocate(nb_dof_loc);
    for (int num_dof = 0; num_dof < res.GetM(); num_dof++)
      {
	res(num_dof).Init(BasisDx_Function(num_dof).Evaluate(point_loc),
			  BasisDy_Function(num_dof).Evaluate(point_loc),
			  BasisDz_Function(num_dof).Evaluate(point_loc));
      }
  }
  
  
  //! displays details about class TetrahedronMassLumped
  ostream& operator <<(ostream& out, const TetrahedronMassLumped& e)
  {
    return out;
  }
  
} // end namespace
  
#define MONTJOIE_FILE_TETRAHEDRON_MASSLUMPED_CXX
#endif
