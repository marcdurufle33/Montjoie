#ifndef MONTJOIE_FILE_PYRAMID_REFERENCE_SPLIT_CXX

namespace Montjoie
{
  //! default constructor
  PyramidReferenceSplit::PyramidReferenceSplit() : VolumeReference()
  {
    type_interpolation = LOBATTO_ELEMENT;
    // type_interpolation = REGULAR_ELEMENT;
    quadrature_equal_nodal = false;
    dof_equal_nodal = true;
    dof_equal_quadrature = false;

    // type_basis = WALKER;
    // type_basis = WIENERS;
    type_basis = LIU;
    
    // normales of the five faces
    normale.Reallocate(5);
    normale(0).Init(0,0,-1); normale(1).Init(0,-1,0);
    normale(2).Init(1,0,1); normale(3).Init(0,1,1); normale(4).Init(-1,0,0);
    
    if (type_basis == LIU)
      {
	normale.Reallocate(5);
	normale(0).Init(0,0,-4); normale(1).Init(0,-2,2);
	normale(2).Init(2,0,2); normale(3).Init(0,2,2); normale(4).Init(-2,0,2);
      }
  }
  
  
  /****************************
   * Initialization functions *
   ****************************/
  
  
  //! 5 vertices
  inline int PyramidReferenceSplit::GetNbVertices()
  { 
    return 5;
  }


  //! 5 faces
  inline int PyramidReferenceSplit::GetNbBoundaries()
  { 
    return 5;
  }
  
  
  inline int PyramidReferenceSplit::GetNbDofBoundary(int i) const
  {
    if (i == 0)
      return nb_dof_quad;
    
    return nb_dof_tri;
  }
  
  
  inline int PyramidReferenceSplit::GetNbQuadBoundary(int i) const
  {
    if (i == 0)
      return points2d_quad.GetM();
    
    return points2d_tri.GetM();
  }


  inline int PyramidReferenceSplit::GetNbNodalBoundary(int i) const
  {
    if (i == 0)
      return points_nodal2d_quad.GetM();
    
    return points_nodal2d_tri.GetM();
  }
  

  /*int PyramidReferenceSplit::GetNbPointsNodal(TinyVector<Real_wp, 3>& x) const
  {
    return points_nodal3d.GetM();
    }*/

  
  //!< 2-D integration point
  inline const R2& PyramidReferenceSplit::PointsQuadratureBoundary(int k, int nf) const
  {
    if (nf == 0)
      return points2d_quad(k);
    
    return points2d_tri(k);
  }
  
  
  //!< 2-D integration weight
  inline const Real_wp& PyramidReferenceSplit::WeightsQuadratureBoundary(int k, int nf) const
  {
    if (nf == 0)
      return weights2d_quad(k);
    
    return weights2d_tri(k);
  }
  

  inline const R2& PyramidReferenceSplit::PointsNodalBoundary(int k, int nf) const
  {
    if (nf == 0)
      return points_nodal2d_quad(k);
    
    return points_nodal2d_tri(k);
  }
  
  
  inline int PyramidReferenceSplit::GetNumNodes2D(int i, int j, int num_loc) const
  {
    if (num_loc > 0)
      return NumNodes2D_tri(i, j);
    
    return NumNodes2D_quad(i,j);
  }
  
  
  void PyramidReferenceSplit::ConstructNumberMap(NumberMap& nmap)
  {
    if (type_basis == WALKER)
      nmap.SetNbDofPyramid(order, nb_dof_loc - 3*order*order - 2);
    else
      nmap.SetNbDofPyramid(order, 0);
    
    nmap.SetNbPointsQuadBoundariesPyramid(order, nb_points_quadrature_boundaries);
  }

  
  void PyramidReferenceSplit::ConstructFiniteElement(int r)
  {
    order = r;
    ConstructQuadrature(r);
    ConstructNodalPoints();
    ComputeCoefficientTransformation();
    
    int nb_points_nodal_elt = points_nodal3d.GetM();
    VectReal_wp phi; VectR3 grad_phi;
    // Checking the base
    // phi_i(x_j) = delta_ij
    for (int i = 0; i < nb_points_nodal_elt; i++)
      {
	ComputeValuesPhiNodalRef(points_nodal3d(i), phi);
	for (int j = 0; j < nb_points_nodal_elt; j++)
	  {
	    if ((i==j)&&(abs(phi(j)-1) > 1e5*epsilon_machine))
	      {
		DISP(i); DISP(phi(j));
		abort();
	      }
	    
	    if ((i!=j)&&(abs(phi(j))>1e5*epsilon_machine))
	      {
		DISP(i); DISP(j); DISP(phi(j));
		abort();
	      }
	  }
      }
    
    VectReal_wp phi_xm, phi_xp, phi_ym, phi_yp, phi_zm, phi_zp; R3 pt;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      {
	Real_wp h = pow(epsilon_machine, 1.0/3.0);
	ComputeGradientPhiNodalRef(points3d(j), grad_phi);
	
	pt = points3d(j); pt(0) -= h;
	ComputeValuesPhiNodalRef(pt, phi_xm);

	pt = points3d(j); pt(0) += h;
	ComputeValuesPhiNodalRef(pt, phi_xp);

	pt = points3d(j); pt(1) -= h;
	ComputeValuesPhiNodalRef(pt, phi_ym);

	pt = points3d(j); pt(1) += h;
	ComputeValuesPhiNodalRef(pt, phi_yp);

	pt = points3d(j); pt(2) -= h;
	ComputeValuesPhiNodalRef(pt, phi_zm);

	pt = points3d(j); pt(2) += h;
	ComputeValuesPhiNodalRef(pt, phi_zp);
	
	for (int i = 0; i < nb_points_nodal_elt; i++)
	  {
	    R3 grad_num;
	    grad_num(0) = (phi_xp(i) - phi_xm(i))/(2.0*h);
	    grad_num(1) = (phi_yp(i) - phi_ym(i))/(2.0*h);
	    grad_num(2) = (phi_zp(i) - phi_zm(i))/(2.0*h);
	    
	    if (grad_num.Distance(grad_phi(i)) > 10.0*sqrt(epsilon_machine))
	      abort();
	  }
	
      }
  }
  
  
  //! construction of quadrature points and nodal points
  void PyramidReferenceSplit::ConstructQuadrature(int order_)
  {
    // quadrature points on the volume
    Globatto<Dimension3> gauss_tet;
    // gauss_r.ConstructQuadrature(order, gauss_r.QUADRATURE_PYRAMID_RADAU_JACOBI);
    gauss_tet.ConstructQuadrature(2*order, gauss_tet.QUADRATURE_TETRAHEDRON_GAUSS);
    if (type_basis == LIU)
      {
	// four tetrahedra
	int N = gauss_tet.GetNbPointsQuad();
	points3d.Reallocate(4*N); weights3d.Reallocate(4*N);
	Real_wp x, y, z; 
	for (int i = 0; i < N; i++)
	  {
	    x = gauss_tet.Points3D(i)(0);
            y = gauss_tet.Points3D(i)(1); z = gauss_tet.Points3D(i)(2);
	    points3d(i).Init(-x+y, -x-y, z); weights3d(i) = 2.0*gauss_tet.Weights3D(i);
	    points3d(i+N).Init(x+y, -x+y, z); weights3d(i+N) = 2.0*gauss_tet.Weights3D(i);
	    points3d(i+2*N).Init(x-y, x+y, z); weights3d(i+N) = 2.0*gauss_tet.Weights3D(i);
	    points3d(i+3*N).Init(-x-y, x-y, z); weights3d(i+N) = 2.0*gauss_tet.Weights3D(i);
	  }
      }
    else
      {
	int N = gauss_tet.GetNbPointsQuad();
	points3d.Reallocate(2*N); weights3d.Reallocate(2*N);
	Real_wp x, y, z;
	for (int i = 0; i < N; i++)
	  {
	    x = gauss_tet.Points3D(i)(0);
            y = gauss_tet.Points3D(i)(1); z = gauss_tet.Points3D(i)(2);
	    points3d(i).Init(x+y, y, z); weights3d(i) = gauss_tet.Weights3D(i);
	    points3d(i+N).Init(y, x+y, z); weights3d(i+N) = gauss_tet.Weights3D(i);
	  }
      }
    
    Matrix<int> coor;
    NumberedMesh<Dimension2>::ConstructQuadrilateralNumbering(order, NumNodes2D_quad, coor);
    NumberedMesh<Dimension2>::ConstructTriangularNumbering(order, NumNodes2D_tri, coor);
    
    nb_points_quadrature_inside = points3d.GetM();
    
    // quadrature points on the triangle
    Globatto<Dimension2> face_gauss;
    face_gauss.ConstructQuadrature(order, face_gauss.QUADRATURE_TRIANGLE_GAUSS);
    points2d_tri = face_gauss.Points2D();
    weights2d_tri = face_gauss.Weights2D();
        
    // quadrature points on the base
    face_gauss.ConstructQuadrature(order, face_gauss.QUADRATURE_GAUSS);
    points2d_quad.Reallocate((order+1)*(order+1));
    weights2d_quad.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  points2d_quad(NumNodes2D_quad(i, j)).Init(face_gauss.Points(i), face_gauss.Points(j));
	  weights2d_quad(NumNodes2D_quad(i,j)) = face_gauss.Weights(i)*face_gauss.Weights(j);
	}
    
    // quadrature points on all the faces
    nb_points_quadrature_boundaries = points2d_quad.GetM() + 4*points2d_tri.GetM();
    
    // we add quadrature points of the face at the end of Points3D
    // Points of integration for the first face z = 0 (square)
    Real_wp one(1);
    if (type_basis == LIU)
      for (int i = 0; i < points2d_quad.GetM(); i++)
	points3d.PushBack(R3(2.0*points2d_quad(i)(0)-one,
                             2.0*points2d_quad(i)(1)-one, Real_wp(0)));
    else
      for (int i = 0; i < points2d_quad.GetM(); i++)
	points3d.PushBack(R3(points2d_quad(i)(0), points2d_quad(i)(1), Real_wp(0)));
    
    // Points of integration for the second face y = 0 (triangle)
    if (type_basis == LIU)
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(2.0*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			     points2d_tri(i)(1) - one, points2d_tri(i)(1)));
    else
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(points2d_tri(i)(0) , Real_wp(0), points2d_tri(i)(1)));
    
    // Points of integration for the third face z+x = 1 (triangle)
    if (type_basis == LIU)
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(one - points2d_tri(i)(1),
			     2.0*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			     points2d_tri(i)(1)));
    else
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(one - points2d_tri(i)(1), points2d_tri(i)(0), points2d_tri(i)(1)));
    
    // Points of integration for the fourth face z+y = 1 (triangle)
    if (type_basis == LIU)
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(2.0*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			     one - points2d_tri(i)(1), points2d_tri(i)(1)));
    else
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(points2d_tri(i)(0), one - points2d_tri(i)(1), points2d_tri(i)(1)));
    
    // Points of integration for the fifth face x = 0 (triangle)
    if (type_basis == LIU)
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(points2d_tri(i)(1) - Real_wp(1),
			     2.0*points2d_tri(i)(0) + points2d_tri(i)(1) - one,
			     points2d_tri(i)(1)));
    else
      for (int i = 0; i < points2d_tri.GetM(); i++)
	points3d.PushBack(R3(Real_wp(0), points2d_tri(i)(0), points2d_tri(i)(1)));
    
    // constructing array this->num_quad_points_surf
    // to get numbers of quadrature points of the faces
    int ind = nb_points_quadrature_inside;
    this->num_quad_points_surf.Reallocate(5);
    this->num_quad_points_surf(0).Reallocate(points2d_quad.GetM());
    for (int n = 1; n < 5; n++)
      this->num_quad_points_surf(n).Reallocate(points2d_tri.GetM());
    
    // first face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      this->num_quad_points_surf(0)(i) = ind++;
    
    // second face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(1)(i) = ind++;
    
    // third face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(2)(i) = ind++;
    
    // fourth face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(3)(i) = ind++;
      
    // fifth face
    for (int i = 0; i < points2d_tri.GetM(); i++)
      this->num_quad_points_surf(4)(i) = ind++;
    
    offset_faceSh.Reallocate(6);
    offset_faceSh(0) = 0;
    for (int i = 1; i < 5; i++)
      offset_faceSh(i) = points2d_quad.GetM() + (i-1)*points2d_tri.GetM();
    
    offset_faceSh(5) = nb_points_quadrature_boundaries;
  }
  
  
  //! construction of nodal points, which are used for interpolation
  void PyramidReferenceSplit::ConstructNodalPoints()
  {
    // number of regular points
    Mesh<Dimension3>::ConstructPyramidalNumbering(order, NumNodes3D,
						  CoordinateNodes);
    
    EdgesNodal.Reallocate(order+1, 8);
    for (int i = 0; i <= order; i++)
      {
	EdgesNodal(i, 0) = NumNodes3D(i, 0, 0);
	EdgesNodal(i, 1) = NumNodes3D(order, i, 0);
	EdgesNodal(i, 2) = NumNodes3D(i, order, 0);
	EdgesNodal(i, 3) = NumNodes3D(0, i, 0);
	EdgesNodal(i, 4) = NumNodes3D(0, 0, i);
	EdgesNodal(i, 5) = NumNodes3D(order-i, 0, i);
	EdgesNodal(i, 6) = NumNodes3D(order-i, order-i, i);
	EdgesNodal(i, 7) = NumNodes3D(0, order-i, i);
      }
    
    FacesDof.Reallocate((order+1)*(order+1), 5);
    FacesDof.Fill(-1);
    // quadrangular face
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  int node = NumNodes2D_quad(i,j);
	  FacesDof(node, 0) = NumNodes3D(i, j, 0);
	}
    
    // triangular faces
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	{
	  int node = NumNodes2D_tri(i,j);
	  FacesDof(node,1) = NumNodes3D(i, 0, j);
	  FacesDof(node,2) = NumNodes3D(order-j, i, j);
	  FacesDof(node,3) = NumNodes3D(i, order-j, j);
	  FacesDof(node,4) = NumNodes3D(0, i, j);
	}
    
    FacesNodal = FacesDof;

    if (type_interpolation == REGULAR_ELEMENT)
      ConstructRegularNodalPoints();
    else
      ConstructLobattoPoints();
    
    ComputeLagrangianFunctions();
    
    if (order > 2)
      NumberedMesh<Dimension3>::
	GetRotationTriangularFace(points_nodal2d_inside_tri, FacesNodal_Rotation_Tri);
    
    if (order > 1)
      {
	NumberedMesh<Dimension3>::
	  GetRotationQuadrilateralFace(NumNodes2D_InsideQuad, FacesNodal_Rotation_Quad);
      }
    
    points_dof3d = points_nodal3d;
    
   }


  //! regularly spaced nodal points
  void PyramidReferenceSplit::ConstructRegularNodalPoints()
  {    
    
  }
  
  
  //! computation of coefFicients to compute quickly Fi and DFi
  void PyramidReferenceSplit::ComputeCoefficientTransformation()
  {
    VectR3 all_points = points_nodal3d;
    all_points.PushBack(points3d);
    all_points.PushBack(points_dof3d);
    int nb_points_all = all_points.GetM();
    
    // straight pyramid
    coefFi.Reallocate(nb_points_all, 5);
    coefDFi_dx.Reallocate(nb_points_all, 5);
    coefDFi_dy.Reallocate(nb_points_all, 5);
    coefDFi_dz.Reallocate(nb_points_all, 5);
    VectReal_wp phi(5);
    VectR3 dphi(5);
	
    Globatto<Real_wp> lob1D;
    lob1D.ConstructQuadrature(order, lob1D.QUADRATURE_LOBATTO);
    for (int i = 0; i < nb_points_all; i++) 
      {
        ComputeValuesPhiFirstOrder(all_points(i), phi);   
        ComputeGradientPhiFirstOrder(all_points(i), dphi);  
        for (int j = 0; j < 5; j++)
	  {
	    coefFi(i,j) = phi(j);
	    coefDFi_dx(i,j) = dphi(j)(0);
	    coefDFi_dy(i,j) = dphi(j)(1);
	    coefDFi_dz(i,j) = dphi(j)(2);
	  }
      }
    
    int nb_points_nodal_elt = points_nodal3d.GetM();
    coefFi_curve.Reallocate(nb_points_all, nb_points_nodal_elt);
    coefFi_curve.Fill(0);
    Real_wp x, y, z, xs, ys, zeta, coef, coefs;
    Real_wp one(1), half(0.5), quarter(0.25);
    for (int i = 3*order*order+2; i < nb_points_nodal_elt; i++)
      {
	x = points_nodal3d(i)(0);
	y = points_nodal3d(i)(1);
	z = points_nodal3d(i)(2);
	
	// linear part
	ComputeValuesPhiFirstOrder(points_nodal3d(i), phi);
	for (int k = 0; k < 5; k++)
	  coefFi_curve(i, k) = phi(k);
	
	// first edge
	zeta = half + half*x/(one - z);
	coef = half*(one - y - z);
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 0)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) -= (one-zeta)*coef;
	coefFi_curve(i, 1) -= zeta*coef;
	
	// second edge
	zeta = half + half*y/(one - z);
	coef = half*(one + x - z);
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 1)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) -= (one-zeta)*coef;
	coefFi_curve(i, 2) -= zeta*coef;

	// third edge
	zeta = half + half*x/(one-z);
	coef = half*(one + y - z);
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 2)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) -= (one-zeta)*coef;
	coefFi_curve(i, 2) -= zeta*coef;
	
	// fourth edge
	zeta = half + half*y/(one-z);
	coef = half*(one - x - z);
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 3)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) -= (one-zeta)*coef;
	coefFi_curve(i, 3) -= zeta*coef;
	
	// fifth edge
	zeta = z;
	coef = quarter*(one-x-z)*(one-y-z)/((one-z)*(one-z));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 4)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// sixth edge
	zeta = z;
	coef = quarter*(one+x-z)*(one-y-z)/((one-z)*(one-z));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 5)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// seventh edge
	zeta = z;
	coef = quarter*(one+x-z)*(one+y-z)/((one-z)*(one-z));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 6)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 2) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// eigth edge
	zeta = z;
	coef = quarter*(one-x-z)*(one+y-z)/((one-z)*(one-z));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 7)) += coef*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) -= (one-zeta)*coef;
	coefFi_curve(i, 4) -= zeta*coef;
	
	// quadrangular base
	xs = half + half*x/(one-z);
	ys = half + half*y/(one-z);
	coef = one-z;
	for (int j = 0; j <= order; j++)
	  for (int k = 0; k <= order; k++)
	    coefFi_curve(i, NumNodes3D(j, k, 0)) += coef*lob1D.EvaluatePhi(j, xs)*lob1D.EvaluatePhi(k, ys);

	for (int k = 0; k <= order; k++)
	  {
	    coefFi_curve(i, EdgesNodal(k, 0)) -= coef*(one-ys)*lob1D.EvaluatePhi(k, xs);
	    coefFi_curve(i, EdgesNodal(k, 1)) -= coef*xs*lob1D.EvaluatePhi(k, ys);
	    coefFi_curve(i, EdgesNodal(k, 2)) -= coef*ys*lob1D.EvaluatePhi(k, xs);
	    coefFi_curve(i, EdgesNodal(k, 3)) -= coef*(one-xs)*lob1D.EvaluatePhi(k, ys);
	  }
	
	coefFi_curve(i, 0) += coef*(one-xs)*(one-ys);
	coefFi_curve(i, 1) += coef*xs*(one-ys);
	coefFi_curve(i, 2) += coef*xs*ys;
	coefFi_curve(i, 3) += coef*(one-xs)*ys;

	// first triangular face
	xs = half*(one+x-z);
	ys = z;
	coef = half*(one-y-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, FacesNodal(j, 1)) += coef*phi(j);
	
	coefFi_curve(i, 0) -= (one-xs-ys)*coef;
	coefFi_curve(i, 1) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 0)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 1) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 5)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 4)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	// second triangular face
	xs = half*(one+y-z);
	ys = z;
	coef = half*(one+x-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, FacesNodal(j, 2)) += coef*phi(j);
	
	coefFi_curve(i, 1) -= (one-xs-ys)*coef;
	coefFi_curve(i, 2) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 1)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 2) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 6)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 2) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 5)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 1) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	// third triangular face
	xs = half*(one+x-z);
	ys = z;
	coef = half*(one+y-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, FacesNodal(j, 3)) += coef*phi(j);
	
	coefFi_curve(i, 3) -= (one-xs-ys)*coef;
	coefFi_curve(i, 2) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 2)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 2) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 7)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 6)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 2) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	// fourth triangular face
	xs = half*(one+y-z);
	ys = z;
	coef = half*(one-x-z)/(one-z);
	function_tri.ComputeValuesPhiNodalRef(R2(xs, ys), phi);
	for (int j = 0; j < points_nodal2d_tri.GetM(); j++)
	  coefFi_curve(i, FacesNodal(j, 4)) += coef*phi(j);
	
	coefFi_curve(i, 0) -= (one-xs-ys)*coef;
	coefFi_curve(i, 3) -= xs*coef;
	coefFi_curve(i, 4) -= ys*coef;
	
	zeta = xs + half*ys;
	coefs = (one-xs-ys)*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 3)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 3) += zeta*coef*coefs;
	
	zeta = half*ys - half*xs +half;
	coefs = ys*xs/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 7)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 3) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
	zeta = one - ys - half*xs;
	coefs = (one-xs-ys)*ys/(zeta*(one-zeta));
	for (int k = 0; k <= order; k++)
	  coefFi_curve(i, EdgesNodal(k, 4)) -= coef*coefs*lob1D.EvaluatePhi(k, zeta);
	
	coefFi_curve(i, 0) += (one-zeta)*coef*coefs;
	coefFi_curve(i, 4) += zeta*coef*coefs;
	
      }
    
    // for other points, we use lagrangian interpolation
    for (int i = nb_points_nodal_elt; i < nb_points_all; i++)
      {
	ComputeValuesPhiNodalRef(all_points(i), phi);
	for (int j = 0; j < nb_points_nodal_elt; j++)
	  coefFi_curve(i, j) = phi(j);
      }
    
    coefDFi_dx_curve.Reallocate(nb_points_all, nb_points_nodal_elt);
    coefDFi_dy_curve.Reallocate(nb_points_all, nb_points_nodal_elt);
    coefDFi_dz_curve.Reallocate(nb_points_all, nb_points_nodal_elt);
    for (int i = 0; i < nb_points_all; i++)
      {
	ComputeGradientPhiNodalRef(all_points(i), dphi);
	for (int j = 0; j < nb_points_nodal_elt; j++)
	  {
	    coefDFi_dx_curve(i, j) = dphi(j)(0);
	    coefDFi_dy_curve(i, j) = dphi(j)(1);
	    coefDFi_dz_curve(i, j) = dphi(j)(2);
	  }
      }

  }
  
  
  //! constructing points on the pyramid 
  void PyramidReferenceSplit::ConstructLobattoPoints()
  {
    points_nodal1d.Reallocate(order+1);
    for (int k = 0; k <= order; k++)
      points_nodal1d(k) = Real_wp(k)/order;
    
    points_nodal2d_quad.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	points_nodal2d_quad(NumNodes2D_quad(i, j)).Init(Real_wp(i)/order, Real_wp(j)/order);
    
    points_nodal2d_tri.Reallocate((order+1)*(order+2)/2);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order-i; j++)
	points_nodal2d_tri(NumNodes2D_tri(i, j)).Init(Real_wp(i)/order, Real_wp(j)/order);
    
    NumNodes2D_InsideQuad.Reallocate(order-1, order-1);
    for (int i = 1;  i < order; i++)
      for (int j = 1;  j < order; j++)
	NumNodes2D_InsideQuad(i-1, j-1) = (i-1)*(order-1) + j-1;
    
    int nb_points_nodal_inside_tri = 0;
    NumNodes2D_InsideTri.Reallocate(order-1, order-1);
    NumNodes2D_InsideTri.Fill(-1);
    points_nodal2d_inside_tri.Reallocate((order-1)*(order-2)/2);
    for (int j = 1; j < order; j++)
      for (int k = 1; k < order-j; k++)
	{
	  int node = NumNodes2D_tri(j, k);
	  points_nodal2d_inside_tri(nb_points_nodal_inside_tri) = points_nodal2d_tri(node);
	  NumNodes2D_InsideTri(j-1, k-1) = nb_points_nodal_inside_tri;
	  nb_points_nodal_inside_tri++;
	}
      
    int nb_points_nodal_elt = (order+1)*(order+2)*(order+3)/6 + order*order;
    if (type_basis != WALKER)
      nb_points_nodal_elt -= order-1;
    
    points_nodal3d.Reallocate(nb_points_nodal_elt);
    for (int k = 0; k <= order; k++)
      for (int i = 0; i <= order-k; i++)
	for (int j = 0; j <= order-k; j++)
	  if (NumNodes3D(i, j, k) < nb_points_nodal_elt)
	    {
	      int node = NumNodes3D(i, j, k);
	      Real_wp x(0), y(0), z;
	      z = Real_wp(k)/order;
	      if (k != order)
		{
		  x = Real_wp(i)/(order-k)*(1.0-z);
		  y = Real_wp(j)/(order-k)*(1.0-z);
		}
	      
	      points_nodal3d(node).Init(x, y, z);
	    }
    
    // for Liu stuff, we have to change for symmetric pyramids
    if (type_basis == LIU)
      {
	Real_wp x, y, z;
	for (int i = 0; i < nb_points_nodal_elt; i++)
	  {
	    x = points_nodal3d(i)(0);
	    y = points_nodal3d(i)(1);
	    z = points_nodal3d(i)(2);
	    points_nodal3d(i).Init(2.0*x+z-1.0, 2.0*y+z-1.0, z);
	  }
      } 
    // points_nodal3d.WriteText("pts.dat");
  }
  
  
  //! computation of basis functions
  void PyramidReferenceSplit::ComputeLagrangianFunctions()
  {
    // 3-D nodal functions
    int nb_points_nodal_elt = points_nodal3d.GetM();
    nb_dof_loc = nb_points_nodal_elt;
    
    int ind = 0;
    if (type_basis == WALKER)
      {
	// vandermonde computation to force delta_ij
	Matrix<Real_wp> VDM(nb_points_nodal_elt, nb_points_nodal_elt);
	VDM.Fill(0);
	// part with Pr
	for (int m = 0; m <= order; m++)
	  for (int n = 0; n <= order-m; n++)
	    for (int p = 0; p <= order-m-n; p++)
	      {
		// monomial x^m y^n z^p
		for (int j = 0; j < nb_points_nodal_elt; j++)
		  VDM(ind, j) = pow(points_nodal3d(j)(0), m)*pow(points_nodal3d(j)(1), n)
                    *pow(points_nodal3d(j)(2), p);
		
		ind++;
	      }
	
	// teta part
	int nb_other = nb_dof_loc - (order+1)*(order+2)*(order+3)/6;
	power_other_coef.Reallocate(nb_other);
	other_coef_lagrange.Reallocate(nb_dof_loc); // DISP(nb_other);
	nb_other = 0;
	for (int n = 0; n <= order; n++)
	  for (int m = n; m <= order; m++)
	    for (int l = n; l <= order; l++)
	      if (l+m-n > order)
		{
		  int q = min(l-n, m-n);
		  DISP(n); DISP(m); DISP(l); DISP(q);
		  // storing powers for monomials in T1
		  power_other_coef(nb_other)(0) = q;
		  power_other_coef(nb_other)(1) = l-n-q;
		  power_other_coef(nb_other)(2) = m-n;
		  power_other_coef(nb_other)(3) = n;
		  Real_wp x, y, z, val(1);
		  
		  for (int j = 0; j < nb_points_nodal_elt; j++)
		    {
		      x = points_nodal3d(j)(0); y = points_nodal3d(j)(1); z = points_nodal3d(j)(2);
		      if (x > y)
			val = pow(x+z, q)*pow(x, l-n-q)*pow(y, m-n)*pow(z, n);
		      else
			val = pow(y+z, q)*pow(x, l-n)*pow(y, m-n-q)*pow(z, n);
		      
		      VDM(ind, j) = val;
		    }
		  
		  ind++; nb_other++;
		}
	
	// DISP(nb_other); DISP(power_other_coef);
	// DISP(points_nodal3d); DISP(VDM); VDM.WriteText("Ah.dat");
	// filling functions with VDM
	/* if (order == 2)
	  {
	    LagrangePr.Reallocate(nb_dof_loc);
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		other_coef_lagrange(j).Reallocate(nb_other);
		other_coef_lagrange(j).Fill(0);       
	      }
	    
	    MultivariatePolynomial<Real_wp> x(3,1), y(3,1), z(3,1), one(3,0);
	    one(0,0,0) = 1.0; x(1,0,0) = 1.0; y(0,1,0) = 1.0; z(0,0,1) = 1.0;
	    LagrangePr(0) = one - 3*x -3*y - 3*z + 2*x*x + 9*x*y + 2*y*y + 4*y*z + 2*z*z + 4*x*z;
	    other_coef_lagrange(0)(0) = -6.0; other_coef_lagrange(0)(1) = 4.0;
            other_coef_lagrange(0)(2) = -1.0; other_coef_lagrange(0)(3) = -6.0;
	    
	    LagrangePr(1) = -1*x + 2*x*x + 3*x*y;
	    other_coef_lagrange(1)(0) = -2.0; other_coef_lagrange(1)(1) = 4.0;
            other_coef_lagrange(1)(2) = 1.0; other_coef_lagrange(1)(3) = -6.0;
	    
	    LagrangePr(2) = x*y;
	    other_coef_lagrange(2)(0) = -2.0; other_coef_lagrange(2)(1) = 4.0;
            other_coef_lagrange(2)(2) = -1.0; other_coef_lagrange(2)(3) = -2.0;
	    
	    LagrangePr(3) = -1*y + 3*x*y + 2*y*y;
	    other_coef_lagrange(3)(0) = -6.0; other_coef_lagrange(3)(1) = 4.0;
            other_coef_lagrange(3)(2) = 1.0; other_coef_lagrange(3)(3) = -2.0;
	    
	    LagrangePr(4) = -1*z + 2*z*z;
	    
	    LagrangePr(5) = 4*x -4*x*x -12*x*y - 4*x*z;
	    other_coef_lagrange(5)(0) = 8.0; other_coef_lagrange(5)(1) = -8.0;
            other_coef_lagrange(5)(2) = 0.0; other_coef_lagrange(5)(3) = -12.0;
	    
	    LagrangePr(6) = 4*y -12*x*y - 4*y*y - 4*y*z;
	    other_coef_lagrange(6)(0) = 12.0; other_coef_lagrange(6)(1) = -8.0;
            other_coef_lagrange(6)(2) = 0.0; other_coef_lagrange(6)(3) = 8.0;
	    
	    LagrangePr(7) = 16*x*y;
	    other_coef_lagrange(7)(0) = -16.0; other_coef_lagrange(7)(1) = 16.0;
            other_coef_lagrange(7)(2) = 0.0; other_coef_lagrange(7)(3) = -16.0;
	    
	    LagrangePr(8) = -4*x*y;
	    other_coef_lagrange(8)(0) = 4.0; other_coef_lagrange(8)(1) = -8.0;
            other_coef_lagrange(8)(2) = 0.0; other_coef_lagrange(8)(3) = 8.0;
	    
	    LagrangePr(9) = -4*x*y;
	    other_coef_lagrange(9)(0) = 8.0; other_coef_lagrange(9)(1) = -8.0;
            other_coef_lagrange(9)(2) = 0.0; other_coef_lagrange(9)(3) = 4.0;
	    
	    LagrangePr(10) = 4*z - 4*y*z - 4*z*z - 4*x*z;
	    other_coef_lagrange(10)(2) = 4.0;
	    
	    LagrangePr(11) = 4*x*z;
	    other_coef_lagrange(11)(2) = -4.0;
	    
	    LagrangePr(12) = 0*one;
	    other_coef_lagrange(12)(2) = 4.0;
	    
	    LagrangePr(13) = 4*y*y;
	    other_coef_lagrange(13)(2) = -4.0;
	  }
	else
	{*/
	    GetInverse(VDM);
	    
	    ind = 0;
	    LagrangePr.Reallocate(nb_dof_loc);
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		LagrangePr(j).SetOrder(3, order);
		LagrangePr(j).Fill(0);
		other_coef_lagrange(j).Reallocate(nb_other);
		other_coef_lagrange(j).Fill(0);       
	      }
	    
	    // filling coefficients
	    for (int m = 0; m <= order; m++)
	      for (int n = 0; n <= order-m; n++)
		for (int p = 0; p <= order-m-n; p++)
		  {
		    for (int j = 0; j < nb_dof_loc; j++)
		      LagrangePr(j)(m, n, p) = VDM(j, ind);
		    
		    ind++;
		  }
	    
	    for (int k = 0; k < nb_other; k++)
	      {
		for (int j = 0; j < nb_dof_loc; j++)
		  other_coef_lagrange(j)(k) = VDM(j, ind);
		
		ind++;
	      }
	    // }
	
	DxLagrangePr.Reallocate(nb_dof_loc); 
	DyLagrangePr.Reallocate(nb_dof_loc); 
	DzLagrangePr.Reallocate(nb_dof_loc);
	
	for (int k = 0; k < nb_dof_loc; k++)
	  {
	    DxLagrangePr(k) = Derivate(LagrangePr(k), 0);
	    DyLagrangePr(k) = Derivate(LagrangePr(k), 1);
	    DzLagrangePr(k) = Derivate(LagrangePr(k), 2);	
	  }
      }
  }
  
  
  /****************
   * Fj transform *
   ****************/
  
  
  //! computes res = Fi(point)
  /*!
    \param[in] s list of vertices of the pyramid
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit pyramid
    \param[out] res the result of transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void PyramidReferenceSplit::Fj(const VectR3& s,const SetPoints<Dimension3>& PTReel,
				   const R3& point, R3& res,
				   const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjCurve(PTReel, point, res);
    else
      FjLinear(s, point, res);
  }
  
  
  //! computes res = DFi(point)
  /*!
    \param[in] s list of vertices of the pyramid
    \param[in] PTReel "reference points" after transformation Fi
    \param[in] point local coordinates on the unit pyramid
    \param[out] res the jacobian matrix
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void PyramidReferenceSplit::DFj(const VectR3& s, const SetPoints<Dimension3>& PTReel,
				    const R3& point, Matrix3_3& res,
				    const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      DFjCurve(PTReel, point, res);
    else
      DFjLinear(s, point, res);
  }
  
  
  //! transformation Fi in the case of straight pyramid
  void PyramidReferenceSplit::FjLinear(const VectR3& s, const R3& point, R3& res) const
  {
    // Fi = phi1*S0 + phi2*S1 + phi3*S2 + phi4*S3 + z*S4
    res.Zero();
    VectReal_wp phi; ComputeValuesPhiFirstOrder(point, phi);

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 3; j++)
        res(j) += phi(i)*s(i)(j);	
  }


  //! transformation DFi in the case of straight pyramid
  inline void PyramidReferenceSplit
  ::DFjLinear(const VectR3& s,const R3& point,Matrix3_3& res) const
  {
    res.Zero();
    VectR3 dphi; ComputeGradientPhiFirstOrder(point, dphi);

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 3; j++)
	{
	  res(j,0) += dphi(i)(0)*s(i)(j);
	  res(j,1) += dphi(i)(1)*s(i)(j);
	  res(j,2) += dphi(i)(2)*s(i)(j);
	}
  }
  
  
  //! transformation Fi in the case of curved pyramid
  void PyramidReferenceSplit
  ::FjCurve(const SetPoints<Dimension3>& PTReel,const R3& pointloc,R3& res) const
  {
    res.Fill(0); VectReal_wp phi;
    ComputeValuesPhiNodalRef(pointloc, phi);
    for (int node = 0; node < points_nodal3d.GetM(); node++)
      Add(phi(node), PTReel.GetPointNodal(node), res);
  }
  
  
  //! transformation DFi in the case of curved pyramid
  void PyramidReferenceSplit
  ::DFjCurve(const SetPoints<Dimension3>& PTReel,const R3& pointloc,Matrix3_3& res) const
  {
    R3 grad; VectR3 grad_phi;
    res.Fill(0); ComputeGradientPhiNodalRef(pointloc, grad_phi);
    for (int node = 0; node < points_nodal3d.GetM(); node++)
      {
	grad = grad_phi(node);
	res(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	res(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
	res(2,0) += PTReel.GetPointNodal(node)(2)*grad(0);
	res(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	res(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
	res(2,1) += PTReel.GetPointNodal(node)(2)*grad(1);
	res(0,2) += PTReel.GetPointNodal(node)(0)*grad(2);
	res(1,2) += PTReel.GetPointNodal(node)(1)*grad(2);
	res(2,2) += PTReel.GetPointNodal(node)(2)*grad(2);
      }
  }
  
  
  //! returns true if the point is outside the unit pyramid
  /*!
    \param[in] Xn local coordinates of the point
    \param[in] epsilon threshold 
    \return true if the point is outside
  */
  inline bool PyramidReferenceSplit
  ::OutsideReferenceElement(const VectR3& s, const R3& Xn, const Real_wp& epsilon) const
  {
    if (type_basis == LIU)
      {
	if (Xn(2) < -epsilon)
	  return true;
	if ((Xn(2)+Xn(0))>(Real_wp(1)+epsilon)||(Xn(2)+Xn(1))>(Real_wp(1)+epsilon))
	  return true;
	if ((Xn(0) < -epsilon)||(Xn(1) < -epsilon))
	  return true;
      }	
    else
      {
	if (Xn(2)<0)
	  return true;
	if ((Xn(2)-Xn(0))>(Real_wp(1)+epsilon)||(Xn(2)-Xn(1))>(Real_wp(1)+epsilon))
	  return true;
	if ((Xn(2)+Xn(0))>(Real_wp(1)+epsilon)||(Xn(2)+Xn(1))>(Real_wp(1)+epsilon))
	  return true;
      }
    
    return false;
  }
  
  
  //! returns distance of the local point to the boundary of the element
  /*!
    \param[in] pointloc local coordinates of the point
    \return distance to the boundary (negative if the point is outside)
   */
  inline Real_wp PyramidReferenceSplit::GetDistanceToBoundary(const R3& pointloc) const
  {
    Real_wp distance = pointloc(2); 
    if (type_basis == LIU)
      {
	distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)+pointloc(0)));
	distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)-pointloc(0))); 
	distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)+pointloc(1)));
	distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)-pointloc(1)));
      }
    else
      {
	distance = min(distance, Real_wp(pointloc(0)));
	distance = min(distance, Real_wp(pointloc(1)));
	distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)-pointloc(0))); 
	distance = min(distance, Real_wp(Real_wp(1)-pointloc(2)-pointloc(1)));
      }
    return distance;
  }
  
  
  //! project the local point on the boundary, if outside the unit pyramid
  inline void PyramidReferenceSplit::ProjectPointOnBoundary(R3& pointloc) const
  {
    pointloc(2) = max( Real_wp(0), pointloc(2) );
    Real_wp one(1);
    if (type_basis == LIU)
      {
	if (pointloc(2)-pointloc(0)<one)
	  pointloc(0) = Real_wp(pointloc(2)-one);
	else if (pointloc(2)+pointloc(0)<one)
	  pointloc(0) = Real_wp(one-pointloc(2));
	if (pointloc(2)-pointloc(1)<one)
	  pointloc(1) = Real_wp(pointloc(2)-one);
	else if (pointloc(2)+pointloc(1)<one)
	  pointloc(1) = Real_wp(one-pointloc(2));
      }
    else
      {
	if (pointloc(2) + pointloc(0) > one)
	  pointloc(2) = Real_wp(one - pointloc(0));
	
	if (pointloc(2) + pointloc(1) > one)
	  pointloc(2) = Real_wp(one - pointloc(1));
	
	pointloc(0) = max(Real_wp(0), pointloc(0));
	pointloc(1) = max(Real_wp(0), pointloc(1));
      }
  }
  
  
  /**********************
   * FjElem and DFjElem *
   **********************/
  
  
  //! computes res = Fi(point) for all reference points
  /*!
    \param[in] s list of vertices of the element
    \param[out] res "reference points" after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void PyramidReferenceSplit::FjElem(const VectR3& s, SetPoints<Dimension3>& res,
				       const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad);
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureCurve(s, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      {
	FjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  FjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  FjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = Fi(point) for all nodal points
  inline void PyramidReferenceSplit::FjElemNodal(const VectR3& s, SetPoints<Dimension3>& res,
					    const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      FjElemNodalCurve(s, res, mesh, nquad);
    else
      FjElemNodalLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all quadrature points
  inline void PyramidReferenceSplit::FjElemQuadrature(const VectR3& s, SetPoints<Dimension3>& res,
						 const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  FjElemQuadratureCurve(s, res, mesh, nquad);
      }
    else
      FjElemQuadratureLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for all dof points
  inline void PyramidReferenceSplit::FjElemDof(const VectR3& s, SetPoints<Dimension3>& res,
					  const Mesh<Dimension3>& mesh, int nquad) const
  {
    if (mesh.Element(nquad).IsCurved())
      {
	FjElemNodalCurve(s, res, mesh, nquad);
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  FjElemDofCurve(s, res, mesh, nquad);
      }
    else
      FjElemDofLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all reference points
  /*!
    \param[in] s list of vertices of the element
    \param[in] PTReel "reference points" after transformation Fi
    \param[out] res jacobian matrices DFi
    \param[in] mesh given mesh
    \param[in] nquad element number
   */
  inline void PyramidReferenceSplit::DFjElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
				    SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh,
				    int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      {
	DFjElemNodalLinear(s, res);
	if (!this->quadrature_equal_nodal)
	  DFjElemQuadratureLinear(s, res);
	else
	  res.CopyNodalToQuadrature();
	
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else if (this->dof_equal_quadrature)
	  res.CopyQuadratureToDof();
	else
	  DFjElemDofLinear(s, res);
      }
  }
  
  
  //! computes res = DFi(point) for all nodal points
  inline void PyramidReferenceSplit
  ::DFjElemNodal(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                 SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      DFjElemNodalCurve(s,PTReel,res,mesh,nquad);
    else
      DFjElemNodalLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all quadrature points
  inline void PyramidReferenceSplit
  ::DFjElemQuadrature(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                      SetMatrices<Dimension3>& res,
                      const Mesh<Dimension3>& mesh,int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->quadrature_equal_nodal)
	  res.CopyNodalToQuadrature();
	else
	  DFjElemQuadratureCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemQuadratureLinear(s,res);
  }
  
  
  //! computes res = DFi(point) for all dof points
  inline void PyramidReferenceSplit
  ::DFjElemDof(const VectR3& s, const SetPoints<Dimension3>& PTReel,
               SetMatrices<Dimension3>& res,
               const Mesh<Dimension3>& mesh, int nquad) const
  {
    res.InitSetPoints(PTReel);
    if (mesh.Element(nquad).IsCurved())
      {
	DFjElemNodalCurve(s, PTReel, res, mesh, nquad);
	if (this->dof_equal_nodal)
	  res.CopyNodalToDof();
	else
	  DFjElemDofCurve(s, PTReel, res, mesh, nquad);
      }
    else
      DFjElemDofLinear(s, res);
  }
  
  
  //! computes res = Fi(point) for points located on the faces
  /*!
    \param[in] s list of vertices of the quad
    \param res references points after transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
    \param[in] ne local face number
   */
  inline void PyramidReferenceSplit::
  FjSurfaceElem(const VectR3& s, SetPoints<Dimension3>& res,
		const Mesh<Dimension3>& mesh, int nquad, int ne) const
  {
    int N;
    if (ne==0)
      N = points2d_quad.GetM();
    else
      N = points2d_tri.GetM();
     
    for (int k = 0; k < N; k++)
      res.SetPointQuadratureBoundary(k, res.GetPointQuadrature(this->num_quad_points_surf(ne)(k)));
  }
  
  
  //! computes res = DFi(point) for points located on the faces
  /*!
    \param[in] s list of vertices of the quad
    \param[in] PTReel references points after transformation Fi
    \param res jacobian matrices
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
    \param[in] ne local face number
    this method computes the normal and surfacic element (ds) on each quadrature point
   */
  inline void PyramidReferenceSplit
  ::DFjSurfaceElem(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                   SetMatrices<Dimension3>& res,
                   const Mesh<Dimension3>& mesh, int nquad, int ne) const
  {
    R3 normale_fj;
    Real_wp dsj(0); 
    Matrix3_3 mat_dfj,dfjm1;
    int N;
    if (ne==0)
      N = points2d_quad.GetM();
    else
      N = points2d_tri.GetM();
    
    for (int k = 0; k < N; k++)
      {
	// on récupère la matrice jacobienne sur le point de quadrature de la face
	res.SetPointQuadratureBoundary(k,res.GetPointQuadrature(this->num_quad_points_surf(ne)(k)));
	// on utilise la formule n = J_i / ds  DF_i^{-t} \hat{n} 
	mat_dfj = res.GetPointQuadrature(this->num_quad_points_surf(ne)(k));
	Real_wp deter = Det(mat_dfj);
	if (deter != Real_wp(0))
	  {
	    GetInverse(mat_dfj, dfjm1);
	    MltTrans(dfjm1, normale(ne), normale_fj);
	    Mlt(deter, normale_fj);
	    dsj = Norm2(normale_fj);
	    Mlt(Real_wp(Real_wp(1)/dsj), normale_fj);
	  }
	else
	  {
	    dsj = 0.0;
	    normale_fj.Zero();
	  }

	res.SetNormaleQuadratureBoundary(k,normale_fj);
	res.SetDsQuadratureBoundary(k,dsj);
      }
  }
  
  
  //! nodal points in the case of straight pyramid
  void PyramidReferenceSplit::FjElemNodalLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    for (int n = 0; n < points_nodal3d.GetM(); n++)
      {
        res_n.Fill(0);
        for (int i = 0; i < 5; i++)
          for (int j = 0; j < 3; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointNodal(n,res_n);
      }
  }
  
  
  //! quadrature points in the case of straight pyramid
  void PyramidReferenceSplit
  ::FjElemQuadratureLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {  
    R3 res_n;
    for (int n = points_nodal3d.GetM(); n < points_nodal3d.GetM() + points3d.GetM(); n++)
      {
        res_n.Fill(0);
        for (int i = 0; i < 5; i++)
          for (int j = 0; j< 3 ; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointQuadrature(n-points_nodal3d.GetM(),res_n);
      }
  }
  
  
  //! dof points in the case of straight pyramid
  void PyramidReferenceSplit::FjElemDofLinear(const VectR3& s, SetPoints<Dimension3>& res) const
  {
    R3 res_n;
    for (int n = points_nodal3d.GetM()+points3d.GetM();
         n < points_nodal3d.GetM()+points3d.GetM()+points_dof3d.GetM(); n++)
      {
        res_n.Fill(0);
	for (int i = 0; i < 5; i++)
          for (int j = 0; j < 3; j++)
            res_n(j) += coefFi(n,i)*s(i)(j);
	
	res.SetPointDof(n-points_nodal3d.GetM()-points3d.GetM(),res_n);
      }

  }
  
  
  //! transformation DFi for nodal points in the case of straight pyramid
  void PyramidReferenceSplit
  ::DFjElemNodalLinear(const VectR3& s, SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    for (int n = 0; n < points_nodal3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 5; i++)
         for (int j = 0; j < 3; j++)
          {
	    mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	    mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	    mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	  }
        res.SetPointNodal(n,mat);	  
      }
  }
  
  
  //! transformation DFi for quadrature points in the case of straight pyramid
  void PyramidReferenceSplit
  ::DFjElemQuadratureLinear(const VectR3& s,SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    for (int n = points_nodal3d.GetM(); n < points_nodal3d.GetM()+points3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 5; i++)
         for (int j = 0; j < 3; j++)
          {
	    mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	    mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	    mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	  }
        res.SetPointQuadrature(n-points_nodal3d.GetM(),mat);	  
      }
  }
  
  
  //! transformation DFi for dof points in the case of straight pyramid
  void PyramidReferenceSplit::DFjElemDofLinear(const VectR3& s,SetMatrices<Dimension3>& res) const
  {
    Matrix3_3 mat;
    for (int n = points_nodal3d.GetM()+points3d.GetM();
         n < points_nodal3d.GetM()+points3d.GetM()+points_dof3d.GetM(); n++)
      {
	mat.Fill(0);
	for (int i = 0; i < 5; i++)
         for (int j = 0; j < 3; j++)
          {
	    mat(j,0) += coefDFi_dx(n,i)*s(i)(j);
	    mat(j,1) += coefDFi_dy(n,i)*s(i)(j);
	    mat(j,2) += coefDFi_dz(n,i)*s(i)(j);
	  }
        res.SetPointDof(n-points_nodal3d.GetM()-points3d.GetM(),mat);	  
      }
  }
  
    
  //! transformation Fi for nodal points in the case of curved pyramid
  /*!
    \param[in] s four vertices of the element
    \param[out] res references points after the transformation Fi
    \param[in] mesh given mesh
    \param[in] nquad element number in the mesh
   */
  void PyramidReferenceSplit
  ::FjElemNodalCurve(const VectR3& s, SetPoints<Dimension3>& res,
                     const Mesh<Dimension3>& mesh,int nquad) const
  {
    VolumeReference::FjElemNodalCurve(s, res, mesh, nquad);
    
    R3 res_n;
    // for other nodal points, we use coefFi_curve
    for (int node = 3*order*order+2; node < points_nodal3d.GetM(); node++)
      {
	res_n.Zero();
	for (int k = 0; k < 3*order*order+2; k++)
	  Add(coefFi_curve(node, k), res.GetPointNodal(k), res_n);
	
	res.SetPointNodal(node, res_n);
      }
  }
  
  
  //! transformation Fi for quadrature points in the case of curved pyramid
  void PyramidReferenceSplit::FjElemQuadratureCurve(const VectR3& s, SetPoints<Dimension3>& res,
					   const Mesh<Dimension3>& mesh, int nquad) const
  {
    R3 res_n; int offset = points_nodal3d.GetM();
    for (int i = 0; i < points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < points_nodal3d.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }
  
  
  //! transformation Fi for dof points in the case of curved pyramid
  void PyramidReferenceSplit
  ::FjElemDofCurve(const VectR3& s, SetPoints<Dimension3>& res,
                   const Mesh<Dimension3>& mesh,int nquad) const
  {
    R3 res_n; int offset = points_nodal3d.GetM() + points3d.GetM();
    for (int i = 0; i < points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < points_nodal3d.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  //! transformation DFi for nodal points in the case of curved pyramid
  void PyramidReferenceSplit
  ::DFjElemNodalCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                      SetMatrices<Dimension3>& res,const Mesh<Dimension3>& mesh,int nquad) const
  {
    Matrix3_3 tmp;
    R3 grad; grad.Fill(0);
    for (int i = 0; i < points_nodal3d.GetM(); i++)
      {
	tmp.Fill(0);
	for (int node = 0; node < points_nodal3d.GetM(); node++)
	  {
	    grad(0) = coefDFi_dx_curve(i, node);
	    grad(1) = coefDFi_dy_curve(i, node);
	    grad(2) = coefDFi_dz_curve(i, node);
	    tmp(0,0) += PTReel.GetPointNodal(node)(0)*grad(0);
	    tmp(1,0) += PTReel.GetPointNodal(node)(1)*grad(0);
	    tmp(2,0) += PTReel.GetPointNodal(node)(2)*grad(0);
	    tmp(0,1) += PTReel.GetPointNodal(node)(0)*grad(1);
	    tmp(1,1) += PTReel.GetPointNodal(node)(1)*grad(1);
	    tmp(2,1) += PTReel.GetPointNodal(node)(2)*grad(1);
	    tmp(0,2) += PTReel.GetPointNodal(node)(0)*grad(2);
	    tmp(1,2) += PTReel.GetPointNodal(node)(1)*grad(2);
	    tmp(2,2) += PTReel.GetPointNodal(node)(2)*grad(2);
	  }
	res.SetPointNodal(i,tmp);
      } 
  }
  
  
  //! transformation DFi for quadrature points in the case of curved pyramid
  void PyramidReferenceSplit
  ::DFjElemQuadratureCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                           SetMatrices<Dimension3>& res,
                           const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n; int offset = points_nodal3d.GetM();
    for (int i = 0; i < points3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < points_nodal3d.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointQuadrature(i,res_n);
      }
  }
  
  
  //! transformation DFi for dof points in the case of curved pyramid
  void PyramidReferenceSplit
  ::DFjElemDofCurve(const VectR3& s, const SetPoints<Dimension3>& PTReel,
                    SetMatrices<Dimension3>& res, const Mesh<Dimension3>& mesh, int nquad) const
  {
    Matrix3_3 res_n; int offset = points_nodal3d.GetM() + points3d.GetM();
    for (int i = 0; i < points_dof3d.GetM(); i++)
      {
	res_n.Fill(0);
	for (int node = 0; node < points_nodal3d.GetM(); node++)
	  Add(coefFi_curve(offset+i, node), res.GetPointNodal(node), res_n);
	
	res.SetPointDof(i,res_n);
      }
  }
  
  
  /*******************
   * Other functions *
   *******************/
  
  
  //! returns \f$ \varphi_{node}(pointloc) \f$
  /*!
    where \f$ \varphi_{node} \f$ is the basis function associated
    with the nodal point \f$ \xi_{node} \f$
  */
  void PyramidReferenceSplit
  ::ComputeValuesPhiNodalRef(const R3& pointloc, VectReal_wp& res) const  
  {   
    res.Reallocate(points_nodal3d.GetM());
    if (type_basis == WALKER)
      {
	Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
	
	VectReal_wp val_monome(power_other_coef.GetM());
	for (int k = 0; k < power_other_coef.GetM(); k++)
	  {
	    if (x > y)
	      val_monome(k) = pow(x+z, power_other_coef(k)(0))*pow(x, power_other_coef(k)(1))
		*pow(y, power_other_coef(k)(2))*pow(z, power_other_coef(k)(3));
	    else
	      val_monome(k) = pow(y+z, power_other_coef(k)(0))
                *pow(x, power_other_coef(k)(1)+power_other_coef(k)(0))
		*pow(y, power_other_coef(k)(2)-power_other_coef(k)(0))
                *pow(z, power_other_coef(k)(3));
	  }
	
	for (int node = 0; node < res.GetM(); node++)
	  {
	    res(node) = LagrangePr(node).Evaluate(pointloc);
	    for (int k = 0; k < power_other_coef.GetM(); k++)
	      res(node) += val_monome(k)*other_coef_lagrange(node)(k);
	  }	
      }
    else if (type_basis == WIENERS)
      {
	Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
	if (order == 1)
	  {
	    if (x > y)
	      {
		res(0) = (1.0-x)*(1.0-y) + z*(y-1.0);
		res(1) = x*(1.0-y) - z*y;
		res(2) = x*y + z*y;
		res(3) = (1.0-x)*y - z*y;
	      }
	    else
	      {
		res(0) = (1.0-x)*(1.0-y) + z*(x-1.0);
		res(1) = x*(1.0-y) - z*x;
		res(2) = x*y + z*x;
		res(3) = (1.0-x)*y - z*x;
	      }
	    
	    res(4) = z;
	  }
	else
	  {
	    if (x > y)
	      {
		res(0) = ((1.0-x)*(1.0-y)+z*(y-1.0))*(1.0-2*x-2*y-2*z);
		res(1) = (x*(1.0-y)-z*y)*(2*x-2*y-1.0);
		res(2) = (x*y+z*y)*(2*x+2*y+2*z-3.0);
		res(3) = ((1.0-x)*y-z*y)*(2*y-2*x-1.0);
		res(5) = 4.0*x*((1.0-x)*(1.0-y) + z*(y-1.0)) - 2.0*y*z*(1.0-x-z);
		res(6) = 4.0*y*(x*(1.0-y) - z*y) + 2.0*y*z*(1.0-x-z);
		res(7) = 4.0*x*((1.0-x)*y - z*y) + 2.0*y*z*(1.0-x-z);
		res(8) = 4.0*y*((1.0-x)*(1.0-y) + z*(y-1.0)) - 2.0*y*z*(1.0-x-z);
		res(9) = 4.0*z*((1.0-x)*(1.0-y) + z*(y-1.0));
		res(10) = 4.0*z*(x*(1.0-y) - z*y);
		res(11) = 4.0*z*(x*y + z*y);
		res(12) = 4.0*z*((1.0-x)*y - z*y);
	      }
	    else
	      {
		res(0) = ((1.0-x)*(1.0-y)+z*(x-1.0))*(1.0-2*x-2*y-2*z);
		res(1) = (x*(1.0-y)-z*x)*(2*x-2*y-1.0);
		res(2) = (x*y+z*x)*(2*x+2*y+2*z-3.0);
		res(3) = ((1.0-x)*y-z*x)*(2*y-2*x-1.0);
		res(5) = 4.0*x*((1.0-x)*(1.0-y) + z*(x-1.0)) - 2.0*x*z*(1.0-y-z);
		res(6) = 4.0*y*(x*(1.0-y) - z*x) + 2.0*x*z*(1.0-y-z);
		res(7) = 4.0*x*((1.0-x)*y - z*x) + 2.0*x*z*(1.0-y-z);
		res(8) = 4.0*y*((1.0-x)*(1.0-y) + z*(x-1.0)) - 2.0*x*z*(1.0-y-z);
		res(9) = 4.0*z*((1.0-x)*(1.0-y) + z*(x-1.0));
		res(10) = 4.0*z*(x*(1.0-y) - z*x);
		res(11) = 4.0*z*(x*y + z*x);
		res(12) = 4.0*z*((1.0-x)*y - z*x);
	      }
	    
	    res(4) = z*(2*z-1.0);
	  }
      }
    else if (type_basis == LIU)
      {
	Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
	if (order == 1)
	  {
	    if (y < x)
	      {
		if (y < -x)
		  {
		    res(0) = 0.25*(x-1.0)*(y-z-1.0) - 0.5*z;
		    res(1) = 0.25*(x+1.0)*(-y+z+1.0) - 0.5*z;
		    res(2) = 0.25*(x+1.0)*(y-z+1.0);
		    res(3) = 0.25*(-x+1.0)*(y-z+1.0);
		    res(4) = z;
		  }
		else
		  {
		    res(0) = 0.25*(x-1.0)*(y-z-1.0) + 0.25*z*(x+y-2.0);
		    res(1) = 0.25*(x+1.0)*(-y+z+1.0) - 0.25*z*(x+y+2.0);
		    res(2) = 0.25*(x+1.0)*(y-z+1.0) + 0.25*z*(x+y); 
		    res(3) = 0.25*(-x+1.0)*(y-z+1.0) - 0.25*z*(x+y);
		    res(4) = z;
		  }
	      }
	    else
	      {
		if (y > -x)
		  {
		    res(0) = 0.25*(x-1.0)*(y+z-1.0);
		    res(1) = 0.25*(x+1.0)*(-y-z+1.0);
		    res(2) = 0.25*(x+1.0)*(y-z+1.0) + 0.5*x*z;
		    res(3) = 0.25*(-x+1.0)*(y-z+1.0) - 0.5*x*z;
		    res(4) = z;
		  }
		else
		  {
		    res(0) = 0.25*(x-1.0)*(y-z-1.0) + 0.25*z*(x-y-2.0);
		    res(1) = 0.25*(x+1.0)*(-y+z+1.0) - 0.25*z*(x-y+2.0);
		    res(2) = 0.25*(x+1.0)*(y-z+1.0) + 0.25*z*(x-y);
		    res(3) = 0.25*(-x+1.0)*(y-z+1.0) - 0.25*z*(x-y);
		    res(4) = z;
		  }
	      }
	  }	
      }
  }
  
  
  //! computation of \f$ \nabla \varphi_{node}(pointloc) \f$
  /*!
    where \f$ \varphi_{node} \f$ is the basis function
    associated with the nodal point \f$ \xi_{node} \f$
  */
  void PyramidReferenceSplit::ComputeGradientPhiNodalRef(const R3& pointloc, VectR3& dphi) const
  {
    dphi.Reallocate(points_nodal3d.GetM());
    if (type_basis == WALKER)
      {
	Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
	
	int N = power_other_coef.GetM();
	VectReal_wp val_dx(N), val_dy(N), val_dz(N);
	for (int k = 0; k < power_other_coef.GetM(); k++)
	  {
	    if (x > y)
	      {
		if (power_other_coef(k)(0) > 0)
		  {
		    val_dx(k) = power_other_coef(k)(0)*pow(x+z, power_other_coef(k)(0)-1)
                      *pow(x, power_other_coef(k)(1))
		      *pow(y, power_other_coef(k)(2))*pow(z, power_other_coef(k)(3));
		    
		    val_dz(k) = val_dx(k);
		  }
		
		if (power_other_coef(k)(1) > 0)
		  val_dx(k) += power_other_coef(k)(1)*pow(x+z, power_other_coef(k)(0))
                    *pow(x, power_other_coef(k)(1)-1)
		    *pow(y, power_other_coef(k)(2))*pow(z, power_other_coef(k)(3));
		
		if (power_other_coef(k)(2) > 0)
		  val_dy(k) += power_other_coef(k)(2)*pow(x+z, power_other_coef(k)(0))
                    *pow(x, power_other_coef(k)(1))
		    *pow(y, power_other_coef(k)(2)-1)*pow(z, power_other_coef(k)(3));
		
		if (power_other_coef(k)(3) > 0)
		  val_dz(k) += power_other_coef(k)(3)*pow(x+z, power_other_coef(k)(0))
                    *pow(x, power_other_coef(k)(1))
		    *pow(y, power_other_coef(k)(2))*pow(z, power_other_coef(k)(3)-1);
	      }
	    else
	      {
		if (power_other_coef(k)(0) > 0)
		  {
		    val_dy(k) = power_other_coef(k)(0)*pow(y+z, power_other_coef(k)(0)-1)
                      *pow(x, power_other_coef(k)(1)+power_other_coef(k)(0))
		      *pow(y, power_other_coef(k)(2)-power_other_coef(k)(0))
                      *pow(z, power_other_coef(k)(3));
		    
		    val_dz(k) = val_dy(k);
		  }
		
		if (power_other_coef(k)(1)+power_other_coef(k)(0) > 0)
		  val_dx(k) += (power_other_coef(k)(1)+power_other_coef(k)(0))
                    *pow(y+z, power_other_coef(k)(0))
		    *pow(x, power_other_coef(k)(1)+power_other_coef(k)(0)-1)
		    *pow(y, power_other_coef(k)(2)-power_other_coef(k)(0))
                    *pow(z, power_other_coef(k)(3));
		
		if (power_other_coef(k)(2)-power_other_coef(k)(0) > 0)
		  val_dy(k) += (power_other_coef(k)(2)-power_other_coef(k)(0))
                    *pow(y+z, power_other_coef(k)(0))
		    *pow(x, power_other_coef(k)(1)+power_other_coef(k)(0))
		    *pow(y, power_other_coef(k)(2)-power_other_coef(k)(0)-1)
                    *pow(z, power_other_coef(k)(3));
		
		if (power_other_coef(k)(3) > 0)
		  val_dz(k) += power_other_coef(k)(3)*pow(y+z, power_other_coef(k)(0))
		    *pow(x, power_other_coef(k)(1)+power_other_coef(k)(0))
		    *pow(y, power_other_coef(k)(2)-power_other_coef(k)(0))
                    *pow(z, power_other_coef(k)(3)-1);
	      }
	  }
	   
	for (int node = 0; node < dphi.GetM(); node++)
	  {
	    dphi(node)(0) = DxLagrangePr(node).Evaluate(pointloc);
	    dphi(node)(1) = DyLagrangePr(node).Evaluate(pointloc);
	    dphi(node)(2) = DzLagrangePr(node).Evaluate(pointloc);
	    
	    for (int k = 0; k < N; k++)
	      {
		dphi(node)(0) += val_dx(k)*other_coef_lagrange(node)(k);
		dphi(node)(1) += val_dy(k)*other_coef_lagrange(node)(k);
		dphi(node)(2) += val_dz(k)*other_coef_lagrange(node)(k);
	      }
	  }
      }
    else if (type_basis == WIENERS)
      {
	Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
	if (order == 1)
	  {
	    if (x > y)
	      {
		dphi(0).Init(-(1.0-y), -(1.0-x) + z, y-1.0);
		dphi(1).Init(1.0-y, -x - z, -y);
		dphi(2).Init(y, x+z, y);
		dphi(3).Init(-y, 1.0-x-z, -y);
	      }
	    else
	      {
		dphi(0).Init(-(1.0-y) +z, -(1.0-x), x-1.0);
		dphi(1).Init(1.0-y-z, -x, -x);
		dphi(2).Init(y+z, x, x);
		dphi(3).Init(-y-z, 1.0-x, -x);
	      }
		
	    dphi(4).Init(0, 0, 1.0);
	  }
	else
	  {
	    Real_wp v1, v2;
	    if (x > y)
	      {
		v1 = ((1.0-x)*(1.0-y)+z*(y-1.0));
		v2 = 1.0-2*x-2*y-2*z;
		dphi(0)(0) = -(1.0-y)*v2 -2.0*v1;
		dphi(0)(1) = (z+x-1.0)*v2 -2.0*v1;
		dphi(0)(2) = (y-1.0)*v2 -2.0*v1;

		v1 = x*(1.0-y)-z*y;
		v2 = 2*x-2*y-1.0;
		dphi(1)(0) = (1.0-y)*v2 + 2.0*v1;
		dphi(1)(1) = (-z-x)*v2 - 2.0*v1;
		dphi(1)(2) = -y*v2;
		
		v1 = x*y+z*y;
		v2 = 2*x+2*y+2*z-3.0;
		dphi(2)(0) = y*v2 + 2.0*v1;
		dphi(2)(1) = (x+z)*v2 + 2.0*v1;
		dphi(2)(2) = y*v2 + 2.0*v1;
		
		v1 = (1.0-x)*y-z*y;
		v2 = 2*y-2*x-1.0;
		dphi(3)(0) = -y*v2 - 2.0*v1;
		dphi(3)(1) = (1.0-x-z)*v2 + 2.0*v1;
		dphi(3)(2) = -y*v2;
		
		dphi(5)(0) = 4.0*((1.0-x)*(1.0-y) + z*(y-1.0)) - 4.0*x*(1.0-y) + 2.0*y*z;
		dphi(5)(1) = -4.0*x*(1.0-x-z) - 2.0*z*(1.0-x-z);
		dphi(5)(2) = 4.0*x*(y-1.0) - 2.0*y*(1.0-x-z) + 2.0*y*z;
		
		dphi(6)(0) = 4.0*y*(1.0-y) - 2.0*y*z;
		dphi(6)(1) = 4.0*(x*(1.0-y) - z*y) - 4.0*y*(x+z) + 2.0*z*(1.0-x-z);
		dphi(6)(2) = -4.0*y*y + 2.0*y*(1.0-x-z) - 2.0*y*z;
		
		dphi(7)(0) = 4.0*((1.0-x)*y - z*y) -4.0*x*y - 2.0*y*z;
		dphi(7)(1) = 4.0*x*(1.0-x-z) + 2.0*z*(1.0-x-z);
		dphi(7)(2) = -4.0*x*y + 2.0*y*(1.0-x-z) - 2.0*y*z;
		
		dphi(8)(0) = -4.0*y*(1.0-y) + 2.0*y*z;
		dphi(8)(1) = 4.0*((1.0-x)*(1.0-y) + z*(y-1.0)) + 4.0*y*(x-1.0+z) - 2.0*z*(1.0-x-z);
		dphi(8)(2) = 4.0*y*(y-1.0) - 2.0*y*(1.0-x-z) + 2.0*y*z;
		
		dphi(9)(0) = 4.0*z*(y-1.0);
		dphi(9)(1) = 4.0*z*(x-1.0+z);
		dphi(9)(2) = 4.0*((1.0-x)*(1.0-y) + z*(y-1.0)) + 4.0*z*(y-1.0);
		
		dphi(10)(0) = 4.0*z*(1.0-y);
		dphi(10)(1) = 4.0*z*(-x-z);
		dphi(10)(2) = 4.0*(x*(1.0-y) - z*y) - 4.0*z*y;
		
		dphi(11)(0) = 4.0*z*y;
		dphi(11)(1) = 4.0*z*(x+z);
		dphi(11)(2) = 4.0*(x*y + z*y) + 4.0*z*y;
		
		dphi(12)(0) = -4.0*z*y;
		dphi(12)(1) = 4.0*z*(1.0-x-z);
		dphi(12)(2) = 4.0*((1.0-x)*y - z*y) - 4.0*z*y;
	      }
	    else
	      {
		v1 = ((1.0-x)*(1.0-y)+z*(x-1.0));
		v2 = (1.0-2*x-2*y-2*z);
		dphi(0)(0) = (z+y-1.0)*v2 -2.0*v1;
		dphi(0)(1) = (x-1.0)*v2 -2.0*v1;
		dphi(0)(2) = (x-1.0)*v2 -2.0*v1;
		    
		v1 = x*(1.0-y)-z*x;
		v2 = 2*x-2*y-1.0;
		dphi(1)(0) = (1.0-y-z)*v2 + 2.0*v1;
		dphi(1)(1) = -x*v2 - 2.0*v1;
		dphi(1)(2) = -x*v2;

		v1 = x*y+z*x;
		v2 = 2*x+2*y+2*z-3.0;
		dphi(2)(0) = (y+z)*v2 + 2.0*v1;
		dphi(2)(1) = x*v2 + 2.0*v1;
		dphi(2)(2) = x*v2 + 2.0*v1;
		
		v1 = (1.0-x)*y-z*x;
		v2 = 2*y-2*x-1.0;
		dphi(3)(0) = -(y+z)*v2 - 2.0*v1;
		dphi(3)(1) = (1.0-x)*v2 + 2.0*v1;
		dphi(3)(2) = -x*v2;

		dphi(5)(0) = 4.0*((1.0-x)*(1.0-y) + z*(x-1.0)) + 4.0*x*(y+z-1.0) - 2.0*z*(1.0-y-z);
		dphi(5)(1) = 4.0*x*(x-1.0) + 2.0*x*z;
		dphi(5)(2) = 4.0*x*(x-1.0) - 2.0*x*(1.0-y-z) + 2.0*x*z;
		
		dphi(6)(0) = 4.0*y*(1.0-y-z) + 2.0*z*(1.0-y-z);
		dphi(6)(1) = 4.0*(x*(1.0-y) - z*x) - 4.0*y*x - 2.0*x*z;
		dphi(6)(2) = -4.0*y*x + 2.0*x*(1.0-y-z) - 2.0*x*z;
		
		dphi(7)(0) = 4.0*((1.0-x)*y - z*x) + 4.0*x*(-y-z) + 2.0*z*(1.0-y-z);
		dphi(7)(1) = 4.0*x*(1.0-x) - 2.0*x*z;
		dphi(7)(2) = -4.0*x*x + 2.0*x*(1.0-y-z) -2.0*x*z;
		
		dphi(8)(0) = -4.0*y*(1.0-y-z) -2.0*z*(1.0-y-z);
		dphi(8)(1) = 4.0*((1.0-x)*(1.0-y) + z*(x-1.0)) + 4.0*y*(x-1.0) + 2.0*x*z;
		dphi(8)(2) = 4.0*y*(x-1.0) - 2.0*x*(1.0-y-z) + 2.0*x*z;
		
		dphi(9)(0) = 4.0*z*(y-1.0+z);
		dphi(9)(1) = 4.0*z*(x-1.0);
		dphi(9)(2) = 4.0*((1.0-x)*(1.0-y) + z*(x-1.0)) + 4.0*z*(x-1.0);
		
		dphi(10)(0) = 4.0*z*(1.0-y-z);
		dphi(10)(1) = -4.0*z*x;
		dphi(10)(2) = 4.0*(x*(1.0-y) - z*x) - 4.0*z*x;
		
		dphi(11)(0) = 4.0*z*(y+z);
		dphi(11)(1) = 4.0*z*x;
		dphi(11)(2) = 4.0*(x*y + z*x) + 4.0*z*x;
		
		dphi(12)(0) = 4.0*z*(-y-z);
		dphi(12)(1) = 4.0*z*(1.0-x);
		dphi(12)(2) = 4.0*((1.0-x)*y - z*x) - 4.0*z*x;
	      }
	    
	    dphi(4).Init(0, 0, 4*z-1.0);
	    
	  }
      }
    else if (type_basis == LIU)
      {
	Real_wp x = pointloc(0), y = pointloc(1), z = pointloc(2);
	dphi(4).Init(0, 0, 1.0);
	if (y < x)
	  {
	    if (y < -x)
	      {
		dphi(0).Init(0.25*(y-z-1.0), 0.25*(x-1.0), -0.25*(x+1.0)); 
		dphi(1).Init(0.25*(-y+z+1.0), -0.25*(x+1.0), 0.25*(x-1.0)); 
		dphi(2).Init(0.25*(y-z+1.0), 0.25*(x+1.0), -0.25*(x+1.0)); 
		dphi(3).Init(-0.25*(y-z+1.0), 0.25*(-x+1.0), -0.25*(-x+1.0)); 
	      }
	    else
	      {
		dphi(0).Init(0.25*(y-1.0), 0.25*(x+z-1.0), 0.25*(y-1.0)); 
		dphi(1).Init(0.25*(-y+1.0), 0.25*(-x-1.0-z), 0.25*(-y-1.0)); 
		dphi(2).Init(0.25*(y+1.0), 0.25*(x+1.0+z), 0.25*(y-1.0)); 
		dphi(3).Init(0.25*(-y-1.0), 0.25*(-x+1.0-z), 0.25*(-1.0-y)); 
	      }
	  }
	else
	  {
	    if (y > -x)
	      {
		dphi(0).Init(0.25*(y+z-1.0), 0.25*(x-1.0), 0.25*(x-1.0)); 
		dphi(1).Init(0.25*(-y-z+1.0), -0.25*(x+1.0), -0.25*(x+1.0)); 
		dphi(2).Init(0.25*(y+z+1.0), 0.25*(x+1.0), 0.25*(x-1.0)); 
		dphi(3).Init(-0.25*(y+z+1.0), 0.25*(-x+1.0), 0.25*(-x-1.0)); 
	      }
	    else
	      {
		dphi(0).Init(0.25*(y-1.0), 0.25*(x-1.0-z), 0.25*(-1.0-y)); 
		dphi(1).Init(0.25*(-y+1.0), 0.25*(-x-1.0+z), 0.25*(y-1.0)); 
		dphi(2).Init(0.25*(y+1.0), 0.25*(x+1.0-z), 0.25*(-y-1.0)); 
		dphi(3).Init(0.25*(-y-1.0), 0.25*(-x+1.0+z), 0.25*(-1.0+y)); 
	      }
	  }
      }
  }    
  
    
  //! returns 3-D coordinates of a point on a face
  /*!
    \param[in] num_loc local face number in the quad
    \param[in] point_loc coordinate on the face
    \param[out] res local coordinates of the point in the unit square
  */
  inline void PyramidReferenceSplit
  ::GetLocalCoordOnBoundary(int num_loc, const R2& point_loc, R3& res) const
  {
    if (type_basis == LIU)
      {
	switch(num_loc)
	  {
	  case 0:
	    res.Init(2*point_loc(0)-1.0, 2*point_loc(1)-1.0, 0);
	    break;
	  case 1:
	    res.Init(2*point_loc(0)+point_loc(1)-Real_wp(1),point_loc(1)-Real_wp(1),point_loc(1));
	    break;
	  case 2:
	    res.Init(Real_wp(1)-point_loc(1),2*point_loc(0)+point_loc(1)-Real_wp(1),point_loc(1));
            break;
	  case 3:
	    res.Init(2*point_loc(0)+point_loc(1)-Real_wp(1),Real_wp(1)-point_loc(1),point_loc(1));
            break;
	  case 4:
	    res.Init(point_loc(1)-Real_wp(1),2*point_loc(0)+point_loc(1)-Real_wp(1),point_loc(1));
            break;
	  }
      }
    
    switch(num_loc)
      {
      case 0:
	res.Init(point_loc(0), point_loc(1), 0);break;
      case 1:
	res.Init(point_loc(0), 0, point_loc(1)); break;
      case 2:
	res.Init(1-point_loc(1), point_loc(0), point_loc(1)); break;
      case 3:
	res.Init(point_loc(0), 1-point_loc(1), point_loc(1)); break;
      case 4:
	res.Init(0, point_loc(0), point_loc(1)); break;
      }
  }
  
  
  void PyramidReferenceSplit::ComputeValuesPhiFirstOrder(const R3& point, VectReal_wp & phi) const
  {
    Real_wp x = point(0), y = point(1), z = point(2);
    phi.Reallocate(5);
    if (type_basis != LIU)
      {
	Real_wp teta(0);
	if (x > y)
	  teta = y*(x+z);
	else
	  teta = x*(y+z);
	
	phi(0) = 1.0 - x - y - z + teta;
	phi(1) = x - teta; phi(2) = teta; phi(3) = y-teta; phi(4) = z;
      }
    else
      {
	phi(4) = z;
	if (y < x)
	  {
	    if (y < -x)
	      {
		phi(0) = 0.25*(x-1.0)*(y-z-1.0) - 0.5*z;
		phi(1) = 0.25*(x+1.0)*(-y+z+1.0) - 0.5*z;
		phi(2) = 0.25*(x+1.0)*(y-z+1.0);
		phi(3) = 0.25*(-x+1.0)*(y-z+1.0);
	      }
	    else
	      {
		phi(0) = 0.25*(x-1.0)*(y-z-1.0) + 0.25*z*(x+y-2.0);
		phi(1) = 0.25*(x+1.0)*(-y+z+1.0) - 0.25*z*(x+y+2.0);
		phi(2) = 0.25*(x+1.0)*(y-z+1.0) + 0.25*z*(x+y);
		phi(3) = 0.25*(-x+1.0)*(y-z+1.0) - 0.25*z*(x+y);
	      }
	  }
	else
	  {
	    if (y > -x)
	      {
		phi(0) = 0.25*(x-1.0)*(y+z-1.0);
		phi(1) = 0.25*(x+1.0)*(-y-z+1.0);
		phi(2) = 0.25*(x+1.0)*(y-z+1.0) + 0.5*x*z;
		phi(3) = 0.25*(-x+1.0)*(y-z+1.0) - 0.5*x*z;
	      }
	    else
	      {
		phi(0) = 0.25*(x-1.0)*(y-z-1.0) + 0.25*z*(x-y-2.0);
		phi(1) = 0.25*(x+1.0)*(-y+z+1.0) - 0.25*z*(x-y+2.0);
		phi(2) = 0.25*(x+1.0)*(y-z+1.0) + 0.25*z*(x-y);
		phi(3) = 0.25*(-x+1.0)*(y-z+1.0) - 0.25*z*(x-y);
	      }
	  }	  
      }
  }
  
  
  void PyramidReferenceSplit::ComputeGradientPhiFirstOrder(const R3& point, VectR3& dphi) const
  {
    Real_wp x = point(0), y = point(1), z = point(2);
    dphi.Reallocate(5);
    if (type_basis != LIU)
      {
	Real_wp dteta_dx(0), dteta_dy(0), dteta_dz(0);
	if (x > y)
	  {
	    dteta_dx = y; dteta_dy = x+z; dteta_dz = y;
	  }
	else
	  {
	    dteta_dx = y+z; dteta_dy = x; dteta_dz = x;
	  }
	
	dphi(0)(0) = -1.0 + dteta_dx; dphi(0)(1) = -1.0 + dteta_dy; dphi(0)(2) = -1.0 + dteta_dz;
	dphi(1)(0) = 1.0 - dteta_dx; dphi(1)(1) = - dteta_dy; dphi(1)(2) = - dteta_dz;
	dphi(2)(0) = dteta_dx; dphi(2)(1) = dteta_dy; dphi(2)(2) = dteta_dz;
	dphi(3)(0) = - dteta_dx; dphi(3)(1) = 1.0 - dteta_dy; dphi(3)(2) = - dteta_dz;
	dphi(4)(0) = 0; dphi(4)(1) = 0; dphi(4)(2) = 1.0;
      }
    else
      {
	dphi(4)(0) = 0; dphi(4)(1) = 0; dphi(4)(2) = 1.0;
	if (y < x)
	  {
	    if (y < -x)
	      {
		dphi(0)(0) = 0.25*(y-z-1.0); dphi(0)(1) = 0.25*(x-1.0); dphi(0)(2) = -0.25*(x+1.0);
		dphi(1)(0) = 0.25*(-y+z+1.0); dphi(1)(1) = -0.25*(x+1.0);
                dphi(1)(2) = 0.25*(x-1.0);
		dphi(2)(0) = 0.25*(y-z+1.0); dphi(2)(1) = 0.25*(x+1.0); dphi(2)(2) = -0.25*(x+1.0);
		dphi(3)(0) = -0.25*(y-z+1.0); dphi(3)(1) = 0.25*(-x+1.0);
                dphi(3)(2) = -0.25*(-x+1.0);
	      }
	    else
	      {
		dphi(0)(0) = 0.25*(y-1.0); dphi(0)(1) = 0.25*(x+z-1.0); dphi(0)(2) = 0.25*(y-1.0);
		dphi(1)(0) = 0.25*(-y+1.0); dphi(1)(1) = 0.25*(-x-1.0-z);
                dphi(1)(2) = 0.25*(-y-1.0);
		dphi(2)(0) = 0.25*(y+1.0); dphi(2)(1) = 0.25*(x+1.0+z); dphi(2)(2) = 0.25*(y-1.0);
		dphi(3)(0) = 0.25*(-y-1.0); dphi(3)(1) = 0.25*(-x+1.0-z);
                dphi(3)(2) = 0.25*(-1.0-y);
	      }
	  }
	else
	  {
	    if (y > -x)
	      {
		dphi(0)(0) = 0.25*(y+z-1.0); dphi(0)(1) = 0.25*(x-1.0); dphi(0)(2) = 0.25*(x-1.0);
		dphi(1)(0) = 0.25*(-y-z+1.0); dphi(1)(1) = -0.25*(x+1.0);
                dphi(1)(2) = -0.25*(x+1.0);
		dphi(2)(0) = 0.25*(y+z+1.0); dphi(2)(1) = 0.25*(x+1.0); dphi(2)(2) = 0.25*(x-1.0);
		dphi(3)(0) = -0.25*(y+z+1.0); dphi(3)(1) = 0.25*(-x+1.0);
                dphi(3)(2) = 0.25*(-x-1.0);
	      }
	    else
	      {
		dphi(0)(0) = 0.25*(y-1.0); dphi(0)(1) = 0.25*(x-1.0-z); dphi(0)(2) = 0.25*(-1.0-y);
		dphi(1)(0) = 0.25*(-y+1.0); dphi(1)(1) = 0.25*(-x-1.0+z);
                dphi(1)(2) = 0.25*(y-1.0);
		dphi(2)(0) = 0.25*(y+1.0); dphi(2)(1) = 0.25*(x+1.0-z); dphi(2)(2) = 0.25*(-y-1.0);
		dphi(3)(0) = 0.25*(-y-1.0); dphi(3)(1) = 0.25*(-x+1.0+z);
                dphi(3)(2) = 0.25*(-1.0+y);
	      }
	  }
      }
  }
  
  
  //! displays informations class PyramidReferenceSplit
  ostream& operator <<(ostream& out, const PyramidReferenceSplit& e)
  {
    out<<static_cast<const VolumeReference&>(e);
    out<<"Type interpolation "<<e.type_interpolation<<endl;
    return out;
  }
  
} // end namespace

#define MONTJOIE_FILE_PYRAMID_REFERENCE_SPLIT_CXX
#endif
