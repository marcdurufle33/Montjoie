#ifndef MONTJOIE_FILE_HEXAHEDRON_DG_JACOBI_CXX

namespace Montjoie
{
  
  //! default constructor
  HexahedronDgJacobi::HexahedronDgJacobi() : HexahedronDgGauss()
  {
    quadrature_equal_nodal = false;
    dof_equal_nodal = false;
    dof_equal_quadrature = true;
    
    mass_lumping = true;
  }
  
  
  //! construction of basis functions
  void HexahedronDgJacobi::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						  int rsurf_tri, int rsurf_quad,
						  int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    HexahedronGauss::ConstructFiniteElement(r);
    
    // we choose points for integration and dofs of the base
    base1D.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    base1D.ComputeGradPhi(1e3*epsilon_machine);
    Weights1D_B = base1D.Weights(); Points1D_B = base1D.Points();
    
    // we choose points for integration and dofs in the direction z
    z1D.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_RADAU_JACOBI);
    z1D.ComputeGradPhi(1e3*epsilon_machine);
    Weights1D_Z = z1D.Weights(); Points1D_Z = z1D.Points();

    int nb_points_quadrature = (order+1)*(order+1)*(order+1);
    nb_points_quadrature_inside = nb_points_quadrature;
    weights3d.Reallocate(nb_points_quadrature); points3d.Reallocate(nb_points_quadrature);
    for (int i3 = 0; i3 <= order; i3++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i1 = 0; i1 <= order; i1++)
	  {
	    weights3d(NumNodes3D(i1,i2,i3))
              = Weights1D_B(i1)*Weights1D_B(i2)*Weights1D_Z(i3)/((1-Points1D_Z(i3)));
	    
            points3d(NumNodes3D(i1,i2,i3))(0) = Points1D_B(i1);
	    points3d(NumNodes3D(i1,i2,i3))(1) = Points1D_B(i2);
	    points3d(NumNodes3D(i1,i2,i3))(2) = Points1D_Z(i3);
	  }

    // we choose points of integration for the faces
    face1D.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_GAUSS);
    face1D.ComputeGradPhi(1e3*epsilon_machine);
    Weights1D_F = face1D.Weights(); Points1D_F = face1D.Points();

    weights2d_quad.Reallocate((order+1)*(order+1));
    points2d_quad.Reallocate((order+1)*(order+1));
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  weights2d_quad(NumNodes2D_quad(i, j)) = Weights1D_F(i)*Weights1D_F(j);
	  points2d_quad(NumNodes2D_quad(i, j)).Init(Points1D_F(i), Points1D_F(j));
	}
	  
    points_dof3d = points3d;
    
    points3d.Resize(nb_points_quadrature+6*points2d_quad.GetM());
    int nb = nb_points_quadrature;
    // Points of integration for the first face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(nb++).Init(0.0, points2d_quad(i)(0), points2d_quad(i)(1));
    
    // Points of integration for the second face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(nb++).Init(points2d_quad(i)(0), 0.0, points2d_quad(i)(1));
    
    // Points of integration for the third face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(nb++).Init(points2d_quad(i)(0), points2d_quad(i)(1), 0.0);
    
    // Points of integration for the fourth face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(nb++).Init(points2d_quad(i)(0), points2d_quad(i)(1), 1.0);
    
    // Points of integration for the fifth face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(nb++).Init(points2d_quad(i)(0), 1.0, points2d_quad(i)(1));
    
    // Points of integration for the sixth face
    for (int i = 0; i < points2d_quad.GetM(); i++)
      points3d(nb++).Init(1.0, points2d_quad(i)(0), points2d_quad(i)(1));
    
    nb = nb_points_quadrature;
    this->num_quad_points_surf.Reallocate(6);
    for (int num_loc = 0; num_loc < 6; num_loc++)
      {
	this->num_quad_points_surf(num_loc).Reallocate(points2d_quad.GetM());
	for (int i = 0; i < points2d_quad.GetM(); i++)
	  this->num_quad_points_surf(num_loc)(i) = nb++;
      }
    
    nb_points_quadrature += 6*points2d_quad.GetM();
    
    ValGauss_Extremity.Reallocate(order+1,2);
    for (int i = 0; i <= order; i++)
      {
	ValGauss_Extremity(i, 0) = base1D.EvaluatePhi(i, Real_wp(0));
	ValGauss_Extremity(i, 1) = base1D.EvaluatePhi(i, Real_wp(1));
      }
    
    ValJacobi_Extremity.Reallocate(order+1,2);
    for (int i = 0; i <= order; i++)
      {
	ValJacobi_Extremity(i, 0) = z1D.EvaluatePhi(i, Real_wp(0));
	ValJacobi_Extremity(i, 1) = z1D.EvaluatePhi(i, Real_wp(1));
      }

    // we construct I_GL
    lobatto1D.ConstructQuadrature(order, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    Points1D_L = lobatto1D.Points();
    J_GL.Reallocate(order+1,order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
        J_GL(i,j) = z1D.EvaluatePhi(i,Points1D_L(j));
    
    
    this->ComputeCoefficientTransformation();
    Value_Phi.Reallocate(nb_dof_loc, points3d.GetM()); VectReal_wp phi;
    Gradient_Phi.Reallocate(nb_dof_loc, points3d.GetM()); VectR3 dphi;
    for (int i = 0; i < points3d.GetM(); i++)
      {
	ComputeValuesPhiRef(points3d(i), phi);
	ComputeGradientPhiRef(points3d(i), dphi);
	for (int j = 0; j < nb_dof_loc; j++)
	  {
	    Value_Phi(j, i) = phi(j);
	    Gradient_Phi(j, i) = dphi(j);
	  }
      }
    
    this->ConstructStiffnessMatrixDG(*this);
  }
  
  
    //! returns \f$ \varphi_{num\_dof}(point\_loc) \f$
  /*!
    \param[in] num_dof dof number
    \param[in] point_loc coordinates of point in unit cube
   */
  void HexahedronDgJacobi::ComputeValuesPhiRef(const R3& point_loc, VectReal_wp& phi) const
  {
    VectReal_wp phix(order+1), phiy(order+1), phiz(order+1);
    for (int i = 0; i <= order; i++)
      {
	phix(i) = base1D.EvaluatePhi(i, point_loc(0));
	phiy(i) = base1D.EvaluatePhi(i, point_loc(1));
	phiz(i) = z1D.EvaluatePhi(i, point_loc(2));
      }
    
    phi.Reallocate(nb_dof_loc);
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  {
	    int node = NumNodes3D(i1, i2, i3);
	    phi(node) = phix(i1)*phiy(i2)*phiz(i3);
	  }
  }

  
  //! returns \f$ \varphi_{num\_dof}(point\_loc) \f$
  /*!
    \param[in] num_dof dof number
    \param[in] point_loc coordinates of point in unit cube
   */
  void HexahedronDgGauss::ComputeGradientPhiRef(const R3& point_loc, VectR3& dphi) const
  {
    VectReal_wp phix(order+1), phiy(order+1), phiz(order+1);
    VectReal_wp dphix(order+1), dphiy(order+1), dphiz(order+1);
    for (int i = 0; i <= order; i++)
      {
	phix(i) = base1D.EvaluatePhi(i, point_loc(0));
	phiy(i) = base1D.EvaluatePhi(i, point_loc(1));
	phiz(i) = z1D.EvaluatePhi(i, point_loc(2));
	
	dphix(i) = base1D.EvaluatePhiGrad(i, point_loc(0));
	dphiy(i) = base1D.EvaluatePhiGrad(i, point_loc(1));
	dphiz(i) = z1D.EvaluatePhiGrad(i, point_loc(2));
      }
    
    dphi.Reallocate(nb_dof_loc);
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  {
	    int node = NumNodes3D(i1, i2, i3);
	    dphi(node).Init(dphix(i1)*phiy(i2)*phiz(i3),
                            phix(i1)*dphiy(i2)*phiz(i3), phix(i1)*phiy(i2)*dphiz(i3));
	  }
  }
  
   
  //! computation of U on nodal points
  /*!
    \param[in] MatricesElem jacobian matrices
    \param[in] Uloc components of U on dofs
    \param[out] Uloc_node values of U on nodal points
    \param[in] mesh considered mesh
    \param[in] i element number in the mesh
   */
  template<class Vector1, class Vector2>
  void HexahedronDgJacobi::
  ComputeNodalValues(const SetMatrices<Dimension>& MatricesElem,
		     const Vector1& Uloc, Vector2& Uloc_node,
		     const Mesh<Dimension3>& mesh, int i) const
  {
    Vector2 U1(nb_dof_loc), U2(nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumNodes3D(i, j, k);
	    U1(node).Zero();
	    for (int m = 0; m <= order; m++)
	      Add(this->G_GL(m, k), Uloc(NumNodes3D(i, j, m)), U1(node));
	  }

    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumNodes3D(i, j, k);
	    U2(node).Zero();
	    for (int m = 0; m <= order; m++)
	      Add(this->G_GL(m, j), U1(NumNodes3D(i, m, k)), U2(node));
	  }
    
    Uloc_node.Reallocate(nb_dof_loc);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    int node = NumNodes3D(i, j, k);
	    Uloc_node(node).Zero();
	    for (int m = 0; m <= order; m++)
	      Add(J_GL(m, i), U2(NumNodes3D(m, j, k)), Uloc_node(node));
	  }
  }
  
  
  //! computation of u on quadrature points of a face
  /*!
    \param[in] feval dof components of u
    \param[out] res values of u on quadrature points of the face
    \param[in] num_loc local position of the face in the element
  */
  template<class Vector1, class Vector2>
  void HexahedronDgJacobi::
  ComputeValueBoundaryRef(const Vector1& feval, Vector2& res, int num_loc) const
  {
    // to do
    abort();
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \hat{\nabla} \hat{\varphi} \f$
  /*!
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral agains gradient of each basis function
    \warning feval(i) is assumed to be already multiplied by the weight of integration 
   */
  template<class Vector1,class Vector2>
  void HexahedronDgJacobi::
  ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(0);
    for (int i1 = 0; i1 <= order; i1++)
      for (int i2 = 0; i2 <= order; i2++)
	for (int i3 = 0; i3 <= order; i3++)
	  for (int k = 0; k <= order; k++)
	    res(NumNodes3D(i1,i2,i3)) += base1D.GradPhi(i1, k)*feval(NumNodes3D(k,i2,i3))(0)
	      + base1D.GradPhi(i2, k)*feval(NumNodes3D(i1,k,i3))(1)
	      + z1D.GradPhi(i3,k)*feval(NumNodes3D(i1,i2,k))(2);
    
  }
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains each basis function
    \param[in] num_loc local edge number
    \warning feval(i) is assumed to be already multiplied by the weight of integration 
   */
  template<class Vector1, class Vector2>
  void HexahedronDgJacobi::
  ComputeIntegralSurfaceRef(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    // to do
    abort();
  }
  
  
  //! computation of \f$ \int_\Sigma f \hat{\varphi}_i' ds \f$
  /*!
    \param[in] feval values of the function f on quadrature points of the edge
    \param[out] res integral agains surfacic derivative of each basis function
    \param[in] num_loc local edge number
    \warning it is supposed that feval(i) is already multiplied by the weight of integration 
   */
  template<class Vector1,class Vector2>
  void HexahedronDgJacobi::
  ComputeIntegralSurfaceGradientRef(const Vector1 & feval,
				    Vector2& res, int num_loc) const
  {
    res.Fill(0.0);
  }
  
} // namespace Montjoie
  
#define MONTJOIE_FILE_HEXAHEDRON_DG_JACOBI_CXX
#endif
