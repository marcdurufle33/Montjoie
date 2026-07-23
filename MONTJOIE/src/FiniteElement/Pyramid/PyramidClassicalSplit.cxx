#ifndef MONTJOIE_FILE_PYRAMID_CLASSICAL_SPLIT_CXX

namespace Montjoie
{
  
  //! default constructor
  PyramidClassicalSplit::PyramidClassicalSplit() : PyramidReferenceSplit()
  {
    quadrature_equal_nodal = false;
    dof_equal_nodal = true;
    dof_equal_quadrature = false;
  }
  
  
  void PyramidClassicalSplit::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						     int rsurf_tri, int rsurf_quad,
						     int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    PyramidReferenceSplit::ConstructFiniteElement(r);
    ConstructFunctions();
    ConstructMassMatrix();
    ConstructStiffnessMatrix();
    ComputeCoefficientTransformation();
    
  }

  
  //! construction of basis functions
  void PyramidClassicalSplit::ConstructFunctions()
  {
    nb_dof_loc = points_nodal3d.GetM();
        
    nb_dof_quad = (order+1)*(order+1);
    nb_dof_tri = (order+1)*(order+2)/2;
    
    if (type_basis != WALKER)
      nb_dof_quad = 4*order;
    
    // DISP(Points3D); DISP(Weights3D); DISP(PointsDof3D);
  }
  
  
  //! construction of mass matrix
  void PyramidClassicalSplit::ConstructMassMatrix()
  {
    Value_Phi.Reallocate(nb_dof_loc, points3d.GetM());
    VectReal_wp phi;
    for (int i = 0; i < points3d.GetM(); i++)
      {
	ComputeValuesPhiRef(points3d(i), phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Value_Phi(j,i) = phi(j);
      }
    
    mass_matrix.Reallocate(nb_dof_loc, nb_dof_loc);
    Real_wp vloc;
    
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = i; j < nb_dof_loc; j++)
	{
	  vloc = 0.0;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    vloc += weights3d(k)*Value_Phi(i,k)*Value_Phi(j,k);
	  
	  mass_matrix(i,j) = vloc;
	}
    
    mass_matrix_chol = mass_matrix;
    GetCholesky(mass_matrix_chol);

  }
  
  
  //! construction of stiffness matrix
  void PyramidClassicalSplit::ConstructStiffnessMatrix()
  {
    Gradient_Phi.Reallocate(nb_dof_loc, points3d.GetM());
    VectR3 grad_phi;
    for (int i = 0; i < points3d.GetM(); i++)
      {
	ComputeGradientPhiRef(points3d(i), grad_phi);
	for (int j = 0; j < nb_dof_loc; j++)
	  Gradient_Phi(j, i) = grad_phi(j);
      }
     
    stiffness_matrix00.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix11.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix01.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix10.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix22.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix02.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix20.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix12.Reallocate(nb_dof_loc,nb_dof_loc);
    stiffness_matrix21.Reallocate(nb_dof_loc,nb_dof_loc);
    
    gradient_matrix0.Reallocate(nb_dof_loc,nb_dof_loc);
    gradient_matrix1.Reallocate(nb_dof_loc,nb_dof_loc);
    gradient_matrix2.Reallocate(nb_dof_loc,nb_dof_loc);
    Real_wp vloc1,vloc2,vloc3;
    
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  vloc1 = 0.0;
	  vloc2 = 0.0;
	  vloc3 = 0.0;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    {
	      vloc1 += weights3d(k)*Gradient_Phi(i,k)(0)*Gradient_Phi(j,k)(0);
	      vloc2 += weights3d(k)*Gradient_Phi(i,k)(1)*Gradient_Phi(j,k)(1);
	      vloc3 += weights3d(k)*Gradient_Phi(i,k)(0)*Gradient_Phi(j,k)(1);
	    }
	  stiffness_matrix00(i,j) = vloc1;
	  stiffness_matrix11(i,j) = vloc2;
	  stiffness_matrix01(i,j) = vloc3;
	  
	  vloc1 = 0.0;
	  vloc2 = 0.0;
	  vloc3 = 0.0;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    {
	      vloc1 += weights3d(k)*Gradient_Phi(i,k)(2)*Gradient_Phi(j,k)(2);
	      vloc2 += weights3d(k)*Gradient_Phi(i,k)(2)*Gradient_Phi(j,k)(0);
	      vloc3 += weights3d(k)*Gradient_Phi(i,k)(1)*Gradient_Phi(j,k)(2);
	    }
	  stiffness_matrix22(i,j) = vloc1;
	  stiffness_matrix20(i,j) = vloc2;
	  stiffness_matrix12(i,j) = vloc3;
	  
	  vloc1 = 0.0;
	  vloc2 = 0.0;
	  vloc3 = 0.0;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    {
	      vloc1 += weights3d(k)*Gradient_Phi(i,k)(1)*Gradient_Phi(j,k)(0);
	      vloc2 += weights3d(k)*Gradient_Phi(i,k)(0)*Gradient_Phi(j,k)(2);
	      vloc3 += weights3d(k)*Gradient_Phi(i,k)(2)*Gradient_Phi(j,k)(1);
	    }
	  stiffness_matrix10(i,j) = vloc1;
	  stiffness_matrix02(i,j) = vloc2;
	  stiffness_matrix21(i,j) = vloc3;
	  
	  vloc1 = 0.0;
	  vloc2 = 0.0;
	  vloc3 = 0.0;
	  for (int k = 0; k < nb_points_quadrature_inside; k++)
	    {
	      vloc1 += weights3d(k)*Gradient_Phi(i,k)(0)*Value_Phi(j,k);
	      vloc2 += weights3d(k)*Gradient_Phi(i,k)(1)*Value_Phi(j,k);
	      vloc3 += weights3d(k)*Gradient_Phi(i,k)(2)*Value_Phi(j,k);
	    }
	  gradient_matrix0(i,j) = vloc1;
	  gradient_matrix1(i,j) = vloc2;
	  gradient_matrix2(i,j) = vloc3;
	}
    // DISP(Gradient_Phi); DISP(Value_Phi);
  }
  
  
  //! computation of elementary interaction 
  /*!
    \param[in] i row number
    \param[in] j column number
    \param[in] mass mass coefficient
    \param[in] C stiffness tensor C
    \param[in] D stiffness tensor D
    \param[in] E stiffness tensor E
    \param[out] vloc result
    vloc = \f$ \int_{\hat{T}} mass \varphi_i \varphi_j + C*\nabla \varphi_i \nabla \varphi_j
    + D*\nabla \varphi_j \varphi_i + E \varphi_j \nabla \varphi_i \f$
  */
  template<class T, class Matrix1, class Vector1>
  void PyramidClassicalSplit
  ::ComputeElementaryInteraction(int i, int j, const T& mass, const Matrix1& C,
                                 const Vector1& D, const Vector1& E, T& vloc) const
  {
    // stiffness part
    vloc = stiffness_matrix00(i,j)*C(0,0) + stiffness_matrix11(i,j)*C(1,1)
      + stiffness_matrix22(i,j) * C(2,2) + stiffness_matrix01(i,j)*C(0,1)
      + stiffness_matrix02(i,j) * C(0,2) + stiffness_matrix12(i,j)*C(1,2) +
      stiffness_matrix10(i,j)*C(1,0) + stiffness_matrix20(i,j) * C(2,0)
      + stiffness_matrix21(i,j)*C(2,1);
    
    // gradient part
    vloc += gradient_matrix0(j,i)*D(0) + gradient_matrix1(j,i)*D(1)
      + gradient_matrix2(j,i)*D(2) + gradient_matrix0(i,j)*E(0)
      + gradient_matrix1(i,j)*E(1) + gradient_matrix2(i,j)*E(2);
    
    // mass part
    vloc += mass_matrix(i,j)*mass;
  }
  
  
  void PyramidClassicalSplit::ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const
  {
    ComputeValuesPhiNodalRef(x, phi);
  }
  
  
  void PyramidClassicalSplit::ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const
  {
    ComputeGradientPhiNodalRef(x, grad_phi);
  }
  
  
  //! returns \f$ \varphi_i(\xi_j) \f$
  /*! 
    num_dof = i : degree of freedom number
    num_points = j : quadrature point number
  */
  Real_wp PyramidClassicalSplit::GetValuePhi(int num_dof, int num_point) const
  {
    return Value_Phi(num_dof, num_point);
  }
  
  
  //! returns \f$ \varphi^{2D}_k(\xi_j) \f$
  Real_wp PyramidClassicalSplit::GetValuePhiOnBoundary(int k, int j, int num_loc) const
  {
    return Value_Phi(k, this->num_quad_points_surf(num_loc)(j));
  }

  
  //! computes \f$ \nabla \varphi_{numdof} (\xi_{numpoint}) \f$
  const R3& PyramidClassicalSplit::
  GetGradientPhi(int num_dof, int num_point) const
  {
    return Gradient_Phi(num_dof, num_point);
  }
    
  
  //! returns \f$ \varphi_k(\xi_j) \f$
  /*!
    where \f$ \varphi \mbox{ and } \xi \f$
    are restriction of functions on boundary and quadrature point on the boundary
  */
  const R3& PyramidClassicalSplit::
    GetGradientPhiOnBoundary(int node, int num_point, int num_loc) const
  {
    return GetGradientPhi(node, this->num_quad_points_surf(num_loc)(num_point));
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \hat{\varphi}_i \f$
  /*!
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral agains each basis function
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1,class Vector2>
  void PyramidClassicalSplit::ComputeIntegralRef(const Vector1 & feval, Vector2& res) const
  {
    res.Reallocate(nb_dof_loc); res.Fill(0);
    for (int j = 0; j < nb_dof_loc; j++)
      {
	for (int i = 0; i < nb_points_quadrature_inside; i++)
	  {
	    res(j) += feval(i)*GetValuePhi(j,i);
	  }
      }
  }
  
  
  //! computation of \f$ \int_{\hat{K}} f \hat{\nabla} \hat{\varphi} \f$
  /*!
    \param[in] feval values of the function f on quadrature points
    \param[out] res integral against gradient of each basis function
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1,class Vector2>
  void PyramidClassicalSplit::ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const
  {
    res.Fill(0);
    for (int j = 0; j < nb_dof_loc; j++)
      for (int i = 0; i < nb_points_quadrature_inside; i++)
	res(j) += DotProd(Gradient_Phi(j,i), feval(i));
  }
  
  
  //! computation of \f$ \int_\Sigma f \phi_i ds \f$
  /*!
    \param[in] feval values of the function f on quadrature points of the face
    \param[out] res integral against each basis function
    \param[in] num_loc local face number
    \warning feval(i) is assumed to be already multiplied by the integration weight
   */
  template<class Vector1,class Vector2>
  void PyramidClassicalSplit
  ::ComputeIntegralSurfaceRef(const Vector1 & feval, Vector2& res, int num_loc) const
  {
    res.Fill(0);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < feval.GetM(); j++)
	res(i) += feval(j)*Value_Phi(i, this->num_quad_points_surf(num_loc)(j));
  }

  
  //! not impemented
  template<class Vector1,class Vector2>
  void PyramidClassicalSplit
  ::ComputeIntegralSurfaceGradientRef(const Vector1 & feval, Vector2& res, int num_loc) const
  {
  }
  
  
  //! displays details about class PyramidClassicalSplit
  ostream& operator <<(ostream& out, const PyramidClassicalSplit& e)
  {
    out<<static_cast<const PyramidReferenceSplit&>(e);
    return out;
  }
  
} // end namespace
  
#define MONTJOIE_FILE_PYRAMID_CLASSICAL_SPLIT_CXX
#endif
