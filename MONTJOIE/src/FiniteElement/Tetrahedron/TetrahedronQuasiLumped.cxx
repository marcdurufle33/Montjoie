#ifndef MONTJOIE_FILE_TETRAHEDRON_QUASI_LUMPED_CXX

namespace Montjoie
{
  //! default constructor
  TetrahedronQuasiLumped::TetrahedronQuasiLumped() : TetrahedronClassical()
  {
    type_interpolation = LOBATTO_ELEMENT;
  }
    
  // construction of finite element
  void TetrahedronQuasiLumped::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						      int rsurf_tri, int rsurf_quad,
						      int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    TetrahedronClassical::ConstructFiniteElement(r);
    
    // quadrature points on tetrahedron
    Globatto<Real_wp> lob;
    lob.ConstructQuadrature(order, lob.QUADRATURE_TETRAHEDRON_QUASILUMPED);
    points3d = lob.Points3D();
    weights3d = lob.Weights3D();
    points1d = lob.Points();
    weights1d = lob.Weights();
    
    nb_points_quadrature_inside = lob.GetNbPointsQuad();
    int nb_points_quadrature = nb_points_quadrature_inside + 4*points2d_tri.GetM();

    points3d.Resize(nb_points_quadrature);
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
    for (int i = 0; i < points2d_tri.GetM(); i++)
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
    
    // renumbering Points3D, Weights3D 
    IVect new_number(nb_points_quadrature_inside); new_number.Fill(-1);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_points_quadrature_inside; j++)
	if (points3d(j) == points_dof3d(i))
	  new_number(j) = i;
    
    int ind = nb_dof_loc;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      if (new_number(j) == -1)
	new_number(j) = ind++;
    
    VectR3 OldPoints = points3d; VectReal_wp OldWeights = weights3d;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	points3d(new_number(i)) = OldPoints(i);
	weights3d(new_number(i)) = OldWeights(i);
      }
    
    // on teste la formule de quadrature
    for (int p = 0; p <= 2*order; p++)
      {
	bool integration_ok = true;
	for (int i = 0; i <= p; i++)
	  for (int j = 0; j <= p-i; i++)
	    {
	      int k = p-i-j;
	      // on integre x^i y^j z^k
	      Real_wp val_exact = tgamma(i+1)*tgamma(j+1)*tgamma(k+1)/Real_wp(tgamma(p+4));
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
    
    this->ConstructMassMatrix();
    this->ConstructStiffnessMatrix();
    
    this->ComputeCoefficientTransformation();
  }


  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */    
  template<class Vector1, class Vector2>
  void TetrahedronQuasiLumped::ApplyCh(const Vector1& Vh, Vector2& Uh) const
  {
    for (int i = 0; i < nb_dof_loc; i++)
      {
	// identity part
	Uh(i) = Vh(i);
	
	// part with matrix P
	for (int j = nb_dof_loc; j < nb_points_quadrature_inside; j++)
	  Uh(i) += Value_Phi(i, j)*Vh(j);
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
  void TetrahedronQuasiLumped::ApplyChTranspose(const Vector1& Uh, Vector2& Vh) const
  {
    // identity part
    for (int i = 0; i < nb_dof_loc; i++)
      Vh(i) = Uh(i);
    
    // part with matrix P
    for (int i = nb_dof_loc; i < nb_points_quadrature_inside; i++)
      {
	Vh(i) = 0;
	for (int j = 0; j < nb_dof_loc; j++)
	  Vh(i) += Value_Phi(j, i)*Uh(j);
      }
  }
    
  
  //! displaying informations about TetrahedronQuasiLumped
  ostream& operator <<(ostream& out, const TetrahedronQuasiLumped& e)
  {
    return out;
  }

} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_QUASI_LUMPED_HXX
#endif
