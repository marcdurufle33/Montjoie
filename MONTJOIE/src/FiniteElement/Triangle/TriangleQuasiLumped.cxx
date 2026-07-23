#ifndef MONTJOIE_FILE_TRIANGLE_QUASI_LUMPED_CXX

namespace Montjoie
{
  TriangleQuasiLumped::TriangleQuasiLumped() : TriangleClassical()
  {
  }
    
  // construction of finite element
  void TriangleQuasiLumped::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						   int rsurf, int type_surf)
  {
    /*order = r;
    type_interpolation = this->LOBATTO_ELEMENT;
    
    // quadrature points on triangles
    Globatto<Dimension2> lob;
    lob.ConstructQuadrature(order, lob.QUADRATURE_TRIANGLE_QUASILUMPED);
    points2d = lob.Points2D();
    weights2d = lob.Weights2D();
    points1d = lob.Points();
    weights1d = lob.Weights();
            
    int ind = 3;
    this->num_quad_points_surf.Reallocate(3);
    // first edge
    this->num_quad_points_surf(0).Reallocate(order+1);
    this->num_quad_points_surf(0)(0) = 0;
    this->num_quad_points_surf(0)(order) = 1;
    for (int i = 1; i < order; i++)
      this->num_quad_points_surf(0)(i) = ind++;
    
    // second edge
    this->num_quad_points_surf(1).Reallocate(order+1);
    this->num_quad_points_surf(1)(0) = 1;
    this->num_quad_points_surf(1)(order) = 2;
    for (int i = 1; i < order; i++)
      this->num_quad_points_surf(1)(i) = ind++;
    
    // third edge
    this->num_quad_points_surf(2).Reallocate(order+1);
    this->num_quad_points_surf(2)(0) = 2;
    this->num_quad_points_surf(2)(order) = 0;
    for (int i = 1; i < order; i++)
      this->num_quad_points_surf(2)(i) = ind++;
    
    nb_points_quadrature_inside = points2d.GetM();
    nb_points_quadrature_boundaries = 4*(order+1);

    this->ConstructNodalPoints(r);
    this->ConstructFunctions();
    
    // renumbering Points2D, Weights2D 
    IVect new_number(nb_points_quadrature_inside); new_number.Fill(-1);
    for (int i = 0; i < nb_dof_loc; i++)
      for (int j = 0; j < nb_points_quadrature_inside; j++)
	if (points2d(j) == points_dof2d(i))
	  new_number(j) = i;
    
    ind = nb_dof_loc;
    for (int j = 0; j < nb_points_quadrature_inside; j++)
      if (new_number(j) == -1)
	new_number(j) = ind++;
    
    VectR2 OldPoints = points2d; VectReal_wp OldWeights = weights2d;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	points2d(new_number(i)) = OldPoints(i);
	weights2d(new_number(i)) = OldWeights(i);
      }
    
    DISP(points2d); DISP(weights2d);
    // on teste la formule de quadrature
    for (int k = 0; k <= 2*order; k++)
      {
	bool integration_ok = true;
	for (int i = 0; i <= k; i++)
	  {
	    int j = k-i;
	    // on integre x^i y^j
	    Real_wp val_exact = tgamma(i+1)*tgamma(j+1)/tgamma(k+3);
	    Real_wp val_app = 0;
	    for (int m = 0; m < nb_points_quadrature_inside; m++)
	      val_app += weights2d(m)*pow(points2d(m)(0), i)*pow(points2d(m)(1), j);
	    
	    if (abs(val_app-val_exact) > 1e5*epsilon_machine)
	      {
		DISP(val_app); DISP(val_exact); DISP(i); DISP(j);
		integration_ok = false;
	      }
	    }
	
	if (!integration_ok)
	  {
	    cout<<"L'ordre "<<k<<" n'est pas integre exactement"<<endl;
	    break;
	  }
      }

    this->ConstructMassMatrix();
    this->ConstructStiffnessMatrix();
    
    this->ComputeCoefficientTransformation();

    EdgesDof.Reallocate(order+1, 3);    
    for (int i = 0; i <= order; i++)
    for (int num_loc = 0; num_loc < 3; num_loc++)
    EdgesDof(i, num_loc) = this->num_quad_points_surf(num_loc)(i);
    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);    
    DISP(Value_Phi);
    */
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
  void TriangleQuasiLumped::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
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
  void TriangleQuasiLumped::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
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
  
  
  void TriangleQuasiLumped::ApplyCh(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }

  
  void TriangleQuasiLumped::ApplyCh(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChGen(Vh, Uh);
  }


  void TriangleQuasiLumped::ApplyChTranspose(const VectReal_wp& Vh, VectReal_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }

  
  void TriangleQuasiLumped::ApplyChTranspose(const VectComplex_wp& Vh, VectComplex_wp& Uh) const
  {
    ApplyChTransposeGen(Vh, Uh);
  }
    
  
  ostream& operator <<(ostream& out, const TriangleQuasiLumped& e)
  {
    return out;
  }

} // namespace Montjoie

#define MONTJOIE_FILE_TRIANGLE_QUASI_LUMPED_CXX
#endif
