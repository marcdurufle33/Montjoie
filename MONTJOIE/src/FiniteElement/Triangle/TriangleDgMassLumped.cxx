#ifndef MONTJOIE_FILE_TRIANGLE_DG_MASSLUMPED_CXX

namespace Montjoie
{
  
  //! default constructor
  TriangleMassLumped::TriangleMassLumped() : TriangleClassical()
  {    
    this->Fb_geom.quadrature_equal_nodal = false;
    this->Fb_geom.dof_equal_nodal = false;
    this->Fb_geom.dof_equal_quadrature = false;
  }
  
  
  //! how to number the mesh
  void TriangleMassLumped::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    TriangleClassical::ConstructNumberMap(nmap, dg);
    
    if (dg == ElementReference_Base::CONTINUOUS)
      nmap.SetNbDofTriangle(order, nb_dof_loc-3*order);
  }
  
  
  //! constructing finite element
  void TriangleMassLumped::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
						  int rsurf, int type_surf)
  {
    type_interpolation = TriangleGeomReference::MASS_LUMPED_BASIS;
    TriangleReference<1>::
      ConstructFiniteElement(r, rgeom, r, TriangleQuadrature::QUADRATURE_MASS_LUMPED); 
    
    ConstructFunctions();
    
    ComputeLagrangianBubbleFunctions(r);
    
    ConstructElementaryMatrix(*this);
  }
  
  
  //! Computation of Lagragian functions
  void TriangleMassLumped::ComputeLagrangianBubbleFunctions(int r)
  {        
    // 2-D basis functions
    InverseBasisVDM.Reallocate(nb_dof_loc, nb_dof_loc);
    InverseBasisVDM.Fill(Real_wp(0));
    
    // P_order
    int nb_coef_Pk = (r+1)*(r+2)/2;
    int nb_polynome = 0;
    // we use canonical basis
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r-i; j++)
	{
	  for (int k = 0; k < nb_dof_loc; k++)
	    InverseBasisVDM(nb_polynome, k) = pow(this->PointsDofND(k)(0),i)*pow(this->PointsDofND(k)(1),j);
	  
	  nb_polynome++;
	}
    
    // Base completion
    // allocation of additional polynoms
    Vector<MultivariatePolynomial<Real_wp> > additional_polynoms(nb_dof_loc - nb_coef_Pk);
    MultivariatePolynomial<Real_wp> x(2,1), y(2,1), one(2, 0);
    one(0, 0) = 1.0;
    x(1, 0) = 1.0; y(0, 1) = 1.0;
    
    if (r == 2)
      {
	// Bubble   (1-x-y)*x*y
	additional_polynoms(0) = (one-x-y)*x*y;
	for (int k = 0; k < nb_dof_loc; k++)
	  InverseBasisVDM(nb_dof_loc-1, k) = (1-this->PointsDofND(k)(0)-this->PointsDofND(k)(1))
	    *this->PointsDofND(k)(0)*this->PointsDofND(k)(1);
      }
    else if (r==3)
      {
	// DISP(nb_points_nodal_elt);
	if (nb_dof_loc == 13)
	  {
	    // 13 dof
	    additional_polynoms(0) = (one-x-y)*x*y*x;
	    additional_polynoms(1) = (one-x-y)*x*y*y;
	    additional_polynoms(2) = (one-x-y)*x*y*x*y;
	    for (int m = 0; m < additional_polynoms.GetM(); m++)
	      for (int k = 0; k < nb_dof_loc; k++)
		InverseBasisVDM(nb_coef_Pk+m, k)
                  = additional_polynoms(m).Evaluate(this->PointsNodalND(k));
	  }
	else
	  {
	    // 15 dof
	    additional_polynoms(0) = x*y*(one-x-y)*x;
	    additional_polynoms(1) = x*y*(one-x-y)*y; 
	    additional_polynoms(2) = x*y*(one-x-y)*x*y;
	    additional_polynoms(3) = x*y*(one-x-y)*x*x;
	    additional_polynoms(4) = x*y*(one-x-y)*y*y;
	    for (int m = 0; m < additional_polynoms.GetM(); m++)
	      for (int k = 0; k < nb_dof_loc; k++)
		InverseBasisVDM(nb_coef_Pk+m, k)
                  = additional_polynoms(m).Evaluate(this->PointsNodalND(k));
	  }

      }
    
    Seldon::GetInverse(InverseBasisVDM); 

    LagrangeFunction.Reallocate(nb_dof_loc); 
    for (int k = 0;  k < nb_dof_loc; k++)
      {
	LagrangeFunction(k).SetOrder(2, r);
	LagrangeFunction(k).Fill(0);
	// Pk part
	nb_polynome = 0;
	for (int i = 0; i <= r; i++)
	  for (int j = 0; j <= r-i; j++)
	    {
	      LagrangeFunction(k)(i,j) = InverseBasisVDM(k, nb_polynome);
	      nb_polynome++;
	    }
	
	// base completion with additional_polynoms
	for (int m = 0; m < additional_polynoms.GetM(); m++)
	  Add(InverseBasisVDM(k, m+nb_coef_Pk), additional_polynoms(m), LagrangeFunction(k));
      }  
    
    // we check that phi_i(\xi_j) = delta_{i,j}
    for (int j = 0; j < nb_dof_loc; j++)
      for (int i = 0; i < nb_dof_loc; i++)
	{
	  Real_wp val = LagrangeFunction(i).Evaluate(this->PointsDofND(j));
	  Real_wp val_exact = 0;
	  if (i == j)
	    val_exact = 1.0;
	  
	  if (abs(val-val_exact) > 1e5*epsilon_machine)
	    {
	      DISP(i); DISP(j); DISP(val);
	      abort();
	    }
	}
	    
    DxLagrangeFunction.Reallocate(nb_dof_loc); 
    DyLagrangeFunction.Reallocate(nb_dof_loc);
    
    for (int k = 0; k < nb_dof_loc; k++)
      {
	DerivatePolynomial(LagrangeFunction(k), DxLagrangeFunction(k), 0);
	DerivatePolynomial(LagrangeFunction(k), DyLagrangeFunction(k), 1);
      }
    
  }

  
  //! Evaluating basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] phi values of basis functions on pointloc
  */
  void TriangleMassLumped::ComputeValuesPhiRef(const R2& pointloc, VectReal_wp& phi) const
  {
    phi.Reallocate(this->GetNbPointsDof());
    for (int i = 0; i < this->GetNbPointsDof(); i++)
      phi(i) = LagrangeFunction(i).Evaluate(pointloc);
  }
  
  
  //! Evaluating gradient of basis functions on a point of the element
  /*!
    \param[in] pointloc local point where functions are evaluated
    \param[out] grad_phi gradient of basis functions on pointloc
  */
  void TriangleMassLumped::ComputeGradientPhiRef(const R2& pointloc, VectR2& grad_phi) const
  {
    grad_phi.Reallocate(nb_dof_loc);
    for (int i = 0; i < nb_dof_loc; i++)
      grad_phi(i).Init(DxLagrangeFunction(i).Evaluate(pointloc),
                       DyLagrangeFunction(i).Evaluate(pointloc));
  }
  
}

#define MONTJOIE_FILE_TRIANGLE_DG_MASSLUMPED_CXX
#endif
