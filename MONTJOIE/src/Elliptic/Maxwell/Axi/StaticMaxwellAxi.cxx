#ifndef MONTJOIE_FILE_STATIC_MAXWELL_AXI_CXX

namespace Montjoie
{

  //! marks dofs related to E_theta
  /*!
    \param[out] IsDofTeta IsDofTeta(i) true if dof i is related to E_theta
  */
  void EllipticProblem<StaticMaxwellEquation_HcurlAxi>
  ::GetDof_Eteta(VectBool& IsDofTeta) const
  {
    // ddl E^teta
    IsDofTeta.Reallocate(this->nodl); IsDofTeta.Fill(false);
    for (int i = this->offset_dof_unknown(1); i < this->offset_dof_unknown(2); i++)
      IsDofTeta(i) = true;
  }
  
  
  //! no need to mark dof related to E_theta
  template<class TypeEquation>
  void StaticMaxwell_Axi_DG<TypeEquation>::GetDof_Eteta(VectBool& IsDofTeta) const
  {
    // ddl E^teta
    IsDofTeta.Reallocate(this->nodl); IsDofTeta.Fill(false);
  }
  
  
  template<class TypeEquation>
  void StaticMaxwell_Axi_DG<TypeEquation>
  ::EvaluateH_MixedFormulation(const Vector<Complexe>& U0,
			       Vector<Complexe>& EvalH_Nodal) const
  {
  }
  

  //! \f$ E_\theta \; H_\theta \f$ are derived
  template<class TypeEquation, class T0, class Vector1>
  void StaticMaxwellEquation_Axi_DG
  ::GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
			const GlobalGenericMatrix<T0>& nat_mat,
			Vector1& unknown_to_derive, Vector1& fct_test_to_derive) 
  {
    unknown_to_derive.Fill(false);
    fct_test_to_derive.Fill(false);
    // on a besoin de ne deriver Eteta et Hteta
    fct_test_to_derive(1) = true;
    fct_test_to_derive(4) = true;
    unknown_to_derive(1) = true;
    unknown_to_derive(4) = true;
  }
  
  
  //! \f$ Un \, = \, \sum A^i Vn(i) \f$
  template<class Vector1, class Matrix1, class T0, class TypeEquation>
  void StaticMaxwellEquation_Axi_DG
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& vars,
			 int i,int  k, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, const Matrix1& Vn, Vector1& Un)
  {
    Un.Fill(0);
    Real_wp r = vars.Glob_rtilde(i)(k); // DISP(i); DISP(k); DISP(r);
    Un(0) = r*Vn(4)(1); Un(2) = -r*Vn(4)(0);
    Un(3) = r*Vn(1)(1); Un(5) = -r*Vn(1)(0);
    Mlt(nat_mat.GetCoefStiffness(), Un);
  }
  
  
  //! Vn(i) = A^i Un
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void StaticMaxwellEquation_Axi_DG
  ::ApplyGradientFctTest(const EllipticProblem<TypeEquation>& vars,
			 int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			 int ref, Vector1& Un, Vector2& Vn)
  {
    Vn.Fill(0);
    Real_wp r = vars.Glob_rtilde(i)(j); // DISP(i); DISP(k); DISP(r);
    Vn(1)(0) = r*Un(5); Vn(4)(0) = r*Un(2);
    Vn(1)(1) = -r*Un(3); Vn(4)(1) = -r*Un(0);
    Mlt(-nat_mat.GetCoefStiffness(), Vn);
  }
  
  
  //! matrices Ai
  template<class TypeEquation, class T0, class MatStiff>
  void StaticMaxwellEquation_Axi_DG
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    Real_wp r = vars.Glob_rtilde(num_elem)(jloc)*nat_mat.GetCoefStiffness();
    FillZero(Dgrad_phi); FillZero(Ephi_grad);
    
    Dgrad_phi(2,4)(0) = -r; Dgrad_phi(5,1)(0) = -r;
    Dgrad_phi(0,4)(1) = r; Dgrad_phi(3,1)(1) = r;
    
    Ephi_grad(4,2)(0) = -r; Ephi_grad(1,5)(0) = -r;
    Ephi_grad(4,0)(1) = r; Ephi_grad(1,3)(1) = r;
  }
  
  
  //! matrix N (boundary condition
  template<class Matrix1, class GenericPb, class T0>
  void StaticMaxwellEquation_Axi_DG
  ::GetNabc(Matrix1& Nabc, R2& normale, int ref, int iquad, int k,
	    const GlobalGenericMatrix<T0>& nat_mat, int ref2,
	    const GenericPb& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Nabc.Fill(0);
    Real_wp r = vars.Glob_rtilde(iquad)(k);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	Nabc(0,4) = r*normale(1); Nabc(1,3) = -r*normale(1); Nabc(1,5) = r*normale(0);
	Nabc(2,4) = -r*normale(0); Nabc(3,1) = -r*normale(1); Nabc(4,0) = r*normale(1);
	Nabc(4,2) = -r*normale(0); Nabc(5,1) = r*normale(0);
	Mlt(nat_mat.GetCoefStiffness(), Nabc);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	Nabc(0,4) = -r*normale(1); Nabc(1,3) = r*normale(1); Nabc(1,5) = -r*normale(0);
	Nabc(2,4) = r*normale(0); Nabc(3,1) = r*normale(1); Nabc(4,0) = -r*normale(1);
	Nabc(4,2) = r*normale(0); Nabc(5,1) = -r*normale(0);
	Mlt(nat_mat.GetCoefStiffness(), Nabc);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	Real_wp n00 = r*normale(0)*normale(0), n11 = r*normale(1)*normale(1),
	  n01 = r*normale(0)*normale(1);
	
	Nabc(0,0) = n11; Nabc(0,2) = -n01; Nabc(1,1) = r;
	Nabc(2,0) = -n01; Nabc(2,2) = n00; Nabc(3,3) = -n11;
	Nabc(3,5) = n01; Nabc(4,4) = -r; Nabc(5,3) = n01; Nabc(5,5) = -n00;
	Mlt(nat_mat.GetCoefStiffness(), Nabc);
      }
  }
  
  
  //! multiplication by matrix N
  template<class Vector1, class TypeEquation, class T0>
  void StaticMaxwellEquation_Axi_DG
  ::MltNabc(R2& normale, int ref, const Vector1& Vn, Vector1& Un, int num_elem1, int k,
	    const GlobalGenericMatrix<T0>& nat_mat, int ref2, 
	    const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension2, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref); 
    Un.Fill(0);
    Real_wp r = vars.Glob_rtilde(num_elem1)(k);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	Un(0) = r*Vn(4)*normale(1);
	Un(1) = -r*(Vn(3)*normale(1)-Vn(5)*normale(0));
	Un(2) = -r*Vn(4)*normale(0);
	Un(3) = -r*Vn(1)*normale(1);
	Un(4) = r*(Vn(0)*normale(1)-Vn(2)*normale(0));
	Un(5) = r*Vn(1)*normale(0);
	Mlt(nat_mat.GetCoefStiffness(), Un);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	Un(0) = -r*Vn(4)*normale(1);
	Un(1) = r*(Vn(3)*normale(1)-Vn(5)*normale(0));
	Un(2) = r*Vn(4)*normale(0);
	Un(3) = r*Vn(1)*normale(1);
	Un(4) = -r*(Vn(0)*normale(1)-Vn(2)*normale(0));
	Un(5) = -r*Vn(1)*normale(0);
	Mlt(nat_mat.GetCoefStiffness(), Un);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	Real_wp En = Vn(0)*normale(1)-Vn(2)*normale(0);
	Real_wp Hn = Vn(3)*normale(1)-Vn(5)*normale(0);
	Un(0) = r*En*normale(1);
	Un(1) = r*Vn(1);
	Un(2) = -r*En*normale(0);
	Un(3) = -r*Hn*normale(1);
	Un(4) = -r*Vn(4);
	Un(5) = r*Hn*normale(0);
	Mlt(nat_mat.GetCoefStiffness(), Un);
      }
    // DISP(Un); DISP(Vn);
  }    
  
  
  //! multiplication by matrix C
  template<class TypeEquation, class T0, class Vector1>
  void StaticMaxwellEquation_Axi_DG
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
		    const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& Un, Vector1& Vn)
  {
    Vn.Fill(0);
    //Real_wp r = var.Glob_rtilde(i)(j);
    int m = var.GetCurrentModeNumber();
    Vn(0) = -m*Un(5); Vn(1) = -Un(5); Vn(2) = m*Un(3)-Un(4);
    Vn(3) = m*Un(2); Vn(4) = -Un(2); Vn(5) = -m*Un(0)-Un(1);
    Mlt(nat_mat.GetCoefStiffness(), Vn);
  }
  
  
  //! matrix C
  template<class TypeEquation, class T0, class MatMass>
  void StaticMaxwellEquation_Axi_DG
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		  int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& Cj)
  {
    FillZero(Cj);
    
    int m = vars.GetCurrentModeNumber();
    Cj(0,5) = -m; Cj(1,5) = -1.0; Cj(2,3) = m; Cj(2,4) = -1.0;
    Cj(3,2) = m; Cj(4,2) = -1.0; Cj(5,0) = -m; Cj(5,1) = -1.0;
    Mlt(nat_mat.GetCoefStiffness(), Cj);
  }
  
} 

#define MONTJOIE_FILE_STATIC_MAXWELL_AXI_CXX
#endif
