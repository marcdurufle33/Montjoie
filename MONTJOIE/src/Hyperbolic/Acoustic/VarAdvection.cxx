#ifndef MONTJOIE_FILE_VAR_ADVECTION_CXX

namespace Montjoie
{

  template<class T, class Dim>  template<class TypeEquation, class T0, class MatStiff>
  void AdvectionEquation_Base<T, Dim>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    Dgrad_phi.Zero();
    
    Ephi_grad(0, 0) = vars.ref_v0(ref).GetCoefficient(vars, num_elem, jloc);
    Ephi_grad *= 0.5*vars.ref_rho(ref).GetCoefficient(vars, num_elem, jloc)*nat_mat.GetCoefStiffness();
    
    Dgrad_phi(0, 0) = -Ephi_grad(0, 0);
  }
  
  
  template<class T, class Dim> template<class TypeEquation, class T0, class MatMass>
  void AdvectionEquation_Base<T, Dim>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
		  int ref, MatMass& mass)
  {
    mass(0, 0) = vars.ref_rho(ref).GetCoefficient(vars, num_elem, jloc);
    T sigma = vars.ref_sigma(ref).GetCoefficient(vars, num_elem, jloc);
    T m_iomega; vars.GetMiomega(m_iomega);
    mass(0, 0) *= (m_iomega*nat_mat.GetCoefMass() + sigma*nat_mat.GetCoefDamping());      
  }
  
  
  // computation of impedance
  template<class Complexe, class Dimension>
  void ImpedanceABC<Complexe, AdvectionEquation<Dimension> >
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			 const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
			 const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    if (k == 0)
      {
	int nb_points = Mat.GetNbPointsQuadratureBoundary();
	this->stored_coef_phi(0).Reallocate(nb_points);
	this->stored_coef_phi(0).Fill(0);
      }
    
    Complexe feval = DotProd(Mat.GetNormaleQuadratureBoundary(k),
			     var_helm.ref_v0(ref).GetCoefficient(var_problem, num_elem, 0));
    
    feval = 0.5*abs(feval)*var_helm.ref_rho(ref).GetCoefficient(var_problem, num_elem, 0)
      *nat_mat.GetCoefStiffness();
    
    this->stored_coef_phi(0)(k) = feval;
  }

  // computation of impedance
  template<class Dimension> 
  void ImpedanceABC<Complex_wp, HarmonicAdvectionEquation<Dimension> > 
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			 const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
			 const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    if (k == 0)
      {
	int nb_points = Mat.GetNbPointsQuadratureBoundary();
	this->stored_coef_phi(0).Reallocate(nb_points);
	this->stored_coef_phi(0).Fill(0);
      }
    
    Complex_wp feval = DotProd(Mat.GetNormaleQuadratureBoundary(k),
			       var_helm.ref_v0(ref).GetCoefficient(var_problem, num_elem, 0));
      
    feval = 0.5*abs(feval)*var_helm.ref_rho(ref).GetCoefficient(var_problem, num_elem, 0)
      *nat_mat.GetCoefStiffness();
    
    this->stored_coef_phi(0)(k) = feval;
  }
  

  template<class Dimension>    
  void EllipticProblem<AdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class Dimension>    
  void EllipticProblem<AdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicAdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  //! class to solve advection equation with H1 elements in time-harmonic domain
  template<class Dimension>
  void EllipticProblem<HarmonicAdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class T, class Dim> template<class Matrix1, class T0, class GenericPb>
  void AdvectionEquationDG_Base<T, Dim>
  ::GetNabc(Matrix1& Nabc, typename Dim::R_N& normale,
	    int ref, int iquad, int k,const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
	    const GenericPb& vars, const ElementReference<Dim, 1>& Fb)
  {
    Nabc.Fill(0);
    int cond = vars.mesh.GetBoundaryCondition(ref);
    if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	T0 v0_dot_n = DotProd(normale, vars.ref_v0(ref_d).GetCoefficient(vars, iquad, k));
        Nabc(0, 0) = abs(v0_dot_n)*vars.ref_rho(ref_d).GetCoefficient(vars, iquad, k);
        Nabc(0, 0) *= nat_mat.GetCoefStiffness();
      }
  }
  
  
  template<class T, class Dim>  template<class Matrix1, class T0, class GenericPb>
  void AdvectionEquationDG_Base<T, Dim>
  ::GetPenalDG(Matrix1& Nabc, typename Dim::R_N& normale,
	       int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	       const GenericPb& vars, const ElementReference<Dim, 1>& Fb)
  {
    T0 v0_dot_n = DotProd(normale, vars.ref_v0(ref).GetCoefficient(vars, iquad, k));      
    Nabc(0, 0) = abs(v0_dot_n)*nat_mat.GetCoefDamping()*vars.alpha_penalization;
    Nabc(0, 0) *= vars.ref_rho(ref).GetCoefficient(vars, iquad, k);
  }
  
  
  template<class Dimension>
  void EllipticProblem<AdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class Dimension>
  void EllipticProblem<AdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  

  template<class Dimension>
  void EllipticProblem<AdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class Dimension>
  void EllipticProblem<AdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicAdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicAdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }
  

  template<class Dimension>
  void EllipticProblem<HarmonicAdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicAdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class T, class Dim> template<class TypeEquation, class T0, class MatStiff>
  void DoubleAdvectionEquation_Base<T, Dim>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
		     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi); FillZero(Ephi_grad);
    
    Ephi_grad(0, 0) = vars.ref_v0(ref).GetCoefficient(vars, num_elem, jloc);
    Ephi_grad(0, 0) *= 0.5*vars.ref_rho(ref).GetCoefficient(vars, num_elem, jloc)*nat_mat.GetCoefStiffness();
    Ephi_grad(1, 1) = Ephi_grad(0, 0);
    
    Dgrad_phi(0, 0) = -Ephi_grad(0, 0);
    Dgrad_phi(1, 1) = Dgrad_phi(0, 0);
  }
  
  
  template<class T, class Dim> template<class TypeEquation, class T0, class MatMass>
  void DoubleAdvectionEquation_Base<T, Dim>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
		  int ref, MatMass& mass)
  {
    mass(0, 0) = vars.ref_rho(ref).GetCoefficient(vars, num_elem, jloc);
    T sigma = vars.ref_sigma(ref).GetCoefficient(vars, num_elem, jloc);
    T m_iomega; vars.GetMiomega(m_iomega);
    mass(0, 0) *= (m_iomega*nat_mat.GetCoefMass() + sigma*nat_mat.GetCoefDamping());
    mass(1, 1) = mass(0, 0);
    mass(0, 1) = -vars.ref_rho(ref).GetCoefficient(vars, num_elem, jloc)*nat_mat.GetCoefStiffness();
  }
  
  
  template<class T, class Dim> template<class Matrix1, class T0, class GenericPb>
  void DoubleAdvectionEquation_Base<T, Dim>
  ::GetNabc(Matrix1& Nabc, typename Dim::R_N& normale,
	    int ref, int iquad, int k,const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
	    const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Fill(0);
    int cond = vars.mesh.GetBoundaryCondition(ref);
    if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	T0 v0_dot_n = DotProd(normale, vars.ref_v0(ref_d).GetCoefficient(vars, iquad, k));
	if (real(v0_dot_n) > 0)
	  {              
	    Nabc(0, 0) = v0_dot_n*vars.ref_rho(ref_d).GetCoefficient(vars, iquad, k);
	    Nabc(0, 0) *= nat_mat.GetCoefStiffness();
	    Nabc(1, 1) = Nabc(0, 0);
	  }
      }
  }
  
  
  // computation of impedance
  template<class Complexe, class Dimension>
  void ImpedanceABC_DoubleAdvec<Complexe, Dimension>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			 const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
			 const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    if (k == 0)
      {
	int nb_points = Mat.GetNbPointsQuadratureBoundary();
	this->stored_coef_phi(0).Reallocate(nb_points);
	this->stored_coef_phi(0).Fill(0);
      }
    
    Complexe feval = DotProd(Mat.GetNormaleQuadratureBoundary(k),
			     var_helm.ref_v0(ref).GetCoefficient(var_problem, num_elem, 0));
    
    feval *= var_helm.ref_rho(ref).GetCoefficient(var_helm.GetVarProblem(), num_elem, 0);
    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    feval = abs(feval)*nat_mat.GetCoefDamping()*m_iomega;
    this->stored_coef_phi(0)(k) = feval;
  }


  template<class Dimension>
  void EllipticProblem<DoubleAdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class Dimension>
  void EllipticProblem<DoubleAdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }
  
   
  template<class Dimension>
  void EllipticProblem<HarmonicDoubleAdvectionEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class T, class Dim> template<class Matrix1, class T0, class GenericPb>
  void DoubleAdvectionEquationDG_Base<T, Dim>
  ::GetPenalDG(Matrix1& Nabc, typename Dim::R_N& normale,
	       int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	       const GenericPb& vars, const ElementReference<Dim, 1>& Fb)
  {
    T0 v0_dot_n = DotProd(normale, vars.ref_v0(ref).GetCoefficient(vars, iquad, k));
    Nabc(0, 0) = abs(v0_dot_n)*nat_mat.GetCoefDamping()*vars.alpha_penalization;
    Nabc(0, 0) *= vars.ref_rho(ref).GetCoefficient(vars, iquad, k);
    Nabc(1, 1) = Nabc(0, 0);
  }

  
  template<class Dimension>
  void EllipticProblem<DoubleAdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class Dimension>
  void EllipticProblem<DoubleAdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>  
  void EllipticProblem<DoubleAdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class Dimension>  
  void EllipticProblem<DoubleAdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }
  

  template<class Dimension>
  void EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicDoubleAdvectionEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }
  
}

#define MONTJOIE_FILE_VAR_ADVECTION_CXX
#endif
