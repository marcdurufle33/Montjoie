#ifndef MONTJOIE_FILE_VAR_ELASTIC_CXX

namespace Montjoie
{
  
  /*******************
   * ElasticEquation *
   *******************/
  
  
  //! not used
  template<class T, class Dimension> 
  template<class TypeEquation>
  void ElasticEquation_Base<T, Dimension>::
  ComputeMassMatrix(const EllipticProblem<TypeEquation>& var,
		    int i, const ElementReference_Dim<Dimension>& Fb)
  {
  }
    
  
  //! multiplication by elastic tensor C
  template<class T, class Dimension> 
  template<class TypeEquation, class T0, class Vector1>
  void ElasticEquation_Base<T, Dimension>::
  ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var, int i, int j,
		       const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& dU, Vector1& dV)
  {
    var.ref_tensorC(ref).MltMatrix(var, i, j, dU, dV); 
    Mlt(nat_mat.GetCoefStiffness(), dV);
  }
  
  
  //! elastic tensor C is retrieved
  template<class T, class Dimension> 
  template<class TypeEquation, class T0, class MatStiff>
  void ElasticEquation_Base<T, Dimension>::
  GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
		    int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
                    int ref, MatStiff& Cgrad_grad)
  {
    vars.ref_tensorC(ref).GetCoefficient(vars, num_elem, jloc, Cgrad_grad);
    Mlt(nat_mat.GetCoefStiffness(), Cgrad_grad);
  }
  

  template<class T, class Dimension> 
  template<class Matrix1, class T0, class GenericPb>
  void ElasticEquation_Base<T, Dimension>
  ::GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale,
	       int iquad, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
	       int ref, int ref2,
	       const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.SetIdentity();
    Nabc *= vars.alpha_penalization*nat_mat.GetCoefStiffness();
  }
  
  
  //! mass tensor M is computed
  template<class Dimension>
  template<class TypeEquation, class T0, class MatMass>
  void ElasticEquation<Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    T0 coef = vars.ref_rho(ref).GetCoefficient(vars, i, j)*nat_mat.GetCoefMass();
    coef += vars.ref_sigma(ref).GetCoefficient(vars, i, j)*nat_mat.GetCoefDamping();
    
    mass.SetIdentity();
    mass *= coef;
  }
  
  
  //! application of mass tensor, V = M U
  template<class Dimension>
  template<class TypeEquation, class T0, class Vector1>
  void ElasticEquation<Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    T0 coef = var.ref_rho(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefMass();
    coef += var.ref_sigma(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefDamping();
    
    V = U;
    V *= coef;
  }
  

  //! mass tensor M is computed  
  template<class Dimension>
  template<class TypeEquation, class T0, class MatMass>
  void HarmonicElasticEquation<Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& var,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    Complex_wp coef = -var.GetSquareOmega()
      *var.ref_rho(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefMass();
    
    coef -= Iwp*var.GetOmega()
      *var.ref_sigma(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefDamping();
    
    mass.SetIdentity();
    mass *= coef;
  }
  

  //! application of mass tensor, V = M U  
  template<class Dimension>
  template<class TypeEquation, class T0, class Vector1>
  void HarmonicElasticEquation<Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
  {
    Complex_wp coef = -var.GetSquareOmega()
      *var.ref_rho(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefMass();
    
    coef -= Iwp*var.GetOmega()
      *var.ref_sigma(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefDamping();
    
    V = U;
    V *= coef;
  }
  
  
  /*******************
   * VarElastic_Base *
   *******************/
  
  
  //! reads a line of the data file
  void VarElastic_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "DisplayStress")
      {
	if (parameters(0) == "NO")
	  display_stress = false;
	else
	  display_stress = true;
      }
  }
  
  
  //! modification of shift due to adimensionalization
  template<class Complexe>
  void VarElastic_Base::UpdateShiftAdimensionalization(Complexe& shift, Complexe& shift_imag)
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        if (var_problem_base.FirstOrderFormulation())
          {
            shift /= omega_bar;
            shift_imag /= omega_bar;
          }
        else
          {
            shift /= omega_bar*omega_bar;
            shift_imag /= omega_bar*omega_bar;
          }
      }
  }
    
  
  //! modification of eigenvalues and eigenvectors due to adimensionalization
  template<class Complexe>
  void VarElastic_Base::
  UpdateEigenvaluesAdimensionalization(Vector<Complexe>& lambda, Vector<Complexe>& lambda_imag,
                                       Matrix<Complexe, General, ColMajor>& eigen_vec)
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        Complexe coef; SetComplexOne(coef);
        if (var_problem_base.FirstOrderFormulation())
          coef = omega_bar;
        else
          coef = omega_bar*omega_bar;
        
        Mlt(coef, lambda);
        Mlt(coef, lambda_imag);
      }
  }


  template<class T0, class Complexe>
  void VarElastic_Base
  ::ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension2, 2, T0> >& ref_tensorC,
			const VarProblem<Dimension2>& var_problem,
			Vector<Complexe>& val_u,
			Vector<Complexe>& grad_u, int i,
			const GridInterpolation<Dimension2>& var_interp,
			int iquad, bool compute_grad) const
  {
    cout << "Not implemented" << endl;
    abort();
  }

  
  template<class Complexe>
  void VarElastic_Base
  ::ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension2, 2, Complexe> >& ref_tensorC,
			const VarProblem<Dimension2>& var_problem,
			Vector<Complexe>& val_u,
			Vector<Complexe>& grad_u, int i,
			const GridInterpolation<Dimension2>& var_interp,
			int iquad, bool compute_grad) const
  {
    if (!display_stress)
      return;
    
    if (!compute_grad)
      return;
    
    int ref = var_problem.mesh.Element(iquad).GetReference();
    TinyVector<Complexe, 4> dU, dV;
    dU(0) = grad_u(0); dU(1) = grad_u(1);
    dU(2) = grad_u(2); dU(3) = grad_u(3);
    
    if (ref_tensorC(ref).IsVarying())
      {
	VectReal_wp phi;
	var_problem.mesh.ComputeValuesPhiNodalRef(iquad, var_interp.GetLocalCoordinate(i), phi);
	ref_tensorC(ref).MltVectorPoint(var_problem.ElementRho(iquad), phi,
					var_interp.GetGlobalCoordinate(i), dU, dV);
      }
    else
      ref_tensorC(ref).MltVector(var_problem, iquad, 0, dU, dV);

    grad_u(0) = dV(0); grad_u(1) = dV(1);
    grad_u(2) = dV(2); grad_u(3) = dV(3);
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  template<class T0, class Complexe>
  void VarElastic_Base
  ::ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension3, 3, T0> >& ref_tensorC,
			const VarProblem<Dimension3>& var_problem,
			Vector<Complexe>& val_u,
			Vector<Complexe>& grad_u, int i,
			const GridInterpolation<Dimension3>& var_interp,
			int iquad, bool compute_grad) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Complexe>
  void VarElastic_Base
  ::ModifyOutputUnknown(const Vector<ElasticPhysicalIndice<Dimension3, 3, Complexe> >& ref_tensorC,
			const VarProblem<Dimension3>& var_problem,
			Vector<Complexe>& val_u,
			Vector<Complexe>& grad_u, int i,
			const GridInterpolation<Dimension3>& var_interp,
			int iquad, bool compute_grad) const
  {
    if (!display_stress)
      return;

    if (!compute_grad)
      return;
    
    int ref = var_problem.mesh.Element(iquad).GetReference();
    TinyVector<Complexe, 9> dU, dV;
    dU(0) = grad_u(0); dU(1) = grad_u(1); dU(2) = grad_u(2);
    dU(3) = grad_u(3); dU(4) = grad_u(4); dU(5) = grad_u(5);
    dU(6) = grad_u(6); dU(7) = grad_u(7); dU(8) = grad_u(8);
    
    if (ref_tensorC(ref).IsVarying())
      {
	VectReal_wp phi;
	var_problem.mesh.ComputeValuesPhiNodalRef(iquad, var_interp.GetLocalCoordinate(i), phi);
	ref_tensorC(ref).MltVectorPoint(var_problem.ElementRho(iquad), phi,
					var_interp.GetGlobalCoordinate(i), dU, dV);
      }
    else
      ref_tensorC(ref).MltVector(var_problem, iquad, 0, dU, dV);

    grad_u(0) = dV(0); grad_u(1) = dV(1); grad_u(2) = dV(2);
    grad_u(3) = dV(3); grad_u(4) = dV(4); grad_u(5) = dV(5);
    grad_u(6) = dV(6); grad_u(7) = dV(7); grad_u(8) = dV(8);
  }
#endif
  
  
  //! computes u and C nabla u n from u and C nabla u
  template<class Complexe>
  void VarElastic_Base
  ::ComputeEnHnQuadrature(Vector<Vector<Complexe> >& u_quadrature,
			  Vector<Vector<Complexe> >& grad_quadrature,
			  int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
			  bool compute_H, Vector<Complexe>& En_quad, Vector<Complexe>& Hn_quad) const
  {
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(2*nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(2*nb_points_quad);
    
    if (var_problem_base.FirstOrderFormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
	cout << "not implemented" << endl;
	abort();
      }
    else
      {
        for (int j = 0; j < nb_points_quad; j++)
          {
	    En_quad(2*j) = u_quadrature(0)(j);
	    En_quad(2*j+1) = u_quadrature(1)(j);
	    
	    if (compute_H)
	      {
		cout << "not implemented" << endl;
		abort();
		// compute dV = C grad(u)
		//Hn_quad(2*j) = DotProd(dV, normale(j));
		//Hn_quad(2*j+1) = DotProd(dV, normale(j));
	      }
          }
      }
  }  


#ifdef MONTJOIE_WITH_THREE_DIM
  //! computes u and C nabla u n from u and C nabla u
  template<class Complexe>
  void VarElastic_Base
  ::ComputeEnHnQuadrature(Vector<Vector<Complexe> >& u_quadrature,
			  Vector<Vector<Complexe> >& grad_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, Vector<Complexe>& En_quad, Vector<Complexe>& Hn_quad) const
  {
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(3*nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(3*nb_points_quad);
    
    if (var_problem_base.FirstOrderFormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
	cout << "not implemented" << endl;
	abort();
      }
    else
      {
        for (int j = 0; j < nb_points_quad; j++)
          {
	    En_quad(3*j) = u_quadrature(0)(j);
	    En_quad(3*j+1) = u_quadrature(1)(j);
	    En_quad(3*j+2) = u_quadrature(2)(j);
	    
	    if (compute_H)
	      {
		cout << "not implemented" << endl;
		abort();
		// compute dV = C grad(u)
		//Hn_quad(3*j) = DotProd(dV, normale(j));
		//Hn_quad(3*j+1) = DotProd(dV, normale(j));
		//Hn_quad(3*j+2) = DotProd(dV, normale(j));
	      }
          }
      }
  }  
#endif
  

  template<class Complexe>
  void VarElastic_Base
  ::ComputeEnHnNodal(Vector<Vector<Complexe> >& u_nodal, Vector<Vector<Complexe> >& grad_nodal,
		     int num_elem, const Vector<R2>& pts, const Vector<R2>& normale,
		     Vector<Vector<Complexe> >& En_nodal, Vector<Vector<Complexe> >& Hn_nodal) const
  {
    abort();
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  template<class Complexe>
  void VarElastic_Base
  ::ComputeEnHnNodal(Vector<Vector<Complexe> >& u_nodal, Vector<Vector<Complexe> >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<Vector<Complexe> >& En_nodal, Vector<Vector<Complexe> >& Hn_nodal) const
  {
    abort();
  }
#endif
  

  //! returns the coefficient |mu|
  Real_wp VarElastic_Base::GetCoefficientPenaltyStiffness(int ref) const
  {
    return coefficient_mu_penalty(ref);
  }
  

  //! computation of elementary matrix for elastodynamic equation
  /*!
    \param[in] iquad element number
    \param[out] num_dof numbers of degrees of freedom for the element
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] var_problem considered problem
    \param[in] Fb finite element associated with the element
  */
  template<class T0, class Complexe>
  void VarElastic_Base
  ::ComputeElementaryMatrix(const Vector<ElasticPhysicalIndice<Dimension2, 2, T0> >& ref_tensorC,
			    const Vector<ScalarPhysicalIndice<Dimension2, T0> >& ref_rho,
			    const Vector<ScalarPhysicalIndice<Dimension2, T0> >& ref_sigma,
			    int iquad, IVect& num_dof, VirtualMatrix<Complexe>& mat_interac,
			    const GlobalGenericMatrix<Complexe>& nat_mat,
			    const VarProblem<Dimension2>& var_problem,
			    const VarBoundaryCondition_Dim<T0, Dimension2>& var_boundary,
			    const ElementReference<Dimension2, 1>& Fb)
  {
    // for first-order (in time) formulation
    // specific expression with split pml
    
    // without PML, the considered equations are (in 2-D) :
    // du_x/dt - d/dx(sigma_xx) - d/dy(sigma_xy) = 0
    // du_y/dt - d/dx(sigma_xy) - d/dy(sigma_yy) = 0
    // dsigma/dt - C nabla u = 0
    
    // it should be noticed that the symmetry of sigma is exploited
    // such that only sigma_xx, sigma_xy, sigma_yy are computed in 2-D
    // sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz in 3-D
    
    // With PML, u_x, u_y, sigma_xx, sigma_xy, etc are splitted in d unknowns
    //      (u_x^1, u_x^2, u_x^3 for u_x)
    // In order to use continuous approximations for u, we use as unknowns :
    // u = u_x^1 + u_x^2 + u_x^3
    // u^* = u_x^1 - u_x^2
    // u^d = u_x^1 - u_x^3
    // as it is done for acoustics equations
    // For sigma, since they are discontinuous, we keep sigma_xx^1, sigma_xx^2, etc

    // number of integration points
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    
    // number of degrees of freedom
    int nb_dof_elt = Fb.GetNbDof();
    int nb_dof_u = 2*nb_dof_elt;
    int nb_comp = 3;
    int nb_dof_sigma = nb_comp*nb_points_quad;
    int nb_dof_all = nb_dof_u + nb_dof_sigma;
    int offset_sigma = nb_dof_u;
    int offset_sigma2(0);
    
    if (var_problem.InsidePML(iquad))
      {
        nb_dof_all += nb_dof_u + nb_dof_sigma;
        offset_sigma += nb_dof_u;
        offset_sigma2 = offset_sigma + nb_dof_sigma;
      }
    
    // meshes
    const Mesh<Dimension2>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = var_problem.offset_dof_unknown(1);
    int offset_sigma_glob = Nvol*2;
    offset_sigma_glob += var_problem.GetOffsetDofV(iquad);
    
    // reference of the domain
    int ref_domain = var_problem.mesh.Element(iquad).GetReference();
    
    // dof numbers
    num_dof.Reallocate(nb_dof_all); num_dof.Fill(-1);
    IVect Nodle = var_problem.GetDofNumberOnElement(iquad);
    for (int j = 0; j < nb_dof_elt; j++)
      for (int m = 0; m < 2; m++)
        num_dof(j+m*nb_dof_elt) = Nodle(j) + m*Nvol;
    
    for (int j = 0; j < nb_points_quad; j++)
      for (int m = 0; m < nb_comp; m++)
        num_dof(offset_sigma + j*nb_comp + m) = offset_sigma_glob + j*nb_comp + m; 
    
    int i1 = -1;
    // additional dofs for pml
    if (var_problem.InsidePML(iquad))
      {
        // additional dofs for u
        for (int j = 0; j < nb_dof_elt; j++)
          {
            int num_pml = mesh_num.GetDofPML(Nodle(j));
            for (int m = 0; m < 2; m++)
	      num_dof(nb_dof_u + m*nb_dof_elt + j) = Nvol*m + mesh_num.GetNbDof() + num_pml;
          }
	
        // additional dofs for sigma
        for (int j = 0; j < nb_points_quad; j++)
          {
            for (int m = 0; m < nb_comp; m++)
              num_dof(offset_sigma2 + j*nb_comp + m)
                = offset_sigma_glob + (nb_points_quad+j)*nb_comp + m;
	  }
	
        i1 = iquad - mesh.GetNbElt() + var_boundary.GetNbEltPML();
      }
    
    // initialisation of the matrix
    mat_interac.Clear();
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();

    // if true, stiffness matrix must be computed
    bool stiff = false;
    if (nat_mat.GetCoefStiffness() != Real_wp(0))
      stiff = true;

    bool affine = var_problem.mesh.IsElementAffine(iquad);
    bool pml = var_problem.InsidePML(iquad);

    // loop over dofs 
    Real_wp jacobian; Complexe poids;
    VectReal_wp ones(nb_dof_elt), val_phi(nb_points_quad), grad_chap_phi(2*nb_points_quad);
    
    Vector<Complexe> feval(nb_points_quad), fevalTauX(nb_points_quad), fevalTauY(nb_points_quad),
      contrib(nb_dof_elt);
    
    Vector<R2> grad_phi(nb_points_quad); 
    TinyVector<Real_wp, 2> vec_u;
    TinyVector<Complexe, 2> tau;
    val_phi.Fill(0); grad_chap_phi.Fill(0);
    feval.Fill(0); contrib.Fill(0);
    
    T0 iomega, c0000, c0001, c0011,
      c0101, c0111, c1111, rho, rho_damping;
    
    Complexe coef, coef_damping;
    TinyMatrix<TinyMatrix<T0, General, 2, 2>, General, 2, 2> C;
    var_problem.GetMiomega(iomega);
    
    // mass matrix for sigma
    for (int k = 0; k < nb_points_quad; k++)
      {
	if (affine)
	  jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	else
	  jacobian = var_problem.Glob_jacobian(iquad)(k);
	
	coef = iomega*jacobian*nat_mat.GetCoefMass();              
	mat_interac.SetEntry(offset_sigma + 3*k, offset_sigma + 3*k, coef);
	mat_interac.SetEntry(offset_sigma + 3*k+1, offset_sigma + 3*k+1, coef);
        mat_interac.SetEntry(offset_sigma + 3*k+2, offset_sigma + 3*k+2, coef);
	if (pml)
	  {
	    tau = var_boundary.GetTauPML(i1, k);
	    tau *= jacobian*nat_mat.GetCoefDamping();
	    
	    mat_interac.AddInteraction(offset_sigma + 3*k, offset_sigma + 3*k, tau(0));
	    mat_interac.AddInteraction(offset_sigma + 3*k+1, offset_sigma + 3*k+1, tau(0));
	    mat_interac.AddInteraction(offset_sigma + 3*k+2, offset_sigma + 3*k+2, tau(0)); 
            
	    mat_interac.SetEntry(offset_sigma2 + 3*k, offset_sigma2 + 3*k, coef + tau(1));
            mat_interac.SetEntry(offset_sigma2 + 3*k+1, offset_sigma2 + 3*k+1, coef + tau(1));
	    mat_interac.SetEntry(offset_sigma2 + 3*k+2, offset_sigma2 + 3*k+2, coef + tau(1)); 
	  }
      }
    
    // stiffness matrix and mass matrix for u
    Real_wp max_dfjm1 = 1.0;
    for (int i = 0; i < nb_dof_elt; i++)
      {
        ones.Fill(0); ones(i) = 1.0;
        // computation of phi_i(xi_k) and grad phi_i(xi_k)
        Fb.ApplyChTranspose(ones, val_phi);
        
        if (stiff)
          {
            if (Fb.UseQuadraturePointsForRh())
              Fb.ApplyRhQuadratureTranspose(val_phi, grad_chap_phi);
            else
              Fb.ApplyRhTranspose(ones, grad_chap_phi);
            
            if (affine)
              {
                Matrix2_2 dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
                max_dfjm1 = MaxAbs(dfjm1);
                for (int k = 0; k < nb_points_quad; k++)
                  {
                    CopyVector(grad_chap_phi, k, vec_u);
                    MltTrans(dfjm1, vec_u, grad_phi(k));
                  }
              }
            else
              {
                max_dfjm1 = MaxAbs(var_problem.Glob_DFjm1(iquad)(0));
                for (int k = 0; k < nb_points_quad; k++)
                  {
                    CopyVector(grad_chap_phi, k, vec_u);
                    MltTrans(var_problem.Glob_DFjm1(iquad)(k), vec_u, grad_phi(k));
                  }
              }
          }
	
	
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    if ((abs(val_phi(k)) > epsilon_machine)
		|| (Norm2(grad_phi(k)) > max_dfjm1*epsilon_machine))
	      {
		if (affine)
		  jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
		else
		  jacobian = var_problem.Glob_jacobian(iquad)(k);
		
		if (stiff)
		  {
		    poids = Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
		    ref_tensorC(ref_domain).GetCoefficient(var_problem, iquad, k, C);
		    
		    c0000 = C(0, 0)(0, 0);
		    c0001 = C(0, 0)(0, 1);
		    c0011 = C(0, 1)(0, 1);
		    c0101 = C(0, 0)(1, 1);
		    c0111 = C(0, 1)(1, 1);
		    c1111 = C(1, 1)(1, 1);
                    
		    // part div (sigma)
		    mat_interac.SetEntry(i, offset_sigma + 3*k, poids*grad_phi(k)(0));
		    mat_interac.SetEntry(i, offset_sigma + 3*k+1, poids*grad_phi(k)(1));
                    
		    mat_interac.SetEntry(nb_dof_elt+i, offset_sigma + 3*k+1, poids*grad_phi(k)(0));
		    mat_interac.SetEntry(nb_dof_elt+i, offset_sigma + 3*k+2, poids*grad_phi(k)(1));
                    
		    if (pml)
		      {
			// part div(sigma^2)
			mat_interac.SetEntry(i, offset_sigma2 + 3*k, poids*grad_phi(k)(0));
			mat_interac.SetEntry(i, offset_sigma2 + 3*k+1, poids*grad_phi(k)(1));
                        
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma2 + 3*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma2 + 3*k+2,
                                             poids*grad_phi(k)(1));
			
			// part for u^*
			mat_interac.SetEntry(2*nb_dof_elt + i, offset_sigma + 3*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(2*nb_dof_elt + i, offset_sigma + 3*k+1,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(3*nb_dof_elt+i, offset_sigma + 3*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(3*nb_dof_elt+i, offset_sigma + 3*k+2,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(2*nb_dof_elt + i, offset_sigma2 + 3*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(2*nb_dof_elt + i, offset_sigma2 + 3*k+1,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(3*nb_dof_elt+i, offset_sigma2 + 3*k+1,
                                             poids*grad_phi(k)(0));
                        
			mat_interac.SetEntry(3*nb_dof_elt+i, offset_sigma2 + 3*k+2,
                                             -poids*grad_phi(k)(1));
			
			// part for sigma^1, sigma^2
			mat_interac.SetEntry(offset_sigma + 3*k, i, -poids*(c0000*grad_phi(k)(0)));
			mat_interac.SetEntry(offset_sigma + 3*k, nb_dof_elt+i,
                                             -poids*(c0001*grad_phi(k)(0)));
			
			mat_interac.SetEntry(offset_sigma + 3*k+1, i, -poids*(c0001*grad_phi(k)(0)));
			mat_interac.SetEntry(offset_sigma + 3*k+1, nb_dof_elt+i,
                                             -poids*(c0101*grad_phi(k)(0)));
			
			mat_interac.SetEntry(offset_sigma + 3*k+2, i, -poids*(c0011*grad_phi(k)(0)));
			mat_interac.SetEntry(offset_sigma + 3*k+2, nb_dof_elt+i,
                                             -poids*(c0111*grad_phi(k)(0)));
			
			mat_interac.SetEntry(offset_sigma2 + 3*k, i, -poids*(c0001*grad_phi(k)(1)));
			mat_interac.SetEntry(offset_sigma2 + 3*k, nb_dof_elt+i,
                                             -poids*(c0011*grad_phi(k)(1)));
			
			mat_interac.SetEntry(offset_sigma2 + 3*k+1, i, -poids*(c0101*grad_phi(k)(1)));
			mat_interac.SetEntry(offset_sigma2 + 3*k+1, nb_dof_elt+i,
                                             -poids*(c0111*grad_phi(k)(1)));
			
			mat_interac.SetEntry(offset_sigma2 + 3*k+2, i, -poids*(c0111*grad_phi(k)(1)));
			mat_interac.SetEntry(offset_sigma2 + 3*k+2, nb_dof_elt+i,
                                             -poids*(c1111*grad_phi(k)(1)));
		      }
		    else
		      {
			// part C grad(u)
			mat_interac.SetEntry(offset_sigma + 3*k, i,
                                             -poids*(c0000*grad_phi(k)(0) + c0001*grad_phi(k)(1)));
			mat_interac.SetEntry(offset_sigma + 3*k, nb_dof_elt+i,
                                             -poids*(c0001*grad_phi(k)(0) + c0011*grad_phi(k)(1)));
			
			mat_interac.SetEntry(offset_sigma + 3*k+1, i,
                                             -poids*(c0001*grad_phi(k)(0) + c0101*grad_phi(k)(1)));
			mat_interac.SetEntry(offset_sigma + 3*k+1, nb_dof_elt+i,
                                             -poids*(c0101*grad_phi(k)(0) + c0111*grad_phi(k)(1)));
			
			mat_interac.SetEntry(offset_sigma + 3*k+2, i,
                                             -poids*(c0011*grad_phi(k)(0) + c0111*grad_phi(k)(1)));
			mat_interac.SetEntry(offset_sigma + 3*k+2, nb_dof_elt+i,
                                             -poids*(c0111*grad_phi(k)(0) + c1111*grad_phi(k)(1)));
		      }
		  }
		
		// mass matrix for u
		coef = iomega*jacobian*nat_mat.GetCoefMass();              
		rho = ref_rho(ref_domain).GetCoefficient(var_problem, iquad, k);
		coef_damping = jacobian*nat_mat.GetCoefDamping();              
		rho_damping = ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, k);
		
		if (pml)
		  {
		    tau = var_boundary.GetTauPML(i1, k);
		    tau *= jacobian*nat_mat.GetCoefDamping();
		    
		    // part used for damping matrix of u
		    fevalTauX(k) = 0.5*(tau(0) + tau(1))*rho*val_phi(k);
		    fevalTauY(k) = 0.5*(tau(0) - tau(1))*rho*val_phi(k);
		  }
		
		// for u
		coef *= rho; coef_damping *= rho_damping;
		feval(k) = (coef + coef_damping)*val_phi(k);
	      }
	    else
	      {
		feval(k) = 0.0;
		fevalTauX(k) = 0.0;
		fevalTauY(k) = 0.0;
	      }
	  }
	
	Fb.ApplyCh(feval, contrib);
        
	for (int j = 0; j < nb_dof_elt; j++)
          if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
            {
              mat_interac.SetEntry(i, j, contrib(j));
              mat_interac.SetEntry(nb_dof_elt+i, nb_dof_elt+j, contrib(j));
            }
	
	if (pml)
	  {
	    for (int j = 0; j < nb_dof_elt; j++)
              if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                {
                  mat_interac.SetEntry(2*nb_dof_elt+i, 2*nb_dof_elt+j, contrib(j));
                  mat_interac.SetEntry(3*nb_dof_elt+i, 3*nb_dof_elt+j, contrib(j));
                }
	    
	    Fb.ApplyCh(fevalTauX, contrib);
	    for (int j = 0; j < nb_dof_elt; j++)
              if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                {
                  mat_interac.AddInteraction(i, j, contrib(j));
                  mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt+j, contrib(j));
                  mat_interac.AddInteraction(2*nb_dof_elt+i, 2*nb_dof_elt+j, contrib(j));
                  mat_interac.AddInteraction(3*nb_dof_elt+i, 3*nb_dof_elt+j, contrib(j));
                }
	    
	    Fb.ApplyCh(fevalTauY, contrib);
	    for (int j = 0; j < nb_dof_elt; j++)
              if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                {
                  mat_interac.AddInteraction(i, 2*nb_dof_elt+j, contrib(j));
                  mat_interac.AddInteraction(2*nb_dof_elt+i, j, contrib(j));
                  mat_interac.AddInteraction(nb_dof_elt+i, 3*nb_dof_elt+j, contrib(j));
                  mat_interac.AddInteraction(3*nb_dof_elt+i, nb_dof_elt+j, contrib(j));
                }
	  }
      }    
  }
  

#ifdef MONTJOIE_WITH_THREE_DIM
  //! computation of elementary matrix for elastodynamic equation
  /*!
    \param[in] iquad element number
    \param[out] num_dof numbers of degrees of freedom for the element
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] var_problem considered problem
    \param[in] Fb finite element associated with the element
  */
  template<class T0, class Complexe>
  void VarElastic_Base
  ::ComputeElementaryMatrix(const Vector<ElasticPhysicalIndice<Dimension3, 3, T0> >& ref_tensorC,
			    const Vector<ScalarPhysicalIndice<Dimension3, T0> >& ref_rho,
			    const Vector<ScalarPhysicalIndice<Dimension3, T0> >& ref_sigma,
			    int iquad, IVect& num_dof, VirtualMatrix<Complexe>& mat_interac,
			    const GlobalGenericMatrix<Complexe>& nat_mat,
			    const VarProblem<Dimension3>& var_problem,
			    const VarBoundaryCondition_Dim<T0, Dimension3>& var_boundary,
			    const ElementReference<Dimension3, 1>& Fb)
  {
    // for first-order (in time) formulation
    // specific expression with split pml
    
    // without PML, the considered equations are (in 2-D) :
    // du_x/dt - d/dx(sigma_xx) - d/dy(sigma_xy) = 0
    // du_y/dt - d/dx(sigma_xy) - d/dy(sigma_yy) = 0
    // dsigma/dt - C nabla u = 0
    
    // it should be noticed that the symmetry of sigma is exploited
    // such that only sigma_xx, sigma_xy, sigma_yy are computed in 2-D
    // sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz in 3-D
    
    // With PML, u_x, u_y, sigma_xx, sigma_xy, etc are splitted in d unknowns
    //      (u_x^1, u_x^2, u_x^3 for u_x)
    // In order to use continuous approximations for u, we use as unknowns :
    // u = u_x^1 + u_x^2 + u_x^3
    // u^* = u_x^1 - u_x^2
    // u^d = u_x^1 - u_x^3
    // as it is done for acoustics equations
    // For sigma, since they are discontinuous, we keep sigma_xx^1, sigma_xx^2, etc

    // number of integration points
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    
    // number of degrees of freedom
    int nb_dof_elt = Fb.GetNbDof();
    int nb_dof_u = 3*nb_dof_elt;
    int nb_comp = 6;
    int nb_dof_sigma = nb_comp*nb_points_quad;
    int nb_dof_all = nb_dof_u + nb_dof_sigma;
    int offset_sigma = nb_dof_u;
    int offset_sigma2(0), offset_sigma3(0);
    
    if (var_problem.InsidePML(iquad))
      {
        nb_dof_all += 2*(nb_dof_u + nb_dof_sigma);
        offset_sigma += 2*nb_dof_u;
        offset_sigma2 = offset_sigma + nb_dof_sigma;
        offset_sigma3 = offset_sigma2 + nb_dof_sigma;
      }
    
    // meshes
    const Mesh<Dimension3>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = var_problem.offset_dof_unknown(1);
    int Nvol_pml = mesh_num.GetNbDofPML();
    int offset_sigma_glob = Nvol*3;
    offset_sigma_glob += var_problem.GetOffsetDofV(iquad);
    
    // reference of the domain
    int ref_domain = var_problem.mesh.Element(iquad).GetReference();
    
    // dof numbers
    num_dof.Reallocate(nb_dof_all); num_dof.Fill(-1);
    IVect Nodle = var_problem.GetDofNumberOnElement(iquad);
    for (int j = 0; j < nb_dof_elt; j++)
      for (int m = 0; m < 3; m++)
        num_dof(j+m*nb_dof_elt) = Nodle(j) + m*Nvol;
    
    for (int j = 0; j < nb_points_quad; j++)
      for (int m = 0; m < nb_comp; m++)
        num_dof(offset_sigma + j*nb_comp + m) = offset_sigma_glob + j*nb_comp + m; 
    
    int i1 = -1;
    // additional dofs for pml
    if (var_problem.InsidePML(iquad))
      {
        // additional dofs for u
        for (int j = 0; j < nb_dof_elt; j++)
          {
            int num_pml = mesh_num.GetDofPML(Nodle(j));
            for (int m = 0; m < 3; m++)
              {
                num_dof(nb_dof_u + m*nb_dof_elt + j) = Nvol*m + mesh_num.GetNbDof() + num_pml;
		num_dof(nb_dof_u + (m+3)*nb_dof_elt + j)
		  = Nvol*m + Nvol_pml + mesh_num.GetNbDof() + num_pml;
              }
          }
        
        // additional dofs for sigma
        for (int j = 0; j < nb_points_quad; j++)
          {
            for (int m = 0; m < nb_comp; m++)
              num_dof(offset_sigma2 + j*nb_comp + m)
                = offset_sigma_glob + (nb_points_quad+j)*nb_comp + m;
            
	    for (int m = 0; m < nb_comp; m++)
	      num_dof(offset_sigma3 + j*nb_comp + m)
		= offset_sigma_glob + (2*nb_points_quad+j)*nb_comp + m;
          }
        
        i1 = iquad - mesh.GetNbElt() + var_boundary.GetNbEltPML();
      }
    
    // initialisation of the matrix
    mat_interac.Clear();
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();

    // if true, stiffness matrix must be computed
    bool stiff = false;
    if (nat_mat.GetCoefStiffness() != Real_wp(0))
      stiff = true;

    bool affine = var_problem.mesh.IsElementAffine(iquad);
    bool pml = var_problem.InsidePML(iquad);

    // loop over dofs 
    Real_wp jacobian; Complexe poids;
    VectReal_wp ones(nb_dof_elt), val_phi(nb_points_quad), grad_chap_phi(3*nb_points_quad);
    
    Vector<Complexe> feval(nb_points_quad), fevalTauX(nb_points_quad), fevalTauY(nb_points_quad),
      fevalTauZ(nb_points_quad), contrib(nb_dof_elt);
    
    Vector<TinyVector<Real_wp, 3> > grad_phi(nb_points_quad); 
    TinyVector<Real_wp, 3> vec_u;
    TinyVector<Complexe, 3> tau;
    val_phi.Fill(0); grad_chap_phi.Fill(0);
    feval.Fill(0); contrib.Fill(0);
    
    T0 iomega, c0000, c0001, c0011,
      c0101, c0111, c1111, rho, rho_damping;
    
    T0 c0002, c0012, c0022, c0102, c0112,
      c0122, c0202, c0211, c0212, c0222, c1112, c1122, c1212, c1222, c2222;
    
    Complexe coef, coef_damping;
    TinyMatrix<TinyMatrix<T0, General, 3, 3>, General, 3, 3> C;
    
    var_problem.GetMiomega(iomega);
    
    // mass matrix for sigma
    for (int k = 0; k < nb_points_quad; k++)
      {
	if (affine)
	  jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
	else
	  jacobian = var_problem.Glob_jacobian(iquad)(k);
	
	coef = iomega*jacobian*nat_mat.GetCoefMass();              
	mat_interac.SetEntry(offset_sigma + 6*k, offset_sigma + 6*k, coef);
	mat_interac.SetEntry(offset_sigma + 6*k+1, offset_sigma + 6*k+1, coef);
	mat_interac.SetEntry(offset_sigma + 6*k+2, offset_sigma + 6*k+2, coef);
	mat_interac.SetEntry(offset_sigma + 6*k+3, offset_sigma + 6*k+3, coef);
        mat_interac.SetEntry(offset_sigma + 6*k+4, offset_sigma + 6*k+4, coef);
	mat_interac.SetEntry(offset_sigma + 6*k+5, offset_sigma + 6*k+5, coef);
        
	if (pml)
	  {
	    tau = var_boundary.GetTauPML(i1, k);
	    tau *= jacobian*nat_mat.GetCoefDamping();
	    
	    mat_interac.AddInteraction(offset_sigma + 6*k, offset_sigma + 6*k, tau(0));
	    mat_interac.AddInteraction(offset_sigma + 6*k+1, offset_sigma + 6*k+1, tau(0));
	    mat_interac.AddInteraction(offset_sigma + 6*k+2, offset_sigma + 6*k+2, tau(0)); 
	    mat_interac.AddInteraction(offset_sigma + 6*k+3, offset_sigma + 6*k+3, tau(0)); 
	    mat_interac.AddInteraction(offset_sigma + 6*k+4, offset_sigma + 6*k+4, tau(0)); 
	    mat_interac.AddInteraction(offset_sigma + 6*k+5, offset_sigma + 6*k+5, tau(0)); 
            
	    mat_interac.SetEntry(offset_sigma2 + 6*k, offset_sigma2 + 6*k, coef + tau(1));
	    mat_interac.SetEntry(offset_sigma2 + 6*k+1, offset_sigma2 + 6*k+1, coef + tau(1));
	    mat_interac.SetEntry(offset_sigma2 + 6*k+2, offset_sigma2 + 6*k+2, coef + tau(1)); 
	    mat_interac.SetEntry(offset_sigma2 + 6*k+3, offset_sigma2 + 6*k+3, coef + tau(1)); 
	    mat_interac.SetEntry(offset_sigma2 + 6*k+4, offset_sigma2 + 6*k+4, coef + tau(1)); 
	    mat_interac.SetEntry(offset_sigma2 + 6*k+5, offset_sigma2 + 6*k+5, coef + tau(1)); 
            
	    mat_interac.SetEntry(offset_sigma3 + 6*k, offset_sigma3 + 6*k, coef + tau(2));
	    mat_interac.SetEntry(offset_sigma3 + 6*k+1, offset_sigma3 + 6*k+1, coef + tau(2));
	    mat_interac.SetEntry(offset_sigma3 + 6*k+2, offset_sigma3 + 6*k+2, coef + tau(2)); 
	    mat_interac.SetEntry(offset_sigma3 + 6*k+3, offset_sigma3 + 6*k+3, coef + tau(2)); 
	    mat_interac.SetEntry(offset_sigma3 + 6*k+4, offset_sigma3 + 6*k+4, coef + tau(2)); 
	    mat_interac.SetEntry(offset_sigma3 + 6*k+5, offset_sigma3 + 6*k+5, coef + tau(2)); 
	  }
      }
    
    if (Fb.LumpedMassMatrix() && (!stiff))
      {
        // condensation de masse, on evite les doubles boucles
        for (int k = 0; k < nb_dof_elt; k++)
          {
            if (affine)
              jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
            else
              jacobian = var_problem.Glob_jacobian(iquad)(k);
            
            int i = k, j = k;
            coef = iomega*jacobian*nat_mat.GetCoefMass();              
            rho = ref_rho(ref_domain).GetCoefficient(var_problem, iquad, k);
            coef_damping = jacobian*nat_mat.GetCoefDamping();
            rho_damping = ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, k);
            
            coef *= rho; coef_damping *= rho_damping;
            Complexe vloc = coef+coef_damping;
            mat_interac.SetEntry(k, k, vloc);
            mat_interac.SetEntry(nb_dof_elt + k, nb_dof_elt + k, vloc);
            mat_interac.SetEntry(2*nb_dof_elt + k, 2*nb_dof_elt + k, vloc);
            
            if (pml)
              {
                tau = var_boundary.GetTauPML(i1, k);
                tau *= jacobian*nat_mat.GetCoefDamping();
                
                Complexe coef_x = tau(0)/Real_wp(3)*rho;
                Complexe coef_y = tau(1)/Real_wp(3)*rho;
		Complexe coef_z = tau(2)/Real_wp(3)*rho;
                
                if (abs(vloc) > var_problem.GetThresholdMatrix())
                  {
                    mat_interac.SetEntry(3*nb_dof_elt+k, 3*nb_dof_elt+k, vloc);
                    mat_interac.SetEntry(4*nb_dof_elt+k, 4*nb_dof_elt+k, vloc);
                    mat_interac.SetEntry(5*nb_dof_elt+k, 5*nb_dof_elt+k, vloc);
                    mat_interac.SetEntry(6*nb_dof_elt+k, 6*nb_dof_elt+k, vloc);
                    mat_interac.SetEntry(7*nb_dof_elt+k, 7*nb_dof_elt+k, vloc);
                    mat_interac.SetEntry(8*nb_dof_elt+k, 8*nb_dof_elt+k, vloc);
                  }
                
                if (abs(coef_x) > var_problem.GetThresholdMatrix())
                  {
                    mat_interac.AddInteraction(i, j, coef_x);
                    mat_interac.AddInteraction(i, 3*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(i, 6*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(3*nb_dof_elt + i, j, coef_x);
                    mat_interac.AddInteraction(3*nb_dof_elt + i, 3*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(3*nb_dof_elt + i, 6*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(6*nb_dof_elt + i, j, coef_x);
                    mat_interac.AddInteraction(6*nb_dof_elt + i, 3*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(6*nb_dof_elt + i, 6*nb_dof_elt + j, coef_x);
                    
                    mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(nb_dof_elt+i, 4*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(nb_dof_elt+i, 7*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(4*nb_dof_elt + i, nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(4*nb_dof_elt + i, 4*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(4*nb_dof_elt + i, 7*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(7*nb_dof_elt + i, nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(7*nb_dof_elt + i, 4*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(7*nb_dof_elt + i, 7*nb_dof_elt + j, coef_x);
                    
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 2*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 5*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 8*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(5*nb_dof_elt + i, 2*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(5*nb_dof_elt + i, 5*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(5*nb_dof_elt + i, 8*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(8*nb_dof_elt + i, 2*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(8*nb_dof_elt + i, 5*nb_dof_elt + j, coef_x);
                    mat_interac.AddInteraction(8*nb_dof_elt + i, 8*nb_dof_elt + j, coef_x);
                  }

                if (abs(coef_y) > var_problem.GetThresholdMatrix())
                  {
                    mat_interac.AddInteraction(i, j, coef_y);
                    mat_interac.AddInteraction(i, 3*nb_dof_elt + j, -2.0*coef_y);
                    mat_interac.AddInteraction(i, 6*nb_dof_elt + j, coef_y);
                    mat_interac.AddInteraction(3*nb_dof_elt + i, j, -coef_y);
                    mat_interac.AddInteraction(3*nb_dof_elt + i, 3*nb_dof_elt + j, 2.0*coef_y);
                    mat_interac.AddInteraction(3*nb_dof_elt + i, 6*nb_dof_elt + j, -coef_y);
                    
                    mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt + j, coef_y);
                    mat_interac.AddInteraction(nb_dof_elt+i, 4*nb_dof_elt + j, -2.0*coef_y);
                    mat_interac.AddInteraction(nb_dof_elt+i, 7*nb_dof_elt + j, coef_y);
                    mat_interac.AddInteraction(4*nb_dof_elt + i, nb_dof_elt + j, -coef_y);
                    mat_interac.AddInteraction(4*nb_dof_elt + i, 4*nb_dof_elt + j, 2.0*coef_y);
                    mat_interac.AddInteraction(4*nb_dof_elt + i, 7*nb_dof_elt + j, -coef_y);
                    
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 2*nb_dof_elt + j, coef_y);
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 5*nb_dof_elt + j, -2.0*coef_y);
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 8*nb_dof_elt + j, coef_y);
                    mat_interac.AddInteraction(5*nb_dof_elt + i, 2*nb_dof_elt + j, -coef_y);
                    mat_interac.AddInteraction(5*nb_dof_elt + i, 5*nb_dof_elt + j, 2.0*coef_y);
                    mat_interac.AddInteraction(5*nb_dof_elt + i, 8*nb_dof_elt + j, -coef_y);
                  }
                
                if (abs(coef_z) > var_problem.GetThresholdMatrix())
                  {
                    mat_interac.AddInteraction(i, j, coef_z);
                    mat_interac.AddInteraction(i, 3*nb_dof_elt + j, coef_z);
                    mat_interac.AddInteraction(i, 6*nb_dof_elt + j, -2.0*coef_z);
                    mat_interac.AddInteraction(6*nb_dof_elt + i, j, -coef_z);
                    mat_interac.AddInteraction(6*nb_dof_elt + i, 3*nb_dof_elt + j, -coef_z);
                    mat_interac.AddInteraction(6*nb_dof_elt + i, 6*nb_dof_elt + j, 2.0*coef_z);
                    
                    mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt + j, coef_z);
                    mat_interac.AddInteraction(nb_dof_elt+i, 4*nb_dof_elt + j, coef_z);
                    mat_interac.AddInteraction(nb_dof_elt+i, 7*nb_dof_elt + j, -2.0*coef_z);
                    mat_interac.AddInteraction(7*nb_dof_elt + i, nb_dof_elt + j, -coef_z);
                    mat_interac.AddInteraction(7*nb_dof_elt + i, 4*nb_dof_elt + j, -coef_z);
                    mat_interac.AddInteraction(7*nb_dof_elt + i, 7*nb_dof_elt + j, 2.0*coef_z);
                    
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 2*nb_dof_elt + j, coef_z);
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 5*nb_dof_elt + j, coef_z);
                    mat_interac.AddInteraction(2*nb_dof_elt + i, 8*nb_dof_elt + j, -2.0*coef_z);
                    mat_interac.AddInteraction(8*nb_dof_elt + i, 2*nb_dof_elt + j, -coef_z);
                    mat_interac.AddInteraction(8*nb_dof_elt + i, 5*nb_dof_elt + j, -coef_z);
                    mat_interac.AddInteraction(8*nb_dof_elt + i, 8*nb_dof_elt + j, 2.0*coef_z);
                  }            
              }
          } 
        
        return;
      }

    // stiffness matrix and mass matrix for u
    Real_wp max_dfjm1 = 1.0;
    for (int i = 0; i < nb_dof_elt; i++)
      {
        ones.Fill(0); ones(i) = 1.0;
        // computation of phi_i(xi_k) and grad phi_i(xi_k)
        Fb.ApplyChTranspose(ones, val_phi);
        if (stiff)
          {
            if (Fb.UseQuadraturePointsForRh())
              Fb.ApplyRhQuadratureTranspose(val_phi, grad_chap_phi);
            else
              Fb.ApplyRhTranspose(ones, grad_chap_phi);
            
            if (affine)
              {
                Matrix3_3 dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
                max_dfjm1 = MaxAbs(dfjm1);
                for (int k = 0; k < nb_points_quad; k++)
                  {
                    CopyVector(grad_chap_phi, k, vec_u);
                    MltTrans(dfjm1, vec_u, grad_phi(k));
                  }
              }
            else
              {
                max_dfjm1 = MaxAbs(var_problem.Glob_DFjm1(iquad)(0));
                for (int k = 0; k < nb_points_quad; k++)
                  {
                    CopyVector(grad_chap_phi, k, vec_u);
                    MltTrans(var_problem.Glob_DFjm1(iquad)(k), vec_u, grad_phi(k));
                  }
              }
          }
	
	// 3-D case
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    if ((abs(val_phi(k)) > epsilon_machine)
		|| (Norm2(grad_phi(k)) > max_dfjm1*epsilon_machine))
	      {
		if (affine)
		  jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
		else
		  jacobian = var_problem.Glob_jacobian(iquad)(k);
		
		if (stiff)
		  {
		    poids = Fb.WeightsND(k)*nat_mat.GetCoefStiffness();
		    ref_tensorC(ref_domain).GetCoefficient(var_problem, iquad, k, C);
		    
		    c0000 = C(0, 0)(0, 0);
		    c0001 = C(0, 0)(0, 1);
		    c0002 = C(0, 0)(0, 2);
		    c0011 = C(0, 1)(0, 1);
		    c0012 = C(0, 1)(0, 2);
		    c0022 = C(0, 2)(0, 2);
		    c0101 = C(0, 0)(1, 1);
		    c0102 = C(0, 0)(1, 2);
		    c0111 = C(0, 1)(1, 1);
		    c0112 = C(0, 1)(1, 2);
		    c0122 = C(0, 2)(1, 2);
		    c0202 = C(0, 0)(2, 2);
		    c0211 = C(0, 1)(2, 1);
		    c0212 = C(0, 1)(2, 2);
		    c0222 = C(0, 2)(2, 2);
		    c1111 = C(1, 1)(1, 1);
		    c1112 = C(1, 1)(1, 2);
		    c1122 = C(1, 2)(1, 2);
		    c1212 = C(1, 1)(2, 2);
		    c1222 = C(1, 2)(2, 2);
		    c2222 = C(2, 2)(2, 2);
		    
		    // part div(sigma)
		    mat_interac.SetEntry(i, offset_sigma + 6*k, poids*grad_phi(k)(0));
		    mat_interac.SetEntry(i, offset_sigma + 6*k+1, poids*grad_phi(k)(1));
		    mat_interac.SetEntry(i, offset_sigma + 6*k+2, poids*grad_phi(k)(2));
                    
		    mat_interac.SetEntry(nb_dof_elt+i, offset_sigma + 6*k+1, poids*grad_phi(k)(0));
		    mat_interac.SetEntry(nb_dof_elt+i, offset_sigma + 6*k+3, poids*grad_phi(k)(1));
		    mat_interac.SetEntry(nb_dof_elt+i, offset_sigma + 6*k+4, poids*grad_phi(k)(2));
                    
		    mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma + 6*k+2, poids*grad_phi(k)(0));
		    mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma + 6*k+4, poids*grad_phi(k)(1));
		    mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma + 6*k+5, poids*grad_phi(k)(2));
                    
		    if (pml)
		      {
			// part div(sigma^2) and div(sigma^3)
			mat_interac.SetEntry(i, offset_sigma2 + 6*k, poids*grad_phi(k)(0));
			mat_interac.SetEntry(i, offset_sigma2 + 6*k+1, poids*grad_phi(k)(1));
			mat_interac.SetEntry(i, offset_sigma2 + 6*k+2, poids*grad_phi(k)(2));
                        
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma2 + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma2 + 6*k+3,
                                             poids*grad_phi(k)(1));
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma2 + 6*k+4,
                                             poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma2 + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma2 + 6*k+4,
                                             poids*grad_phi(k)(1));
			mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma2 + 6*k+5,
                                             poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(i, offset_sigma3 + 6*k, poids*grad_phi(k)(0));
			mat_interac.SetEntry(i, offset_sigma3 + 6*k+1, poids*grad_phi(k)(1));
			mat_interac.SetEntry(i, offset_sigma3 + 6*k+2, poids*grad_phi(k)(2));
                        
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma3 + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma3 + 6*k+3,
                                             poids*grad_phi(k)(1));
			mat_interac.SetEntry(nb_dof_elt+i, offset_sigma3 + 6*k+4,
                                             poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma3 + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma3 + 6*k+4,
                                             poids*grad_phi(k)(1));
			mat_interac.SetEntry(2*nb_dof_elt+i, offset_sigma3 + 6*k+5,
                                             poids*grad_phi(k)(2));
			
			// part for u^* (u_1 - u_2)
			mat_interac.SetEntry(3*nb_dof_elt + i, offset_sigma + 6*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(3*nb_dof_elt + i, offset_sigma + 6*k+1,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(4*nb_dof_elt+i, offset_sigma + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(4*nb_dof_elt+i, offset_sigma + 6*k+3,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(5*nb_dof_elt + i, offset_sigma + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(5*nb_dof_elt + i, offset_sigma + 6*k+4,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(3*nb_dof_elt + i, offset_sigma2 + 6*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(3*nb_dof_elt + i, offset_sigma2 + 6*k+1,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(4*nb_dof_elt+i, offset_sigma2 + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(4*nb_dof_elt+i, offset_sigma2 + 6*k+3,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(5*nb_dof_elt + i, offset_sigma2 + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(5*nb_dof_elt + i, offset_sigma2 + 6*k+4,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(3*nb_dof_elt + i, offset_sigma3 + 6*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(3*nb_dof_elt + i, offset_sigma3 + 6*k+1,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(4*nb_dof_elt+i, offset_sigma3 + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(4*nb_dof_elt+i, offset_sigma3 + 6*k+3,
                                             -poids*grad_phi(k)(1));
			
			mat_interac.SetEntry(5*nb_dof_elt + i, offset_sigma3 + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(5*nb_dof_elt + i, offset_sigma3 + 6*k+4,
                                             -poids*grad_phi(k)(1));
			
			// part for u^d (u_1 - u_3)
			mat_interac.SetEntry(6*nb_dof_elt + i, offset_sigma + 6*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(6*nb_dof_elt + i, offset_sigma + 6*k+2,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(7*nb_dof_elt+i, offset_sigma + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(7*nb_dof_elt+i, offset_sigma + 6*k+4,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(8*nb_dof_elt + i, offset_sigma + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(8*nb_dof_elt + i, offset_sigma + 6*k+5,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(6*nb_dof_elt + i, offset_sigma2 + 6*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(6*nb_dof_elt + i, offset_sigma2 + 6*k+2,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(7*nb_dof_elt+i, offset_sigma2 + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(7*nb_dof_elt+i, offset_sigma2 + 6*k+4,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(8*nb_dof_elt + i, offset_sigma2 + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(8*nb_dof_elt + i, offset_sigma2 + 6*k+5,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(6*nb_dof_elt + i, offset_sigma3 + 6*k,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(6*nb_dof_elt + i, offset_sigma3 + 6*k+2,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(7*nb_dof_elt+i, offset_sigma3 + 6*k+1,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(7*nb_dof_elt+i, offset_sigma3 + 6*k+4,
                                             -poids*grad_phi(k)(2));
			
			mat_interac.SetEntry(8*nb_dof_elt + i, offset_sigma3 + 6*k+2,
                                             poids*grad_phi(k)(0));
			mat_interac.SetEntry(8*nb_dof_elt + i, offset_sigma3 + 6*k+5,
                                             -poids*grad_phi(k)(2));
			
			// part for sigma_xx^1, sigma_xx^2, sigma_xx^3
			mat_interac.SetEntry(offset_sigma + 6*k, i, -poids*c0000*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k, nb_dof_elt+i,
                                             -poids*c0001*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k, 2*nb_dof_elt+i,
                                             -poids*c0002*grad_phi(k)(0));
			
			mat_interac.SetEntry(offset_sigma2 + 6*k, i, -poids*c0001*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k, nb_dof_elt+i,
                                             -poids*c0011*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k, 2*nb_dof_elt+i,
                                             -poids*c0012*grad_phi(k)(1));
			
			mat_interac.SetEntry(offset_sigma3 + 6*k, i, -poids*c0002*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k, nb_dof_elt+i,
                                             -poids*c0012*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k, 2*nb_dof_elt+i,
                                             -poids*c0022*grad_phi(k)(2));
			
			// part for sigma_xy^1, sigma_xy^2, sigma_xy^3
			mat_interac.SetEntry(offset_sigma + 6*k+1, i, -poids*c0001*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+1, nb_dof_elt+i,
                                             -poids*c0101*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+1, 2*nb_dof_elt+i,
                                             -poids*c0102*grad_phi(k)(0));
                        
			mat_interac.SetEntry(offset_sigma2 + 6*k+1, i, -poids*c0101*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+1, nb_dof_elt+i,
                                             -poids*c0111*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+1, 2*nb_dof_elt+i,
                                             -poids*c0112*grad_phi(k)(1));
			
			mat_interac.SetEntry(offset_sigma3 + 6*k+1, i, -poids*c0102*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+1, nb_dof_elt+i,
                                             -poids*c0112*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+1, 2*nb_dof_elt+i,
                                             -poids*c0122*grad_phi(k)(2));
			
			// part for sigma_xz^1, sigma_xz^2, sigma_xz^3
			mat_interac.SetEntry(offset_sigma + 6*k+2, i, -poids*c0002*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+2, nb_dof_elt+i,
                                             -poids*c0102*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+2, 2*nb_dof_elt+i,
                                             -poids*c0202*grad_phi(k)(0));
			
			mat_interac.SetEntry(offset_sigma2 + 6*k+2, i, -poids*c0102*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+2, nb_dof_elt+i,
                                             -poids*c0211*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+2, 2*nb_dof_elt+i,
                                             -poids*c0212*grad_phi(k)(1));
			
			mat_interac.SetEntry(offset_sigma3 + 6*k+2, i, -poids*c0202*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+2, nb_dof_elt+i,
                                             -poids*c0212*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+2, 2*nb_dof_elt+i,
                                             -poids*c0222*grad_phi(k)(2));
			
			// part for sigma_yy^1, sigma_yy^2, sigma_yy^3
			mat_interac.SetEntry(offset_sigma + 6*k+3, i, -poids*c0011*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+3, nb_dof_elt+i,
                                             -poids*c0111*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+3, 2*nb_dof_elt+i,
                                             -poids*c0211*grad_phi(k)(0));
			
			mat_interac.SetEntry(offset_sigma2 + 6*k+3, i, -poids*c0111*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+3, nb_dof_elt+i,
                                             -poids*c1111*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+3, 2*nb_dof_elt+i,
                                             -poids*c1112*grad_phi(k)(1));
			
			mat_interac.SetEntry(offset_sigma3 + 6*k+3, i, -poids*c0211*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+3, nb_dof_elt+i,
                                             -poids*c1112*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+3, 2*nb_dof_elt+i,
                                             -poids*c1122*grad_phi(k)(2));
                        
			// part for sigma_yz^1, sigma_yz^2, sigma_yz^3
			mat_interac.SetEntry(offset_sigma + 6*k+4, i, -poids*c0012*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+4, nb_dof_elt+i,
                                             -poids*c0112*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+4, 2*nb_dof_elt+i,
                                             -poids*c0212*grad_phi(k)(0));
			
			mat_interac.SetEntry(offset_sigma2 + 6*k+4, i, -poids*c0112*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+4, nb_dof_elt+i,
                                             -poids*c1112*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+4, 2*nb_dof_elt+i,
                                             -poids*c1212*grad_phi(k)(1));
			
			mat_interac.SetEntry(offset_sigma3 + 6*k+4, i, -poids*c0212*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+4, nb_dof_elt+i,
                                             -poids*c1212*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+4, 2*nb_dof_elt+i,
                                             -poids*c1222*grad_phi(k)(2));
			
			// part for sigma_zz^1, sigma_zz^2, sigma_zz^3
			mat_interac.SetEntry(offset_sigma + 6*k+5, i, -poids*c0022*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+5, nb_dof_elt+i,
                                             -poids*c0122*grad_phi(k)(0));
			mat_interac.SetEntry(offset_sigma + 6*k+5, 2*nb_dof_elt+i,
                                             -poids*c0222*grad_phi(k)(0));
			
			mat_interac.SetEntry(offset_sigma2 + 6*k+5, i, -poids*c0122*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+5, nb_dof_elt+i, 
                                             -poids*c1122*grad_phi(k)(1));
			mat_interac.SetEntry(offset_sigma2 + 6*k+5, 2*nb_dof_elt+i,
                                             -poids*c1222*grad_phi(k)(1));
			
			mat_interac.SetEntry(offset_sigma3 + 6*k+5, i, -poids*c0222*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+5, nb_dof_elt+i,
                                             -poids*c1222*grad_phi(k)(2));
			mat_interac.SetEntry(offset_sigma3 + 6*k+5, 2*nb_dof_elt+i,
                                             -poids*c2222*grad_phi(k)(2));         
		      }
		    else
		      {
			// part C grad u
			mat_interac.SetEntry(offset_sigma + 6*k, i,
			  -poids*(c0000*grad_phi(k)(0) + c0001*grad_phi(k)(1)
                                  + c0002*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k, nb_dof_elt+i,
                                             -poids*(c0001*grad_phi(k)(0) + c0011*grad_phi(k)(1)
                                                     + c0012*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k, 2*nb_dof_elt+i,
                                             -poids*(c0002*grad_phi(k)(0) + c0012*grad_phi(k)(1)
                                                     + c0022*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+1, i,
                                             -poids*(c0001*grad_phi(k)(0)
                                                     + c0101*grad_phi(k)(1) + c0102*grad_phi(k)(2)));
                        
			mat_interac.SetEntry(offset_sigma + 6*k+1, nb_dof_elt+i,
                                             -poids*(c0101*grad_phi(k)(0) + c0111*grad_phi(k)(1)
                                                     + c0112*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+1, 2*nb_dof_elt+i,
                                             -poids*(c0102*grad_phi(k)(0) + c0112*grad_phi(k)(1)
                                                     + c0122*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+2, i,
                                             -poids*(c0002*grad_phi(k)(0) + c0102*grad_phi(k)(1)
                                                     + c0202*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+2, nb_dof_elt+i,
                                             -poids*(c0102*grad_phi(k)(0) + c0211*grad_phi(k)(1)
                                                     + c0212*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+2, 2*nb_dof_elt+i,
                                             -poids*(c0202*grad_phi(k)(0) + c0212*grad_phi(k)(1)
                                                     + c0222*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+3, i,
                                             -poids*(c0011*grad_phi(k)(0) + c0111*grad_phi(k)(1)
                                                     + c0211*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+3, nb_dof_elt+i,
                                             -poids*(c0111*grad_phi(k)(0) + c1111*grad_phi(k)(1)
                                                     + c1112*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+3, 2*nb_dof_elt+i,
                                             -poids*(c0211*grad_phi(k)(0) + c1112*grad_phi(k)(1)
                                                     + c1122*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+4, i,
                                             -poids*(c0012*grad_phi(k)(0) + c0112*grad_phi(k)(1)
                                                     + c0212*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+4, nb_dof_elt+i,
                                             -poids*(c0112*grad_phi(k)(0) + c1112*grad_phi(k)(1)
                                                     + c1212*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+4, 2*nb_dof_elt+i,
                                             -poids*(c0212*grad_phi(k)(0) + c1212*grad_phi(k)(1)
                                                     + c1222*grad_phi(k)(2)));
                        
			mat_interac.SetEntry(offset_sigma + 6*k+5, i,
                                             -poids*(c0022*grad_phi(k)(0) + c0122*grad_phi(k)(1)
                                                     + c0222*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+5, nb_dof_elt+i,
                                             -poids*(c0122*grad_phi(k)(0) + c1122*grad_phi(k)(1)
                                                     + c1222*grad_phi(k)(2)));
			
			mat_interac.SetEntry(offset_sigma + 6*k+5, 2*nb_dof_elt+i,
                                             -poids*(c0222*grad_phi(k)(0) + c1222*grad_phi(k)(1)
                                                     + c2222*grad_phi(k)(2)));
		      }
		  }
		
		// mass matrix for sigma
		coef = iomega*jacobian*nat_mat.GetCoefMass();              
		rho = ref_rho(ref_domain).GetCoefficient(var_problem, iquad, k);
		coef_damping = jacobian*nat_mat.GetCoefDamping();
		rho_damping = ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, k);
		
		if (pml)
		  {
		    tau = var_boundary.GetTauPML(i1, k);
		    tau *= jacobian*nat_mat.GetCoefDamping();
                    
		    fevalTauX(k) = tau(0)/Real_wp(3)*rho*val_phi(k);
		    fevalTauY(k) = tau(1)/Real_wp(3)*rho*val_phi(k);
		    fevalTauZ(k) = tau(2)/Real_wp(3)*rho*val_phi(k);
		  }
		
		// for u
		coef *= rho; coef_damping *= rho_damping;
		feval(k) = (coef+coef_damping)*val_phi(k);
	      }
	    else
	      {
		feval(k) = 0.0;
		fevalTauX(k) = 0.0;
		fevalTauY(k) = 0.0;
		fevalTauZ(k) = 0.0;
	      }
	  }
	
	Fb.ApplyCh(feval, contrib);
        
	for (int j = 0; j < nb_dof_elt; j++)
          if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
            {
              mat_interac.SetEntry(i, j, contrib(j));
              mat_interac.SetEntry(nb_dof_elt+i, nb_dof_elt+j, contrib(j));
              mat_interac.SetEntry(2*nb_dof_elt+i, 2*nb_dof_elt+j, contrib(j));
            }
	
	if (pml)
	  {
	    for (int j = 0; j < nb_dof_elt; j++)
	      if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
		{
		  mat_interac.SetEntry(3*nb_dof_elt+i, 3*nb_dof_elt+j, contrib(j));
		  mat_interac.SetEntry(4*nb_dof_elt+i, 4*nb_dof_elt+j, contrib(j));
		  mat_interac.SetEntry(5*nb_dof_elt+i, 5*nb_dof_elt+j, contrib(j));
		  mat_interac.SetEntry(6*nb_dof_elt+i, 6*nb_dof_elt+j, contrib(j));
		  mat_interac.SetEntry(7*nb_dof_elt+i, 7*nb_dof_elt+j, contrib(j));
		  mat_interac.SetEntry(8*nb_dof_elt+i, 8*nb_dof_elt+j, contrib(j));
		}
	    
	    Fb.ApplyCh(fevalTauX, contrib);
	    for (int j = 0; j < nb_dof_elt; j++)
	      if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
		{
		  mat_interac.AddInteraction(i, j, contrib(j));
		  mat_interac.AddInteraction(i, 3*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(i, 6*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(3*nb_dof_elt + i, j, contrib(j));
		  mat_interac.AddInteraction(3*nb_dof_elt + i, 3*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(3*nb_dof_elt + i, 6*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(6*nb_dof_elt + i, j, contrib(j));
		  mat_interac.AddInteraction(6*nb_dof_elt + i, 3*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(6*nb_dof_elt + i, 6*nb_dof_elt + j, contrib(j));
                  
		  mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(nb_dof_elt+i, 4*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(nb_dof_elt+i, 7*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(4*nb_dof_elt + i, nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(4*nb_dof_elt + i, 4*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(4*nb_dof_elt + i, 7*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(7*nb_dof_elt + i, nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(7*nb_dof_elt + i, 4*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(7*nb_dof_elt + i, 7*nb_dof_elt + j, contrib(j));
		  
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 2*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 5*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 8*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(5*nb_dof_elt + i, 2*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(5*nb_dof_elt + i, 5*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(5*nb_dof_elt + i, 8*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(8*nb_dof_elt + i, 2*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(8*nb_dof_elt + i, 5*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(8*nb_dof_elt + i, 8*nb_dof_elt + j, contrib(j));
		}
	    
	    Fb.ApplyCh(fevalTauY, contrib);
	    for (int j = 0; j < nb_dof_elt; j++)
	      if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
		{
		  mat_interac.AddInteraction(i, j, contrib(j));
		  mat_interac.AddInteraction(i, 3*nb_dof_elt + j, -2.0*contrib(j));
		  mat_interac.AddInteraction(i, 6*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(3*nb_dof_elt + i, j, -contrib(j));
		  mat_interac.AddInteraction(3*nb_dof_elt + i, 3*nb_dof_elt + j, 2.0*contrib(j));
		  mat_interac.AddInteraction(3*nb_dof_elt + i, 6*nb_dof_elt + j, -contrib(j));
                  
		  mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(nb_dof_elt+i, 4*nb_dof_elt + j, -2.0*contrib(j));
		  mat_interac.AddInteraction(nb_dof_elt+i, 7*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(4*nb_dof_elt + i, nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(4*nb_dof_elt + i, 4*nb_dof_elt + j, 2.0*contrib(j));
		  mat_interac.AddInteraction(4*nb_dof_elt + i, 7*nb_dof_elt + j, -contrib(j));
                  
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 2*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 5*nb_dof_elt + j, -2.0*contrib(j));
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 8*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(5*nb_dof_elt + i, 2*nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(5*nb_dof_elt + i, 5*nb_dof_elt + j, 2.0*contrib(j));
		  mat_interac.AddInteraction(5*nb_dof_elt + i, 8*nb_dof_elt + j, -contrib(j));
		}
	    
	    Fb.ApplyCh(fevalTauZ, contrib);
	    for (int j = 0; j < nb_dof_elt; j++)
	      if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
		{
		  mat_interac.AddInteraction(i, j, contrib(j));
		  mat_interac.AddInteraction(i, 3*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(i, 6*nb_dof_elt + j, -2.0*contrib(j));
		  mat_interac.AddInteraction(6*nb_dof_elt + i, j, -contrib(j));
		  mat_interac.AddInteraction(6*nb_dof_elt + i, 3*nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(6*nb_dof_elt + i, 6*nb_dof_elt + j, 2.0*contrib(j));
                  
		  mat_interac.AddInteraction(nb_dof_elt+i, nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(nb_dof_elt+i, 4*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(nb_dof_elt+i, 7*nb_dof_elt + j, -2.0*contrib(j));
                  mat_interac.AddInteraction(7*nb_dof_elt + i, nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(7*nb_dof_elt + i, 4*nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(7*nb_dof_elt + i, 7*nb_dof_elt + j, 2.0*contrib(j));
		  
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 2*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 5*nb_dof_elt + j, contrib(j));
		  mat_interac.AddInteraction(2*nb_dof_elt + i, 8*nb_dof_elt + j, -2.0*contrib(j));
		  mat_interac.AddInteraction(8*nb_dof_elt + i, 2*nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(8*nb_dof_elt + i, 5*nb_dof_elt + j, -contrib(j));
		  mat_interac.AddInteraction(8*nb_dof_elt + i, 8*nb_dof_elt + j, 2.0*contrib(j));
		}
	  }
      }        
  }
#endif
  
  
  /******************
   * VarElastic_Dim *
   ******************/  
  
  
  //! allocation of arrays containing physical properties
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::InitIndices(int n)
  {
    ref_rho.Reallocate(n);
    ref_sigma.Reallocate(n);
    ref_tensorC.Reallocate(n);
    ref_tensorA.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
	ref_rho(i).SetIdentity();
	ref_sigma(i).Zero();
	ref_tensorC(i).SetIdentity();
	ref_tensorA(i).SetIdentity();
      }
  }
  

  //! returns the number of different physical media
  template<class Complexe ,class Dimension>
  int VarElastic_Dim<Complexe, Dimension>::GetNbPhysicalIndices() const
  {
    return ref_rho.GetM();
  }

  
  //! copies physical properties of an another given problem
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::
  CopyInputData(const VarElastic_Dim<Complexe, Dimension>& var)
  {
    ref_rho = var.ref_rho;
    ref_sigma = var.ref_sigma;
    ref_tensorC = var.ref_tensorC;
    ref_tensorA = var.ref_tensorA;
  }
  
  
  //! modification of physical indexes according to parameters of data file
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }
    
    ref_rho(i).SetInputData(nb, parameters, parameters(0));
    ref_tensorC(i).SetInputData(nb, parameters, parameters(0));
    ref_tensorA(i) = ref_tensorC(i);
    ref_tensorA(i).GetInverse();
    ref_sigma(i).SetInputData(nb, parameters, parameters(0));
  }
  

  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }
    
    if (name_media == "rho")
      ref_rho(i).SetInputData(nb, parameters, parameters(0));
    else if (name_media == "C")
      {
        ref_tensorC(i).SetInputData(nb, parameters, parameters(0));
        ref_tensorA(i) = ref_tensorC(i);
        ref_tensorA(i).GetInverse();
      }
    else if (name_media == "sigma")
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));    
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }


  //! returns the name associated with the physical index num
  template<class Complexe, class Dimension>
  string VarElastic_Dim<Complexe, Dimension>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("rho");
      case 1: return string("C");
      case 2: return string("sigma");
      }

    return string();
  }
  
  
  //! returns true if physical indexes are varying inside element i
  template<class Complexe, class Dimension>
  bool VarElastic_Dim<Complexe, Dimension>::IsVaryingMedia(int i) const
  {
    if (ref_rho(i).IsVarying()||ref_tensorC(i).IsVarying()||ref_sigma(i).IsVarying())
      return true;
    
    return false;
  }
  

  //! returns the maximal velocity of waves in a media
  template<class Complexe, class Dimension>
  Real_wp VarElastic_Dim<Complexe, Dimension>::GetVelocityOfMedia(int ref) const
  {
    TinyVector<Complexe, Dimension::dim_N> lambda;
    
    Real_wp Lmax = 0.0;
    TinyVector<Real_wp, Dimension::dim_N> kwave;
    int nb_points1 = 30;
    int nb_points2 = 30;
    if (Dimension::dim_N == 2)
      nb_points2 = 1;
    
    for (int m = 0; m < nb_points1; m++)
      for (int n = 0; n < nb_points2; n++)
        {
          Real_wp theta = -0.5*pi_wp + pi_wp*Real_wp(m) / nb_points1;
          Real_wp phi = 2.0*pi_wp*Real_wp(n) / nb_points2;
          SetIncidentAngle(Real_wp(1), kwave, theta, phi);
          GetEigenvalues(ref_tensorC(ref).GetConstant(), kwave, lambda);
          for (int i = 0; i < lambda.GetM(); i++)
            Lmax = max(Lmax, abs(lambda(i)));
        }
    
    return sqrt(Lmax/abs(ref_rho(ref).GetConstant()));
  }

  
  //! returns the velocity of waves at infinity
  template<class Complexe, class Dimension>
  Real_wp VarElastic_Dim<Complexe, Dimension>::GetVelocityOfInfinity() const
  {
    return 1.0;
  }
  

  //! allocation of arrays needed to compute finite element matrix
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::AllocateMassMatrices()
  {    
    this->coefficient_mu_penalty.Reallocate(this->GetNbPhysicalIndices());
    for (int i = 0; i < this->coefficient_mu_penalty.GetM(); i++)
      {
        this->coefficient_mu_penalty(i) = MaxAbs(ref_tensorC(i).GetConstant());
      }
  }
  
  
  //! fills the varying indices that need to be computed  
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::
  GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complexe>* >& rho_complex,
		    IVect& num_ref, IVect& num_index, IVect& num_component,
		    Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    int nb = 0;
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        nb += ref_rho(i).GetNbVaryingMedia();
        nb += ref_sigma(i).GetNbVaryingMedia();
        nb += ref_tensorC(i).GetNbVaryingMedia();
      }
    
    rho_complex.Reallocate(nb);
    num_ref.Reallocate(nb);
    num_index.Reallocate(nb);
    num_component.Reallocate(nb);
    compute_grad.Reallocate(nb);
    compute_hess.Reallocate(nb);
    compute_grad.Fill(false);
    compute_hess.Fill(false);
    nb = 0;
    for (int i = 0; i < ref_rho.GetM(); i++)
      {
        int nb0 = nb;
        ref_rho(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        ref_tensorC(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 1;
            num_ref(j) = i;
          }
	
	nb0 = nb;
        ref_sigma(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 2;
            num_ref(j) = i;
          }
      }
  }  
  
  
  //! adimensionalization of equations if required
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::FinalizeComputationVaryingIndices()
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        Vector<bool> IsRefUsed(ref_rho.GetM());
        IsRefUsed.Fill(false);
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          {
            int ref = var_problem.mesh.Element(i).GetReference();
            IsRefUsed(ref) = true;
          }
        
        int nb_ref_used = 0;
        for (int i = 0; i < ref_rho.GetM(); i++)
          if (IsRefUsed(i))
            nb_ref_used++;
        
        IVect RefUsed(nb_ref_used);
        nb_ref_used = 0;
        for (int i = 0; i < ref_rho.GetM(); i++)
          if (IsRefUsed(i))
            RefUsed(nb_ref_used++) = i;
                
        // retrieving the maximal values for rho and C
        Real_wp rho_max(0), Cmax(0);
        TinySymmetricTensor<Complexe, Dimension::dim_N> cteC;
        for (int n = 0; n < RefUsed.GetM(); n++)
          {
            int ref = RefUsed(n);
            rho_max = max(rho_max, abs(this->ref_rho(ref).GetConstant()));
            cteC = this->ref_tensorC(ref).GetConstant();
            Cmax = max(Cmax, abs(cteC(0, 0, 0, 0)));
          }

#ifdef SELDON_WITH_MPI
        int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
        
        R2 vec_u(rho_max, Cmax), vec_v;
        Vector<int64_t> vec_tmp;
        MpiAllreduce(var_problem.comm_group_mode, &vec_u(0), vec_tmp, &vec_v(0), 2, MPI_MAX);
        rho_max = vec_v(0); Cmax = vec_v(1);
#else
        int rank_proc(0);
#endif
        
        u_bar = 1.0/Cmax;
        omega_bar = sqrt(Cmax/rho_max);
        Real_wp omegap = var_problem.GetOmega()/omega_bar;
        if ((rank_proc == 0) && (var_problem.print_level >= 2))
          {
            cout << "Value of displacement u used in adimensionalization = " << u_bar << endl;
            cout << "Value of pulsation omega used in adimensionalization = " << omega_bar << endl;
            cout << "New value of pulsation = " << omegap << endl;
          }
        
        Complexe coef_rho = square(omega_bar)*u_bar;
        Complexe coef_C = u_bar, coef_A = 1.0/coef_C;
        for (int n = 0; n < RefUsed.GetM(); n++)
          {
            int ref = RefUsed(n);
            this->ref_rho(ref).Mlt(coef_rho);
            this->ref_sigma(ref).Mlt(coef_rho);
            this->ref_tensorC(ref).Mlt(coef_C);
            this->ref_tensorA(ref).Mlt(coef_A);
          }
        
        var_problem.SetOmega(omegap);
      }
  }  
  
  
  //! computation of the number of degrees of freedom
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>
  ::ComputeNumberOfDofs(int& nodl, IVect& offset_dof_unknown,
			IVect& OffsetDofV)
  {

    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    if (var_problem.FirstOrderFormulation())
      {
        nodl += mesh_num.GetNbDofPML()*(Dimension::dim_N-1)*Dimension::dim_N;    
        int nodl_scalar = mesh_num.GetNbDof();
        nodl_scalar += mesh_num.GetNbDofPML()*(Dimension::dim_N-1);
        for (int m = 0; m <= var_problem.nb_unknowns_scal; m++)
	  offset_dof_unknown(m) = m*nodl_scalar;
	
	int Nvol = mesh_num.GetNbDof();
	for (int m = var_problem.nb_unknowns_scal+1; m <= var_problem.nb_unknowns; m++)
	  offset_dof_unknown(m) = offset_dof_unknown(m-1) + Nvol;
	
        int nb_dof_pml = 0, offset_v = 0;
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          {
            int nb_dof_v = 0;
            if (var_problem.InsidePML(i))
              nb_dof_v = (OffsetDofV(i+1) - offset_v)*(Dimension::dim_N-1);
            
            offset_v = OffsetDofV(i+1);
            OffsetDofV(i+1) += nb_dof_v + nb_dof_pml;
            nb_dof_pml += nb_dof_v;
          }
        
        nodl += nb_dof_pml;
      }
  }

  
  //! Adds dofs due to pml for parallel computations
  template<class Complexe, class Dimension>
  void VarElastic_Dim<Complexe, Dimension>::PutOtherGlobalDofs()
  {
    MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    if (var_problem.FirstOrderFormulation())
      {
        int offset_loc = mesh_num.GetNbDof();
        int offset_glob = var_problem.GetNbGlobalMeshDof();
        int offset_loc2 = offset_loc + mesh_num.GetNbDofPML();
        int offset_glob2 = offset_glob + var_problem.GetNbGlobalDofPML();
        int d = Dimension::dim_N;
        
        if (mesh_num.GetNbDofPML() > 0)
          {    
            mesh_num.GlobDofNumber_Subdomain.Resize(var_problem.offset_dof_unknown(1));
            for (int i = 0; i < mesh_num.GetNbDof(); i++)
              {
                int npml = mesh_num.GetDofPML(i);
                if (npml >= 0)
                  {
                    mesh_num.GlobDofNumber_Subdomain(offset_loc + npml)
                      = offset_glob + mesh_num.GlobDofPML_Subdomain(npml);
                    
                    if (Dimension::dim_N == 3)
                      mesh_num.GlobDofNumber_Subdomain(offset_loc2 + npml)
                        = offset_glob2 + mesh_num.GlobDofPML_Subdomain(npml);
                  }
              }
            
	    Vector<IVect>& MatchingDof = var_problem.GetOriginalMatchingDofNeighbor();
            for (int num = 0; num < MatchingDof.GetM(); num++)
              {
                int nb_dof = MatchingDof(num).GetM();
                int nb_dof_pml = 0;
                for (int j = 0; j < nb_dof; j++)
                  if (mesh_num.GetDofPML(MatchingDof(num)(j)) >= 0)
                    nb_dof_pml++;
                
                if (nb_dof_pml > 0)
                  {
                    int nb_old = nb_dof;
                    MatchingDof(num).Resize(nb_old + nb_dof_pml*(d-1));
                    for (int j = 0; j < nb_old; j++)
                      {
                        int npml
                          = mesh_num.GetDofPML(MatchingDof(num)(j));
                        
                        if (npml >= 0)
                          {
                            MatchingDof(num)(nb_dof++) = offset_loc + npml;
                            if (d == 3)
                              MatchingDof(num)(nb_dof++) = offset_loc2 + npml;
                          }
                      }
                  }
              }
          }
      }    
  }
  
  
  /*****************
   * VarElastic_Eq *
   *****************/
  

  template<class TypeEquation>
  VarElastic_Eq<TypeEquation>::VarElastic_Eq()
    : VarElastic_Dim<Complexe, Dimension>(this->GetLeafClass())
  {
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarElastic_Base::SetInputData(description_field, parameters);
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
  }


  //! switching to the first-order formulation (in time)
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::SetFirstOrderFormulation(bool mix_)
  {
    this->mixed_formulation = mix_;
    if (this->mixed_formulation)
      this->mesh_num.compute_dof_pml = true;
  }
  
  
  //! returns true if the associated finite element matrix is symmetric
  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsSymmetricProblem(bool eigen) const
  {
    if (this->FirstOrderFormulation())
      return false;
    
    return VarHarmonic<TypeEquation>::IsSymmetricProblem(eigen);
  }

  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation() && nat_mat.GetCoefStiffness() == Real_wp(0))
      {
        Vector<bool> diag_elt;
        int type_mat = this->GetMassMatrixType(diag_elt);
        if (type_mat == FemMassMatrix::DIAGONAL)
          return true;
      }

    return VarHarmonic<TypeEquation>::IsSparseElementaryMatrix(nat_mat);
  }


  template<class TypeEquation>  
  bool VarElastic_Eq<TypeEquation>::IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation() && nat_mat.GetCoefStiffness() == Complex_wp(0, 0))
      {
        Vector<bool> diag_elt;
        int type_mat = this->GetMassMatrixType(diag_elt);
        if (type_mat == FemMassMatrix::DIAGONAL)
          return true;
      }

    return VarHarmonic<TypeEquation>::IsSparseElementaryMatrix(nat_mat);
  }

  
  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation())
      return false;

    return VarHarmonic<TypeEquation>::IsSymmetricElementaryMatrix(nat_mat);
  }


  template<class TypeEquation>  
  bool VarElastic_Eq<TypeEquation>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>& nat_mat) const
  {
    if (this->FirstOrderFormulation())
      return false;

    return VarHarmonic<TypeEquation>::IsSymmetricElementaryMatrix(nat_mat);
  }


  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {    
    Vector<bool> diag_elt;
    int type_mat = this->GetMassMatrixType(diag_elt);
    if (type_mat == FemMassMatrix::DIAGONAL)
      {
        if (this->FirstOrderFormulation() && (nat_mat.GetCoefDamping() != Real_wp(0)))
          return false;

        if (nat_mat.GetCoefStiffness() == Real_wp(0))
          return true;
      }

    return false;
  }
  

  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsDiagonalElementaryMatrix(const GlobalGenericMatrix<Complex_wp>& nat_mat) const
  {    
    Vector<bool> diag_elt;
    int type_mat = this->GetMassMatrixType(diag_elt);
    if (type_mat == FemMassMatrix::DIAGONAL)
      {
        if (this->FirstOrderFormulation() && (nat_mat.GetCoefDamping() != Complex_wp(0, 0)))
          return false;
                
        if (nat_mat.GetCoefStiffness() == Complex_wp(0, 0))
          return true;
      }

    return false;

  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::ComputeNumberOfDofs()
  {
    VarHarmonic<TypeEquation>::ComputeNumberOfDofs();
    VarElastic_Dim<Complexe, Dimension>::
      ComputeNumberOfDofs(this->nodl, this->offset_dof_unknown,
			  this->OffsetDofV);    
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::PutOtherGlobalDofs()
  {
    VarElastic_Dim<Complexe, Dimension>::PutOtherGlobalDofs();
  }
  

  //! treating Dirichlet condition (retrieving dof numbers)
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::TreatDirichletCondition()
  {
    VarHarmonic<TypeEquation>::TreatDirichletCondition();      
    if (this->FirstOrderFormulation())
      {
	int n = this->is_dof_dirichlet.GetM();
	if (this->nodl > n)
	  {
	    this->is_dof_dirichlet.Resize(this->nodl);
	    for (int i = n; i < this->nodl; i++)
	      this->is_dof_dirichlet(i) = false;
	  }
	
	int d = Dimension::dim_N;
	int nb_dof = 0;
	for (int i = 0; i < this->Dirichlet_dof.GetM(); i++)
	  {
	    int j = this->Dirichlet_dof(i);
	    if (j < this->mesh_num.GetNbDof())
              {
                if (this->mesh_num.GetDofPML(j) >= 0)
                  nb_dof += d*(d-1);
              }
	  }
	
	if (nb_dof > 0)
	  {
	    int nb_old = this->Dirichlet_dof.GetM();
	    int Nvol = this->mesh_num.GetNbDof();
	    int Nvol_pml = this->mesh_num.GetNbDofPML();
	    this->Dirichlet_dof.Resize(nb_old + nb_dof);
            this->nb_dof_dirichlet = nb_old + nb_dof;
	    nb_dof = nb_old;
	    for (int i = 0; i < nb_old; i++)
	      {
		int j = this->Dirichlet_dof(i);
		if (j < Nvol)
                  if (this->mesh_num.GetDofPML(j) >= 0)
                    {
                      int k = this->mesh_num.GetDofPML(j);
                      if (d == 2)
                        {
                          this->Dirichlet_dof(nb_dof) = k + Nvol;
                          this->Dirichlet_dof(nb_dof+1) = k + 2*Nvol + Nvol_pml;
                        }
                      else
                        {
                          this->Dirichlet_dof(nb_dof) = k + Nvol;
                          this->Dirichlet_dof(nb_dof+1) = k + Nvol + Nvol_pml;
                          this->Dirichlet_dof(nb_dof+2) = k + 2*Nvol + 2*Nvol_pml;
                          this->Dirichlet_dof(nb_dof+3) = k + 2*Nvol + 3*Nvol_pml;
                          this->Dirichlet_dof(nb_dof+4) = k + 3*Nvol + 4*Nvol_pml;
                          this->Dirichlet_dof(nb_dof+5) = k + 3*Nvol + 5*Nvol_pml;
                        }
                      
                      nb_dof += d*(d-1);
                    }
	      }
	    
            Sort(this->Dirichlet_dof);
	    for (int i = nb_old; i < this->Dirichlet_dof.GetM(); i++)
	      this->is_dof_dirichlet(this->Dirichlet_dof(i)) = true;
	  }
      }
  }

  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::UpdateShiftAdimensionalization(Real_wp& sr, Real_wp& si)
  {
    VarElastic_Base::UpdateShiftAdimensionalization(sr, si);
  }


  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::UpdateShiftAdimensionalization(Complex_wp& sr, Complex_wp& si)
  {
    VarElastic_Base::UpdateShiftAdimensionalization(sr, si);
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::UpdateEigenvaluesAdimensionalization(Vector<Real_wp>&Lr, Vector<Real_wp>& Li, 
					 Matrix<Real_wp, General, ColMajor>& V)
  {
    VarElastic_Base::UpdateEigenvaluesAdimensionalization(Lr, Li, V);
  }


  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::UpdateEigenvaluesAdimensionalization(Vector<Complex_wp>&Lr, Vector<Complex_wp>& Li, 
					 Matrix<Complex_wp, General, ColMajor>& V)
  {
    VarElastic_Base::UpdateEigenvaluesAdimensionalization(Lr, Li, V);
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<Real_wp>& val_u,
			Vector<Real_wp>& grad_u, int i,
			const GridInterpolation<Dimension>& var_interp,
			int iquad, bool compute_grad) const
  {
    VarElastic_Base::ModifyOutputUnknown(this->ref_tensorC, *this, val_u, grad_u,
					 i, var_interp, iquad, compute_grad);
  }


  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<Complex_wp>& val_u,
			Vector<Complex_wp>& grad_u, int i,
			const GridInterpolation<Dimension>& var_interp,
			int iquad, bool compute_grad) const
  {
    VarElastic_Base::ModifyOutputUnknown(this->ref_tensorC, *this, val_u, grad_u,
					 i, var_interp, iquad, compute_grad);
  }
  

  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                        int, bool, bool) const
  {
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                        int, bool, bool) const
  {
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, Vector<Real_wp>& En_quad, Vector<Real_wp>& Hn_quad) const
  {
    VarElastic_Base::ComputeEnHnQuadrature(u_quadrature, grad_quadrature, num_elem,
					   pts, normale, compute_H, En_quad, Hn_quad);
  }


  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, Vector<Complex_wp>& En_quad, Vector<Complex_wp>& Hn_quad) const
  {
    VarElastic_Base::ComputeEnHnQuadrature(u_quadrature, grad_quadrature, num_elem,
					   pts, normale, compute_H, En_quad, Hn_quad);
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ComputeEnHnNodal(Vector<Vector<Real_wp> >& u_nodal, Vector<Vector<Real_wp> >& grad_nodal,
		     int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
		     Vector<Vector<Real_wp> >& En_nodal, Vector<Vector<Real_wp> >& Hn_nodal) const
  {
    VarElastic_Base::ComputeEnHnNodal(u_nodal, grad_nodal, num_elem, pts, normale,
				      En_nodal, Hn_nodal);
  }


  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::ComputeEnHnNodal(Vector<Vector<Complex_wp> >& u_nodal, Vector<Vector<Complex_wp> >& grad_nodal,
		     int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
		     Vector<Vector<Complex_wp> >& En_nodal, Vector<Vector<Complex_wp> >& Hn_nodal) const
  {
    VarElastic_Base::ComputeEnHnNodal(u_nodal, grad_nodal, num_elem, pts, normale,
				      En_nodal, Hn_nodal);
  }
  
  
  template<class TypeEquation>
  Real_wp VarElastic_Eq<TypeEquation>::GetCoefficientPenaltyStiffness(int ref) const
  {
    return VarElastic_Base::GetCoefficientPenaltyStiffness(ref);
  }

  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::InitIndices(int n)
  {
    VarElastic_Dim<Complexe, Dimension>::InitIndices(n);
  }
  
  
  template<class TypeEquation>
  int VarElastic_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarElastic_Dim<Complexe, Dimension>::GetNbPhysicalIndices();
  }
    
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::CopyInputData(const VarProblem_Base& var_base)
  {
    const VarElastic_Eq<TypeEquation>& var = static_cast<const VarElastic_Eq<TypeEquation>& >(var_base);
    VarHarmonic<TypeEquation>::CopyInputData(var);
    VarElastic_Dim<Complexe, Dimension>::CopyInputData(var);
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    VarElastic_Dim<Complexe, Dimension>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarElastic_Dim<Complexe, Dimension>::SetPhysicalIndex(name_media, i, parameters);
  }
  

  template<class TypeEquation>
  string VarElastic_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarElastic_Dim<Complexe, Dimension>::GetPhysicalIndexName(m);
  }

  
  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarElastic_Dim<Complexe, Dimension>::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool VarElastic_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarElastic_Dim<Complexe, Dimension>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  Real_wp VarElastic_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return VarElastic_Dim<Complexe, Dimension>::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  Real_wp VarElastic_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return VarElastic_Dim<Complexe, Dimension>::GetVelocityOfInfinity();
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::AllocateMassMatrices()
  {
    VarElastic_Dim<Complexe, Dimension>::AllocateMassMatrices();
  }
  
  
  template<class TypeEquation>
  void VarElastic_Eq<TypeEquation>::FinalizeComputationVaryingIndices()
  {
    VarElastic_Dim<Complexe, Dimension>::FinalizeComputationVaryingIndices();
  }

  
  template<class TypeEquation>
  void VarElastic_Cplx<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    VarElastic_Dim<Complexe, Dimension>::GetVaryingIndices(rho_complex, num_ref, num_index, num_component,
							   compute_grad, compute_hess);
  }

  
  template<class TypeEquation>
  void VarElastic_Real<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* >& rho_complex,
		      Vector<PhysicalVaryingMedia<Dimension, Real_wp>* >& rho_real, IVect& num_ref,
		      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		      Vector<bool>& compute_hess)
  {
    VarElastic_Dim<Complexe, Dimension>::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
							   compute_grad, compute_hess);
  }
  
  
  /****************
   * ImpedanceABC *
   ****************/

  
  //! computation of impedance
  template<class Dimension>
  void ImpedanceABC<Complex_wp, HarmonicElasticEquation<Dimension> >
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			 const GlobalGenericMatrix<Complex_wp>& nat_mat, int ref,
			 const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    if (k == 0)
      feval.Reallocate(Mat.GetNbPointsQuadratureBoundary());
    
    Complex_wp c11 = var_elas.ref_tensorC(ref)(0, 0, 0, 0);
    Complex_wp c33 = var_elas.ref_tensorC(ref)(0, 1, 0, 1);
    Complex_wp rho = var_elas.ref_rho(ref);
    rho += Iwp*var_elas.ref_sigma(ref).GetConstant()/var_problem.GetOmega();
    typename Dimension::R_N n = Mat.GetNormaleQuadratureBoundary(k);
    TinyMatrix<Complex_wp, Symmetric, d, d> Nn;
    GetNormalProjector(n, Nn);
    GetTangentialProjector(n, feval(k));
    Mlt(sqrt(abs(c11*rho)), Nn);
    Mlt(sqrt(abs(c33*rho)), feval(k));
    feval(k) += Nn;
    feval(k) *= -Iwp*var_problem.GetOmega();
  }
  
  
  //! application of impedance to phi
  template<class Dimension>
  void ImpedanceABC<Complex_wp, HarmonicElasticEquation<Dimension> >
  ::ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi,
                         const typename Dimension::R_N& dphi, Vector<Complex_wp>& feval_phi)
  {
    TinyVector<Complex_wp, d> phi_, fphi_;
    phi_(n) = phi(0);
    Mlt(feval(j), phi_, fphi_);
    CopyVector(fphi_, 0, feval_phi);
  }


  //! computation of impedance
  template<class T, class Dimension>
  void ImpedanceABC<T, ElasticEquation<Dimension> >
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
			 const GlobalGenericMatrix<T>& nat_mat, int ref,
			 const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    if (k == 0)
      feval.Reallocate(Mat.GetNbPointsQuadratureBoundary());
    
    Real_wp c11 = var_elas.ref_tensorC(ref)(0, 0, 0, 0);
    Real_wp c33 = var_elas.ref_tensorC(ref)(0, 1, 0, 1);
    Real_wp rho = var_elas.ref_rho(ref);
    typename Dimension::R_N n = Mat.GetNormaleQuadratureBoundary(k);
    TinyMatrix<T, Symmetric, d, d> Nn;
    GetNormalProjector(n, Nn);
    GetTangentialProjector(n, feval(k));
    Mlt(sqrt(abs(c11*rho)), Nn);
    Mlt(sqrt(abs(c33*rho)), feval(k));
    feval(k) += Nn;
    feval(k) *= nat_mat.GetCoefDamping();
  }
  
  
  //! application of impedance to phi
  template<class T, class Dimension>
  void ImpedanceABC<T, ElasticEquation<Dimension> >
  ::ApplyImpedancePhi_H1(int n, int j, int, const TinyVector<Real_wp, 1>& phi,
                         const typename Dimension::R_N& dphi, Vector<T>& feval_phi)
  {
    TinyVector<T, d> phi_, fphi_;
    phi_(n) = phi(0);
    Mlt(feval(j), phi_, fphi_);
    CopyVector(fphi_, 0, feval_phi);
  }

  
  /*********************
   * ElasticEquationDG *
   *********************/

  
  //! not used
  template<class T, class Dim> template<class TypeEquation>
  void ElasticEquationDG_Base<T, Dim>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                      int num_elem, const ElementReference_Dim<Dim>& Fb)
  {
  }
  
  
  //! which unknowns and function tests need to be derivated ?
  template<class T, class Dim> template<class GenericPb, class T0, class Matrix1>
  void ElasticEquationDG_Base<T, Dim>
  ::GetNeededDerivative(const GenericPb& vars,
			const GlobalGenericMatrix<T0>& nat_mat, Matrix1& At, Matrix1& A)
  {
    A.Fill(false); At.Fill(false);
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        A(i) = true;
        At(i) = true;
      }
  }
  
  
  //! fills tensors D and E appearing in the variational formulation
  template<class T, class Dim>
  template<class TypeEquation, class T0, class MatStiff>
  void ElasticEquationDG_Base<T, Dim>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                     int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
                     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    T0 s(nat_mat.GetCoefStiffness());
    TinyMatrix<TinyMatrix<T,
      General, Dimension::dim_N, Dimension::dim_N>,
      General, Dimension::dim_N, Dimension::dim_N> C;
    
    vars.ref_tensorC(ref).GetCoefficient(vars, num_elem, jloc, C);
    if (Dimension::dim_N == 2)
      {
        Dgrad_phi(0, 2)(0) = -s;
        Dgrad_phi(0, 3)(1) = -s;
        Dgrad_phi(1, 3)(0) = -s;
        Dgrad_phi(1, 4)(1) = -s;
        
        Ephi_grad(2, 0)(0) = C(0, 0)(0, 0)*s;
        Ephi_grad(2, 1)(1) = C(0, 1)(0, 1)*s;
        Ephi_grad(2, 0)(1) = C(0, 0)(0, 1)*s;
        Ephi_grad(2, 1)(0) = C(0, 0)(0, 1)*s;

        Ephi_grad(3, 0)(0) = C(0, 0)(0, 1)*s;
        Ephi_grad(3, 1)(1) = C(0, 1)(1, 1)*s;
        Ephi_grad(3, 0)(1) = C(0, 0)(1, 1)*s;
        Ephi_grad(3, 1)(0) = C(0, 0)(1, 1)*s;

        Ephi_grad(4, 0)(0) = C(0, 1)(0, 1)*s;
        Ephi_grad(4, 1)(1) = C(1, 1)(1, 1)*s;
        Ephi_grad(4, 0)(1) = C(0, 1)(1, 1)*s;
        Ephi_grad(4, 1)(0) = C(0, 1)(1, 1)*s;
      }
    else
      {
        abort();
      }
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                         int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                         int ref, const TinyVector<TinyVector<T0, 2>, 5>& dU,
                         TinyVector<T0, 5>& V)
  {
    abort();
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                         int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                         int ref, const TinyVector<TinyVector<T0, 3>, 9>& dU,
                         TinyVector<T0, 9>& V)
  {
    abort();
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                         int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                         int ref, const TinyVector<T0, 5>& U,
                         TinyVector<TinyVector<T0, 2>, 5>& dV)
  {
    abort();
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                         int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                         int ref, const TinyVector<T0, 9>& U,
                         TinyVector<TinyVector<T0, 3>, 9>& dV)
  {
    abort();
  }
    
  
  template<class T, class Dim>
  template<class TypeEquation, class T0, class MatMass>
  void ElasticEquationDG_Base<T, Dim>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                  int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
  {
    T rho, m_iomega, sigma;
    vars.GetMiomega(m_iomega);
    rho = vars.ref_rho(ref).GetCoefficient(vars, i, j);
    sigma = vars.ref_sigma(ref).GetCoefficient(vars, i, j);
    T0 coef;
    coef = m_iomega*rho*nat_mat.GetCoefMass() + sigma*nat_mat.GetCoefDamping();
    for (int i = 0; i < Dimension::dim_N; i++)
      mass(i, i) = coef;
    
    coef = m_iomega*nat_mat.GetCoefMass();
    for (int i = Dimension::dim_N; i < mass.GetM(); i++)
      mass(i, i) = coef;
  }
  
  
  template<class T, class Dim>
  template<class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                    const GlobalGenericMatrix<T0>& nat_mat, int ref,
                    const TinyVector<T0, 5>& U, TinyVector<T0, 5>& V)
  {
    abort();
  }
  
  
  template<class T, class Dim> template<class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                    const GlobalGenericMatrix<T0>& nat_mat, int ref,
                    const TinyVector<T0, 9>& U, TinyVector<T0, 9>& V)
  {
    abort();
  }
  
  
  template<class T, class Dim>
  template<class Matrix1, class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale,
            int ref, int num_elem, int jloc,
            const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
            const EllipticProblem<TypeEquation>& vars,
	    const ElementReference<Dim, 1>& Fb)
  {
    T0 s = nat_mat.GetCoefStiffness();
    Nabc.Fill(0);
    int cond = vars.mesh.GetBoundaryCondition(ref);

    TinyMatrix<TinyMatrix<T,
      General, Dimension::dim_N, Dimension::dim_N>,
      General, Dimension::dim_N, Dimension::dim_N> C;
    
    vars.ref_tensorC(ref_d).GetCoefficient(vars, num_elem, jloc, C);
    
    if (Dimension::dim_N == 2)
      {        
        if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
          {
            Nabc(0, 2) = -normale(0)*s;
            Nabc(0, 3) = -normale(1)*s;
            Nabc(1, 3) = -normale(0)*s;
            Nabc(1, 4) = -normale(1)*s;
            
            Nabc(2, 0) = (C(0, 0)(0, 0)*normale(0) + C(0, 0)(0, 1)*normale(1))*s;
            Nabc(2, 1) = (C(0, 1)(0, 1)*normale(1) + C(0, 0)(0, 1)*normale(0))*s;
            Nabc(3, 0) = (C(0, 0)(0, 1)*normale(0) + C(0, 0)(1, 1)*normale(1))*s;
            Nabc(3, 1) = (C(0, 1)(1, 1)*normale(1) + C(0, 0)(1, 1)*normale(0))*s;
            Nabc(4, 0) = (C(0, 1)(0, 1)*normale(0) + C(0, 1)(1, 1)*normale(1))*s;
            Nabc(4, 1) = (C(1, 1)(1, 1)*normale(1) + C(0, 1)(1, 1)*normale(0))*s;
          }
        else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
          {
            Nabc(0, 2) = normale(0)*s;
            Nabc(0, 3) = normale(1)*s;
            Nabc(1, 3) = normale(0)*s;
            Nabc(1, 4) = normale(1)*s;
            
            Nabc(2, 0) = -(C(0, 0)(0, 0)*normale(0) + C(0, 0)(0, 1)*normale(1))*s;
            Nabc(2, 1) = -(C(0, 1)(0, 1)*normale(1) + C(0, 0)(0, 1)*normale(0))*s;
            Nabc(3, 0) = -(C(0, 0)(0, 1)*normale(0) + C(0, 0)(1, 1)*normale(1))*s;
            Nabc(3, 1) = -(C(0, 1)(1, 1)*normale(1) + C(0, 0)(1, 1)*normale(0))*s;
            Nabc(4, 0) = -(C(0, 1)(0, 1)*normale(0) + C(0, 1)(1, 1)*normale(1))*s;
            Nabc(4, 1) = -(C(1, 1)(1, 1)*normale(1) + C(0, 1)(1, 1)*normale(0))*s;
          }
      }
    else
      {
        abort();
      }
  }
  
  
  template<class T, class Dim>
  template<class Vector1, class TypeEquation, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
            int num_elem1, int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d,
            const EllipticProblem<TypeEquation>& vars, const ElementReference<Dim, 1>& Fb)
  {
    abort();
  }

  
  template<class T, class Dim>
  template<class Matrix1, class GenericPb, class T0>
  void ElasticEquationDG_Base<T, Dim>
  ::GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad, int k,
               int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
               const GenericPb& vars, const ElementReference<Dim, 1>& Fb)
  {
    Nabc.Fill(0);
    /* typedef typename NatureMatrix::value_type T0;
    TinyMatrix<TinyVector<T0, 2>, General, 5, 5> Dgrad_phi, Ephi_grad;
    GetGradPhiTensor(vars, iquad, k, nat_mat, ref, Dgrad_phi, Ephi_grad);
    Complex_wp rho = vars.ref_rho(1).GetCoefficient(vars, iquad, k);
    for (int n = 0; n < Dimension::dim_N; n++)
      for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
          Nabc(i, j) += normale(n)*(Ephi_grad(i, j)(n) - Dgrad_phi(i, j)(n));

    for (int j = 0; j < 5; j++)
      {
        Nabc(0, j) /= rho;
        Nabc(1, j) /= rho;
      }
    
    GetAbsoluteValue(Nabc);

    for (int j = 0; j < 5; j++)
      {
        Nabc(0, j) *= rho;
        Nabc(1, j) *= rho;
      }
    
      Mlt(vars.alpha_penalization*nat_mat.GetCoefDamping(), Nabc); */
    for (int i = 0; i < Dimension::dim_N; i++)
      Nabc(i, i) = vars.alpha_penalization*nat_mat.GetCoefDamping();
    
    if (Dimension::dim_N == 2)
      { 
        //typename NatureMatrix::value_type coef = vars.delta_penalization*nat_mat.GetCoefDamping();
        //Nabc(2, 2) = normale(0)*coef;
        //Nabc(2, 3) = normale(1)*coef;
        //Nabc(4, 3) = normale(0)*coef;
        //Nabc(4, 4) = normale(1)*coef;
      }
    else
      {
        abort();
      }
  }
  
  
  template<class T, class Dim>
  template<class T0, class GenericPb>
  void ElasticEquationDG_Base<T, Dim>
  ::MltPenalDG(const typename Dimension::R_N& normale,
               const TinyVector<T0, 5>& Vn, TinyVector<T0, 5>& Un,
               int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
               int ref, int ref2, const GenericPb& vars, const ElementReference<Dim, 1>& Fb)
  {
    abort();
  }
  
  
  template<class T, class Dim>
  template<class T0, class GenericPb>
  void ElasticEquationDG_Base<T, Dim>
  ::MltPenalDG(const typename Dimension::R_N& normale,
               const TinyVector<T0, 9>& Vn, TinyVector<T0, 9>& Un,
               int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
               int ref, int ref2, const GenericPb& vars, const ElementReference<Dim, 1>& Fb)
  {
    abort();
  }

  
  /*******************
   * EllipticProblem *
   *******************/
  

  template<class Dimension>
  void EllipticProblem<ElasticEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    if (!this->FirstOrderFormulation())
      {
        // using the generic function
	Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
					  this->GetReferenceElementH1(i));
        return;
      }
    
    VarElastic_Base::ComputeElementaryMatrix(this->ref_tensorC, this->ref_rho, this->ref_sigma,
					     i, num_dof, mat_elem, nat_mat, *this, *this,
					     this->GetReferenceElementH1(i));
  }


  template<class Dimension>
  void EllipticProblem<ElasticEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  

  template<class Dimension>
  void EllipticProblem<ElasticEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class Dimension>
  void EllipticProblem<ElasticEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquation<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    if (!this->FirstOrderFormulation())
      {
        // using the generic function
	Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
					  this->GetReferenceElementH1(i));
        return;
      }
    
    VarElastic_Base::ComputeElementaryMatrix(this->ref_tensorC, this->ref_rho, this->ref_sigma,
					     i, num_dof, mat_elem, nat_mat, *this, *this,
					     this->GetReferenceElementH1(i));
  }
  

  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquation<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class Dimension>
  void EllipticProblem<ElasticEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  template<class Dimension>
  void EllipticProblem<ElasticEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }  

  template<class Dimension>
  void EllipticProblem<ElasticEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }


  template<class Dimension>
  void EllipticProblem<ElasticEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquationDG<Dimension> >
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }
  

  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void EllipticProblem<HarmonicElasticEquationDG<Dimension> >
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);    
  }

  
  /**********************
   * FemMatrixFreeClass *
   **********************/
  
  
  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddVectorH1(Real_wp(1), nat_mat, trans, lvl, *this, 
		   X, Real_wp(1), Y, false);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, ElasticEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, HarmonicElasticEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Not possible" << endl;
    abort();
  }


  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, HarmonicElasticEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddVectorH1(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		   X, Complex_wp(1, 0), Y, false);
  }

}

#define MONTJOIE_FILE_VAR_ELASTIC_CXX
#endif
