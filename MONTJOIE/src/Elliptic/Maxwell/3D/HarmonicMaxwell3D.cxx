#ifndef MONTJOIE_FILE_HARMONIC_MAXWELL_3D_CXX

namespace Montjoie
{
  
  /******************************
   * HarmonicMaxwellEquation_3D *
   ******************************/
  

  template<class T>
  bool MaxwellEquation3D_Base<T>::store_dfjm1(false);
  
  //! computation of mass coefficients
  /*!
    \param[in,out] var given problem
    \param[in] num_elem element number
    \param[in] N number of quadrature points
    \param[in] Points points after transformation F_i
    \param[in] mat jacobian matrices
  */
  template<class T> template<class TypeEquation>
  void MaxwellEquation3D_Base<T>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                      int num_elem, const ElementReference_Dim<Dimension3>& Fb0)
  {
    const ElementReference<Dimension3, 2>& Fb = dynamic_cast<const ElementReference<Dimension3, 2>& >(Fb0);
    
    int ref = var.mesh.Element(num_elem).GetReference();
    bool variable = var.UseNumericalIntegration(num_elem);
    bool affine = var.mesh.IsElementAffine(num_elem);
    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;

    typedef typename TypeEquation::Complexe Complexe;
    HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass = var.Glob_matMass_elem(num_elem);
    
    int i1 = num_elem - var.mesh.GetNbElt() + var.GetNbEltPML();    
    TinyMatrix<T, General, 3, 3> A_tmp;
    TinyMatrix<T, Symmetric, 3, 3> mu, invMu, epsilon, sigma, Bmass, Astiff;
    epsilon = var.ref_epsilon(ref).GetConstant();
    sigma = var.ref_sigma(ref).GetConstant();
    mu = var.ref_mu(ref).GetConstant();
    
    T momega2, m_iomega;
    var.GetMiomega(m_iomega);
    var.GetMomega2(momega2);
    
    bool diag_mass = false;
    if ((var.OrthogonalElement(num_elem) == 0) && Fb.MassLumpingOrthogonalElement() && affine)
      diag_mass = true;
    
    if ((var.ref_epsilon(ref).GetAnisotropy() == TensorPhysicalIndice<Dimension3, 3, Complexe>::ANISOTROPE) ||
	(var.ref_sigma(ref).GetAnisotropy() == TensorPhysicalIndice<Dimension3, 3, Complexe>::ANISOTROPE) )
      diag_mass = false;

    if (!diag_mass)
      mass.ClearDiagonalDh();
    
    if (!var.FirstOrderFormulation())
      {
        // second-order formulation
        // solved equation : -omega^2 epsilon u - i omega sigma u + rot(1/mu rot u) = f
        
        // we store Bh = -omega^2 Ji(xi_k) DF_i^{-1} epsilon DF_i^{*-1} \omega_k
        // where (xi_k, omega_k) is the quadrature formula used in the element
        // BhSigma = -i omega Ji(xi_k) DF_i^{-1} sigma DF_i^{*-1} \omega_k
        // Ah = 1/Ji(xi_k) DF_i^{*} \mu^{-1} DF_i omega_k
        // all these expressions are computed without omega_k 
        // for constant coefficients (uniform epsilon, mu and affine element)
        mass.ReallocateAhBh(N, N, N);
        
        // loop on quadrature points
        for (int j = 0; j < N; j++)
          {
	    // evaluation of epsilon and mu if variables
	    epsilon = var.ref_epsilon(ref).GetCoefficient(var, num_elem, j);
 	    sigma = var.ref_sigma(ref).GetCoefficient(var, num_elem, j);            
	    mu = var.ref_mu(ref).GetCoefficient(var, num_elem, j);	    

            if (!var.linearize_drude)
              var.ref_drude(ref).ModifyCoefficientMaxwell(var.GetOmega(), epsilon, sigma, mu);
            
            if (var.InsidePML(num_elem))
	      var.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, j);
	    
	    if (mu.IsZero())
	      invMu.Zero();
	    else
	      GetInverse(mu, invMu); 
	    
            // jacob = Ji(xi_k), dfjm1 = DF_i^{-1}
            Real_wp jacob; Matrix3_3 dfjm1;
            if (variable)
              {
                if (affine)
                  {
                    jacob = var.Glob_jacobian(num_elem)(0);
                    dfjm1 = var.Glob_DFjm1(num_elem)(0);                
                    Mlt(1.0/jacob, dfjm1);
                  }
                else
                  {
                    jacob = var.Glob_jacobian(num_elem)(j)/Fb.WeightsND(j);
                    dfjm1 = var.Glob_DFjm1(num_elem)(j);
                    Mlt(1.0/jacob, dfjm1);
                  }
              }
            else
              {
                jacob = var.Glob_jacobian(num_elem)(0);
                dfjm1 = var.Glob_DFjm1(num_elem)(0);
                Mlt(1.0/jacob, dfjm1);
              }
            
            // dfj = DF_i, dfj_trans = DF_i^{*}
            Matrix3_3 dfj, dfj_trans;
            GetInverse(dfjm1, dfj); Transpose(dfj, dfj_trans);
            
            Real_wp poids_mass = Fb.WeightsND(j);
            if (N == 1)
              poids_mass = 1.0;
            
            // mass coefficient
            MltTrans(epsilon, dfjm1, A_tmp);
            Mlt(dfjm1, A_tmp, Bmass); 
            Mlt(momega2*jacob*poids_mass, Bmass);

	    mass.SetBh(j, Bmass);
            
            // stiffness coefficient
            Mlt(invMu, dfj, A_tmp);
            Mlt(dfj_trans, A_tmp, Astiff);
            Mlt(1.0/jacob*poids_mass, Astiff);

	    mass.SetAh(j, Astiff);
            
            // damping coefficient
            MltTrans(sigma, dfjm1, A_tmp);
            Mlt(dfjm1, A_tmp, Bmass);
            
            Bmass *= m_iomega*jacob*poids_mass;
	    mass.SetBhSigma(j, Bmass);
          }

	if (diag_mass)
	  {
	    mass.ClearBhSigma();

	    int Ndof = Fb.GetNbDof();
	    mass.ReallocateDhDiag(Ndof);
	    
	    epsilon = var.ref_epsilon(ref).GetConstant();
	    sigma = var.ref_sigma(ref).GetConstant();
	    if (var.ref_epsilon(ref).IsVarying() || var.ref_sigma(ref).IsVarying())
	      {
		cout << "Not implemented for variables indexes" << endl;
		abort();
	      }

            if (!var.linearize_drude)
              var.ref_drude(ref).ModifyCoefficientMaxwell(var.GetOmega(), epsilon, sigma, mu);
	    
	    Real_wp jacob; Matrix3_3 dfjm1;
	    jacob = var.Glob_jacobian(num_elem)(0);
	    dfjm1 = var.Glob_DFjm1(num_elem)(0);                
	    Mlt(1.0/jacob, dfjm1);
	    
	    MltTrans(epsilon, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Bmass);
	    Mlt(momega2*jacob, Bmass);

	    MltTrans(sigma, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Astiff);
	    Mlt(m_iomega*jacob, Astiff);

	    const VectReal_wp& weights_dof = Fb.WeightsDofND();
	    TinyVector<Complexe, 3> tau, tau_prime;
	    VectR3 s; SetPoints<Dimension3> PointsElem;
	    int num_pml = -1, type_pml = 0;
	    if (var.InsidePML(num_elem))
	      {
		num_pml = var.mesh.Element(num_elem).GetNumberPML();
                type_pml = var.mesh.Element(num_elem).GetTypePML();
		var.mesh.GetVerticesElement(num_elem, s);
		Fb.FjElemDof(s, PointsElem, var.mesh, num_elem);
	      }

	    for (int i = 0; i < Ndof; i++)
	      {
		int n = Fb.GetCoordinateDof(i);
		Complexe coef_pml; SetComplexOne(coef_pml);
		if (var.InsidePML(num_elem))
		  {
		    R3 point_glob = PointsElem.GetPointDof(Fb.GetPointNode_FromDofNumber(i));
		    var.GetDampingFactorPML(point_glob, num_pml, type_pml,
					    tau, tau_prime);
		    
		    Complexe dx = Real_wp(1) / tau(0);
		    Complexe dy = Real_wp(1) / tau(1);
		    Complexe dz = Real_wp(1) / tau(2);
		    switch(n)
		      {
		      case 0 : coef_pml = dy*dz / dx; break;
		      case 1 : coef_pml = dx*dz / dy; break;
		      case 2 : coef_pml = dx*dy / dz; break;
		      }			
		  }
		
		mass.SetDh(i, Bmass(n, n)*weights_dof(i)*coef_pml);
		mass.SetDhSigma(i, Astiff(n, n)*weights_dof(i)*coef_pml);
	      }
	  }
      }
    else
      {
	// only Bh and BhSigma are stored and jacobian matrices
	if (affine)
	  var.Glob_DFj(num_elem).Reallocate(1);
	else
	  var.Glob_DFj(num_elem).Reallocate(N);
		    
        if (diag_mass)
	  {
	    epsilon = var.ref_epsilon(ref).GetConstant();
	    sigma = var.ref_sigma(ref).GetConstant();
	    if (var.ref_epsilon(ref).IsVarying() || var.ref_sigma(ref).IsVarying())
	      {
		cout << "Not implemented for variables indexes" << endl;
		abort();
	      }
            
	    // diagonal case, only diagonal coefficients are stored
	    Real_wp jacob; Matrix3_3 dfjm1, dfj;
	    jacob = var.Glob_jacobian(num_elem)(0);
	    dfjm1 = var.Glob_DFjm1(num_elem)(0);                
	    Mlt(1.0/jacob, dfjm1);
	    
	    GetInverse(dfjm1, dfj);
	    var.Glob_DFj(num_elem)(0) = dfj;
	    
	    int Ndof = Fb.GetNbDof();
	    mass.ReallocateDhDiag(Ndof);

	    MltTrans(epsilon, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Bmass);
	    Mlt(m_iomega*jacob, Bmass);

	    MltTrans(sigma, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Astiff);
	    Mlt(jacob, Astiff);

	    const VectReal_wp& weights_dof = Fb.WeightsDofND();
	    TinyVector<Complexe, 3> tau, tau_prime; R3 tau_r;
	    VectR3 s; SetPoints<Dimension3> PointsElem;
	    int num_pml = -1, type_pml = 0;
	    if (var.InsidePML(num_elem))
	      {
		num_pml = var.mesh.Element(num_elem).GetNumberPML();
		type_pml = var.mesh.Element(num_elem).GetTypePML();
		var.mesh.GetVerticesElement(num_elem, s);
		Fb.FjElemDof(s, PointsElem, var.mesh, num_elem);
		mass.ReallocateTauPML(Ndof);
	      }

	    for (int i = 0; i < Ndof; i++)
	      {
		int n = Fb.GetCoordinateDof(i);
		if (var.InsidePML(num_elem))
		  {
		    R3 point_glob = PointsElem.GetPointDof(Fb.GetPointNode_FromDofNumber(i));
		    var.GetDampingFactorPML(point_glob, num_pml, type_pml,
					    tau, tau_prime);

		    tau_r.Init(realpart(tau(0)), realpart(tau(1)), realpart(tau(2)));
		    mass.SetTauPML(i, tau_r);
		  }
		
		mass.SetDh(i, Bmass(n, n)*weights_dof(i));
		mass.SetDhSigma(i, Astiff(n, n)*weights_dof(i));
	      }
	  }
	else
	  {
	    mass.ReallocateAhBh(0, N, N);
	    // loop on quadrature points
	    for (int j = 0; j < N; j++)
	      {            
		// physical coefficients
		epsilon = var.ref_epsilon(ref).GetCoefficient(var, num_elem, j);
		sigma = var.ref_sigma(ref).GetCoefficient(var, num_elem, j);
		
		mu = var.ref_mu(ref).GetCoefficient(var, num_elem, j);	    
		if (mu.IsZero())
		  invMu.Zero();
		else
		  GetInverse(mu, invMu); 
		
		// jacob = Ji(xi_k), dfjm1 = DF_i^{-1}
		Real_wp jacob; Matrix3_3 dfjm1;
		if (variable)
		  {
		    if (affine)
		      {
			jacob = var.Glob_jacobian(num_elem)(0);
			dfjm1 = var.Glob_DFjm1(num_elem)(0);                
			Mlt(1.0/jacob, dfjm1);
		      }
		    else
		      {
			jacob = var.Glob_jacobian(num_elem)(j)/Fb.WeightsND(j);
			dfjm1 = var.Glob_DFjm1(num_elem)(j);
			Mlt(1.0/jacob, dfjm1);
		      }
		  }
		else
		  {
		    jacob = var.Glob_jacobian(num_elem)(0);
		    dfjm1 = var.Glob_DFjm1(num_elem)(0);
		    Mlt(1.0/jacob, dfjm1);
		  }
		
		// dfj = DF_i, dfj_trans = DF_i^{*}
		Matrix3_3 dfj, dfj_trans;
		GetInverse(dfjm1, dfj); Transpose(dfj, dfj_trans);
		if ((j == 0) || (!affine))
		  var.Glob_DFj(num_elem)(j) = dfj;
		
		Real_wp poids_mass = Fb.WeightsND(j);
		if (N == 1)
		  poids_mass = 1.0;
		
		// mass coefficient
		MltTrans(epsilon, dfjm1, A_tmp);
		Mlt(dfjm1, A_tmp, Bmass); 
		Mlt(m_iomega*jacob*poids_mass, Bmass);
		
		mass.SetBh(j, Bmass);
		
		// damping coefficient
		MltTrans(sigma, dfjm1, A_tmp);
		Mlt(dfjm1, A_tmp, Bmass);
		
		Bmass *= jacob*poids_mass;
		mass.SetBhSigma(j, Bmass);
	      }
	  }
      }
  }
  
  
  //! penalization coefficient
  template<class T> template<class TypeEquation, class T0>
  void MaxwellEquation3D_Base<T>::
  GetPenalizationCoef(const EllipticProblem<TypeEquation>& var,
		      int i, const GlobalGenericMatrix<T0>& nat_mat, int ref, T0& coef)
  {
    coef = var.delta_penalization*nat_mat.GetCoefStiffness()
      /var.coefficient_impedance_absorbing(ref);
  }
  
  
  /***************************
   * HarmonicMaxwell_3D_Base *
   ***************************/


  //! reading of a line of the data file
  void HarmonicMaxwell_3D_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    
    if (description_field == "FileCoefficientsQ")
      {
        if (parameters.GetM() <= 0)
          {
            cout << "FileCoefficientsQ needs one parameter" <<endl;
            cout << "For example, you could write :" << endl;
            cout << "FileCoefficientsQ = fileQ.dat" << endl;
            cout << "Current parameters are " << parameters << endl;
            abort();
          }
        
        file_coefficient_q = parameters(0);
      }
    else if (description_field == "SymmetrizationDrude")
      {
        if (parameters(0) == "YES")
          this->use_symm_drude = true;
        else
          this->use_symm_drude = false;
      }
    else if (description_field == "LinearizationDrude")
      {
        if (parameters(0) == "YES")
          this->linearize_drude = true;
        else
          this->linearize_drude = false;
      }
  }
   

  void HarmonicMaxwell_3D_Base::ComputeMassMatrix(bool sipg)
  {    
    if (!sipg)
      return;
    
    Mesh<Dimension3>& mesh = var_problem.mesh;

    // storing jacobian matrix DFi
    this->Glob_DFj.Reallocate(mesh.GetNbElt());
    this->Glob_invJacobian.Reallocate(mesh.GetNbElt());
    Matrix3_3 dfjm1;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int nb_pts_quad = var_problem.GetNbPointsQuadratureInside(i);
        bool affine = mesh.IsElementAffine(i);    
        if (affine)
          {
            this->Glob_DFj(i).Reallocate(1);
            this->Glob_invJacobian(i).Reallocate(1);
            dfjm1 = var_problem.Glob_DFjm1(i)(0);
            Mlt(1.0/var_problem.Glob_jacobian(i)(0), dfjm1);
            GetInverse(dfjm1, this->Glob_DFj(i)(0));
            this->Glob_invJacobian(i)(0) = 1.0/Det(this->Glob_DFj(i)(0));
          }
        else
          {
            int N = var_problem.Glob_DFjm1(i).GetM() - nb_pts_quad;
            this->Glob_DFj(i).Reallocate(N);
            this->Glob_invJacobian(i).Reallocate(N);
            for (int j = 0; j < N; j++)
              {
                GetInverse(var_problem.Glob_DFjm1(i)(nb_pts_quad+j), this->Glob_DFj(i)(j));
                this->Glob_invJacobian(i)(j) = 1.0/Det(this->Glob_DFj(i)(j));
              }
          }
      }
    
#ifdef SELDON_WITH_MPI
    // for SIPG, we need to know the jacobian matrices from neighboring elements
    // therefore these matrices are exchanged
    // by using ExchangeUfaceDomains/GetUfaceDomains
    MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    int N = var_problem.GetNbPointsQuadratureNeighbor();
    TinyMatrix<Vector<Real_wp>, General, 3, 3> coef_dfjm1;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        coef_dfjm1(j, k).Reallocate(N);
    
    Vector<Vector<Real_wp> > g_send, g_recv;        
    Vector<Vector<int64_t> > g_send_tmp, g_recv_tmp;
    Vector<MPI_Request> request; int tag = 82;
    Real_wp jacob;
    // first step, jacobian matrices DF_i^-1 are collected
    // for all quadrature points of faces shared by different processors
    // these matrices are collected in the array coef_dfjm1
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int num_face = i;
        int ref = mesh.BoundaryRef(i).GetReference();
        int num_elem = mesh.BoundaryRef(i).numElement(0);
        int cond = mesh.GetBoundaryCondition(ref);
        if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            bool affine = mesh.IsElementAffine(num_elem);
            int offset_neighbor = var_problem.GetOffsetNeighboringFace(i);
            int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
            if (affine)
              {
                jacob = var_problem.Glob_jacobian(num_elem)(0);
                dfjm1 = var_problem.Glob_DFjm1(num_elem)(0);
                Mlt(1.0/jacob, dfjm1);
                for (int m = 0; m < nb_points_face; m++)
                  for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                      coef_dfjm1(j, k)(offset_neighbor+m) = dfjm1(j, k);
              }
            else
              {
                int num_loc = mesh.Element(num_elem).GetPositionBoundary(num_face);
                int offset = var_problem.GetNbPointsQuadratureInside(num_elem);
                for (int k = 0; k < num_loc; k++)
                  offset += mesh_num
                    .GetNbPointsQuadratureBoundary(mesh.Element(num_elem)
                                                   .numBoundary(k));
                
                for (int m = 0; m < nb_points_face; m++)
                  for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                      coef_dfjm1(j, k)(offset_neighbor+m)
                        = var_problem.Glob_DFjm1(num_elem)(offset+m)(j, k);
              }
          }
      }
    
    // second step
    // values are exchanged between faces with ExchangeUfaceDomains/GetUfaceDomains
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        {
          var_problem.ExchangeUfaceDomains(coef_dfjm1(j, k), g_send, g_send_tmp,
                                           g_recv, g_recv_tmp, request, tag);
          
          var_problem.GetUfaceDomains(coef_dfjm1(j, k), g_send, g_send_tmp,
                                      g_recv, g_recv_tmp, request, tag);
        }
    
    // last step, jacobian matrices are stored in the array Glob_DFjm1_Neighbor
    this->Glob_DFjm1_Neighbor.Reallocate(N);
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int i = 0; i < N; i++)
          this->Glob_DFjm1_Neighbor(i)(j, k) = coef_dfjm1(j, k)(i);
#endif
  }
  

  //! returns true if the finite element matrix is symmetric
  bool HarmonicMaxwell_3D_Base::IsSymmetricProblem(bool eigen) const
  {
    if (var_problem.FirstOrderFormulation())
      {
	if (eigen)
	  return false;
	
        if ((var_boundary.GetNbGlobalEltPML() > 0) && var_computation.LightStaticCondensation())
          return false;
        
	if (var_computation.GetLeafStaticCondensation())
	  return true;

        if (!this->use_symm_drude)
          if (this->GetNbVectorialDofDrudeAll() > 0)
            return false;        
	
	if ((var_boundary.GetNbGlobalEltPML() > 0) || (!var_computation.GetSymmetrizationUse()))
	  return false;
      }
    else
      {
	if (var_computation.GetLeafStaticCondensation())
	  return true;
        
        if (this->linearize_drude)
          if (this->GetNbVectorialDofDrudeAll() > 0)
            return false;
      }

#ifdef MONTJOIE_WITH_TRANSMISSION
    if (!var_transmission_base.IsSymmetricCondition())
      return false;
#endif
    
    return true;
  }


  void HarmonicMaxwell_3D_Base::ComputeDrudeDofs()
  {
    this->nb_dof_drude_vec = 0;
    if (var_problem.FormulationDG() != ElementReference_Base::CONTINUOUS)
      return;
    
    if (!var_problem.FirstOrderFormulation() && !this->linearize_drude)
      return;
  
    this->OffsetDofDrudeV.Reallocate(var_problem.mesh.GetNbElt()+1);
    this->OffsetDofDrudeV(0) = 0;
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
	int nb_dof_v = 0;
        int ref = var_problem.mesh.Element(i).GetReference();
        if (this->ref_drude(ref).IsEnabled())
          {
            nb_dof_v = 3*var_problem.GetNbPointsQuadratureInside(i);
            nb_dof_v *= this->ref_drude(ref).gamma.GetM();
            this->nb_dof_drude_vec += nb_dof_v;
	  }
        
	this->OffsetDofDrudeV(i+1) = this->OffsetDofDrudeV(i) + nb_dof_v;
      }

    this->nb_dof_drude_vec_all = this->nb_dof_drude_vec;
    //DISP(this->OffsetDofDrudeV); DISP(this->nb_dof_drude_vec);
  }
  
  
  //! case of H(curl) elements
  template<class T>
  void HarmonicMaxwell_3D_Base
  ::ComputeEnHnQuadratureHcurl(Vector<Vector<T> >& u_quadrature,
                               Vector<Vector<T> >& curl_quadrature,
                               int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                               bool compute_H, Vector<T>& En_quad, Vector<T>& Hn_quad) const
  {
    TinyVector<T, 3> Eloc, Hloc, En, Hn;
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(3*nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(3*nb_points_quad);
    
    T coef_miomega; var_problem.GetMiomega(coef_miomega);
    coef_miomega = -1.0/coef_miomega;
    for (int j = 0; j < nb_points_quad; j++)
      {
        Eloc.Init(u_quadrature(0)(j), u_quadrature(1)(j), u_quadrature(2)(j));
        TimesProd(Eloc, normale(j), En);            
        CopyVector(En, j, En_quad);
        
        if (compute_H)
          {
	    // multiplication by -i/omega to recover H
            Hloc.Init(curl_quadrature(0)(j), curl_quadrature(1)(j), curl_quadrature(2)(j));
	    Hloc *= coef_miomega;
            TimesProd(Hloc, normale(j), Hn);
            CopyVector(Hn, j, Hn_quad);
          }
      }
  }


  //! case of Discontinuous Galerkin
  template<class T>
  void HarmonicMaxwell_3D_Base
  ::ComputeEnHnQuadratureDG(Vector<Vector<T> >& u_quadrature,
                            Vector<Vector<T> >& grad_quadrature,
                            int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
                            bool compute_H, Vector<T>& En_quad, Vector<T>& Hn_quad) const
  {
    TinyVector<T, 3> Eloc, Hloc, En, Hn;
    int nb_points_quad = u_quadrature(0).GetM();
    En_quad.Reallocate(3*nb_points_quad);
    if (compute_H)
      Hn_quad.Reallocate(3*nb_points_quad);
    
    if (!var_problem.FirstOrderFormulationDG())
      {
	abort();
      }
    
    for (int j = 0; j < nb_points_quad; j++)
      {
        Eloc(0) = u_quadrature(0)(j);
        Eloc(1) = u_quadrature(1)(j);
        Eloc(2) = u_quadrature(2)(j);

        TimesProd(Eloc, normale(j), En);
        CopyVector(En, j, En_quad);
        
        if (compute_H)
          {
            Hloc(0) = u_quadrature(3)(j);
            Hloc(1) = u_quadrature(4)(j);
            Hloc(2) = u_quadrature(5)(j);
            
            TimesProd(Hloc, normale(j), Hn);
            CopyVector(Hn, j, Hn_quad);
          }
      }
  }

  
  //! computes E \times n and H \times n on nodal points
  template<class T>
  void HarmonicMaxwell_3D_Base
  ::ComputeEnHnNodal(Vector<Vector<T> >& u_nodal, Vector<Vector<T> >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<Vector<T> >& En_nodal, Vector<Vector<T> >& Hn_nodal) const
  {
    TinyVector<T, 3> Eloc, Hloc, En, Hn;
    int nb_points_nodal = u_nodal(0).GetM();
    En_nodal.Reallocate(3); Hn_nodal.Reallocate(3);
    for (int k = 0; k < 3; k++)
      {
	En_nodal(k).Reallocate(nb_points_nodal);
	Hn_nodal(k).Reallocate(nb_points_nodal);
      }
    
    if (var_problem.type_element == 2)
      {
	// H(curl) element
	for (int j = 0; j < nb_points_nodal; j++)
	  {
	    Eloc.Init(u_nodal(0)(j), u_nodal(1)(j), u_nodal(2)(j));
	    TimesProd(Eloc, normale(j), En);
	    En_nodal(0)(j) = En(0);
	    En_nodal(1)(j) = En(1);
	    En_nodal(2)(j) = En(2);
	    
	    Hloc.Init(grad_nodal(0)(j), grad_nodal(1)(j), grad_nodal(2)(j));
	    TimesProd(Hloc, normale(j), Hn);
	    Hn_nodal(0)(j) = Hn(0);
	    Hn_nodal(1)(j) = Hn(1);
	    Hn_nodal(2)(j) = Hn(2);
	  }
      }
    else
      {
	if (!var_problem.FirstOrderFormulationDG())
	  {
	    abort();
	  }
	
	for (int j = 0; j < nb_points_nodal; j++)
	  {
	    Eloc.Init(u_nodal(0)(j), u_nodal(1)(j), u_nodal(2)(j));
	    TimesProd(Eloc, normale(j), En);
	    En_nodal(0)(j) = En(0);
	    En_nodal(1)(j) = En(1);
	    En_nodal(2)(j) = En(2);
	    
	    Hloc.Init(u_nodal(3)(j), u_nodal(4)(j), u_nodal(5)(j));
	    TimesProd(Hloc, normale(j), Hn);
	    Hn_nodal(0)(j) = Hn(0);
	    Hn_nodal(1)(j) = Hn(1);
	    Hn_nodal(2)(j) = Hn(2);
	  }
      }    
  }


  //! computation of the number of degrees of freedom  
  void HarmonicMaxwell_3D_Base
  ::ComputeNumberOfDofs()
  {
    ComputeDrudeDofs();
    
    int nodl = var_problem.GetNbDof();
    int nodl_scalar = var_problem.offset_dof_unknown(1);
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    if (var_problem.FirstOrderFormulation())
      {   
	if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	  {
	    nodl += 6*mesh_num.GetNbDofPML();
	  }
	else
	  {
	    // dofs for unknown E* in PML layers
	    nodl += mesh_num.GetNbDofPML();
	    nodl_scalar += mesh_num.GetNbDofPML();
	    
	    // dofs due to magnetic field H and pml (unknown H* )
	    int nb_dof_pml = 0, offset_v = 0;
	    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
	      {
		int nb_dof_v = 0;
		if (var_problem.InsidePML(i))
		  nb_dof_v = (var_problem.GetOffsetDofV(i+1) - offset_v);
		
		offset_v = var_problem.GetOffsetDofV(i+1);
		var_problem.
		  SetOffsetDofV(i+1, offset_v + nb_dof_v + nb_dof_pml);
		
		nb_dof_pml += nb_dof_v;
	      }
	    
	    nodl += nb_dof_pml;
	  }
        
        if (var_problem.FormulationDG() != ElementReference_Base::HDG)
            nodl += 2*this->nb_dof_drude_vec;
      }
    else
      {
        if (var_problem.FormulationDG() != ElementReference_Base::HDG)
          if (this->linearize_drude)
            nodl += this->nb_dof_drude_vec;
      }

#ifdef MONTJOIE_WITH_TRANSMISSION
    // for transmission conditions, multiplier Lagrange may be added
    var_transmission_base.UpdateNumberOfDofs(nodl, nodl_scalar);
#endif

    //! additional dofs for wires
#ifdef MONTJOIE_WITH_WIRES
    this->IncrementNumberDofs_Wires();
#endif
    
    var_problem.SetNbDof(nodl);
    var_problem.offset_dof_unknown(1) = nodl_scalar;
    //var_problem.SetNbScalarDof(nodl_scalar);
  }
  
  
  //! Adds dofs due to pml for parallel computations
  void HarmonicMaxwell_3D_Base::PutOtherGlobalDofs()
  {
    MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0); 
    if (var_problem.FirstOrderFormulation())
      {
        int offset_loc = mesh_num.GetNbDof();
        int offset_glob = var_problem.GetNbGlobalMeshDof();
        
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
                    MatchingDof(num).Resize(nb_old + nb_dof_pml);
                    for (int j = 0; j < nb_old; j++)
                      {
                        int npml
                          = mesh_num.GetDofPML(MatchingDof(num)(j));
                        if (npml >= 0)
                          MatchingDof(num)(nb_dof++) = offset_loc + npml;
                      }
                  }
              }
          }        
      }

#ifdef MONTJOIE_WITH_TRANSMISSION
    var_transmission_base.UpdateGlobalDofs();
#endif

#ifdef MONTJOIE_WITH_MPI
    MPI_Allreduce(&this->nb_dof_drude_vec, &this->nb_dof_drude_vec_all, 1, MPI_INT, MPI_SUM, var_problem.comm_group_mode);
#else
    this->nb_dof_drude_vec_all = this->nb_dof_drude_vec;
#endif
  }
  

  //! fills all scalar dofs for element i
  IVect HarmonicMaxwell_3D_Base::GetScalarDofNumberOnElement(int i) const
  {
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    IVect Nodle;
    Nodle = mesh_num.number_map.GetDofNumberOnElement(mesh_num, i);
    
    if ((var_problem.FirstOrderFormulation()) && (var_boundary.GetNbGlobalEltPML() > 0))
      {
	int nb_dof = Nodle.GetM();
	bool close_pml = false;
	for (int i = 0; i < nb_dof; i++)
	  if ((Nodle(i) >= 0) && (mesh_num.GetDofPML(Nodle(i)) >= 0))
	    close_pml = true;
	
	if (close_pml)
	  {
	    Nodle.Resize(2*nb_dof);
	    int Nvol = mesh_num.GetNbDof();
	    for (int i = 0; i < nb_dof; i++)
	      {
		if ((Nodle(i) >= 0) && (mesh_num.GetDofPML(Nodle(i)) >= 0))
		  Nodle(nb_dof+i) = Nvol + mesh_num.GetDofPML(Nodle(i));
		else
		  Nodle(nb_dof + i) = -1;
	      }
	  }
      }
    
    return Nodle;
  }
  

  //! treating Dirichlet condition (retrieving dof numbers)  
  void HarmonicMaxwell_3D_Base::TreatDirichletCondition()
  {
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    if ((var_problem.FirstOrderFormulation()) && (var_boundary.GetNbGlobalEltPML() > 0))
      {
	int nb_dof = 0;
	for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	  {
	    int j = var_boundary.GetDirichletDofNumber(i);
	    if (mesh_num.GetDofPML(j) >= 0)
	      nb_dof++;
	  }

	int nb_dof_all = nb_dof;
#ifdef SELDON_WITH_MPI
	MPI_Allreduce(&nb_dof, &nb_dof_all, 1, 
		      MPI_INTEGER, MPI_SUM, var_problem.comm_group_mode);
#endif
	
	if (nb_dof_all > 0)
	  {
	    int Nvol = mesh_num.GetNbDof();
	    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
	      {
		int j = var_boundary.GetDirichletDofNumber(i);
		if (mesh_num.GetDofPML(j) >= 0)
		  {
		    int k = mesh_num.GetDofPML(j);
		    var_boundary.SetDirichletDof(k+Nvol, true);
		  }
	      }
	    
	    var_boundary.UpdateDirichletDofs();
	  }		
      }
  }


  //! writing of the solution on all types of outputs (plane, lines, points ...)
  void HarmonicMaxwell_3D_Base::WriteCoefQ(Vector<VectComplex_wp> & U0)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    //const IVect& ddl_plus = var_transmission_base.GetDofNumberOmegaPlus();
    /*if (ddl_plus.GetM() > 0)
      {
        const IVect& ddl_minus = var_transmission_base.GetDofNumberOmegaMinus();
        int N = ddl_plus.GetM(); int Nvol = var_problem.mesh_num.GetNbDof();
        VectComplex_wp Lambda(Nvol);
        Lambda.Zero();

        //DISP(ddl_minus); DISP(ddl_plus); DISP(Nvol);
        if (U0(0).GetM() > Nvol)
          {
            for (int i = 0; i < N; i++)
              {
                Lambda(ddl_plus(i)) = U0(0)(Nvol + i);
                Lambda(ddl_minus(i)) = U0(0)(Nvol + i);
              }
            
            var_output.WriteOutputFile(Lambda, "Lambda.dat");
          }
      }
    */

    if (file_coefficient_q.size() > 0)
      {
        VectReal_wp coefQ;
        ComputeCoefficientsQ(U0(0), coefQ);
        
        if (nb_proc == 1)
          coefQ.WriteText(file_coefficient_q);
        else
          {
#ifdef SELDON_WITH_MPI
            Vector<int64_t> Qtmp;
            if (rank_proc != 0)
              {
                int nb_elt = coefQ.GetM();
                MPI_Send(&nb_elt, 1, MPI_INTEGER, 0, 44, var_problem.comm_group_mode);
                MPI_Send(var_problem.mesh.GlobElementNumber_Subdomain.GetData(),
                          nb_elt, MPI_INTEGER, 0, 45, var_problem.comm_group_mode);
                
                MpiSend(var_problem.comm_group_mode, coefQ, Qtmp, nb_elt, 0, 46);
              }
            else
              {
                // first receiving the number of elements on each processor
                IVect nb_elt_proc(nb_proc);
                nb_elt_proc(0) = coefQ.GetM();
                int nb_elt_all = nb_elt_proc(0);
                MPI_Status status;
                for (int i = 1; i < nb_elt_proc.GetM(); i++)
                  {
                    MPI_Recv(&nb_elt_proc(i), 1,
                             MPI_INTEGER, i, 44, var_problem.comm_group_mode, &status);
                    
                    nb_elt_all += nb_elt_proc(i);
                  }
                
                VectReal_wp coefQ_glob(nb_elt_all);
                coefQ_glob.Fill(-1);
                for (int i = 0; i < coefQ.GetM(); i++)
                  coefQ_glob(var_problem.mesh.GlobElementNumber_Subdomain(i)) = coefQ(i);
                
                IVect NumElt;
                for (int i = 1; i < nb_elt_proc.GetM(); i++)
                  {
                    NumElt.Reallocate(nb_elt_proc(i));
                    coefQ.Reallocate(nb_elt_proc(i));
                    MPI_Recv(NumElt.GetData(), nb_elt_proc(i),
                             MPI_INTEGER, i, 45, var_problem.comm_group_mode, &status);

                    MpiRecv(var_problem.comm_group_mode, coefQ, Qtmp, nb_elt_proc(i),
                            i, 46, status);
                    
                    for (int j = 0; j < NumElt.GetM(); j++)
                      coefQ_glob(NumElt(j)) = coefQ(j);
                  }
                
                // then writing the coefficients
                coefQ_glob.WriteText(file_coefficient_q);
              }
#endif
          }
      }
  }


  //! computation of elementary matrix for edge finite elements
  /*!
    basic quadrature (loop on all quadrature points)
    \param[in] iquad element number
    \param[in] B mass coefficients 
    \param[in] A stiffness coefficients
    \param[in,out] mat matrix to fill
    \param[in] Fb finite element object
    \param[in] nat_mat additional coefficients
    \param[in] vars given problem
    computation of \f$ -\omega^2 \int_K \varepsilon \varphi_i \varphi_j
    + \int_K \mu^{-1} \nabla \times \varphi_i \nabla \times \varphi_j \f$
    After change of variables from K to \f$ \hat{K} \f$, we get :
    \f$ -\omega^2 \int_{\hat{K}} J_i DF_i^{-1} \varepsilon DF_i^{*-1} \varphi_i \varphi_j 
    + \int_{\hat{K}} J_i DF_i^* \mu^{-1} DF_i \nabla \times \varphi_i \nabla \times \varphi_j \f$
    mass coefficient (3x3 matrix) \f$ B =  J_i DF_i^{-1} \varepsilon DF_i^{*-1} \f$
    stiffness coefficient (3x3 matrix) \f$ A =  J_i DF_i^* \mu^{-1} DF_i \f$ 
   */
  template<class Vect, class Matrix1, class Complexe>
  void HarmonicMaxwell_3D_Base
  ::ComputeVariableElementaryMatrix(int iquad, const Vect& B, const Vect& A, bool compute_mass,
				    Matrix1& mat, const ElementReference<Dimension3, 2>& Fb,
				    const GlobalGenericMatrix<Complexe>& nat_mat)
    
  {
    R3 phi, phi2; int nb_dof_elt = Fb.GetNbDof();
    int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
    
    VectReal_wp val_phi(3*nb_points_quadrature), Ones(nb_dof_elt);
    Vector<Complexe> feval(3*nb_points_quadrature), feval_bis(3*nb_points_quadrature);
    Vector<Complexe> feval_curl(3*nb_points_quadrature);
    VectReal_wp curl_phi(3*nb_points_quadrature);
    Vector<Complexe> contrib(nb_dof_elt), contrib_curl(nb_dof_elt);
    TinyVector<Complexe, 3> vec_u, vec_v;
    val_phi.Fill(0); Ones.Fill(0); feval.Fill(0);
    feval_bis.Fill(0); feval_curl.Fill(0); curl_phi.Fill(0);
    contrib.Fill(0); contrib_curl.Fill(0);
    Complexe one; SetComplexOne(one);
    
    // loop on dofs and quadrature points
    for (int i = 0; i < nb_dof_elt; i++)
      {
        // we compute phi_i(\xi_k) and curl phi_i for all k
        Ones.Fill(0); Ones(i) = 1.0;
        Fb.ApplyChTranspose(Ones, val_phi);
        if (Fb.UseQuadraturePointsForRh())
          Fb.ApplyRhQuadratureTranspose(val_phi, curl_phi);
        else
          Fb.ApplyRhTranspose(Ones, curl_phi);
        
        // application of A and B
        for (int k = 0; k < nb_points_quadrature; k++)
          {
            if (compute_mass)
	      {
		vec_u(0) = val_phi(3*k); vec_u(1) = val_phi(3*k+1); vec_u(2) = val_phi(3*k+2);
		
		Mlt(B(k), vec_u, vec_v);
		
		feval(3*k) = vec_v(0); feval(3*k+1) = vec_v(1); feval(3*k+2) = vec_v(2);
	      }
            
            vec_u(0) = curl_phi(3*k); vec_u(1) = curl_phi(3*k+1); vec_u(2) = curl_phi(3*k+2);
            
            Mlt(A(k), vec_u, vec_v);
            
            feval_curl(3*k) = vec_v(0); feval_curl(3*k+1) = vec_v(1); feval_curl(3*k+2) = vec_v(2);
          }

        // then integration against \nabla phi_j
        if (Fb.UseQuadraturePointsForRh())
          {
            Fb.ApplyRhQuadrature(feval_curl, feval_bis);
            if (compute_mass)
	      Add(one, feval_bis, feval);
	    else
	      feval = feval_bis;
	    
            Fb.ApplyCh(feval, contrib);
          }
        else
          {
            Fb.ApplyRh(feval_curl, contrib_curl);
            if (compute_mass)
	      {
		Fb.ApplyCh(feval, contrib);
		Add(one, contrib_curl, contrib);
	      }
	    else
	      contrib = contrib_curl;
          }
        
        // we got a row of the matrix
        for (int j = 0; j < nb_dof_elt; j++)
          mat.SetEntry(i, j, contrib(j));
      }
  }
  

  void HarmonicMaxwell_3D_Base::FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const
  {
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    i0 = 0; i1 = 0; j0 = 0; j1 = 0;
    if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	i0 = mesh_num.GetNbDof();
	i1 = i0 + var_problem.GetOffsetDofV(var_problem.mesh.GetNbElt());
        bool mode_TE = true;
        for (int ref = 1; ref < this->ref_drude.GetM(); ref++)
          if (ref_drude(ref).IsEnabled())
            {
              if (ref_drude(ref).IsModeTE())
                mode_TE = true;
              else
                mode_TE = false;
            }

        if (this->GetNbVectorialDofDrudeAll() > 0)
          {
            if (!mode_TE)
              {
                j0 = var_problem.GetNbDof() - this->GetNbVectorialDofDrude();
                j1 = var_problem.GetNbDof();
              }
            else
              {
                j0 = var_problem.GetNbDof() - 2*this->GetNbVectorialDofDrude();
                j1 = var_problem.GetNbDof() - this->GetNbVectorialDofDrude();
              }
          }
      }
    else if (var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
	i0 = 3*mesh_num.GetNbDof();
	i1 = i0 + 3*mesh_num.GetNbDof();
      }
    else
      {
	i0 = 0;
	i1 = mesh_num.GetNbDof();
	int Nvol = var_problem.GetOffsetDofV(var_problem.mesh.GetNbElt());
	j0 = i1 + Nvol;
	j1 = j0 + Nvol;
      }
  }
  
  
  void HarmonicMaxwell_3D_Base::GetInternalNodesElement(int i, int nb_dof_loc,
                                                        int& nb_dof_edges, int& nb_dof_int,
                                                        Vector<int>& intern_node) const
  {
    int nb_u = 1;
    if (var_problem.InsidePML(i))
      nb_u = 2;
    
    int nb_dof_elt = var_problem.GetNbLocalDof(i);
    intern_node.Reallocate(nb_dof_loc);
    intern_node.Fill(-1);
    
    nb_dof_edges = nb_u*nb_dof_elt;
    nb_dof_int = nb_dof_loc - nb_dof_edges;

    for (int i = 0; i < nb_dof_edges; i++)
      intern_node(i) = i;
    
    for (int i = nb_dof_edges; i < nb_dof_loc; i++)
      intern_node(i) = -(i-nb_dof_edges)-1;
  }
  
  
  /**********************************
   * HarmonicMaxwell3D_PhysGeomInfo *
   **********************************/

  
  template<class Complexe>
  HarmonicMaxwell3D_PhysGeomInfo<Complexe>::HarmonicMaxwell3D_PhysGeomInfo()
  {
    offset_Ah = 0;
    offset_Bh = 0;
    offset_BhSigma = 0;
    offset_diagDhSigma = 0;
    offset_diagDhStiff = 0;
  }    
  

  template<class Complexe>
  size_t HarmonicMaxwell3D_PhysGeomInfo<Complexe>::GetMemorySize() const
  {
    size_t taille;
    taille = Ah_Bh.GetMemorySize() + diag_Dh.GetMemorySize() + tauPML.GetMemorySize();
    return taille;
  }
  

  template<class Complexe>
  void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::ReallocateAhBh(int N1, int N2, int N3)
  {
    Ah_Bh.Reallocate(N1 + N2 + N3);
    offset_Ah = 0;
    offset_Bh = N1;
    offset_BhSigma = N1+N2;
  }


  template<class Complexe>
  void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::ReallocateBhSurf(int N)
  {
    BhSurf.Reallocate(N);
  }
  
  
  template<class Complexe>
  void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::ReallocateDhDiag(int N, int Ns)
  {
    diag_Dh.Reallocate(2*N + Ns);
    offset_diagDhSigma = N;
    offset_diagDhStiff = 2*N;
  }
  

  template<class Complexe>
  void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::ReallocateTauPML(int N)
  {
    tauPML.Reallocate(N);
  }
  

  template<class Complexe>
  void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::ClearBhSigma()
  {
    Ah_Bh.Resize(offset_Bh); offset_BhSigma = offset_Bh;
  }
    

  template<class Complexe>
  void HarmonicMaxwell3D_PhysGeomInfo<Complexe>::ClearDiagonalDh()
  {
    diag_Dh.Clear(); tauPML.Clear();
  }

  
  /**********************
   * HarmonicMaxwell_3D *
   **********************/
  
  
  //! adds size used by different objects in VarHarmonic
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>
  ::GetMemoryUsed(map<string, size_t>& var) const
  {    
    size_t taille = sizeof(*this) + Seldon::GetMemorySize(this->Glob_DFj)
      + this->Glob_DFjm1_Neighbor.GetMemorySize() + Seldon::GetMemorySize(Glob_invJacobian)
      + coefficient_impedance_absorbing.GetMemorySize();

    for (int i = 0; i < Glob_matMass_elem.GetM(); i++)
      taille += Glob_matMass_elem(i).GetMemorySize();
    
    for (int ref = 0; ref < ref_epsilon.GetM(); ref++)
      taille += ref_epsilon(ref).GetMemorySize();

    for (int ref = 0; ref < ref_invEpsilon.GetM(); ref++)
      taille += ref_invEpsilon(ref).GetMemorySize();

    for (int ref = 0; ref < ref_sigma.GetM(); ref++)
      taille += ref_sigma(ref).GetMemorySize();

    for (int ref = 0; ref < ref_mu.GetM(); ref++)
      taille += ref_mu(ref).GetMemorySize();

    for (int ref = 0; ref < ref_invMu.GetM(); ref++)
      taille += ref_invMu(ref).GetMemorySize();

    var["MaxwellData"] = taille;
        
  }
  

  //! copies input data
  template<class T>
  void HarmonicMaxwell_3D<T>::CopyInputData(const HarmonicMaxwell_3D<T>& var)
  {
    this->ref_drude = var.ref_drude;
    this->linearize_drude = var.linearize_drude;
    
    ref_epsilon = var.ref_epsilon;
    ref_invEpsilon = var.ref_invEpsilon;
    ref_sigma = var.ref_sigma;
    ref_mu = var.ref_mu;
    ref_invMu = var.ref_invMu;
    epsilon0 = var.epsilon0;
    mu0 = var.mu0;
    invMu0 = var.invMu0;
  }
  

  //! returns mass coefficients involved in mass matrix
  template<class T> template<class T0, class Prop>
  void HarmonicMaxwell_3D<T>
  ::GetMassMatrix(int num_elem, int j, const GlobalGenericMatrix<T0>& nat_mat,
		  TinyMatrix<T0, Prop, 3, 3>& Bmass) const
  {
    Bmass = this->Glob_matMass_elem(num_elem).GetBh(j)*nat_mat.GetCoefMass();
    Bmass += this->Glob_matMass_elem(num_elem).GetBhSigma(j)*nat_mat.GetCoefDamping();
  }
  

  //! returns mass coefficients involved in stiffness matrix
  template<class T> template<class T0, class Prop>
  void HarmonicMaxwell_3D<T>
  ::GetStiffMatrix(int num_elem, int j, const GlobalGenericMatrix<T0>& nat_mat,
                   TinyMatrix<T0, Prop, 3, 3>& Astiff) const
  {
    Astiff = this->Glob_matMass_elem(num_elem).GetAh(j)*nat_mat.GetCoefStiffness();
  }


  //! epsilon and mu in case of PML layer
  /*!
    \param[out] epsilon dielectric permittivity
    \param[out] mu magnetic permeability
  */
  template<>
  void HarmonicMaxwell_3D<Complex_wp>
  ::ModifyPhysicalCoefPML(Matrix3_3sym_Complex_wp& epsilon,
			  Matrix3_3sym_Complex_wp& mu, Matrix3_3sym_Complex_wp& sigma,
			  int i1, int i) const
  {
    R3_Complex_wp coef;
    
    Complex_wp dx = 1.0/var_boundary.GetTauPML(i1, i, 0);
    Complex_wp dy = 1.0/var_boundary.GetTauPML(i1, i, 1);
    Complex_wp dz = 1.0/var_boundary.GetTauPML(i1, i, 2);
    
    coef(0) = dy*dz/dx;
    coef(1) = dx*dz/dy;
    coef(2) = dx*dy/dz;

    epsilon(0, 0) *= coef(0);
    epsilon(1, 1) *= coef(1);
    epsilon(2, 2) *= coef(2);

    sigma(0, 0) *= coef(0);
    sigma(1, 1) *= coef(1);
    sigma(2, 2) *= coef(2);

    mu(0, 0) *= coef(0);
    mu(1, 1) *= coef(1);
    mu(2, 2) *= coef(2);
  }
  

  template<>
  void HarmonicMaxwell_3D<Real_wp>
  ::ModifyPhysicalCoefPML(Matrix3_3sym& epsilon,
			  Matrix3_3sym& mu, Matrix3_3sym& sigma,
			  int i1, int i) const
  {
  }
  
      
  //! initialization of indices
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>::InitIndices(int n)
  {
    if (n <= 1)
      {
        cout<<" ATTENTION VOUS N'AVEZ MIS QU'UN SEUL MILIEU 0!!!!!!!!!"<<endl;
        cout<<"Number of domains should be greater than 2"<<endl; abort();
      }
    
    if (PhysicalConstant::adimensionalization != PhysicalConstant::ADIM_NO)
      {
	// non-physical units, epsilon_0 and mu_0 set to 1
	epsilon0 = 1.0; mu0 = 1.0; invMu0 = 1.0;
      }
    else
      {
	// use of physical units
	epsilon0 = PhysicalConstant::epsilon0_permittivity;
	        
        mu0 = PhysicalConstant::mu0_permeability;
        invMu0 = 1.0/mu0;
      }
    
    ref_epsilon.Reallocate(n);
    ref_sigma.Reallocate(n);
    ref_invEpsilon.Reallocate(n);
    ref_mu.Reallocate(n);
    ref_invMu.Reallocate(n);
    this->ref_drude.Reallocate(n);
    
    TinyMatrix<Complexe, Symmetric, 3, 3> zero;
    for (int i = 0; i < n; i++)
      {
	ref_epsilon(i).SetDiagonal(epsilon0);
	ref_invEpsilon(i).SetDiagonal(1.0/epsilon0);
        ref_sigma(i).SetConstant(zero);
	ref_mu(i).SetDiagonal(mu0);
	ref_invMu(i).SetDiagonal(invMu0);
      }
  }
  
  
  //! returns the number of physical media
  template<class Complexe>
  int HarmonicMaxwell_3D<Complexe>::GetNbPhysicalIndices() const
  {
    return ref_epsilon.GetM();
  }

  
  //! reading of data file to modify indices
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>
  ::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_epsilon.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_epsilon.GetM() << endl;
        abort();
      }

    ref_epsilon(i).SetInputData(nb, parameters, parameters(0));
    
    ref_mu(i).SetInputData(nb, parameters, parameters(0));
    
    ref_sigma(i).SetInputData(nb, parameters, parameters(0));
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      ref_sigma(i).Mlt(PhysicalConstant::impedance0);
    
    // multiplication of relative physical indices by epsilon_0, mu_0
    ref_epsilon(i).Mlt(epsilon0);
    ref_mu(i).Mlt(mu0);
  }

  
  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_epsilon.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_epsilon.GetM() << endl;
        abort();
      }

    if (name_media == "epsilon")
      {
        ref_epsilon(i).SetInputData(nb, parameters, parameters(0));
        // multiplication of relative physical indices by epsilon_0, mu_0
        ref_epsilon(i).Mlt(epsilon0);
      }
    else if (name_media == "mu")
      {
        ref_mu(i).SetInputData(nb, parameters, parameters(0));
        ref_mu(i).Mlt(mu0);
      }
    else if (name_media == "sigma")
      {
        ref_sigma(i).SetInputData(nb, parameters, parameters(0));
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          ref_sigma(i).Mlt(PhysicalConstant::impedance0);
      }
    else if (name_media == "Drude")
      {
	this->ref_drude(i).SetInputData(nb, parameters, parameters(0));
	this->ref_drude(i).Adimensionalize(var_problem.GetWaveLengthAdim());
        TinyMatrix<Complexe, Symmetric, 3, 3> eps, mu, sigma;
        eps.SetIdentity(); mu.SetIdentity(); sigma.Zero();

        if (this->ref_drude(i).IsModeTM())
          eps.SetDiagonal(this->ref_drude(i).eps_inf);
        else
          mu.SetDiagonal(this->ref_drude(i).eps_inf);
        
        this->ref_epsilon(i).SetConstant(eps);
        this->ref_mu(i).SetConstant(mu);
        this->ref_sigma(i).SetConstant(sigma);
      }
    else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
  }
  

  //! returns the name associated with the physical index num
  template<class Complexe>
  string HarmonicMaxwell_3D<Complexe>::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("epsilon");
      case 1: return string("mu");
      case 2: return string("sigma");
      }

    return string();
  }  
    
    
  //! fills the arrays rho_complex, rho_real with varying indices of the current problem
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>::
  GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension3, Complexe>* >& rho_complex,
                    IVect& num_ref, IVect& num_index, IVect& num_component,
                    Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    int nb = 0;
    for (int i = 0; i < ref_epsilon.GetM(); i++)
      {
        nb += ref_epsilon(i).GetNbVaryingMedia();
        nb += ref_mu(i).GetNbVaryingMedia();
        nb += ref_sigma(i).GetNbVaryingMedia();
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
    for (int i = 0; i < ref_epsilon.GetM(); i++)
      {
        int nb0 = nb;
        ref_epsilon(i).GetVaryingMedia(nb, rho_complex, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }
        
        nb0 = nb;
        ref_mu(i).GetVaryingMedia(nb, rho_complex, num_component);
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
  
  
  //! returns the velocity of waves in physical media whose reference is ref
  template<class Complexe>
  Real_wp HarmonicMaxwell_3D<Complexe>::GetVelocityOfMedia(int ref) const
  {
    return 1.0/sqrt(abs(this->ref_epsilon(ref)(0, 0)*this->ref_mu(ref)(0, 0)));
  }
  
  
  //! returns the velocity of waves at infinity
  template<class Complexe>
  Real_wp HarmonicMaxwell_3D<Complexe>::GetVelocityOfInfinity() const
  {
    return 1.0/sqrt(abs(this->epsilon0*this->mu0));
  }


  //! computation of invEpsilon and invMu from values of Epsilon and Mu
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>::FinalizeComputationVaryingIndices()
  {
    // computation of inverses    
    for (int ref = 1; ref < ref_mu.GetM(); ref++)
      {
        ref_invEpsilon(ref) = ref_epsilon(ref);
        ref_invEpsilon(ref).GetInverse();

        ref_invMu(ref) = ref_mu(ref);
        ref_invMu(ref).GetInverse();
      }
  }

  
  //! returns true if physical indices are varying inside element i
  template<class Complexe>
  bool HarmonicMaxwell_3D<Complexe>::IsVaryingMedia(int i) const
  {
    if (ref_epsilon(i).IsVarying()||ref_mu(i).IsVarying()||ref_sigma(i).IsVarying())
      return true;
    
    return false;
  }
  

  //! allocation of arrays before computation of mass matrix
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>::AllocateMassMatrices()
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif

    // allocation of arrays needed to compute the finite element matrix
    this->Glob_matMass_elem.Reallocate(var_problem.mesh.GetNbElt());
    
    if (var_problem.FirstOrderFormulation())
      this->Glob_DFj.Reallocate(var_problem.mesh.GetNbElt());
    
    // impedance of physical material
    this->coefficient_impedance_absorbing.Reallocate(this->GetNbPhysicalIndices());
    this->coefficient_impedance_absorbing.Fill(1.0);
    for (int i = 0; i < this->coefficient_impedance_absorbing.GetM(); i++)
      {
	Complexe epsilon = this->ref_epsilon(i)(0,0), mu = this->ref_mu(i)(0,0);
	if ((var_problem.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
	    && var_problem.FirstOrderFormulationDG())
	  this->coefficient_impedance_absorbing(i) = sqrt(mu/epsilon);
	else
	  this->coefficient_impedance_absorbing(i) = sqrt(epsilon/mu);
      }
    
    Mlt(var_boundary.GetImpedanceCoefficientABC(), this->coefficient_impedance_absorbing);

    if ((var_problem.print_level >= 6) && (nb_proc == 1))
      {
        Vector<bool> RefUsed(this->GetNbPhysicalIndices()+1);
        RefUsed.Fill(false);
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          RefUsed(var_problem.mesh.Element(i).GetReference()) = true;
        
        cout << "Frequency = " << var_problem.GetOmega()/(2.0*pi_wp) << endl;
        cout << "Pulsation = " << var_problem.GetOmega() << endl;
        cout << "Wave vector = " << var_problem.GetWaveVector() << endl;
        for (int ref = 1; ref < RefUsed.GetM(); ref++)
          if (RefUsed(ref))
            {
              TinyMatrix<Complexe, Symmetric, 3, 3> epsilon, mu, sigma;
              epsilon = this->ref_epsilon(ref).GetConstant();
              mu = this->ref_mu(ref).GetConstant();
              sigma = this->ref_sigma(ref).GetConstant();
              this->ref_drude(ref).ModifyCoefficientMaxwell(var_problem.GetOmega(), epsilon, sigma, mu);
              
              cout << "Physical Media " << ref << endl;
              cout << "Value of mu = " << mu(0, 0) << endl;
              cout << "Value of epsilon = " << epsilon(0, 0) << endl;
              cout << "Value of sigma = " << sigma(0, 0) << endl;
              cout << "Value of omega^2 epsilon mu = " << 
                var_problem.GetSquareOmega()*mu(0, 0)*epsilon(0, 0) << endl;
            }
      }    
  }


  //! computes coefficients Q_i = \int sigma |E|^2 dx on each element i
  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>
  ::ComputeCoefficientsQ(const VectComplex_wp& U0, VectReal_wp& coefQ) const
  {
    
    // computes values of E on quadrature points
    Vector<Vector<Vector<Complex_wp> > > Equad;
    Vector<Vector<Vector<Complex_wp> > > Hquad;
    Vector<VectComplex_wp> Uvec;
    Uvec.SetData(1, const_cast<VectComplex_wp*>(&U0));
    var_output.ComputeQuadratureUgradU(Uvec, Equad, Hquad, true, false, true);
    
    // computes integrals
    TinyVector<Complex_wp, 3> En;
    coefQ.Reallocate(var_problem.mesh.GetNbElt());
    coefQ.Fill(0);
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
        const VectReal_wp& omega = var_problem.WeightsND(i);
        bool affine = (var_problem.Glob_jacobian(i).GetM() == 1);
        
        int ref = var_problem.mesh.Element(i).GetReference();
        Real_wp val = 0.0, jacob; Complexe sigma;
        for (int j = 0; j < Equad(i)(0).GetM(); j++)
          {
            if (affine)
              jacob = var_problem.Glob_jacobian(i)(0)*omega(j);
            else
              jacob = var_problem.Glob_jacobian(i)(j);
            
	    En.Init(Equad(i)(0)(j), Equad(i)(1)(j), Equad(i)(2)(j));
            sigma = this->ref_sigma(ref).GetCoefficient(var_problem, i, j)(0, 0);
            val += abs(sigma)*jacob*AbsSquare(En);
          }
        
        if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
          val /= PhysicalConstant::impedance0;
        
        coefQ(i) = val;
      }    
    
    Uvec.Nullify();
  }
    

  //! computation of elementary matrix for edge finite element
  /*!
    computation of \f$ -\alpha*\omega^2 \int_K \varepsilon \varphi_i \varphi_j
    + \beta*\int_K \mu^{-1} \nabla \times \varphi_i \nabla \times \varphi_j \f$
    \param[in] iquad element number
    \param[out] num_dof dof numbers
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat coefficients alpha and beta
    \param[in] vars given problem
    \param[in] Fb finite element object
   */
  template<class Complexe> template<class T>
  void HarmonicMaxwell_3D<Complexe>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof,
 			    VirtualMatrix<T>& mat_interac,
			    const GlobalGenericMatrix<T>& nat_mat,
			    const ElementReference<Dimension3, 2>& Fb)
  {
    // int num_loc,order,order_m1,order_p1;
    const Mesh<Dimension3>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    
    // dof numbers for E
    num_dof.Reallocate(Fb.GetNbDof());
    IVect Nodle = var_problem.GetDofNumberOnElement(iquad);
    for (int j = 0; j < Fb.GetNbDof(); j++)
      num_dof(j) = Nodle(j);
    
    // number of integration points
    int Nquad = Fb.GetNbPointsQuadratureInside();
    int ref_domain = var_problem.mesh.Element(iquad).GetReference();
    HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass = this->Glob_matMass_elem(iquad);
	
    // unknowns E, E*, H, H*
    // local offsets of these unknowns
    int offset_Eloc_s(0), offset_Hloc(0), offset_Hloc_s(0);
    int offset_Es = mesh_num.GetNbDof();
    int offset_H = mesh_num.GetNbDof() + mesh_num.GetNbDofPML();
    int offset_Ploc = 0, offset_Qloc = 0;
    
    // dof numbers 
    int nb_unknowns = 1, nb_dof_all = Fb.GetNbDof();    
    if (var_problem.FirstOrderFormulation())
      {
        if (this->ref_drude(ref_domain).IsEnabled())
          {
            if (this->ref_drude(ref_domain).IsModeTM())
              {
                if (abs(this->ref_drude(ref_domain).eps_inf - this->ref_epsilon(ref_domain).GetConstant()(0,0)) > epsilon_machine)
                  {
                    cout << "Epsilon must be equal to eps_infini" << endl;
                    abort();
                  }
              }
            else
              {
                if (abs(this->ref_drude(ref_domain).eps_inf - this->ref_mu(ref_domain).GetConstant()(0,0)) > epsilon_machine)
                  {
                    cout << "Mu must be equal to eps_infini" << endl;
                    abort();
                  }
              }
          }
        
	if (var_computation.GetLeafStaticCondensation())
	  {
            if (var_computation.LightStaticCondensation() && var_problem.InsidePML(iquad))
              {
                nb_unknowns = 2;
                nb_dof_all *= 2;
                offset_Eloc_s = Fb.GetNbDof();
                num_dof.Resize(nb_dof_all);

		// dof numbers for E*
		for (int j = 0; j < Fb.GetNbDof(); j++)
		  {
                    int n = mesh_num.GetDofPML(Nodle(j));
                    if (n < 0)
                      num_dof(j + offset_Eloc_s) = -1;
                    else
                      num_dof(j + offset_Eloc_s) = offset_Es + mesh_num.GetDofPML(Nodle(j));
                  }                
              }
	  }
	else
	  {	    
	    nb_dof_all += 3*Nquad;
	    if (var_problem.InsidePML(iquad))
	      {
		nb_unknowns = 2;
		nb_dof_all *= 2;
		offset_Eloc_s = Fb.GetNbDof();
		offset_Hloc = 2*Fb.GetNbDof();
		offset_Hloc_s = offset_Hloc + 3*Nquad;
	      }
	    else
	      {
		offset_Hloc = Fb.GetNbDof();
		offset_Eloc_s = 0;
		offset_Hloc_s = offset_Hloc;
	      }
	    
	    offset_H += var_problem.GetOffsetDofV(iquad);
	    
	    num_dof.Resize(nb_dof_all);
	    // dof numbers for H
	    for (int j = 0; j < 3*Nquad; j++)
	      num_dof(j + offset_Hloc) = offset_H + j;
	    
	    if (var_problem.InsidePML(iquad))
	      {
		// dof numbers for E* and H*
		for (int j = 0; j < Fb.GetNbDof(); j++)
		  {
                    int n = mesh_num.GetDofPML(Nodle(j));
                    if (n < 0)
                      num_dof(j + offset_Eloc_s) = -1;
                    else
                      num_dof(j + offset_Eloc_s) = offset_Es + mesh_num.GetDofPML(Nodle(j));
                  }
                
		for (int j = 0; j < 3*Nquad; j++)
		  num_dof(j + offset_Hloc_s) = offset_H + 3*Nquad + j;
	      }

            if (this->ref_drude(ref_domain).IsEnabled())
              {
                if (var_problem.InsidePML(iquad))
                  {
                    cout << "Case not implemented" << endl;
                    abort();
                  }

                int nPole = this->ref_drude(ref_domain).gamma.GetM();
                int offset_Q = var_problem.GetNbDof() - this->nb_dof_drude_vec;
                int offset_P = offset_Q - this->nb_dof_drude_vec;
                offset_P += this->OffsetDofDrudeV(iquad);
                offset_Q += this->OffsetDofDrudeV(iquad);
                
                int nb_dof_vec = 3*Nquad;
                nb_dof_all += 2*nb_dof_vec*nPole;
                offset_Ploc = offset_Hloc + nb_dof_vec;
                offset_Qloc = offset_Ploc + nb_dof_vec*nPole;
                
                num_dof.Resize(nb_dof_all);
                for (int j = 0; j < nPole*nb_dof_vec; j++)
                  {
                    num_dof(offset_Ploc + j) = offset_P + j;
                    num_dof(offset_Qloc + j) = offset_Q + j;
                  }
                
                //DISP(iquad); DISP(num_dof);
                //DISP(offset_Ploc); DISP(offset_Qloc); DISP(offset_Hloc);
              }
	  }
      }
    else
      {
        if (this->linearize_drude && this->ref_drude(ref_domain).IsEnabled())
          if (this->ref_drude(ref_domain).IsModeTM() && (!var_computation.GetLeafStaticCondensation()))
            {
              if (var_problem.InsidePML(iquad))
                {
                  cout << "Case not implemented" << endl;
                  abort();
                }

              int nPole = this->ref_drude(ref_domain).gamma.GetM();
              int nb_dof_vec = 3*Nquad;
              nb_dof_all += nPole*nb_dof_vec;
              num_dof.Resize(nb_dof_all);
              int offset_P = var_problem.GetNbDof() - this->nb_dof_drude_vec;
              offset_Ploc = Fb.GetNbDof();
              offset_P += this->OffsetDofDrudeV(iquad);
              for (int j = 0; j < nPole*nb_dof_vec; j++)
                num_dof(offset_Ploc + j) = offset_P + j;
            }
      }

    
    // initialisation of elementary matrix
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();
    
    // variable coefficients, need of a numerical quadrature
    bool variable = var_problem.UseNumericalIntegration(iquad);
    
    int N = Fb.GetNbPointsQuadratureInside();
    bool affine = mesh.IsElementAffine(iquad);
    bool stiff = false;
    T zero; SetComplexZero(zero);
    if (nat_mat.GetCoefStiffness() != zero)
      stiff = true;

    Complexe m_iomega;
    var_problem.GetMiomega(m_iomega);
    Real_wp jacobian; T poids; SetComplexZero(poids);
    bool diag_mass = mass.IsDiagonalMass();
    
    if (var_problem.FirstOrderFormulation() && var_problem.InsidePML(iquad))
      {
	// no damping with pml for first order formulation (not implemented)
	bool presence_damping = false;
	if (!this->ref_sigma(ref_domain).IsZero())
	  presence_damping = true;
	
	if (presence_damping)
	  {
	    cout << "Damping with PML and first order formulation not implemented" << endl;
	    abort();
	  }
      }

    if ((var_problem.FirstOrderFormulation()) && (var_computation.GetLeafStaticCondensation()))
      {
        T coef_drude;
        SetComplexOne(coef_drude);
        
        if (this->ref_drude(ref_domain).IsEnabled())
          {
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            Real_wp eps_inf = this->ref_drude(ref_domain).eps_inf;
            T s = nat_mat.GetCoefStiffness();
            T m = nat_mat.GetCoefMass()*m_iomega;
            
            for (int kp = 0; kp < nPole; kp++)
              {
                Real_wp gamma = this->ref_drude(ref_domain).gamma(kp);
                Real_wp coef_omega2 = this->ref_drude(ref_domain).omega_02(kp);
                Real_wp coef_omegap = this->ref_drude(ref_domain).eps_omega_p2(kp) / eps_inf;
                Real_wp coef_sig = this->ref_drude(ref_domain).eps_sigma(kp) / eps_inf;
                
                coef_drude += (coef_omegap*s*s + coef_sig*s*m) / (coef_omega2*s*s + m*m + gamma*s*m);
              }
          }
        
	if (variable)
	  {
	    Matrix3_3 dfj, dfj_trans, dfjm1, dfjm1_trans;
	    Vector<TinyMatrix<T, General, 3, 3> > Bmass(N), Astiff(N);
            Vector<TinyMatrix<T, General, 3, 3> > Bmass_star, Bmass_cross;
	    TinyMatrix<T, General, 3, 3> A_tmp, B_tmp;
            TinyMatrix<T, Symmetric, 3, 3> eps_tau;

	    int i1 = iquad - mesh.GetNbElt() + var_boundary.GetNbEltPML();
	    bool pml_elt = (var_problem.InsidePML(iquad));
	    bool ortho = (var_problem.OrthogonalElement(iquad) == 0);
	    if ((pml_elt) && (!ortho) && (!var_computation.LightStaticCondensation()))
	      {
		cout << "Non orthogonal PML not implemented for static condensation" << endl;
		abort();
	      }

            if (!diag_mass && pml_elt && var_computation.LightStaticCondensation())
              {
                Bmass_star.Reallocate(N);
                Bmass_cross.Reallocate(N);
              }

	    if (pml_elt && !var_computation.LightStaticCondensation())
	      {
		// only Dirichlet or Neumann condition is accepted
		bool presence_damping = false;
		for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
		  {
		    int num_face = var_problem.mesh.Element(iquad).numBoundary(num_loc);
		    int ref_line = var_problem.mesh.Boundary(num_face).GetReference();
		    if (var_problem.mesh.GetBoundaryCondition(ref_line) == BoundaryConditionEnum::LINE_ABSORBING)
		      presence_damping = true;
		  }
		
		if (presence_damping)
		  {
		    cout << "Only Neumann or Dirichlet condition is accepted with static condensation and PML" << endl;
		    abort();
		  }
	      }
	    
	    for (int k = 0; k < N; k++)
	      {		
		Real_wp jacobian(1);
		if (affine)
		  {
		    jacobian = var_problem.Glob_jacobian(iquad)(0);
		    dfj = this->Glob_DFj(iquad)(0);
		  }
		else
		  {
		    jacobian = var_problem.Glob_jacobian(iquad)(k) / Fb.WeightsND(k);
		    dfj = this->Glob_DFj(iquad)(k);
		  }
		
		Transpose(dfj, dfj_trans);
                GetInverse(dfj, dfjm1);
                Transpose(dfjm1, dfjm1_trans);
                
		T poids_stiff = square(nat_mat.GetCoefStiffness()) / (m_iomega*nat_mat.GetCoefMass()) * Fb.WeightsND(k);
				
		if (!diag_mass)
		  Bmass(k) = this->Glob_matMass_elem(iquad).GetBh(k)*nat_mat.GetCoefMass()
		    + this->Glob_matMass_elem(iquad).GetBhSigma(k)*nat_mat.GetCoefDamping();
		
		if (pml_elt)
		  {
		    TinyVector<Complexe, 3> tau = var_boundary.GetTauPML(i1, k);
		    T m = nat_mat.GetCoefMass()*m_iomega;
		    T dx = m + tau(0)*nat_mat.GetCoefDamping(), dy = m + tau(1)*nat_mat.GetCoefDamping(), dz = m + tau(2)*nat_mat.GetCoefDamping();
		    T coef_x = dy*dz / (dx*m);
		    T coef_y = dx*dz / (dy*m);
		    T coef_z = dx*dy / (dz*m);
		    
		    if (!diag_mass)
                      {
                        if (!var_computation.LightStaticCondensation())
                          {
                            Bmass(k)(0, 0) *= coef_x;
                            Bmass(k)(1, 1) *= coef_y;
                            Bmass(k)(2, 2) *= coef_z;
                          }
                        else
                          {
                            // on rajoute le terme eps tau_{2, 3, 1} dans l'equation de E
                            eps_tau.SetDiagonal(this->ref_epsilon(ref_domain).GetConstant()(0, 0));
                            eps_tau(0, 0) *= tau(1)*nat_mat.GetCoefDamping();
                            eps_tau(1, 1) *= tau(2)*nat_mat.GetCoefDamping();
                            eps_tau(2, 2) *= tau(0)*nat_mat.GetCoefDamping();
                            Mlt(eps_tau, dfjm1_trans, A_tmp);
                            Mlt(dfjm1, A_tmp, B_tmp);
                            Add(T(jacobian*Fb.WeightsND(k)), B_tmp, Bmass(k));

                            // termes pour l'equation dans E*
                            eps_tau.SetIdentity();
                            eps_tau(0, 0) *= dz;
                            eps_tau(1, 1) *= dx;
                            eps_tau(2, 2) *= dy;
                            Mlt(eps_tau, dfjm1_trans, A_tmp);
                            Mlt(dfjm1, A_tmp, Bmass_star(k));
                            Mlt(jacobian*Fb.WeightsND(k), Bmass_star(k));

                            eps_tau.SetIdentity();
                            eps_tau(0, 0) *= dx;
                            eps_tau(1, 1) *= dy;
                            eps_tau(2, 2) *= dz;
                            Mlt(eps_tau, dfjm1_trans, A_tmp);
                            Mlt(dfjm1, A_tmp, Bmass_cross(k));
                            Mlt(-jacobian*Fb.WeightsND(k), Bmass_cross(k));
                          }
                      }
                    
                    eps_tau.SetDiagonal(this->ref_invMu(ref_domain).GetConstant()(0, 0));
                    eps_tau(0, 0) /= coef_x; eps_tau(1, 1) /= coef_y; eps_tau(2, 2) /= coef_z;
                    
                    Mlt(eps_tau, dfj, A_tmp);
                    Mlt(dfj_trans, A_tmp, Astiff(k));
                    Mlt(poids_stiff/jacobian, Astiff(k));
                  }
                else
                  {
                    Mlt(this->ref_invMu(ref_domain).GetConstant(), dfj, A_tmp);
                    Mlt(dfj_trans, A_tmp, Astiff(k));
                    Mlt(poids_stiff/jacobian, Astiff(k));

                  }
              }
            
            if (this->ref_drude(ref_domain).IsEnabled())
              {
                if (this->ref_drude(ref_domain).IsModeTM())
                  for (int k = 0; k < N; k++)
                    Bmass(k) *= coef_drude;
                else
                  {
                    coef_drude = Real_wp(1) / coef_drude;
                    for (int k = 0; k < N; k++)
                      Astiff(k) *= coef_drude;
                  }
	      }
            
	    if (Fb.OptimizedComputationElementaryMatrix())
              {
                if (!diag_mass)
		  Fb.AddVariableMassMatrix(0, 0, Bmass, mat_interac);

                Fb.AddVariableStiffnessMatrix(0, offset_Eloc_s, Astiff, mat_interac);
              }
            else
              ComputeVariableElementaryMatrix(iquad, Bmass, Astiff, !diag_mass,
                                              mat_interac, Fb, nat_mat);
            
            if (mesh_num.drop_interface_pml_dof && pml_elt && var_computation.LightStaticCondensation())
              {
                Vector<T> col(mat_interac.GetM());
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  if (mesh_num.GetDofPML(num_dof(i)) == -1)
                    {
                      mat_interac.GetDenseCol(offset_Eloc_s + i, col);
                      for (int j = 0; j < Fb.GetNbDof(); j++)
                        mat_interac.AddInteraction(j, i, col(j));
                      
                      col.Zero();
                      mat_interac.SetDenseCol(offset_Eloc_s+i, col);
                    }
              }
            
            if (diag_mass)
	      {
		T vloc;
		T m = nat_mat.GetCoefMass()*m_iomega, sig = nat_mat.GetCoefDamping();
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  {
		    if (pml_elt)
		      {
			vloc = mass.GetDh(i)*nat_mat.GetCoefMass();
			TinyVector<Real_wp, 3> tau = mass.GetTauPML(i);
			T dx = m + sig*tau(0), dy = m + sig*tau(1), dz = m + sig*tau(2);
			int n = Fb.GetCoordinateDof(i);
                        if (var_computation.LightStaticCondensation())
                          {
                            T vloc2(0), vloc3(0);
                            switch(n)
                              {
                              case 0 : vloc *= dy/ m; vloc2 = dz; vloc3 = dx; break;
                              case 1 : vloc *= dz / m; vloc2 = dx; vloc3 = dy; break;
                              case 2 : vloc *= dx / m; vloc2 = dy; vloc3 = dz; break;
                              }
                            
                            Complexe invEps = 1.0 / this->ref_epsilon(ref_domain).GetConstant()(0, 0);
                            vloc2 *= mass.GetDh(i)*invEps / m_iomega;
                            vloc3 *= mass.GetDh(i)*invEps / m_iomega;
                            mat_interac.AddInteraction(offset_Eloc_s+i, offset_Eloc_s+i, vloc2);
                            mat_interac.AddInteraction(offset_Eloc_s+i, i, -vloc3);
                          }
                        else
                          {
                            switch(n)
                              {
                              case 0 : vloc *= dy*dz / (dx*m); break;
                              case 1 : vloc *= dx*dz / (dy*m); break;
                              case 2 : vloc *= dx*dy / (dz*m); break;			    
                              }
                          }
		      }
		    else
		      vloc = nat_mat.GetCoefMass()*mass.GetDh(i) + nat_mat.GetCoefDamping()*mass.GetDhSigma(i);
                    
		    mat_interac.AddInteraction(i, i, vloc);
		  }
              }	      
            else
              {
                if (pml_elt && var_computation.LightStaticCondensation())
                  {
                    Fb.AddVariableMassMatrix(offset_Eloc_s, offset_Eloc_s, Bmass_star, mat_interac);
                    Fb.AddVariableMassMatrix(offset_Eloc_s, 0, Bmass_cross, mat_interac);
                  }
              }

            if (mesh_num.drop_interface_pml_dof && pml_elt && var_computation.LightStaticCondensation())
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (mesh_num.GetDofPML(num_dof(i)) == -1)
                  mat_interac.ClearRow(offset_Eloc_s+i);
            
          }
	else
	  {
	    Matrix3_3 dfj = this->Glob_DFj(iquad)(0), dfj_trans;
	    Transpose(dfj, dfj_trans);
	    
	    TinyMatrix<T, General, 3, 3> A_tmp, Astiff, Bmass;
	    
	    Real_wp jacobian = var_problem.Glob_jacobian(iquad)(0);
	    T poids_stiff = square(nat_mat.GetCoefStiffness()) / (m_iomega*nat_mat.GetCoefMass());
	    
	    Mlt(this->ref_invMu(ref_domain).GetConstant(), dfj, A_tmp);
            Mlt(dfj_trans, A_tmp, Astiff);
            Mlt(poids_stiff/jacobian, Astiff);

	    if (!diag_mass)
	      Bmass = this->Glob_matMass_elem(iquad).GetBh(0)*nat_mat.GetCoefMass()
		+ this->Glob_matMass_elem(iquad).GetBhSigma(0)*nat_mat.GetCoefDamping();
	    
	    TinyVector<bool, 4> null_term;
	    null_term.Fill(false);
	    
	    null_term(2) = true;
	    null_term(3) = true;

            if (this->ref_drude(ref_domain).IsEnabled())
              {
                if (this->ref_drude(ref_domain).IsModeTM())
                  Bmass *= coef_drude;
                else
                  {
                    coef_drude = Real_wp(1) / coef_drude;
                    for (int k = 0; k < N; k++)
                      Astiff *= coef_drude;
                  }                
	      }
            
	    if (diag_mass)
	      {
		for (int i = 0; i < Fb.GetNbDof(); i++)
		  mat_interac.AddInteraction(i, i, nat_mat.GetCoefMass()*coef_drude*mass.GetDh(i) +
					     nat_mat.GetCoefDamping()*mass.GetDhSigma(i));		
	      }
	    else
	      Fb.AddConstantMassMatrix(0, 0, Bmass, mat_interac);
	    
	    Fb.AddConstantStiffnessMatrix(0, 0, Astiff, mat_interac);
	  }
      }
    else if (var_problem.FirstOrderFormulation())
      {
        TinyVector<Real_wp, 3> vec_u, vec_v, vec_w;
        TinyVector<T, 3> vec_ur, vec_vr, vec_wr;
        TinyMatrix<T, Symmetric, 3, 3> Mu;
        
        bool sym = var_computation.GetSymmetrizationUse();
	
        // variables for PML
        T mu, epsilon;
        TinyVector<Complexe, 3> tau; R3 tau_r, coef_df;
        int i1 = iquad - mesh.GetNbElt() + var_boundary.GetNbEltPML();
        mu = this->ref_mu(ref_domain).GetConstant()(0, 0);
        epsilon = this->ref_epsilon(ref_domain).GetConstant()(0, 0);
        mu *= nat_mat.GetCoefDamping();
        epsilon *= nat_mat.GetCoefDamping();
        
        // variables for computing elementary matrix
        Vector<T> feval(3*Nquad), contrib(Fb.GetNbDof());
        VectReal_wp val_hat_phi(3*Nquad), curl_hat_phi(3*Nquad);
        VectR3 val_phi(Nquad), curl_phi(Nquad);
        VectReal_wp Ones(Fb.GetNbDof());
        TinyVector<T, 3> vec_uc, vec_vc;
        Vector<Matrix3_3> MatDFjm1;
        contrib.Fill(0);
	TinyVector<int, 3> permut120(1, 2, 0), permut201(2,0,1);
        if (var_problem.InsidePML(iquad))
	  {
	    if (!diag_mass)
	      MatDFjm1.Reallocate(Nquad);
	    else
	      {
		jacobian = var_problem.Glob_jacobian(iquad)(0);
		Matrix3_3 dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
		Mlt(1.0/jacobian, dfjm1);
		coef_df(0) = square(dfjm1(0, 0))*jacobian;
		coef_df(1) = square(dfjm1(1, 1))*jacobian;
		coef_df(2) = square(dfjm1(2, 2))*jacobian;
	      }
	  }

        if (this->ref_drude(ref_domain).IsModeTM())
          MatDFjm1.Reallocate(Nquad);
        
        Ones.Fill(0);
	Real_wp coef_sym(1);
	if (sym)
	  coef_sym = Real_wp(-1);

	Vector<int> dof_Estar(Fb.GetNbDof());
	for (int i = 0; i < Fb.GetNbDof(); i++)
	  dof_Estar(i) = offset_Eloc_s+i;

        if (mesh_num.drop_interface_pml_dof)
          {
            if (var_problem.InsidePML(iquad))
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (mesh_num.GetDofPML(num_dof(i)) == -1)
                  dof_Estar(i) = i;
          }
        else
          {
            int nb_new_dof = 0;
            if ((var_boundary.GetNbGlobalEltPML() > 0) && (!var_problem.InsidePML(iquad)))
              for (int i = 0; i < Fb.GetNbDof(); i++)
                if (mesh_num.GetDofPML(num_dof(i)) >= 0)
                  nb_new_dof++;
            
            if (nb_new_dof > 0)
              {
                nb_dof_all += Fb.GetNbDof();
                offset_Hloc = 2*Fb.GetNbDof();
                offset_Hloc_s = offset_Hloc;
                
                nb_unknowns = 2;
                num_dof.Resize(nb_dof_all);
                // dof numbers for H
                for (int j = 0; j < 3*Nquad; j++)
                  num_dof(j + offset_Hloc) = offset_H + j;
                
                mat_interac.Reallocate(nb_dof_all, nb_dof_all);
                mat_interac.Zero();
                
                for (int i = 0; i < Fb.GetNbDof(); i++)
                  {
                    if (mesh_num.GetDofPML(num_dof(i)) >= 0)
                      {
                        num_dof(Fb.GetNbDof() + i) = offset_Es + mesh_num.GetDofPML(num_dof(i));
                        dof_Estar(i) = Fb.GetNbDof() + i;
                      }
                    else
                      num_dof(Fb.GetNbDof() + i) = -1;
                  }
              }
          }
        
        // loop over rows of the elementary matrix
	const VectReal_wp& weights_dof = Fb.WeightsDofND();
	for (int i = 0; i < Fb.GetNbDof(); i++)
          {
            Ones(i) = 1.0;
            // computation of phi and curl phi on quadrature points
            Fb.ApplyChTranspose(Ones, val_hat_phi);
            if (stiff)
              {
                if (Fb.UseQuadraturePointsForRh())
                  Fb.ApplyRhQuadratureTranspose(val_hat_phi, curl_hat_phi);
                else              
                  Fb.ApplyRhTranspose(Ones, curl_hat_phi);
                
                // curl of phi on the real element
                if (affine)
                  {
                    jacobian = var_problem.Glob_jacobian(iquad)(0);
                    Matrix3_3 dfj = this->Glob_DFj(iquad)(0);
                    for (int k = 0; k < Nquad; k++)
                      {
                        vec_u.Init(curl_hat_phi(3*k), curl_hat_phi(3*k+1), curl_hat_phi(3*k+2)); 
                        Mlt(dfj, vec_u, curl_phi(k));
                        Mlt(1.0/jacobian, curl_phi(k));
                      }
                  }
                else
                  {
                    for (int k = 0; k < Nquad; k++)
                      {
                        jacobian = var_problem.Glob_jacobian(iquad)(k)/Fb.WeightsND(k);
                        vec_u.Init(curl_hat_phi(3*k), curl_hat_phi(3*k+1), curl_hat_phi(3*k+2));
                        Mlt(this->Glob_DFj(iquad)(k), vec_u, curl_phi(k));
                        Mlt(1.0/jacobian, curl_phi(k));
                      }
                  }
              }
            
            // value of phi on the real element
            if ((var_problem.InsidePML(iquad) && !diag_mass) || this->ref_drude(ref_domain).IsModeTM())
              {
                if (affine)
                  {
                    jacobian = var_problem.Glob_jacobian(iquad)(0);
                    Matrix3_3 dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
                    for (int k = 0; k < Nquad; k++)
                      {
                        vec_u.Init(val_hat_phi(3*k), val_hat_phi(3*k+1), val_hat_phi(3*k+2));                       
                        MatDFjm1(k) = dfjm1; Mlt(1.0/jacobian, MatDFjm1(k));
                        MltTrans(MatDFjm1(k), vec_u, val_phi(k));
                      }
                  }
                else
                  {
                    for (int k = 0; k < Nquad; k++)
                      {
                        jacobian = var_problem.Glob_jacobian(iquad)(k)/Fb.WeightsND(k);
                        vec_u.Init(val_hat_phi(3*k), val_hat_phi(3*k+1), val_hat_phi(3*k+2));
                        MatDFjm1(k) = var_problem.Glob_DFjm1(iquad)(k); Mlt(1.0/jacobian, MatDFjm1(k));
                        MltTrans(MatDFjm1(k), vec_u, val_phi(k));
                      }
                  }
              }
                        
	    if (diag_mass)
	      {
		// part epsilon dE/dt + sigma E
		mat_interac.SetEntry(i, i, nat_mat.GetCoefMass()*mass.GetDh(i)
                                     + nat_mat.GetCoefDamping()*mass.GetDhSigma(i));

		if (var_problem.InsidePML(iquad))
		  {
                    tau_r = mass.GetTauPML(i);

		    // part (tau1, tau2, tau0) epsilon E
		    int n = Fb.GetCoordinateDof(i);
		    mat_interac.AddInteraction(i, i, epsilon * coef_df(n) * weights_dof(i) * tau_r(permut120(n)));
                    if (num_dof(offset_Eloc_s+i) >= 0)
                      {
                        // part dE*/dt - dE/dt
                        T vloc = m_iomega * coef_df(n) * nat_mat.GetCoefMass() * weights_dof(i);
                        mat_interac.AddInteraction(offset_Eloc_s+i, offset_Eloc_s+i, vloc);
                        mat_interac.SetEntry(offset_Eloc_s+i, i, -vloc);
                        
                        // part (tau2, tau0, tau1) E*
                        vloc = nat_mat.GetCoefDamping() * coef_df(n) * weights_dof(i) * tau_r(permut201(n));
                        mat_interac.AddInteraction(offset_Eloc_s+i, offset_Eloc_s+i, vloc);
                        
                        // part - (tau0, tau1, tau2) E
                        vloc = nat_mat.GetCoefDamping() * coef_df(n) * weights_dof(i) * tau_r(n);
                        mat_interac.AddInteraction(offset_Eloc_s+i, i, -vloc);
                      }
		  }

	      }
	    else
	      {		
		// part epsilon dE/dt + sigma E
		for (int k = 0; k < Nquad; k++)
		  {
		    vec_ur.Init(val_hat_phi(3*k), val_hat_phi(3*k+1), val_hat_phi(3*k+2));
		    if (variable)
		      {
			Mlt(mass.GetBh(k), vec_ur, vec_vr);
			Mlt(mass.GetBhSigma(k), vec_ur, vec_wr);
		      }
		    else
		      {
			Mlt(mass.GetBh(0), vec_ur, vec_vr);
			Mlt(mass.GetBhSigma(0), vec_ur, vec_wr);
			vec_vr *= Fb.WeightsND(k); vec_wr *= Fb.WeightsND(k);
		      }
		    
		    feval(3*k) = vec_vr(0)*nat_mat.GetCoefMass() + vec_wr(0)*nat_mat.GetCoefDamping();
		    feval(3*k+1)
		      = vec_vr(1)*nat_mat.GetCoefMass() + vec_wr(1)*nat_mat.GetCoefDamping();
		    
		    feval(3*k+2)
		      = vec_vr(2)*nat_mat.GetCoefMass() + vec_wr(2)*nat_mat.GetCoefDamping();
		  }
		
		Fb.ApplyCh(feval, contrib);
		for (int j = 0; j < Fb.GetNbDof(); j++)
		  mat_interac.SetEntry(i, j, contrib(j));
	      }

            if (this->ref_drude(ref_domain).IsModeTM())
              {
                int nPole = this->ref_drude(ref_domain).gamma.GetM();
                T m = m_iomega*nat_mat.GetCoefMass();
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp coef_Q(1), coef_sig = this->ref_drude(ref_domain).eps_sigma(kp);
                    if (!this->use_symm_drude)
                      coef_Q = this->ref_drude(ref_domain).eps_omega_p2(kp);
                    else
                      coef_sig /= this->ref_drude(ref_domain).eps_omega_p2(kp); 
                    
                    // part +Q in E-equation for Drude's model
                    // and +E in Q-equation
                    for (int k = 0; k < Nquad; k++)
                      {
                        int offset_Q2 = offset_Qloc + 3*Nquad*kp + 3*k;
                        if (affine)
                          jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                        else
                          jacobian = var_problem.Glob_jacobian(iquad)(k);
                        
                        poids = jacobian*nat_mat.GetCoefStiffness();
                        mat_interac.SetEntry(i, offset_Q2, poids*val_phi(k)(0));
                        mat_interac.SetEntry(i, offset_Q2+1, poids*val_phi(k)(1));
                        mat_interac.SetEntry(i, offset_Q2+2, poids*val_phi(k)(2));
                        
                        poids = -poids*coef_sym*coef_Q;
                        poids -= jacobian*coef_sym*m*coef_sig;
                        mat_interac.SetEntry(offset_Q2, i, poids*val_phi(k)(0));
                        mat_interac.SetEntry(offset_Q2+1, i, poids*val_phi(k)(1));
                        mat_interac.SetEntry(offset_Q2+2, i, poids*val_phi(k)(2));
                      }
                  }
              }
            
            // part - curl H  /  curl E
            if (stiff)
              for (int k = 0; k < Nquad; k++)
                {
                  if (affine)
                    jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                  else
                    jacobian = var_problem.Glob_jacobian(iquad)(k);
                  
                  poids = jacobian*nat_mat.GetCoefStiffness();
                  mat_interac.SetEntry(i, offset_Hloc_s+3*k, -poids*curl_phi(k)(0));
                  mat_interac.SetEntry(i, offset_Hloc_s+3*k+1, -poids*curl_phi(k)(1));
                  mat_interac.SetEntry(i, offset_Hloc_s+3*k+2, -poids*curl_phi(k)(2));
                  
                  mat_interac.SetEntry(offset_Hloc+3*k, dof_Estar(i), coef_sym*poids*curl_phi(k)(0));
                  mat_interac.SetEntry(offset_Hloc+3*k+1, dof_Estar(i), coef_sym*poids*curl_phi(k)(1));
                  mat_interac.SetEntry(offset_Hloc+3*k+2, dof_Estar(i), coef_sym*poids*curl_phi(k)(2));
                }

            // additional terms due to PML
            if (var_problem.InsidePML(iquad) && (!diag_mass))
              {
                // part (tau1, tau2, tau0) epsilon E
                for (int k = 0; k < Nquad; k++)
                  {                    
                    if (affine)
                      jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                    else
                      jacobian = var_problem.Glob_jacobian(iquad)(k);
                    
                    tau = var_boundary.GetTauPML(i1, k);
                    poids = jacobian*epsilon;
                    vec_uc(0) = tau(1)*val_phi(k)(0)*poids;
                    vec_uc(1) = tau(2)*val_phi(k)(1)*poids;
                    vec_uc(2) = tau(0)*val_phi(k)(2)*poids;
                    Mlt(MatDFjm1(k), vec_uc, vec_vc);
                    
                    feval(3*k) = vec_vc(0);
                    feval(3*k+1) = vec_vc(1);
                    feval(3*k+2) = vec_vc(2);
                  }
                
                Fb.ApplyCh(feval, contrib);
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  mat_interac.AddInteraction(i, j, contrib(j));
                
                // part dE*/dt - dE/dt
                for (int k = 0; k < Nquad; k++)
                  {
                    if (affine)
                      jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                    else
                      jacobian = var_problem.Glob_jacobian(iquad)(k);
                    
                    poids = m_iomega*jacobian*nat_mat.GetCoefMass();
                    vec_uc(0) = val_phi(k)(0)*poids;
                    vec_uc(1) = val_phi(k)(1)*poids;
                    vec_uc(2) = val_phi(k)(2)*poids;
                    Mlt(MatDFjm1(k), vec_uc, vec_vc);
                    
                    feval(3*k) = vec_vc(0);
                    feval(3*k+1) = vec_vc(1);
                    feval(3*k+2) = vec_vc(2);
                  }
                
                Fb.ApplyCh(feval, contrib);
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  {
                    mat_interac.SetEntry(offset_Eloc_s+i, offset_Eloc_s+j, contrib(j));
                    mat_interac.SetEntry(offset_Eloc_s+i, j, -contrib(j));
                  }
                                
                // part (tau2, tau0, tau1) E*
                for (int k = 0; k < Nquad; k++)
                  {                    
                    if (affine)
                      jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                    else
                      jacobian = var_problem.Glob_jacobian(iquad)(k);
                    
                    tau = var_boundary.GetTauPML(i1, k);
                    poids = jacobian*nat_mat.GetCoefDamping();
                    vec_uc(0) = tau(2)*val_phi(k)(0)*poids;
                    vec_uc(1) = tau(0)*val_phi(k)(1)*poids;
                    vec_uc(2) = tau(1)*val_phi(k)(2)*poids;
                    Mlt(MatDFjm1(k), vec_uc, vec_vc);
                    
                    feval(3*k) = vec_vc(0);
                    feval(3*k+1) = vec_vc(1);
                    feval(3*k+2) = vec_vc(2);
                  }
                
                Fb.ApplyCh(feval, contrib);
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  mat_interac.AddInteraction(offset_Eloc_s+i, offset_Eloc_s+j, contrib(j));
                
                // part - (tau0, tau1, tau2) E
                for (int k = 0; k < Nquad; k++)
                  {                    
                    if (affine)
                      jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                    else
                      jacobian = var_problem.Glob_jacobian(iquad)(k);
                    
                    tau = var_boundary.GetTauPML(i1, k);
                    poids = jacobian*nat_mat.GetCoefDamping();
                    vec_uc(0) = tau(0)*val_phi(k)(0)*poids;
                    vec_uc(1) = tau(1)*val_phi(k)(1)*poids;
                    vec_uc(2) = tau(2)*val_phi(k)(2)*poids;
                    Mlt(MatDFjm1(k), vec_uc, vec_vc);
                    
                    feval(3*k) = vec_vc(0);
                    feval(3*k+1) = vec_vc(1);
                    feval(3*k+2) = vec_vc(2);
                  }
                
                Fb.ApplyCh(feval, contrib);
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  mat_interac.AddInteraction(offset_Eloc_s+i, j, -contrib(j));
              }
	    
            Ones(i) = 0.0;
          }

        // part mu dH/dt
        for (int k = 0; k < Nquad; k++)
          {
            if (affine)
              jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
            else
              jacobian = var_problem.Glob_jacobian(iquad)(k);
            
            Mu = this->ref_mu(ref_domain).GetCoefficient(var_problem, iquad, k);
            Mu *= coef_sym*m_iomega*jacobian*nat_mat.GetCoefMass();
            
            mat_interac.SetEntry(offset_Hloc+3*k, offset_Hloc+3*k, Mu(0, 0));
            mat_interac.SetEntry(offset_Hloc+3*k, offset_Hloc+3*k+1, Mu(0, 1));
            mat_interac.SetEntry(offset_Hloc+3*k, offset_Hloc+3*k+2, Mu(0, 2));
            mat_interac.SetEntry(offset_Hloc+3*k+1, offset_Hloc+3*k, Mu(1, 0));
            mat_interac.SetEntry(offset_Hloc+3*k+1, offset_Hloc+3*k+1, Mu(1, 1));
            mat_interac.SetEntry(offset_Hloc+3*k+1, offset_Hloc+3*k+2, Mu(1, 2));
            mat_interac.SetEntry(offset_Hloc+3*k+2, offset_Hloc+3*k, Mu(2, 0));
            mat_interac.SetEntry(offset_Hloc+3*k+2, offset_Hloc+3*k+1, Mu(2, 1));
            mat_interac.SetEntry(offset_Hloc+3*k+2, offset_Hloc+3*k+2, Mu(2, 2));
          }
        
        if (var_problem.InsidePML(iquad))
          for (int k = 0; k < Nquad; k++)
            {
              if (affine)
                jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
              else
                jacobian = var_problem.Glob_jacobian(iquad)(k);
              
              tau = var_boundary.GetTauPML(i1, k);
	      
              // part (tau1, tau2, tau0) mu H
              poids = mu*jacobian;
              mat_interac.AddInteraction(offset_Hloc+3*k, offset_Hloc+3*k, tau(1)*poids);
              mat_interac.AddInteraction(offset_Hloc+3*k+1, offset_Hloc+3*k+1, tau(2)*poids);
              mat_interac.AddInteraction(offset_Hloc+3*k+2, offset_Hloc+3*k+2, tau(0)*poids);

              // part dH*/dt - dH/dt
              poids = m_iomega*jacobian*nat_mat.GetCoefMass();
              mat_interac.SetEntry(offset_Hloc_s+3*k, offset_Hloc_s+3*k, poids);
              mat_interac.SetEntry(offset_Hloc_s+3*k+1, offset_Hloc_s+3*k+1, poids);
              mat_interac.SetEntry(offset_Hloc_s+3*k+2, offset_Hloc_s+3*k+2, poids);
              
              mat_interac.SetEntry(offset_Hloc_s+3*k, offset_Hloc+3*k, -poids);
              mat_interac.SetEntry(offset_Hloc_s+3*k+1, offset_Hloc+3*k+1, -poids);
              mat_interac.SetEntry(offset_Hloc_s+3*k+2, offset_Hloc+3*k+2, -poids);

              // part (tau2, tau0, tau1) H*
              poids = nat_mat.GetCoefDamping()*jacobian;
              mat_interac.AddInteraction(offset_Hloc_s+3*k, offset_Hloc_s+3*k, tau(2)*poids);
              mat_interac.AddInteraction(offset_Hloc_s+3*k+1, offset_Hloc_s+3*k+1, tau(0)*poids);
              mat_interac.AddInteraction(offset_Hloc_s+3*k+2, offset_Hloc_s+3*k+2, tau(1)*poids);
              
              // part -(tau0, tau1, tau2) H
              mat_interac.AddInteraction(offset_Hloc_s+3*k, offset_Hloc+3*k, -tau(0)*poids);
              mat_interac.AddInteraction(offset_Hloc_s+3*k+1, offset_Hloc+3*k+1, -tau(1)*poids);
              mat_interac.AddInteraction(offset_Hloc_s+3*k+2, offset_Hloc+3*k+2, -tau(2)*poids);  
            }

        if (this->ref_drude(ref_domain).IsEnabled())
          {
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            T m = m_iomega*nat_mat.GetCoefMass();
            for (int kp = 0; kp < nPole; kp++)
              {
                Real_wp gamma = this->ref_drude(ref_domain).gamma(kp);
                Real_wp coef_omega2 = this->ref_drude(ref_domain).omega_02(kp);
                Real_wp coef_eps_inf = this->ref_drude(ref_domain).eps_omega_p2(kp);
                
                Real_wp coef1(1), coef2(1);
                if (this->use_symm_drude)
                  {
                    coef1 = Real_wp(1)/coef_eps_inf;
                    coef2 = coef_omega2 * coef1;
                  }
                
                Real_wp c1 = -1, csig = -this->ref_drude(ref_domain).eps_sigma(kp);
                if (this->ref_drude(ref_domain).IsModeTE())
                  {
                    if (!this->use_symm_drude)
                      c1 = -coef_eps_inf;
                    else
                      csig /= coef_eps_inf;
                    
                    coef1 *= coef_sym;
                    coef2 *= coef_sym;
                  }
                
                for (int k = 0; k < Nquad; k++)
                  {
                    if (affine)
                      jacobian = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(k);
                    else
                      jacobian = var_problem.Glob_jacobian(iquad)(k);
                    
                    if (this->ref_drude(ref_domain).IsModeTE())
                      {
                        // part +Q in H-equation
                        poids = jacobian*nat_mat.GetCoefStiffness()*coef_sym;
                        mat_interac.SetEntry(offset_Hloc+3*k, offset_Qloc+3*k, poids);
                        mat_interac.SetEntry(offset_Hloc+3*k+1, offset_Qloc+3*k+1, poids);
                        mat_interac.SetEntry(offset_Hloc+3*k+2, offset_Qloc+3*k+2, poids);
                        
                        // and - eps_inf omega_p^2 H - eps_inf sigma (-i omega H) in Q-equation
                        poids = jacobian*(nat_mat.GetCoefStiffness()*c1 + m*csig);
                        mat_interac.SetEntry(offset_Qloc+3*k, offset_Hloc+3*k, poids);
                        mat_interac.SetEntry(offset_Qloc+3*k+1, offset_Hloc+3*k+1, poids);
                        mat_interac.SetEntry(offset_Qloc+3*k+2, offset_Hloc+3*k+2, poids);
                      }
                    
                    // part -i omega (omega_0^2) / (eps_inf omega_p^2) P in P-equation
                    poids = jacobian*m_iomega*nat_mat.GetCoefMass()*coef2;
                    mat_interac.SetEntry(offset_Ploc+3*k, offset_Ploc+3*k, poids);
                    mat_interac.SetEntry(offset_Ploc+3*k+1, offset_Ploc+3*k+1, poids);
                    mat_interac.SetEntry(offset_Ploc+3*k+2, offset_Ploc+3*k+2, poids);
                    
                    // part -(omega_0^2) / (eps_inf omega_p^2) Q in P-equation
                    poids = -jacobian*nat_mat.GetCoefStiffness()*coef2;
                    mat_interac.SetEntry(offset_Ploc+3*k, offset_Qloc+3*k, poids);
                    mat_interac.SetEntry(offset_Ploc+3*k+1, offset_Qloc+3*k+1, poids);
                    mat_interac.SetEntry(offset_Ploc+3*k+2, offset_Qloc+3*k+2, poids);
                    
                    // symmetric part -(omega_0^2)/ (eps_inf omega_p^2) P in Q-equation
                    poids = jacobian*nat_mat.GetCoefStiffness()*coef_omega2*coef1*coef_sym;
                    mat_interac.SetEntry(offset_Qloc+3*k, offset_Ploc+3*k, poids);
                    mat_interac.SetEntry(offset_Qloc+3*k+1, offset_Ploc+3*k+1, poids);
                    mat_interac.SetEntry(offset_Qloc+3*k+2, offset_Ploc+3*k+2, poids);
                    
                    // part i omega / (eps_inf omega_p^2) Q - gamma / (eps_inf omega_p^2) Q in Q-equation
                    poids = coef_sym*jacobian*coef1*(m_iomega*nat_mat.GetCoefMass() + gamma*nat_mat.GetCoefDamping());
                    mat_interac.SetEntry(offset_Qloc+3*k, offset_Qloc+3*k, poids);
                    mat_interac.SetEntry(offset_Qloc+3*k+1, offset_Qloc+3*k+1, poids);
                    mat_interac.SetEntry(offset_Qloc+3*k+2, offset_Qloc+3*k+2, poids);
                  }

                offset_Ploc += 3*Nquad; offset_Qloc += 3*Nquad;
              }
          }
      }
    else
      {
        T coef_drude; SetComplexOne(coef_drude);
        if (var_computation.GetLeafStaticCondensation() && this->linearize_drude
            && this->ref_drude(ref_domain).IsEnabled() && this->ref_drude(ref_domain).IsModeTM())
          {
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            Real_wp eps_inf = this->ref_drude(ref_domain).eps_inf;
            T s = nat_mat.GetCoefStiffness();
            T m = nat_mat.GetCoefMass()*m_iomega*m_iomega;
            T d = nat_mat.GetCoefDamping()*m_iomega;
            
            for (int kp = 0; kp < nPole; kp++)
              {
                Real_wp gamma = this->ref_drude(ref_domain).gamma(kp);
                Real_wp coef_omega2 = this->ref_drude(ref_domain).omega_02(kp);
                Real_wp coef_omegap = this->ref_drude(ref_domain).eps_omega_p2(kp) / eps_inf;
                Real_wp coef_sig = this->ref_drude(ref_domain).eps_sigma(kp) / eps_inf;
                
                coef_drude += (coef_omegap*s + coef_sig*d) / (coef_omega2*s + m + gamma*d);
              }
          }
        
        // second-order formulation of Maxwell's equations
        if (variable)
          {
            Vector<TinyMatrix<T, General, 3, 3> > Bmass(N), Astiff(N);
            for (int j = 0; j < N; j++)
              {
                if (!diag_mass)
		  {
                    this->GetMassMatrix(iquad, j, nat_mat, Bmass(j));
                    Bmass(j) *= coef_drude;
                  }

		this->GetStiffMatrix(iquad, j, nat_mat, Astiff(j));
              }
            
            if (Fb.OptimizedComputationElementaryMatrix())
              {
                if (!diag_mass)
		  Fb.AddVariableMassMatrix(0, 0, Bmass, mat_interac);

		Fb.AddVariableStiffnessMatrix(0, 0, Astiff, mat_interac);
              }
            else
              ComputeVariableElementaryMatrix(iquad, Bmass, Astiff, !diag_mass,
                                              mat_interac, Fb, nat_mat);
          }
        else
          {
            TinyMatrix<T, General, 3, 3> Bmass_cte, Astiff_cte;
            if (!diag_mass)
	      {
		this->GetMassMatrix(iquad, 0, nat_mat, Bmass_cte);
                Bmass_cte *= coef_drude;
		Fb.AddConstantMassMatrix(0, 0, Bmass_cte, mat_interac);
	      }
	    
	    this->GetStiffMatrix(iquad, 0, nat_mat, Astiff_cte);
	    Fb.AddConstantStiffnessMatrix(0, 0, Astiff_cte, mat_interac);
          }

        if (diag_mass)
          for (int i = 0; i < Fb.GetNbDof(); i++)
            {
              mat_interac.AddInteraction(i, i, nat_mat.GetCoefMass()*mass.GetDh(i) +
                                         nat_mat.GetCoefDamping()*mass.GetDhSigma(i));
            }

        if (this->linearize_drude && this->ref_drude(ref_domain).IsEnabled()
            && this->ref_drude(ref_domain).IsModeTM() && !var_computation.GetLeafStaticCondensation())
          {
            Complexe m_omega2 = m_iomega*m_iomega;
            int nPole = this->ref_drude(ref_domain).gamma.GetM();
            R3 vec_u;  T poidsE = m_omega2*nat_mat.GetCoefMass();
            Vector<T> poidsP(nPole), poidsDiag(nPole);
            for (int k = 0; k < nPole; k++)
              {
                poidsP(k) = -this->ref_drude(ref_domain).eps_omega_p2(k)*nat_mat.GetCoefStiffness()
                  - this->ref_drude(ref_domain).eps_sigma(k)*nat_mat.GetCoefDamping()*m_iomega;
                poidsDiag(k) = m_omega2*nat_mat.GetCoefMass()
                  + m_iomega*nat_mat.GetCoefDamping()*this->ref_drude(ref_domain).gamma(k)
                  + this->ref_drude(ref_domain).omega_02(k)*nat_mat.GetCoefStiffness();
              }
            
            VectR3 val_phi;
            for (int i = 0; i < Nquad; i++)
              {
                Fb.GetValuePhiOnQuadraturePoint(i, val_phi);
                Matrix3_3 dfjm1; Real_wp jacob_weighted;
                if (affine)
                  {
                    dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
                    jacob_weighted = var_problem.Glob_jacobian(iquad)(0)*Fb.WeightsND(i);
                  }
                else
                  {
                    dfjm1 = var_problem.Glob_DFjm1(iquad)(i);
                    jacob_weighted = var_problem.Glob_jacobian(iquad)(i);
                  }
                
                Mlt(Fb.WeightsND(i), dfjm1);
                // part -omega^2 \sum P_k in E-equation
                // and -c_k E - i omega sigma_k E in P-equation
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  {
                    MltTrans(dfjm1, val_phi(j), vec_u);
                    for (int k = 0; k < nPole; k++)
                      {
                        int offset_P2 = offset_Ploc + 3*Nquad*k + 3*i;
                        mat_interac.SetEntry(j, offset_P2, poidsE*vec_u(0));
                        mat_interac.SetEntry(j, offset_P2+1, poidsE*vec_u(1));
                        mat_interac.SetEntry(j, offset_P2+2, poidsE*vec_u(2));

                        mat_interac.SetEntry(offset_P2, j, poidsP(k)*vec_u(0));
                        mat_interac.SetEntry(offset_P2+1, j, poidsP(k)*vec_u(1));
                        mat_interac.SetEntry(offset_P2+2, j, poidsP(k)*vec_u(2));
                      }
                  }
                
                // diagonal part for P
                for (int k = 0; k < nPole; k++)
                  {
                    int offset_P2 = offset_Ploc + 3*Nquad*k + 3*i;
                    mat_interac.SetEntry(offset_P2, offset_P2, poidsDiag(k)*jacob_weighted);
                    mat_interac.SetEntry(offset_P2+1, offset_P2+1, poidsDiag(k)*jacob_weighted);
                    mat_interac.SetEntry(offset_P2+2, offset_P2+2, poidsDiag(k)*jacob_weighted);
                  }
              }
          }
      }
    
    // modification of the matrix if global dofs are obtained by linear combination of local dofs
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad, nb_unknowns);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad, nb_unknowns);

  }


  //! computes the elementary matrix for HDG formulation
  template<class Complexe> template<class T>
  void HarmonicMaxwell_3D<Complexe>
  ::ComputeElementaryMatrixHdg(int iquad, IVect& num_dof, VirtualMatrix<T>& mat_interac,
			       CondensationBlockSolver_Base<T>& solver_c,
			       const GlobalGenericMatrix<T>& nat_mat,
			       const ElementReference<Dimension3, 2>& Fb)
  {
    int nb_dof_elt = Fb.GetNbDof();

    const Mesh<Dimension3>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    
    int i1 = iquad - mesh.GetNbElt() + var_boundary.GetNbEltPML();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    int ref_domain = var_problem.mesh.Element(iquad).GetReference();
    bool affine = mesh.IsElementAffine(iquad);
    bool treat_inside_stiff = solver_c.TreatInsideStiffness();
    
    int offset_E, offset_H, offset_Eloc, offset_Hloc;
    
    bool optim_condensation = false;
    int ic = solver_c.GetCondensedElementNumber();
    if (var_computation.GetLeafStaticCondensation())
      {
	CondensationBlockSolver_MaxwellHdg3D<T, Complexe>& solver
	  = dynamic_cast<CondensationBlockSolver_MaxwellHdg3D<T, Complexe>& >(solver_c);
	
	if ((Fb.LumpedMassMatrix()) && (!Fb.DiscontinuousElement()))
	  optim_condensation = true;
	
	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  {
	    int num_face = mesh.Element(iquad).numBoundary(num_loc);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    if (rf != Fb.GetOrder())
	      optim_condensation = false;
	  }

	solver.SetOptimizedCondensation(ic, optim_condensation);
      }

    // first Lambda, then E and H
    int nb_dof_lambda = mesh_num.Element(iquad).GetNbDof();   

    offset_E = mesh_num.GetNbDof();
    offset_H = var_problem.GetOffsetDofUnknown(2);

    offset_Eloc = nb_dof_lambda;
    offset_Hloc = nb_dof_elt + nb_dof_lambda;
    
    int nb_dof_all = nb_dof_lambda + 2*nb_dof_elt;
    if (optim_condensation)
      nb_dof_all = nb_dof_lambda;
    
    num_dof.Reallocate(nb_dof_all);    
    for (int j = 0; j < mesh_num.Element(iquad).GetNbDof(); j++)
      num_dof(j) = mesh_num.Element(iquad).GetNumberDof(j);

    if (!optim_condensation)
      {
	for (int j = 0; j < nb_dof_elt; j++)
	  num_dof(offset_Eloc + j) = offset_E + var_problem.GetOffsetDofV(iquad) + j;
	
	for (int j = 0; j < nb_dof_elt; j++)
	  num_dof(offset_Hloc + j) = offset_H + var_problem.GetOffsetDofV(iquad) + j;
      }

    //DISP(iquad); DISP(num_dof);
    mat_interac.Reallocate(nb_dof_all, nb_dof_all);
    mat_interac.Zero();

    T zero; SetComplexZero(zero);
    T one; SetComplexOne(one);
    bool stiff = false;
    T s = nat_mat.GetCoefStiffness();
    if (nat_mat.GetCoefStiffness() != zero)
      stiff = true;

    bool mass = false;
    T m = nat_mat.GetCoefMass(), sig = nat_mat.GetCoefDamping();
    if ((nat_mat.GetCoefMass() != zero) || (nat_mat.GetCoefDamping() != zero))
      mass = true;

    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    T coef_tau = -var_problem.alpha_penalization;
    if (var_problem.upwind_fluxes)
      coef_tau = -var_problem.alpha_penalization*this->coefficient_impedance_absorbing(ref_domain);

    coef_tau *= s;

    Complexe coef_H = m_iomega;
    if (var_computation.GetSymmetrizationUse())
      coef_H = -coef_H;
    
    // volume integrals
    if (optim_condensation)
      {
	CondensationBlockSolver_MaxwellHdg3D<T, Complexe>& solver
	  = dynamic_cast<CondensationBlockSolver_MaxwellHdg3D<T, Complexe>& >(solver_c);

	const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	
	// we compute Dh and inverse of Bh
	TinyMatrix<Complexe, Symmetric, 3, 3> epsilon, mu, sigma, Bmass, Cdamp, Astiff;
	TinyMatrix<Complexe, General, 3, 3> A_tmp;
	TinyMatrix<T, Symmetric, 3, 3> As;
	Matrix3_3 dfjm1, dfj_trans; Real_wp jacob, jacob_weighted;
	Vector<TinyMatrix<T, Symmetric, 3, 3> >& invBhVol = solver.GetInverseBh(ic);
	invBhVol.Reallocate(nb_points_quad);
	Matrix<T, General, ArrayRowSparse> mat_invBh(nb_dof_elt, nb_dof_elt);
	Matrix<T, Symmetric, RowSymPacked>& Btilde = solver.GetBtildeMatrix(ic);
	Btilde.Reallocate(nb_dof_elt, nb_dof_elt);
	Btilde.Zero();
	if (treat_inside_stiff)
	  for (int i = 0; i < nb_points_quad; i++)
	    {
	      epsilon = this->ref_epsilon(ref_domain).GetCoefficient(var_problem, iquad, i);
	      sigma = this->ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, i);
	      mu = this->ref_mu(ref_domain).GetCoefficient(var_problem, iquad, i);
	      // indices are modified inside PML
	      if ((var_problem.InsidePML(iquad)) && (!var_problem.FirstOrderFormulation()))
		this->ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, i);
	      
	      var_problem.GetInverseJacobianMatrix(iquad, i, affine, Fb.GetGeometricElement(),
						   dfjm1, jacob, jacob_weighted);

	      MltTrans(epsilon, dfjm1, A_tmp);
	      Mlt(dfjm1, A_tmp, Bmass); 
	      Mlt(m_iomega*jacob_weighted, Bmass);

	      MltTrans(sigma, dfjm1, A_tmp);
	      Mlt(dfjm1, A_tmp, Cdamp); 
	      Mlt(jacob_weighted, Cdamp);

	      int ix, iy, iz;
	      Fb_hex.GetDofNumber_FromPointNode(i, ix, iy, iz);

	      As = Bmass*nat_mat.GetCoefMass() + Cdamp*nat_mat.GetCoefDamping();
	      Btilde(ix, ix) = As(0, 0); Btilde(iy, iy) = As(1, 1); Btilde(iz, iz) = As(2, 2);
	      Btilde(ix, iy) = As(0, 1); Btilde(ix, iz) = As(0, 2); Btilde(iy, iz) = As(1, 2);
	      
	      MltTrans(mu, dfjm1, A_tmp);
	      Mlt(dfjm1, A_tmp, Astiff);

	      GetInverse(Astiff);
	      As = Astiff;
	      
	      As *= square(nat_mat.GetCoefStiffness()) / (jacob_weighted*coef_H*nat_mat.GetCoefMass());
	      invBhVol(i) = As;
	      if (var_computation.GetSymmetrizationUse())
		As = -As; 
	      
	      mat_invBh.AddInteraction(ix, ix, As(0, 0));
	      mat_invBh.AddInteraction(ix, iy, As(0, 1));
	      mat_invBh.AddInteraction(ix, iz, As(0, 2));
	      mat_invBh.AddInteraction(iy, ix, As(1, 0));
	      mat_invBh.AddInteraction(iy, iy, As(1, 1));
	      mat_invBh.AddInteraction(iy, iz, As(1, 2));
	      mat_invBh.AddInteraction(iz, ix, As(2, 0));
	      mat_invBh.AddInteraction(iz, iy, As(2, 1));
	      mat_invBh.AddInteraction(iz, iz, As(2, 2));
	    }
	
	int offset = 0;
	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  offset += Fb.GetNbQuadBoundary(num_loc);
	
	Matrix<T> schur(2*offset, 2*offset);
	schur.Zero();
	
	// loop over boundaries of the element
	Complexe phase, phase_conj;
	TinyVector<Real_wp, 3> vec_ur;
	Matrix<T, General, ArrayRowSparse>& SnD = solver.GetSnD(ic);
	SnD.Reallocate(2*offset, nb_dof_elt);
	Matrix<T, General, ArrayRowSparse>& Mint = solver.GetMint(ic);
	Mint.Reallocate(2*offset, nb_dof_elt);
	offset = 0;
	for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	  {
	    int num_face, ref_boundary, rf, rot; bool new_face;
	    const Matrix<int>& FacesQuadRotation =
	      var_boundary.GetGeometryPhaseData(iquad, num_loc,
					       num_face, ref_boundary, rf, new_face, rot,
					       phase, phase_conj);
	    
 	    if (rf != Fb.GetOrder())
	      {
		cout << "Variable order not implemented with optimized static condensation" << endl;
		abort();
	      }

	    // for hexas only
	    int nx = 0, ny = 1;
	    if ((num_loc == 0) || (num_loc == 5))
	      { nx = 1; ny = 2; }
	    else if ((num_loc == 1) || (num_loc == 4))
	      { nx = 0; ny = 2; }	

	    const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(var_problem.GetSurfaceFiniteElement(num_face));
	    int order = Fb_quad.GetOrder();
            
	    TinyMatrix<Real_wp, Symmetric, 3, 3> DF_DFt;
	    Real_wp dsj(1); TinyMatrix<Real_wp, Symmetric, 2, 2> mat_DF_DFt;
	    for (int k = 0; k < Fb.GetNbQuadBoundary(num_loc); k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  dsj = var_problem.Glob_dsj(num_face)(k);
		else
		  dsj = var_problem.Glob_dsj(num_face)(krot);
		
		int kvol = Fb.GetQuadNumber(num_loc, k);
		var_problem.GetInverseJacobianMatrix(iquad, kvol, affine,
						     Fb.GetGeometricElement(),
						     dfjm1, jacob, jacob_weighted);

		GetInverse(dfjm1, dfj_trans); Transpose(dfj_trans);
		MltTrans(dfj_trans, dfj_trans,  DF_DFt);
		Mlt(Fb_quad.WeightsND(k)/dsj, DF_DFt);
		
		mat_DF_DFt(0, 0) = DF_DFt(ny, ny);
		mat_DF_DFt(0, 1) = -DF_DFt(nx, ny);
		mat_DF_DFt(1, 1) = DF_DFt(nx, nx);

		int kx = Fb_quad.GetXdofNumber(k);
		int ky = Fb_quad.GetYdofNumber(k);
		Real_wp sx = 1, sy = 1;
		if ((kx >= 2*(order+1)) && (kx < 4*(order+1)))
		  sx = -sx;
		
		if ((ky >= 2*(order+1)) && (ky < 4*(order+1)))
		  sy = -sy;

		int kvx, kvy, kvz;
		Fb_hex.GetDofNumber_FromPointNode(kvol, kvx, kvy, kvz);
		int jx = kvx, jy = kvy;
		if ((num_loc == 0) || (num_loc == 5))
		  { jx = kvy; jy = kvz; }
		else if ((num_loc == 1) || (num_loc == 4))
		  { jx = kvx; jy = kvz; }

		T vloc = -Fb.WeightsQuadratureBoundary(k, num_loc);
		if ((num_loc == 0) || (num_loc == 2) || (num_loc == 4))
		  vloc = -vloc;
		
		if (treat_inside_stiff)
		  {
		    if ((ref_boundary == 0) ||
			(mesh.GetBoundaryCondition(ref_boundary) != BoundaryConditionEnum::LINE_DIRICHLET))
		      {			
			SnD.AddInteraction(2*offset + kx, jy, vloc*sx);
			SnD.AddInteraction(2*offset + ky, jx, -vloc*sy);
		      }
		  }

		T coef_abc(coef_tau);
		if (var_computation.GetSymmetrizationUse())
		  coef_abc = -coef_tau;
		
		this->ModifyCoefficientBC_HDG(mesh, ref_boundary,
					      k, nat_mat.GetCoefStiffness(), coef_H,
					      ref_domain, iquad, num_loc, coef_abc);

		if (var_computation.GetSymmetrizationUse())
		  coef_abc = -coef_abc;
		
		schur(2*offset+kx, 2*offset+kx) -= coef_abc*mat_DF_DFt(0, 0);
		schur(2*offset+ky, 2*offset+ky) -= coef_abc*mat_DF_DFt(1, 1);
		schur(2*offset+kx, 2*offset+ky) -= coef_abc*sx*sy*mat_DF_DFt(0, 1);
		schur(2*offset+ky, 2*offset+kx) -= coef_abc*sx*sy*mat_DF_DFt(0, 1);
		
		if (treat_inside_stiff)
		  {
		    vloc = coef_tau;
		    Mint.AddInteraction(2*offset+kx, jx, -vloc*sx*mat_DF_DFt(0, 0));
		    Mint.AddInteraction(2*offset+kx, jy, -vloc*sx*mat_DF_DFt(0, 1));
		    Mint.AddInteraction(2*offset+ky, jx, -vloc*sy*mat_DF_DFt(1, 0));
		    Mint.AddInteraction(2*offset+ky, jy, -vloc*sy*mat_DF_DFt(1, 1));
		    
		    //Ca(offset + k) = coef_abc*s*mat_DF_DFt;
		    //Cl(offset + k) = coef_tau*s*mat_DF_DFt;
		    Btilde.AddInteraction(jx, jx, coef_tau*mat_DF_DFt(0, 0));
		    Btilde.AddInteraction(jx, jy, coef_tau*mat_DF_DFt(0, 1));
		    Btilde.AddInteraction(jy, jx, coef_tau*mat_DF_DFt(1, 0));
		    Btilde.AddInteraction(jy, jy, coef_tau*mat_DF_DFt(1, 1));
		  }
	      }
	    
	    offset += Fb.GetNbQuadBoundary(num_loc);
	  }

	if (!treat_inside_stiff)
	  {
	    mat_interac = schur;
	    return;
	  }
	
	// then we compute matrix Btilde = D_h + (R* - Sn) B_h^{-1} (R - Sn*) + C
	const Matrix<Real_wp, General, ArrayRowSparse>& R = Fb_hex.GetPermutedRh();
	Matrix<T, General, ArrayRowSparse> Rc;
	Copy(R, Rc);
	Matrix<T, General, ArrayRowSparse> Bmat(nb_dof_elt, nb_dof_elt), Btilde_sp(nb_dof_elt, nb_dof_elt);
	Mlt(mat_invBh, Rc, Bmat);
	Transpose(Rc); Mlt(Rc, Bmat, Btilde_sp);

	for (int i = 0; i < Btilde_sp.GetM(); i++)
	  for (int j = 0; j < Btilde_sp.GetRowSize(i); j++)
	    if (Btilde_sp.Index(i, j) >= i)
	      Btilde(i, Btilde_sp.Index(i, j)) += Btilde_sp.Value(i, j);

	// Btilde is replaced by its inverse
	GetInverse(Btilde);
	
	// we compute Mint = -Cl^T + Sn^D B_h^{-1} (R - Sn^T)
	MltAdd(-one, SnD, Bmat, one, Mint);
	
	// then we form schur = -Ca + Mint Btilde^{-1} Mint^T
	Matrix<T> Adense; Matrix<T, General, ArrayRowSparse> Mint_T;
	Adense.Reallocate(2*offset, nb_dof_elt);
	Mint_T.Reallocate(nb_dof_elt, 2*offset);
	Mlt(Mint, Btilde, Adense);

	Mint_T = Mint; Transpose(Mint_T);
	MltAdd(one, Adense, Mint_T, one, schur);	
	
	Adense.Clear(); Mint_T.Clear();
	Btilde_sp.Reallocate(2*offset, 2*offset);	
	Bmat.Reallocate(2*offset, nb_dof_elt);

	// we add -Sn^D B_h^{-1} (Sn^D)^T to schur
	Mlt(SnD, mat_invBh, Bmat);
	
	MltAdd(-one, SeldonNoTrans, Bmat, SeldonTrans, SnD, one, Btilde_sp);
	
	for (int i = 0; i < Btilde_sp.GetM(); i++)
	  for (int j = 0; j < Btilde_sp.GetRowSize(i); j++)
	    schur(i, Btilde_sp.Index(i, j)) += Btilde_sp.Value(i, j);

	if (var_computation.GetSymmetrizationUse())
	  Mlt(Real_wp(-1), schur);
	
	// the elementary matrix is the Schur complement for static condensation
	mat_interac = schur;

	//mat_interac.Write("mat_elem_cond.dat");
	//int test_input; cout << "waiting" << endl; cin >> test_input;

	mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad, 1);
	mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad, 1);
	
	return;
      }
    else if (Fb.LumpedMassMatrix())
      {
	const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	if (stiff)
	  {
	    // stiffness part E \cdot \rot psi and -rot H \cdot phi
	    // we group the term \int_K rot E \cdot \psi - \int n \times E \cdot \psi
	    // in the single term E \cdot rot psi (idem for - rot H \cdot phi
	    VectR3 curl_phi;
	    for (int i = 0; i < nb_points_quad; i++)
	      {
		T poids = Fb.WeightsND(i)*nat_mat.GetCoefStiffness();
		T poids_sym = poids;
		if (var_computation.GetSymmetrizationUse())
		  poids_sym = -poids;

		int ix, iy, iz;
		Fb_hex.GetDofNumber_FromPointNode(i, ix, iy, iz);
		Fb.GetCurlPhiOnQuadraturePoint(i, curl_phi);
		for (int j = 0; j < nb_dof_elt; j++)
		  if (!curl_phi(j).IsZero())
		    {
		      mat_interac.SetEntry(offset_Eloc + ix, offset_Hloc + j, -poids*curl_phi(j)(0));
		      mat_interac.SetEntry(offset_Eloc + iy, offset_Hloc + j, -poids*curl_phi(j)(1));
		      mat_interac.SetEntry(offset_Eloc + iz, offset_Hloc + j, -poids*curl_phi(j)(2));

		      mat_interac.SetEntry(offset_Hloc + j, offset_Eloc + ix, poids_sym*curl_phi(j)(0));
		      mat_interac.SetEntry(offset_Hloc + j, offset_Eloc + iy, poids_sym*curl_phi(j)(1));
		      mat_interac.SetEntry(offset_Hloc + j, offset_Eloc + iz, poids_sym*curl_phi(j)(2));
		    }
	      }
	  }

	if (mass)
	  {
	    // mass part
	    TinyMatrix<Complexe, Symmetric, 3, 3> epsilon, mu, sigma, Bmass, Amass, Cdamp;
	    Matrix3_3 dfjm1; Real_wp jacob, jacob_weighted;
	    TinyMatrix<Complexe, General, 3, 3> A_tmp;
	    for (int i = 0; i < nb_points_quad; i++)
	      {
		epsilon = this->ref_epsilon(ref_domain).GetCoefficient(var_problem, iquad, i);
		mu = this->ref_mu(ref_domain).GetCoefficient(var_problem, iquad, i);
		sigma = this->ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, i);

		if ((var_problem.InsidePML(iquad)) && (!var_problem.FirstOrderFormulation()))
		  this->ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, i);
		
		var_problem.GetInverseJacobianMatrix(iquad, i, affine, Fb.GetGeometricElement(),
						     dfjm1, jacob, jacob_weighted);
		
		MltTrans(epsilon, dfjm1, A_tmp);
		Mlt(dfjm1, A_tmp, Bmass); 
		Mlt(m_iomega*jacob_weighted, Bmass);

		MltTrans(mu, dfjm1, A_tmp);
		Mlt(dfjm1, A_tmp, Amass); 		
		if (var_computation.GetSymmetrizationUse())
		  Mlt(-m_iomega*jacob_weighted, Amass);
		else
		  Mlt(m_iomega*jacob_weighted, Amass);
		
		MltTrans(sigma, dfjm1, A_tmp);
		Mlt(dfjm1, A_tmp, Cdamp); 
		Mlt(jacob_weighted, Cdamp);

		int ix, iy, iz;
		Fb_hex.GetDofNumber_FromPointNode(i, ix, iy, iz);
		
		mat_interac.SetEntry(offset_Eloc+ix, offset_Eloc+ix, Bmass(0, 0)*m + Cdamp(0, 0)*sig);
		mat_interac.SetEntry(offset_Eloc+ix, offset_Eloc+iy, Bmass(0, 1)*m + Cdamp(0, 1)*sig);
		mat_interac.SetEntry(offset_Eloc+ix, offset_Eloc+iz, Bmass(0, 2)*m + Cdamp(0, 2)*sig);
		mat_interac.SetEntry(offset_Eloc+iy, offset_Eloc+ix, Bmass(1, 0)*m + Cdamp(1, 0)*sig);
		mat_interac.SetEntry(offset_Eloc+iy, offset_Eloc+iy, Bmass(1, 1)*m + Cdamp(1, 1)*sig);
		mat_interac.SetEntry(offset_Eloc+iy, offset_Eloc+iz, Bmass(1, 2)*m + Cdamp(1, 2)*sig);
		mat_interac.SetEntry(offset_Eloc+iz, offset_Eloc+ix, Bmass(2, 0)*m + Cdamp(2, 0)*sig);
		mat_interac.SetEntry(offset_Eloc+iz, offset_Eloc+iy, Bmass(2, 1)*m + Cdamp(2, 1)*sig);
		mat_interac.SetEntry(offset_Eloc+iz, offset_Eloc+iz, Bmass(2, 2)*m + Cdamp(2, 2)*sig);
		
		mat_interac.SetEntry(offset_Hloc+ix, offset_Hloc+ix, Amass(0, 0)*m);
		mat_interac.SetEntry(offset_Hloc+ix, offset_Hloc+iy, Amass(0, 1)*m);
		mat_interac.SetEntry(offset_Hloc+ix, offset_Hloc+iz, Amass(0, 2)*m);
		mat_interac.SetEntry(offset_Hloc+iy, offset_Hloc+ix, Amass(1, 0)*m);
		mat_interac.SetEntry(offset_Hloc+iy, offset_Hloc+iy, Amass(1, 1)*m);
		mat_interac.SetEntry(offset_Hloc+iy, offset_Hloc+iz, Amass(1, 2)*m);
		mat_interac.SetEntry(offset_Hloc+iz, offset_Hloc+ix, Amass(2, 0)*m);
		mat_interac.SetEntry(offset_Hloc+iz, offset_Hloc+iy, Amass(2, 1)*m);
		mat_interac.SetEntry(offset_Hloc+iz, offset_Hloc+iz, Amass(2, 2)*m);
	      }
	  }
      }
    else if (!Fb.UsePiolaTransform())
      {
	// case where Piola transform is not used

	// stiffness part E \cdot \rot psi and -rot H \cdot phi
	// we group the term \int_K rot E \cdot \psi - \int n \times E \cdot \psi
	// in the single term E \cdot rot psi (idem for - rot H \cdot phi)
	VectR3 grad_phi; VectReal_wp val_phi; R3 grad_loc;
	const ElementReference<Dimension3, 1>& elt_H1 = Fb.GetScalarElement();
	int Ns = elt_H1.GetNbDof();
	Matrix3_3 dfjm1; Real_wp jacob, jacob_weighted;
	TinyMatrix<Complexe, Symmetric, 3, 3> epsilon, mu, sigma;
	TinyMatrix<T, Symmetric, 3, 3> Bmass, Amass;
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    var_problem.GetInverseJacobianMatrix(iquad, k, affine, Fb.GetGeometricElement(),
						 dfjm1, jacob, jacob_weighted);

	    elt_H1.GetValuePhiOnQuadraturePoint(k, val_phi);
	    
	    if (stiff)
	      {
		elt_H1.GetGradientPhiOnQuadraturePoint(k, grad_phi);
		for (int i = 0; i < Ns; i++)
		  {
		    grad_loc = grad_phi(i);
		    MltTrans(dfjm1, grad_loc, grad_phi(i));
		  }

		T poids = jacob_weighted*nat_mat.GetCoefStiffness();
		T poids_sym = poids;
		if (var_computation.GetSymmetrizationUse())
		  poids_sym = -poids;

		for (int i = 0; i < Ns; i++)
		  for (int j = 0; j < Ns; j++)
		    {
		      int ix = i, iy = Ns + i, iz = 2*Ns + i;
		      int jx = j, jy = Ns + j, jz = 2*Ns + j;
		      mat_interac.AddInteraction(offset_Eloc+ix, offset_Hloc+jz, -poids*grad_phi(j)(1)*val_phi(i));
		      mat_interac.AddInteraction(offset_Eloc+ix, offset_Hloc+jy, poids*grad_phi(j)(2)*val_phi(i));
		      mat_interac.AddInteraction(offset_Eloc+iy, offset_Hloc+jx, -poids*grad_phi(j)(2)*val_phi(i));
		      mat_interac.AddInteraction(offset_Eloc+iy, offset_Hloc+jz, poids*grad_phi(j)(0)*val_phi(i));
		      mat_interac.AddInteraction(offset_Eloc+iz, offset_Hloc+jy, -poids*grad_phi(j)(0)*val_phi(i));
		      mat_interac.AddInteraction(offset_Eloc+iz, offset_Hloc+jx, poids*grad_phi(j)(1)*val_phi(i));

		      mat_interac.AddInteraction(offset_Hloc+jz, offset_Eloc+ix, poids_sym*grad_phi(j)(1)*val_phi(i));
		      mat_interac.AddInteraction(offset_Hloc+jy, offset_Eloc+ix, -poids_sym*grad_phi(j)(2)*val_phi(i));
		      mat_interac.AddInteraction(offset_Hloc+jx, offset_Eloc+iy, poids_sym*grad_phi(j)(2)*val_phi(i));
		      mat_interac.AddInteraction(offset_Hloc+jz, offset_Eloc+iy, -poids_sym*grad_phi(j)(0)*val_phi(i));
		      mat_interac.AddInteraction(offset_Hloc+jy, offset_Eloc+iz, poids_sym*grad_phi(j)(0)*val_phi(i));
		      mat_interac.AddInteraction(offset_Hloc+jx, offset_Eloc+iz, -poids_sym*grad_phi(j)(1)*val_phi(i));
		    }
	      }

	    if (mass)
	      {
		epsilon = this->ref_epsilon(ref_domain).GetCoefficient(var_problem, iquad, k);
		mu = this->ref_mu(ref_domain).GetCoefficient(var_problem, iquad, k);
		sigma = this->ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, k);
		
		if ((var_problem.InsidePML(iquad)) && (!var_problem.FirstOrderFormulation()))
		  this->ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, k);

		T coef1 = m*m_iomega*jacob_weighted;
		T coef2 = sig*jacob_weighted;
		Bmass = coef1*epsilon + coef2*sigma;

		if (var_computation.GetSymmetrizationUse())
		  coef1 = -coef1;
		
		Amass = coef1*mu;

		for (int i = 0; i < Ns; i++)
		  for (int j = 0; j < Ns; j++)
		    {
		      int ix = i, iy = Ns + i, iz = 2*Ns + i;
		      int jx = j, jy = Ns + j, jz = 2*Ns + j;
		      Real_wp prod = val_phi(i)*val_phi(j);
		      mat_interac.AddInteraction(offset_Eloc+ix, offset_Eloc+jx, Bmass(0, 0)*prod);
		      mat_interac.AddInteraction(offset_Eloc+ix, offset_Eloc+jy, Bmass(0, 1)*prod);
		      mat_interac.AddInteraction(offset_Eloc+ix, offset_Eloc+jz, Bmass(0, 2)*prod);
		      mat_interac.AddInteraction(offset_Eloc+iy, offset_Eloc+jx, Bmass(1, 0)*prod);
		      mat_interac.AddInteraction(offset_Eloc+iy, offset_Eloc+jy, Bmass(1, 1)*prod);
		      mat_interac.AddInteraction(offset_Eloc+iy, offset_Eloc+jz, Bmass(1, 2)*prod);
		      mat_interac.AddInteraction(offset_Eloc+iz, offset_Eloc+jx, Bmass(2, 0)*prod);
		      mat_interac.AddInteraction(offset_Eloc+iz, offset_Eloc+jy, Bmass(2, 1)*prod);
		      mat_interac.AddInteraction(offset_Eloc+iz, offset_Eloc+jz, Bmass(2, 2)*prod);

		      mat_interac.AddInteraction(offset_Hloc+ix, offset_Hloc+jx, Amass(0, 0)*prod);
		      mat_interac.AddInteraction(offset_Hloc+ix, offset_Hloc+jy, Amass(0, 1)*prod);
		      mat_interac.AddInteraction(offset_Hloc+ix, offset_Hloc+jz, Amass(0, 2)*prod);
		      mat_interac.AddInteraction(offset_Hloc+iy, offset_Hloc+jx, Amass(1, 0)*prod);
		      mat_interac.AddInteraction(offset_Hloc+iy, offset_Hloc+jy, Amass(1, 1)*prod);
		      mat_interac.AddInteraction(offset_Hloc+iy, offset_Hloc+jz, Amass(1, 2)*prod);
		      mat_interac.AddInteraction(offset_Hloc+iz, offset_Hloc+jx, Amass(2, 0)*prod);
		      mat_interac.AddInteraction(offset_Hloc+iz, offset_Hloc+jy, Amass(2, 1)*prod);
		      mat_interac.AddInteraction(offset_Hloc+iz, offset_Hloc+jz, Amass(2, 2)*prod);
		    }		
	      }
	  }
      }
    else
      {
	HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass = this->Glob_matMass_elem(iquad);
	bool diag_mass = mass.IsDiagonalMass();

	//DISP(diag_mass); DISP(var_problem.GetOmega());
	// mass terms : epsilon dE/dt + sigma E and mu dH/dt
	if (diag_mass)
	  {
	    for (int i = 0; i < Fb.GetNbDof(); i++)
	      {		
		mat_interac.SetEntry(offset_Eloc + i, offset_Eloc + i, nat_mat.GetCoefMass()*mass.GetDh(i)
                                     + nat_mat.GetCoefDamping()*mass.GetDhSigma(i));

		if (var_computation.GetSymmetrizationUse())
		  mat_interac.SetEntry(offset_Hloc + i, offset_Hloc + i, -nat_mat.GetCoefMass()*mass.GetDhStiff(i));
		else
		  mat_interac.SetEntry(offset_Hloc + i, offset_Hloc + i, nat_mat.GetCoefMass()*mass.GetDhStiff(i));
	      }
	  }
	else
	  {
	    Vector<TinyMatrix<T, General, 3, 3> > Amass(nb_points_quad);
	    for (int k = 0; k < nb_points_quad; k++)
	      Amass(k) = nat_mat.GetCoefMass()*mass.GetBh(k) + nat_mat.GetCoefDamping()*mass.GetBhSigma(k);

	    Fb.AddVariableMassMatrix(offset_Eloc, offset_Eloc, Amass, mat_interac);
	    
	    if (var_computation.GetSymmetrizationUse())
	      for (int k = 0; k < nb_points_quad; k++)
		Amass(k) = -nat_mat.GetCoefMass()*mass.GetAh(k);
	    else
	      for (int k = 0; k < nb_points_quad; k++)
		Amass(k) = nat_mat.GetCoefMass()*mass.GetAh(k);

	    Fb.AddVariableMassMatrix(offset_Hloc, offset_Hloc, Amass, mat_interac);
	  }
	
	VectReal_wp Ones(Fb.GetNbDof()), contrib(Fb.GetNbDof());
	VectReal_wp val_phi(3*nb_points_quad);
	val_phi.Zero(); contrib.Zero();
	// stiffness terms : E rot(psi) et - rot(H) phi
	for (int i = 0; i < Fb.GetNbDof(); i++)
	  {
	    //Fb.GetValueSinglePhiQuadrature(i, val_phi);
	    Ones.Zero(); Ones(i) = 1.0;
	    Fb.ApplyChTranspose(Ones, val_phi);
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		val_phi(3*k) *= Fb.WeightsND(k);
		val_phi(3*k+1) *= Fb.WeightsND(k);
		val_phi(3*k+2) *= Fb.WeightsND(k);
	      }
	    
	    Fb.ApplyRh(val_phi, contrib);
	    for (int j = 0; j < Fb.GetNbDof(); j++)
	      {
		mat_interac.AddInteraction(offset_Hloc + j, offset_Eloc + i, -s*contrib(j));
		if (var_computation.GetSymmetrizationUse())
		  mat_interac.AddInteraction(offset_Eloc + i, offset_Hloc + j, -s*contrib(j));
		else
		  mat_interac.AddInteraction(offset_Eloc + i, offset_Hloc + j, s*contrib(j));
	      }
	  }
      }

    if ((!stiff) || (var_problem.FormulationDG() != ElementReference_Base::HDG))
      return;
    
    // surface integrals
    Real_wp dsj; Complexe phase, phase_conj; 
    Vector<TinyMatrix<Real_wp, Symmetric, 2, 2> > mat_DF_DFt;
    Matrix3_3 dfjm1, dfj_trans; Real_wp jacob, jacob_weighted;
    TinyMatrix<Real_wp, Symmetric, 3, 3> DF_DFt;
    TinyMatrix<Real_wp, Symmetric, 2, 2> sDF_DFt;
    
    int offset_L = 0;
    // loop over boundaries of the element K
    for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
      {
	int num_face, ref_boundary, rf, rot; bool new_face;
	const Matrix<int>& FacesQuadRotation =
	  var_boundary.GetGeometryPhaseData(iquad, num_loc,
					    num_face, ref_boundary, rf, new_face, rot,
					    phase, phase_conj);
	
	const ElementReference<Dimension2, 2>& Fb_s = dynamic_cast<const ElementReference<Dimension2, 2>& >(var_problem.GetSurfaceFiniteElement(num_face));
	int Nquad = Fb_s.GetNbPointsQuadratureInside();
	
	if (Fb.LumpedMassMatrix())
	  {
	    // for hexas only
	    int nx = 0, ny = 1;
	    if ((num_loc == 0) || (num_loc == 5))
	      { nx = 1; ny = 2; }
	    else if ((num_loc == 1) || (num_loc == 4))
	      { nx = 0; ny = 2; }	
	    
	    mat_DF_DFt.Reallocate(Nquad);	    
	    for (int k = 0; k < Fb_s.GetNbPointsQuadratureInside(); k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  dsj = var_problem.Glob_dsj(num_face)(k);
		else
		  dsj = var_problem.Glob_dsj(num_face)(krot);
		
		int npoint = Fb.GetQuadNumber(num_loc, k);
		var_problem.GetInverseJacobianMatrix(iquad, npoint, affine,
						     Fb.GetGeometricElement(),
						     dfjm1, jacob, jacob_weighted);
		
		GetInverse(dfjm1, dfj_trans); Transpose(dfj_trans);
		MltTrans(dfj_trans, dfj_trans,  DF_DFt);
		Mlt(Fb_s.WeightsND(k)/dsj, DF_DFt);
		
		mat_DF_DFt(k)(0, 0) = DF_DFt(ny, ny);
		mat_DF_DFt(k)(0, 1) = -DF_DFt(nx, ny);
		mat_DF_DFt(k)(1, 1) = DF_DFt(nx, nx);
	      }
	    
	    const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	    const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);

	    int order = Fb_quad.GetOrder();
	    for (int k = 0; k < Fb_s.GetNbPointsQuadratureInside(); k++)
	      {
		int npoint = Fb.GetQuadNumber(num_loc, k);
		int kx, ky, kz;
		Fb_hex.GetDofNumber_FromPointNode(npoint, kx, ky, kz);
		int jx = kx, jy = ky;
		if ((num_loc == 0) || (num_loc == 5))
		  { jx = ky; jy = kz; }
		else if ((num_loc == 1) || (num_loc == 4))
		  { jx = kx; jy = kz; }
		
		sDF_DFt = mat_DF_DFt(k);
		
		// part coef_tau (n \times E) \cdot (n \times phi) 
		mat_interac.AddInteraction(offset_Eloc+jx, offset_Eloc+jx, coef_tau*sDF_DFt(0, 0));
		mat_interac.AddInteraction(offset_Eloc+jx, offset_Eloc+jy, coef_tau*sDF_DFt(0, 1));
		mat_interac.AddInteraction(offset_Eloc+jy, offset_Eloc+jx, coef_tau*sDF_DFt(1, 0));
		mat_interac.AddInteraction(offset_Eloc+jy, offset_Eloc+jy, coef_tau*sDF_DFt(1, 1));
				
		kx = Fb_quad.GetXdofNumber(k);
		ky = Fb_quad.GetYdofNumber(k);
		Real_wp sx = 1, sy = 1;
		if ((kx >= 2*(order+1)) && (kx < 4*(order+1)))
		  sx = -sx;
		
		if ((ky >= 2*(order+1)) && (ky < 4*(order+1)))
		  sy = -sy;
		
		// part -coef_tau (n \times lambda) \cdot (n \times phi) 
		mat_interac.AddInteraction(offset_Eloc+jx, offset_L + kx, -coef_tau*sx*sDF_DFt(0, 0));
		mat_interac.AddInteraction(offset_Eloc+jx, offset_L + ky, -coef_tau*sy*sDF_DFt(0, 1));
		mat_interac.AddInteraction(offset_Eloc+jy, offset_L + kx, -coef_tau*sx*sDF_DFt(1, 0));
		mat_interac.AddInteraction(offset_Eloc+jy, offset_L + ky, -coef_tau*sy*sDF_DFt(1, 1));
		
		T coef_taus = coef_tau;
		if (var_computation.GetSymmetrizationUse())
		  coef_taus = -coef_taus;
		
		// part coef_tau (n \times q) \cdot (n \times E)
		mat_interac.AddInteraction(offset_L + kx, offset_Eloc+jx, coef_taus*sx*sDF_DFt(0, 0));
		mat_interac.AddInteraction(offset_L + ky, offset_Eloc+jx, coef_taus*sy*sDF_DFt(0, 1));
		mat_interac.AddInteraction(offset_L + kx, offset_Eloc+jy, coef_taus*sx*sDF_DFt(1, 0));
		mat_interac.AddInteraction(offset_L + ky, offset_Eloc+jy, coef_taus*sy*sDF_DFt(1, 1));
		
		if ((ref_boundary == 0) ||
		    (mesh.GetBoundaryCondition(ref_boundary) != BoundaryConditionEnum::LINE_DIRICHLET))
		  {
		    // following contributions are not present for Dirichlet dofs
		    T coef = s*Fb_s.WeightsND(k);
		    if ((num_loc == 0) || (num_loc == 2) || (num_loc == 4))
		      coef = -coef;
		    
		    if (var_computation.GetSymmetrizationUse())
		      coef = -coef;
		    
		    // part -n \times \psi \cdot \lambda
		    mat_interac.AddInteraction(offset_Hloc+jx, offset_L+ky, -sy*coef);
		    mat_interac.AddInteraction(offset_Hloc+jy, offset_L+kx, sx*coef);
		    
		    // part -n \times H \cdot q
		    mat_interac.AddInteraction(offset_L+ky, offset_Hloc+jx, -sy*coef);
		    mat_interac.AddInteraction(offset_L+kx, offset_Hloc+jy, sx*coef);
		  }

		this->ModifyCoefficientBC_HDG(mesh, ref_boundary, k,
					      nat_mat.GetCoefStiffness(), coef_H,
					      ref_domain, iquad, num_loc, coef_taus);
		
		// part - coef_tau n \times lambda \cdot n \times q
		Real_wp sxy = sx*sy;
		mat_interac.AddInteraction(offset_L + kx, offset_L+kx, -coef_taus*sDF_DFt(0, 0));
		mat_interac.AddInteraction(offset_L + ky, offset_L+kx, -coef_taus*sxy*sDF_DFt(0, 1));
		mat_interac.AddInteraction(offset_L + kx, offset_L+ky, -coef_taus*sxy*sDF_DFt(1, 0));
		mat_interac.AddInteraction(offset_L + ky, offset_L+ky, -coef_taus*sDF_DFt(1, 1));
	      }
	  }
	else if (!Fb.UsePiolaTransform())
	  {
	    VectR3 vert;
	    mesh.GetVerticesElement(iquad, vert);

	    SetPoints<Dimension3> PointsElem;
	    SetMatrices<Dimension3> MatricesElem;
	    Fb.FjElem(vert, PointsElem, mesh, iquad);
	    Fb.DFjElem(vert, PointsElem, MatricesElem, mesh, iquad);

	    Fb.FjSurfaceElem(vert, PointsElem, mesh, iquad, num_loc);
	    Fb.DFjSurfaceElem(vert, PointsElem, MatricesElem, mesh, iquad, num_loc);

	    R3 normale; VectR3 val_phi, val_L; VectR2 val_phi2D;
	    for (int k = 0; k < Fb_s.GetNbPointsQuadratureInside(); k++)
	      {
		int krot = FacesQuadRotation(rot, k);

		if (new_face)
		  {
		    normale = var_problem.Glob_normale(num_face)(k);
		    dsj = var_problem.Glob_dsj(num_face)(k);
		  }
		else
		  {
		    normale = var_problem.Glob_normale(num_face)(krot);
		    dsj = var_problem.Glob_dsj(num_face)(krot);
		    normale = -normale;
		  }

		int npoint = Fb.GetQuadNumber(num_loc, k);
		Fb.GetValuePhiOnQuadraturePoint(npoint, val_phi);
		
		dfjm1 = MatricesElem.GetPointQuadratureBoundary(k);
		GetInverse(dfjm1);
		
		Fb_s.GetValuePhiOnQuadraturePoint(k, val_phi2D);
		val_L.Reallocate(Fb_s.GetNbDof());
		for (int j = 0; j < Fb_s.GetNbDof(); j++)
		  {
		    R3 phi_loc = Fb.GetTangentialVector(num_loc, val_phi2D(j));
		    MltTrans(dfjm1, phi_loc, val_L(j));
		  }

		R3 vec_u, vec_v;
		for (int i = 0; i < nb_dof_elt; i++)
		  if (!val_phi(i).IsZero())
		    for (int j = 0; j < nb_dof_elt; j++)
		      if (!val_phi(j).IsZero())
			{
			  // part tau n \times E \cdot n \times \phi
			  TimesProd(normale, val_phi(i), vec_u);
			  TimesProd(normale, val_phi(j), vec_v);
			  T vloc = coef_tau*DotProd(vec_u, vec_v)*dsj*Fb_s.WeightsND(k);
			  mat_interac.AddInteraction(offset_Eloc+i, offset_Eloc+j, vloc);
			}

		for (int i = 0; i < nb_dof_elt; i++)
		  if (!val_phi(i).IsZero())
		    for (int j = 0; j < Fb_s.GetNbDof(); j++)
		      if (!val_phi2D(j).IsZero())
			{
			  // part -tau n \times lambda \cdot n \times \phi
			  TimesProd(normale, val_phi(i), vec_u);
			  TimesProd(normale, val_L(j), vec_v);
			  T vloc = -coef_tau*DotProd(vec_u, vec_v)*dsj*Fb_s.WeightsND(k);
			  mat_interac.AddInteraction(offset_Eloc+i, offset_L+j, vloc);
			  
			  // part tau n \times q \cdot n \times E
			  if (!var_computation.GetSymmetrizationUse())
			    vloc = -vloc;
			  
			  mat_interac.AddInteraction(offset_L+j, offset_Eloc+i, vloc);
			  
			  // part n \times \lambda \cdot psi
			  vloc = DotProd(vec_v, val_phi(i))*dsj*Fb_s.WeightsND(k);
			  if (var_computation.GetSymmetrizationUse())
			    vloc = -vloc;
			  
			  mat_interac.AddInteraction(offset_Hloc+i, offset_L+j, vloc);
			  
			  // part n \times q \cdot H
			  mat_interac.AddInteraction(offset_L+j, offset_Hloc+i, vloc);
			}
		
		for (int i = 0; i < Fb_s.GetNbDof(); i++)
		  if (!val_phi2D(i).IsZero())
		    for (int j = 0; j < Fb_s.GetNbDof(); j++)
		      if (!val_phi2D(j).IsZero())
			{
			  TimesProd(normale, val_L(i), vec_u);
			  TimesProd(normale, val_L(j), vec_v);
			  T vloc = -coef_tau*DotProd(vec_u, vec_v)*dsj*Fb_s.WeightsND(k);
			  if (var_computation.GetSymmetrizationUse())
			    vloc = -vloc;
			  
			  mat_interac.AddInteraction(offset_L+i, offset_L+j, vloc);
			}
	      }
	    
	  }
	else
	  {
	    VectR3 vert;
	    mesh.GetVerticesElement(iquad, vert);
	    
	    SetPoints<Dimension3> PointsElem;
	    SetMatrices<Dimension3> MatricesElem;
	    Fb.FjElem(vert, PointsElem, mesh, iquad);
	    Fb.DFjElem(vert, PointsElem, MatricesElem, mesh, iquad);

	    Fb.FjSurfaceElem(vert, PointsElem, mesh, iquad, num_loc);
	    Fb.DFjSurfaceElem(vert, PointsElem, MatricesElem, mesh, iquad, num_loc);

	    int Nquad = Fb_s.GetNbPointsQuadratureInside();
	    Vector<Matrix3_3sym> mat3_DF_DFt(Nquad);
	    Matrix3_3 invDs_DFi;
	    for (int k = 0; k < Nquad; k++)
	      {
		invDs_DFi = MatricesElem.GetPointQuadratureBoundary(k);
		Mlt(Fb_s.WeightsND(k)/MatricesElem.GetDsQuadratureBoundary(k), invDs_DFi);
		Transpose(invDs_DFi);
		Mlt(invDs_DFi, MatricesElem.GetPointQuadratureBoundary(k), mat3_DF_DFt(k));
	      }
	    
	    // loop on dofs for volume unknowns phi
	    VectReal_wp Ones(Fb.GetNbDof()); Ones.Zero();
	    VectReal_wp phi(3*Fb_s.GetNbPointsQuadratureInside());
	    VectReal_wp feval(3*Fb_s.GetNbPointsQuadratureInside());
	    R3 vec_u, vec_v; R2 vec_u2;
	    VectReal_wp feval2D(2*Fb_s.GetNbPointsQuadratureInside());
	    VectReal_wp fevalHn(2*Fb_s.GetNbPointsQuadratureInside());
	    VectReal_wp contrib2D(Fb_s.GetNbDof()); contrib2D.Zero();
	    for (int j = 0; j < Fb.GetNbDof(); j++)
	      if (Fb.DiscontinuousElement() || Fb.IsTangentialDof(j, num_loc))
		{
		  Ones.Zero(); Ones(j) = 1.0;
		  Fb.ApplyShTranspose(num_loc, Ones, phi, rf);

		  // we compute phi = n \times (1/ds DF_i DF_i^* (E \times n))
		  for (int k = 0; k < Nquad; k++)
		    {
		      vec_u.Init(phi(3*k), phi(3*k+1), phi(3*k+2));
		      TimesProd(vec_u, Fb.NormaleLoc(num_loc), vec_v);
		      vec_u2 = Fb.TransposeTangentialVector(num_loc, vec_v);
		      vec_u2 *= Fb_s.WeightsND(k);
		      fevalHn(2*k) = vec_u2(0); fevalHn(2*k+1) = vec_u2(1);
		      
		      Mlt(mat3_DF_DFt(k), vec_v, vec_u);
		      TimesProd(Fb.NormaleLoc(num_loc), vec_u, vec_v);
		      feval(3*k) = vec_v(0);
		      feval(3*k+1) = vec_v(1); feval(3*k+2) = vec_v(2);
		      vec_u2 = Fb.TransposeTangentialVector(num_loc, vec_v);
		      feval2D(2*k) = vec_u2(0); feval2D(2*k+1) = vec_u2(1);
		    }
		  
		  // part coef_tau (n \times E) \cdot (n \times phi)
		  Ones.Zero();
		  Fb.ApplySh(Real_wp(1.0), num_loc, feval, Ones, rf);
		  
		  for (int k = 0; k < Fb.GetNbDof(); k++)
		    if (abs(Ones(k)) > var_computation.GetThresholdMatrix())
		      mat_interac.AddInteraction(offset_Eloc+j, offset_Eloc+k, coef_tau*s*Ones(k));
		  
		  // part -coef_tau (n \times lambda) \cdot (n \times phi)
		  // part coef_tau (n \times q) \cdot (n \times E)
		  T coef_taus = coef_tau;
		  if (var_computation.GetSymmetrizationUse())
		    coef_taus = -coef_taus;

		  Fb_s.ApplyCh(feval2D, contrib2D);
		  for (int k = 0; k < Fb_s.GetNbDof(); k++)
		    if (abs(contrib2D(k)) > var_computation.GetThresholdMatrix())
		      {
			mat_interac.AddInteraction(offset_Eloc + j, offset_L+k, -coef_tau*contrib2D(k));
			mat_interac.AddInteraction(offset_L+k, offset_Eloc+j, coef_taus*contrib2D(k));
		      }

		  if ((ref_boundary == 0) ||
		      (mesh.GetBoundaryCondition(ref_boundary) != BoundaryConditionEnum::LINE_DIRICHLET))
		    {
		      // following contributions are not present for Dirichlet dofs
		      T coef = s;
		      if (var_computation.GetSymmetrizationUse())
			coef = -coef;
		      
		      Fb_s.ApplyCh(fevalHn, contrib2D);
		      for (int k = 0; k < Fb_s.GetNbDof(); k++)
			if (abs(contrib2D(k)) > var_computation.GetThresholdMatrix())
			  {
			    mat_interac.AddInteraction(offset_Hloc+j, offset_L+k, coef*contrib2D(k));
			    mat_interac.AddInteraction(offset_L+k, offset_Hloc+j, coef*contrib2D(k));
			  }
		    }		  
		}

	    Ones.Reallocate(Fb_s.GetNbDof());
	    Vector<T> feval_cplx(2*Nquad), contrib_cplx(Fb_s.GetNbDof());	    
	    for (int j = 0; j < Fb_s.GetNbDof(); j++)
	      {
		Ones.Zero(); Ones(j) = 1.0;
		Fb_s.ApplyChTranspose(Ones, feval2D);

		// we compute phi = n \times (1/ds DF_i DF_i^* (\lambda \times n))
		for (int k = 0; k < Nquad; k++)
		  {
		    vec_u2.Init(feval2D(2*k), feval2D(2*k+1));
		    vec_u = Fb.GetTangentialVector(num_loc, vec_u2);		    
		    TimesProd(vec_u, Fb.NormaleLoc(num_loc), vec_v);
		    Mlt(mat3_DF_DFt(k), vec_v, vec_u);
		    TimesProd(Fb.NormaleLoc(num_loc), vec_u, vec_v);
		    vec_u2 = Fb.TransposeTangentialVector(num_loc, vec_v);

		    T coef_taus = coef_tau;
		    if (var_computation.GetSymmetrizationUse())
		      coef_taus = -coef_taus;
		    
		    this->ModifyCoefficientBC_HDG(mesh, ref_boundary, k,
						  nat_mat.GetCoefStiffness(), coef_H,
						  ref_domain, iquad, num_loc, coef_taus);

		    feval_cplx(2*k) = vec_u2(0)*coef_taus;
		    feval_cplx(2*k+1) = vec_u2(1)*coef_taus;
		  }
		
		Fb_s.ApplyCh(feval_cplx, contrib_cplx);
		for (int k = 0; k < Fb_s.GetNbDof(); k++)
		  if (abs(contrib_cplx(k)) > var_computation.GetThresholdMatrix())
		    mat_interac.AddInteraction(offset_L + j, offset_L+k, -contrib_cplx(k));
	      }
	    
	  }
	
	offset_L += Fb_s.GetNbDof();
      }
    
    //mat_interac.Write("mat_elem.dat");
    //int test_input; cout << "waiting" << endl; cin >> test_input;
	
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad, 1);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad, 1);
  }


  //! modification of coefficient due to boundary conditions
  template<class Complexe> template<class T>
  void HarmonicMaxwell_3D<Complexe>
  ::ModifyCoefficientBC_HDG(const Mesh<Dimension3>& mesh, int ref_boundary, int k,
			    const T& s, const Complexe& coef_v,
			    int ref_domain, int iquad, int num_loc, T& coef_abc) const
  {
    if (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_ABSORBING)
      coef_abc *= Real_wp(2);
    else if (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
        Complexe coef_imped; SetComplexOne(coef_imped);
	var_boundary.MltParamCondition(ref_boundary, 0, coef_imped);
	var_boundary.MltMuIntegrationByParts(ref_domain, iquad, num_loc, k, coef_imped);
	coef_abc += -coef_imped/coef_v*s;
      }
  }


  template<class Complexe>
  void HarmonicMaxwell_3D<Complexe>::ComputeTauCoefficient()
  {
    if (var_problem.FormulationDG() != ElementReference_Base::HDG)
      return;

    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    var_problem.Glob_CoefPenalDG.Zero();
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
        int ref = var_problem.mesh.Element(i).GetReference();
        if (this->ref_epsilon(ref).IsVarying() || this->ref_mu(ref).IsVarying())
          {
            int offset = var_problem.GetNbPointsQuadratureInside(i);
            for (int num_loc = 0; num_loc < var_problem.mesh.Element(i).GetNbBoundary(); num_loc++)
              {
                int nf = var_problem.mesh.Element(i).numBoundary(num_loc);
                int nb_pts_face = mesh_num.GetNbPointsQuadratureBoundary(nf);
                Real_wp epsilon, mu;
                for (int k = 0; k < nb_pts_face; k++)
                  {
                    epsilon = abs(ref_epsilon(ref).GetCoefficient(var_problem, i, offset+k)(0, 0));
                    mu = abs(ref_mu(ref).GetCoefficient(var_problem, i, offset+k)(0, 0));                    
                    var_problem.Glob_CoefPenalDG(nf) = max(var_problem.Glob_CoefPenalDG(nf), sqrt(epsilon/mu));
                  }

                offset += nb_pts_face;
              }
          }
        else
          for (int num_loc = 0; num_loc < var_problem.mesh.Element(i).GetNbBoundary(); num_loc++)
            {
              int nf = var_problem.mesh.Element(i).numBoundary(num_loc);
              var_problem.Glob_CoefPenalDG(nf) = max(var_problem.Glob_CoefPenalDG(nf),
                                                     abs(this->coefficient_impedance_absorbing(ref)));
            }
      }
  }


  //! Adds numerical fluxes for sipg formulation
  template<class Complexe> template<class T>
  void HarmonicMaxwell_3D<Complexe>
  ::AddElementaryFluxesSipg(VirtualMatrix<T>& mat_sp,
			    const GlobalGenericMatrix<T>& nat_mat,
			    int offset_row, int offset_col)
  {
    const Mesh<Dimension3>& mesh = var_problem.mesh;
    int nb_neighbor = 0;
    
    // loop over all the faces
    // for sake of simplicity, we do a loop on elements and then on local faces
    // Hence, faces may be selected twice (because they can belong to two elements)
    // but only interactions 11 and 12 are added (and not the interactions 21 and 22)
    // if 1 and 2 refers to the element 1 and element 2 adjoining the face
    // so that each value is added only once
    for (int iquad = 0; iquad < mesh.GetNbElt(); iquad++)
      for (int num_pos1_face = 0; num_pos1_face < mesh.Element(iquad).GetNbBoundary();
           num_pos1_face++)
	{
	  int num_face = mesh.Element(iquad).numBoundary(num_pos1_face);
          
          if (var_problem.FaceHasToBeConsideredForBoundaryIntegral(num_face))
            {
              int num_elem2 = mesh.Boundary(num_face).numElement(0);
              int ref = mesh.Boundary(num_face).GetReference();
              bool new_face = var_problem.IsNewFace(iquad)(num_pos1_face);
              if ((num_elem2 == iquad)&&(mesh.Boundary(num_face).GetNbElements()==2))
                num_elem2 = var_problem.mesh.Boundary(num_face).numElement(1);

	      bool neighbor_face = (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR);
	      
	      if (!neighbor_face)
		{
		  // for second order formulation, boundary condition is set with AddMatrixImpedanceBoundary
		  if ( (mesh.GetBoundaryCondition(ref) != 0) && (!var_problem.FirstOrderFormulationDG()))
		    continue;
		}
	      
              // AddElementaryFluxDG is called for the selected face
	      AddElementaryFluxSipg(mat_sp, nat_mat, iquad, num_pos1_face, num_elem2,
				    num_face, ref, new_face, nb_neighbor,
				    offset_row, offset_col, var_problem.GetReferenceElementHcurl(iquad));
            }
        }
  }
  

  //! Adding boundary integrals for Symmetric Interior Penalty Galerkin and Maxwell's equation
  /*!
    adding surface term  \f$ + \displaystyle \sum_{\mbox{e edge or face of the mesh} }
    \int_{Gamma_e} \{ 1/mu curl u \} [ n \times \varphi ]  + \{ 1/mu curl phi \} [ n \times u ] 
    - 1/2 alpha [ n \times u ] [n \times \varphi]  \f$
    where e is an edge shared by elements K- and K+, n is the outward normale of K-,
    and \{ u \} = 1/2 (u+ + u- \} ,  [ u ] = (u+ - u-) 
    \param[in] mat_sp the terms are added to this sparse matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] vars considered problem
    \param[in] iquad number of the first element
    \param[in] num_pos1_face local face number in the first element
    \param[in] num_elem2 number of the second element (if existing)
    \param[in] num_face global face number
    \param[in] ref reference of the face
    \param[in] new_face is the face a new face ? (first time that this face is selected)
    \param[in] nb_neighbor for parallel execution
    \param[in] Fb finite element associated to the first element
    \param[in] Fb2 finite element associated to the second element
  */
  template<class Complexe> template<class T>
  void HarmonicMaxwell_3D<Complexe>
  ::AddElementaryFluxSipg(VirtualMatrix<T>& mat_sp,
			  const GlobalGenericMatrix<T>& nat_mat,
			  int iquad, int num_pos1_face, int num_elem2,
			  int num_face, int ref, bool new_face, int& nb_neighbor,
			  int offset_row, int offset_col,
			  const ElementReference<Dimension3, 2>& Fb)
  {
    const Mesh<Dimension3>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    
    int rf = mesh_num.GetOrderQuadrature(num_face);
    int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
    const VectReal_wp& PoidsFlux = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
    
    bool affine = mesh.IsElementAffine(iquad);
    bool affine_elt2 = mesh.IsElementAffine(num_elem2);
        
    int ref_d = mesh.Element(iquad).GetReference();
    TinyVector<T, 3> vec_u, vec_v, vec_w;
    IVect Nodle = var_problem.GetDofNumberOnElement(iquad);
    IVect Nodle2 = var_problem.GetDofNumberOnElement(num_elem2);
    
    Complexe phase(1), phase_conj(1), cone(1);
    int rot1 = mesh.Element(iquad).GetOrientationBoundary(num_pos1_face), rot2 = 0;

    int cond = mesh.GetBoundaryCondition(ref);
    bool face_on_gamma = false;
    if (cond != BoundaryConditionEnum::LINE_INSIDE)
      {
        if (cond != BoundaryConditionEnum::LINE_NEIGHBOR)
          face_on_gamma = true;
      }

    int num_pos2_face = -1;
    int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0; //int ref2 = -1;
    if (num_elem2 != iquad)
      {
	num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(num_face);
	if (num_pos2_face < 0)
	  {
	    int nf2 = mesh_num.GetPeriodicBoundary(num_face);
            var_boundary.GetPeriodicPhase(num_face, phase);
            var_boundary.GetPeriodicPhase(nf2, phase_conj);
            num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(nf2);
	  }

	//ref2 = mesh.Element(num_elem2).GetReference();
	rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
	rot = mesh_num.GetRotationFace(rot1, rot2, nv);
      }
    
    int krot, num_point = Fb.GetNbPointsQuadratureInside();
    for (int k = 0; k < num_pos1_face; k++)
      num_point += mesh_num.GetNbPointsQuadratureBoundary(mesh.Element(iquad).numBoundary(k));

    const ElementReference<Dimension3, 2>* Fb2_ptr = &var_problem.GetReferenceElementHcurl(num_elem2);
    bool neighbor_face = (cond == BoundaryConditionEnum::LINE_NEIGHBOR);

    IVect NumRotQuad(nb_points_face);
    int offset_neighbor = 0;
#ifdef SELDON_WITH_MPI
    int proc2 = -1;
    if (neighbor_face)
      {
	int type_per = mesh_num.GetPeriodicityTypeForBoundary(num_face);
	if (type_per >= 0)
	  {
	    var_boundary.GetPeriodicPhase(num_face, phase);
	    phase_conj = conjugate(phase);
	  }

	offset_neighbor = var_problem.GetOffsetNeighboringFace(num_face);
        proc2 = var_problem.GetProcessorNeighboringFace(num_face);
        num_pos2_face = var_problem.GetLocalPositionNeighboringFace(num_face);
        int pos_loc_neighbor = mesh_num.GetLocalEdgeNumberNeighborElement(num_face);
        Fb2_ptr = dynamic_cast<const ElementReference<Dimension3, 2>* >(&var_problem.GetNeighborReferenceElement(pos_loc_neighbor));
        rot = var_problem.GetRotationNeighboringFace(num_face);
        Nodle2 = var_problem.GetNodleNeighboringFace(num_face);
        //ref2 = var_problem.GetRefDomainNeighboringFace(num_face);
        NumRotQuad.Reallocate(nb_points_face);
      }
#endif

    const ElementReference<Dimension3, 2>& Fb2 = *Fb2_ptr;
    
    int num_point2 = Fb2.GetNbPointsQuadratureInside();
    for (int k = 0; k < num_pos2_face; k++)
      num_point2 += mesh_num.GetNbPointsQuadratureBoundary(mesh.Element(num_elem2).numBoundary(k));
    
    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
    
    // loop over quadrature points of the face
    Vector<TinyMatrix<T, General, 3, 3> > MatDfCurl;
    Vector<Matrix3_3> MatDfTangentIntra, MatDfTangentExtra;
    MatDfCurl.Reallocate(nb_points_face);
    MatDfTangentIntra.Reallocate(nb_points_face);
    MatDfTangentExtra.Reallocate(nb_points_face);
    Real_wp dsj, jacob; R3 normale; Matrix3_3 dfjm1, mat_dfj, ProjTimesN;
    Vector<T> Poids(nb_points_face); Vector<T> PoidsPenal(nb_points_face);
    TinyMatrix<Complexe, Symmetric, 3, 3> epsilon, mu, invMu, sigma;
    int i1 = iquad - mesh.GetNbElt() + var_boundary.GetNbEltPML();    
    for (int k = 0; k < nb_points_face; k++)
      {
        krot = FacesQuadRotation(rot, k);
        
        if (new_face)
          {
            normale = var_problem.Glob_normale(num_face)(k);
            dsj = var_problem.Glob_dsj(num_face)(k);
          }
        else
          {
            normale = var_problem.Glob_normale(num_face)(krot);
            dsj = var_problem.Glob_dsj(num_face)(krot);
            Mlt(Real_wp(-1), normale);
          }

        // weights for boundary integrals
        Poids(k) = dsj*PoidsFlux(k)*nat_mat.GetCoefStiffness();
        
        // retrieving mu and invMu on the quadrature point
        int k2 = Fb.GetQuadNumber(num_pos1_face, k);
	epsilon = this->ref_epsilon(ref_d).GetCoefficient(var_problem, iquad, k2);
	mu = this->ref_mu(ref_d).GetCoefficient(var_problem, iquad, k2);
	sigma = this->ref_sigma(ref_d).GetCoefficient(var_problem, iquad, k2);
        if (var_problem.InsidePML(iquad))
          {
            this->ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, k2);
            GetInverse(mu, invMu); 
          }
        else
          {
            if (mu.IsZero())
              invMu.Zero();
            else
              GetInverse(mu, invMu); 
          }
	
        // retrieving DFi
        if (affine)
          {
            jacob = var_problem.Glob_jacobian(iquad)(0);
            dfjm1 = var_problem.Glob_DFjm1(iquad)(0);
            Mlt(1.0/jacob, dfjm1);
            GetInverse(dfjm1, mat_dfj);
          }
        else
          {
            dfjm1 = var_problem.Glob_DFjm1(iquad)(num_point);
            GetInverse(dfjm1, mat_dfj);
            jacob = Det(mat_dfj);
          }
        
        // matrix for 1/mu curl => 1/Ji mu^{-1} DF_i
        Mlt(invMu, mat_dfj, MatDfCurl(k));
        Mlt(1.0/jacob, MatDfCurl(k));
        
        PoidsPenal(k) = Poids(k)*var_problem.alpha_penalization;
        if (var_problem.automatic_choice_penalization)
          {
            // interior penalty formulation, the penalty coefficient is well known
            PoidsPenal(k) *= var_problem.Glob_CoefPenalDG(num_face);
          }
        
        // matrix for n x phi => 1/ds DF_i ProjTimesN
        normale = Fb.NormaleLoc(num_pos1_face);
	
        // ProjTimes u = n x u
        ProjTimesN(0, 1) = -normale(2); ProjTimesN(0, 2) = normale(1);
        ProjTimesN(1, 0) = normale(2); ProjTimesN(1, 2) = -normale(0);
        ProjTimesN(2, 0) = -normale(1); ProjTimesN(2, 1) = normale(0);
        
        Mlt(mat_dfj, ProjTimesN, MatDfTangentIntra(k));
        Mlt(1.0/dsj, MatDfTangentIntra(k));
        
        // same matrix for external part
        if (!face_on_gamma)
          {
            if (neighbor_face)
              {
                dfjm1 = this->Glob_DFjm1_Neighbor(offset_neighbor+k);
                GetInverse(dfjm1, mat_dfj);
                jacob = Det(mat_dfj);
              }
            else
              {
                if (affine_elt2)
                  {
                    jacob = var_problem.Glob_jacobian(num_elem2)(0);
                    dfjm1 = var_problem.Glob_DFjm1(num_elem2)(0);
                    Mlt(1.0/jacob, dfjm1);
                    GetInverse(dfjm1, mat_dfj);
                  }
                else
                  {
                    dfjm1 = var_problem.Glob_DFjm1(num_elem2)(num_point2 + krot);
                    GetInverse(dfjm1, mat_dfj);
                    jacob = Det(mat_dfj);
                  }
              }
            
            normale = Fb2.NormaleLoc(num_pos2_face);
	    
            // ProjTimes u = n x u
            ProjTimesN(0, 1) = -normale(2); ProjTimesN(0, 2) = normale(1);
            ProjTimesN(1, 0) = normale(2); ProjTimesN(1, 2) = -normale(0);
            ProjTimesN(2, 0) = -normale(1); ProjTimesN(2, 1) = normale(0);
            
            Mlt(mat_dfj, ProjTimesN, MatDfTangentExtra(k));
            Mlt(-1.0/dsj, MatDfTangentExtra(k));
          }
                
        NumRotQuad(k) = krot;
        num_point++;
      }

    Vector<T> contrib(Fb.GetNbDof()), contrib2(Fb2.GetNbDof());
    contrib.Zero(); contrib2.Zero();
    VectReal_wp val_phi(3*nb_points_face), val_phi_quad(3*Fb.GetNbPointsQuadratureInside());
    Vector<T> feval(3*nb_points_face), feval_quad(3*Fb.GetNbPointsQuadratureInside());
    VectReal_wp Ones(Fb.GetNbDof());
    feval_quad.Zero(); feval.Zero(); val_phi.Zero(); val_phi_quad.Zero();
    
    VectReal_wp curl_phi(3*nb_points_face);
    Vector<T> feval_quad2(3*Fb2.GetNbPointsQuadratureInside());
    feval_quad2.Zero();
    
    // loop over dofs
    for (int i = 0; i < Fb.GetNbDof(); i++)
      {
        Ones.Fill(0); Ones(i) = 1.0;
        if (Fb.UseQuadraturePointsForSh())
          {
            Fb.ApplyChTranspose(Ones, val_phi_quad);
            Fb.ApplyShQuadratureTranspose(num_pos1_face, val_phi_quad, val_phi, rf);
            Fb.ApplyNablaShQuadratureTranspose(num_pos1_face, val_phi_quad, curl_phi, rf);
          }
        else
          {
            Fb.ApplyShTranspose(num_pos1_face, Ones, val_phi, rf);
            Fb.ApplyNablaShTranspose(num_pos1_face, Ones, curl_phi, rf);
          }
       

        if (Norm2(val_phi) > 10.0*epsilon_machine)
          {
            if (!face_on_gamma)
              {
                // internal contribution (n \times u)-  (n \times phi)- 
                for (int k = 0; k < nb_points_face; k++)
                  {
                    CopyVector(val_phi, k, vec_u);
                    Mlt(MatDfTangentIntra(k), vec_u, vec_v);
                    MltTrans(MatDfTangentIntra(k), vec_v, vec_w);
                    feval(3*k) = PoidsPenal(k)*vec_w(0); feval(3*k+1) = PoidsPenal(k)*vec_w(1);
                    feval(3*k+2) = PoidsPenal(k)*vec_w(2);
                  }
                
                if (Fb.UseQuadraturePointsForSh())
                  {
                    feval_quad.Fill(0);
                    Fb.ApplyShQuadrature(cone, num_pos1_face, feval, feval_quad, rf);
                    Fb.ApplyCh(feval_quad, contrib);
                  }
                else
                  {
                    contrib.Fill(0);
                    Fb.ApplySh(cone, num_pos1_face, feval, contrib, rf);
                  }
                
                // adding contributions
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                    {
                      int dofp = offset_row + Nodle(i);
                      int dofq = offset_col + Nodle(j);
                      mat_sp.AddInteraction(dofp, dofq, -contrib(j));
                    }
                
                // external contribution
                for (int k = 0; k < nb_points_face; k++)
                  {
                    CopyVector(val_phi, k, vec_u);
                    Mlt(MatDfTangentIntra(k), vec_u, vec_v);
                    MltTrans(MatDfTangentExtra(k), vec_v, vec_w);
                    feval(3*NumRotQuad(k)) = PoidsPenal(k)*vec_w(0);
                    feval(3*NumRotQuad(k)+1) = PoidsPenal(k)*vec_w(1);
                    feval(3*NumRotQuad(k)+2) = PoidsPenal(k)*vec_w(2);
                  }
                
                if (Fb2.UseQuadraturePointsForSh())
                  {
                    feval_quad2.Fill(0);
                    Fb2.ApplyShQuadrature(cone, num_pos2_face, feval, feval_quad2, rf);
                    Fb2.ApplyCh(feval_quad2, contrib2);
                  }
                else
                  {
                    contrib2.Fill(0);
                    Fb2.ApplySh(cone, num_pos2_face, feval, contrib2, rf);
                  }
                
                // adding contributions
                if (neighbor_face)
                  {                              
#ifdef SELDON_WITH_MPI                          
                    // interactions with another processor
                    for (int j = 0; j < Fb2.GetNbDof(); j++)
                      if (abs(contrib2(j)) > var_problem.GetThresholdMatrix())
                        {
                          int dofp = offset_row + Nodle(i);
                          int dofq = offset_col + Nodle2(j);
                          mat_sp.AddDistantInteraction(dofp, dofq, proc2, phase_conj*contrib2(j));
                        }
#endif
                  }
                else
                  {
                    for (int j = 0; j < Fb2.GetNbDof(); j++)
                      if (abs(contrib2(j)) > var_problem.GetThresholdMatrix())
                        {
                          int dofp = offset_row + Nodle(i);
                          int dofq = offset_col + Nodle2(j);
                          mat_sp.AddInteraction(dofp, dofq, phase_conj*contrib2(j));
                        }
                  }
              }
          }
        
        if (Norm2(curl_phi) > 10.0*epsilon_machine)
          {
            if (!face_on_gamma)
              {
                      
                /**********************
                 * Second order terms *
                 **********************/
                
                // internal contribution -(1/mu rot u)-  n \times phi- 
                //   and -(1/ mu rot phi)- n \times u-
                for (int k = 0; k < nb_points_face; k++)
                  {
                    CopyVector(curl_phi, k, vec_u);
                    Mlt(MatDfCurl(k), vec_u, vec_v);
                    MltTrans(MatDfTangentIntra(k), vec_v, vec_w);
                    feval(3*k) = Poids(k)*vec_w(0);
                    feval(3*k+1) = Poids(k)*vec_w(1); feval(3*k+2) = Poids(k)*vec_w(2);
                  }
                
                if (Fb.UseQuadraturePointsForSh())
                  {
                    feval_quad.Fill(0);
                    Fb.ApplyShQuadrature(cone, num_pos1_face, feval, feval_quad, rf);
                    Fb.ApplyCh(feval_quad, contrib);
                  }
                else
                  {
                    contrib.Fill(0);
                    Fb.ApplySh(cone, num_pos1_face, feval, contrib, rf);
                  }
                
                // adding contributions
                for (int j = 0; j < Fb.GetNbDof(); j++)
                  if (abs(contrib(j)) > var_problem.GetThresholdMatrix())
                    {
                      int dofp = Nodle(i);
                      int dofq = Nodle(j);
                      if (var_problem.sipg_formulation)
                        {
                          // SIPG formulation
                          mat_sp.AddInteraction(offset_row + dofp, offset_col + dofq, -contrib(j));
                          mat_sp.AddInteraction(offset_row + dofq, offset_col + dofp, -contrib(j));
                        }
                      else
                        {
                          // NIPG formulation
                          mat_sp.AddInteraction(offset_row + dofp, offset_col + dofq, contrib(j));
                          mat_sp.AddInteraction(offset_row + dofq, offset_col + dofp, -contrib(j));
                        }
                    }
                          
                // external contribution  (1/mu rot u)-  n \times phi+
                //  and (1/mu rot u)- n \times u+
                for (int k = 0; k < nb_points_face; k++)
                  {
                    CopyVector(curl_phi, k, vec_u);
                    Mlt(MatDfCurl(k), vec_u, vec_v);
                    MltTrans(MatDfTangentExtra(k), vec_v, vec_w);
                    feval(3*NumRotQuad(k)) = Poids(k)*vec_w(0);
                    feval(3*NumRotQuad(k)+1) = Poids(k)*vec_w(1);
                    feval(3*NumRotQuad(k)+2) = Poids(k)*vec_w(2);
                  }
                
                if (Fb2.UseQuadraturePointsForSh())
                  {
                    feval_quad2.Fill(0);
                    Fb2.ApplyShQuadrature(cone, num_pos2_face, feval, feval_quad2, rf);
                    Fb2.ApplyCh(feval_quad2, contrib2);
                  }
                else
                  {
                    contrib2.Fill(0);
                    Fb2.ApplySh(cone, num_pos2_face, feval, contrib2, rf);
                  }
                
                if (neighbor_face)
                  {                              
#ifdef SELDON_WITH_MPI                          
                    // interactions with another processor
                    for (int j = 0; j < Fb2.GetNbDof(); j++)
                      if (abs(contrib2(j)) > var_problem.GetThresholdMatrix())
                        {
                          int dofp = Nodle(i);
                          int dofq = Nodle2(j);
                          if (var_problem.sipg_formulation)
                            {
                              mat_sp.AddDistantInteraction(offset_row + dofp, offset_col + dofq,
                                                           proc2, phase*contrib2(j));
                              mat_sp.AddRowDistantInteraction(offset_row + dofq,
                                                              offset_col + dofp, proc2,
                                                              phase_conj*contrib2(j));
                            }
                          else
                            {
                              mat_sp.AddDistantInteraction(offset_row + dofp, offset_col + dofq,
                                                           proc2, -phase*contrib2(j));
                              mat_sp.AddRowDistantInteraction(offset_row + dofq, offset_col + dofp,
                                                              proc2, phase_conj*contrib2(j));
                            }
                        }
#endif
                  }
                else
                  {                    
                    for (int j = 0; j < Fb2.GetNbDof(); j++)
                      if (abs(contrib2(j)) > var_problem.GetThresholdMatrix())
                        {
                          int dofp = Nodle(i);
                          int dofq = Nodle2(j);
                          if (var_problem.sipg_formulation)
                            {
                              // SIPG formulation
                              mat_sp.AddInteraction(offset_row + dofp, offset_col + dofq,
                                                    phase*contrib2(j));
                              mat_sp.AddInteraction(offset_row + dofq, offset_col + dofp,
                                                    phase_conj*contrib2(j));
                            }
                          else
                            {
                              // NIPG formulation
                              mat_sp.AddInteraction(offset_row + dofp, offset_col + dofq,
                                                    -phase*contrib2(j));
                              mat_sp.AddInteraction(offset_row + dofq, offset_col + dofp,
                                                    phase_conj*contrib2(j));
                            }
                        }
                  }
              }
          }
      }
    
    if (neighbor_face)
      nb_neighbor++;
  
  }

  
  /*************************
   * HarmonicMaxwell_3D_Eq *
   *************************/
  

  //! if mix_ is true, we use a first-order formulation
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::SetFirstOrderFormulation(bool mix_)
  {
    this->mixed_formulation = mix_;
    if (this->mixed_formulation)
      {
        this->mesh_num.compute_dof_pml = true;
        this->compute_dfjm1 = true;
      }
  }

  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::SetTypeEquation(const string& type_equation)
  {
    if ((type_equation == "HARMONIC_MAXWELL_SIPG") || (type_equation == "STATIC_MAXWELL_SIPG")
        || (type_equation == "TIME_MAXWELL_SIPG"))
      {
        this->dg_formulation = ElementReference_Base::DISCONTINUOUS;
	this->compute_dfjm1 = true;
      }

    if (this->dg_formulation != ElementReference_Base::CONTINUOUS)
      {
	this->alpha_penalization = -Real_wp(1);
	this->delta_penalization = -Real_wp(1);
      }
    
    VarHarmonic<TypeEquation>::SetTypeEquation(type_equation);
  }
  

  //! computation of the number of degrees of freedom
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::ComputeNumberOfDofs()
  {
    VarHarmonic<TypeEquation>::ComputeNumberOfDofs();
    HarmonicMaxwell_3D<Complexe>::ComputeNumberOfDofs();
  }


  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::PutOtherGlobalDofs()
  {
    VarHarmonic<TypeEquation>::PutOtherGlobalDofs();
    HarmonicMaxwell_3D<Complexe>::PutOtherGlobalDofs();
  }


  template<class TypeEquation>
  IVect HarmonicMaxwell_3D_Eq<TypeEquation>::GetScalarDofNumberOnElement(int i) const
  {
    return HarmonicMaxwell_3D_Base::GetScalarDofNumberOnElement(i);
  }
  

  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::GetInternalNodesElement(int i, int nb_dof_loc,
                            int& nb_dof_edges, int& nb_dof_int,
                            Vector<int>& intern_node) const
  {
    if (this->light_static_condensation)
      HarmonicMaxwell_3D_Base::GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges,
                                                       nb_dof_int, intern_node);
    else
      VarHarmonic<TypeEquation>::GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges,
                                                         nb_dof_int, intern_node);    
  }


  //! other initializations to perform
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::PerformOtherInitializations()
  { 
    VarHarmonic<TypeEquation>::PerformOtherInitializations();
        
    if (this->FirstOrderFormulation())
      {
        TypeEquation::store_dfjm1 = true;    
        if (this->IsComplexProblem())
          to_complex(Iwp/this->GetOmega(), this->coefficient_volumic_source);
      }

    if (this->linearize_drude)
      {
        this->compute_dfjm1 = true;
        TypeEquation::store_dfjm1 = true;
      }
  }
        

  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::ComputeTauCoefficient()
  {
    HarmonicMaxwell_3D<Complexe>::ComputeTauCoefficient();
  }


  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::TreatDirichletCondition()
  {
    VarHarmonic<TypeEquation>::TreatDirichletCondition();
    HarmonicMaxwell_3D<Complexe>::TreatDirichletCondition();
  }


  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::WriteDatas(Vector<VectReal_wp>& U0)
  {
    VarHarmonic<TypeEquation>::WriteDatas(U0);
  }

  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::WriteDatas(Vector<VectComplex_wp>& U0)
  {
    VarHarmonic<TypeEquation>::WriteDatas(U0);
    HarmonicMaxwell_3D<Complexe>::WriteCoefQ(U0);
  }

  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::GetMemoryUsed(map<string, size_t>& var) const
  {
    VarHarmonic<TypeEquation>::GetMemoryUsed(var);
    HarmonicMaxwell_3D<Complexe>::GetMemoryUsed(var);
  }
  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    HarmonicMaxwell_3D<Complexe>::SetInputData(description_field, parameters);
  }
  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::CopyInputData(const VarProblem_Base& var_base)
  {
    const HarmonicMaxwell_3D_Eq<TypeEquation>& var
      = static_cast<const HarmonicMaxwell_3D_Eq<TypeEquation>& >(var_base);
    
    VarHarmonic<TypeEquation>::CopyInputData(var);
    HarmonicMaxwell_3D<Complexe>::CopyInputData(var);
  }

  template<class TypeEquation>
  bool HarmonicMaxwell_3D_Eq<TypeEquation>
  ::DiagonalScalarMassMatrix(const ElementReference_Dim<Dimension3>&, int i) const
  {
    if (this->InsidePML(i) && this->FirstOrderFormulation())
      return false;
    
    return this->Glob_matMass_elem(i).IsDiagonalMass();
  }

  template<class TypeEquation>
  bool HarmonicMaxwell_3D_Eq<TypeEquation>
  ::BlockDiagonalScalarMassMatrix(const ElementReference_Dim<Dimension3>& Fb, int i) const
  {
    if (this->Glob_matMass_elem(i).IsDiagonalMass())
      return true;

    if (Fb.LumpedMassMatrix())
      return true;
    
    return false;
  }

  template<class TypeEquation>
  bool HarmonicMaxwell_3D_Eq<TypeEquation>
  ::IsSymmetricProblem(bool eigen) const
  {
    bool sym = HarmonicMaxwell_3D<Complexe>::IsSymmetricProblem(eigen);
    if (!sym)
      return false;
    
    return VarHarmonic<TypeEquation>::IsSymmetricProblem(eigen);
  }

  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::InitIndices(int n)
  {
    HarmonicMaxwell_3D<Complexe>::InitIndices(n);
  }
  
  template<class TypeEquation>
  int HarmonicMaxwell_3D_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return HarmonicMaxwell_3D<Complexe>::GetNbPhysicalIndices();
  }
  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    HarmonicMaxwell_3D<Complexe>::SetIndices(i, parameters);
  }
  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    HarmonicMaxwell_3D<Complexe>::SetPhysicalIndex(name_media, i, parameters);
  }
   
  template<class TypeEquation>
  string HarmonicMaxwell_3D_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return HarmonicMaxwell_3D<Complexe>::GetPhysicalIndexName(m);
  }
 
  template<class TypeEquation>
  Real_wp HarmonicMaxwell_3D_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return HarmonicMaxwell_3D<Complexe>::GetVelocityOfMedia(ref);
  }
  
  template<class TypeEquation>
  Real_wp HarmonicMaxwell_3D_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return HarmonicMaxwell_3D<Complexe>::GetVelocityOfInfinity();
  }

  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::FinalizeComputationVaryingIndices()
  {
    HarmonicMaxwell_3D<Complexe>::FinalizeComputationVaryingIndices();
  }
  
  template<class TypeEquation>
  bool HarmonicMaxwell_3D_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return HarmonicMaxwell_3D<Complexe>::IsVaryingMedia(i);
  }

  template<class TypeEquation>
  bool HarmonicMaxwell_3D_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return HarmonicMaxwell_3D<Complexe>::IsVaryingMedia(i);
  }
  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::AllocateMassMatrices()
  {
    HarmonicMaxwell_3D<Complexe>::AllocateMassMatrices();
  }
  
  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>::ComputeMassMatrix(bool compute_rho, bool delete_points)
  {
    bool sipg = false;
    if ((this->FormulationDG() == ElementReference_Base::DISCONTINUOUS) && (!this->FirstOrderFormulationDG()))
      sipg = true;
    
    if (sipg)
        TypeEquation::store_dfjm1 = true;
    
    VarHarmonic<TypeEquation>::ComputeMassMatrix();
    HarmonicMaxwell_3D<Complexe>::ComputeMassMatrix(sipg);
  }


  template<class TypeEquation>
  VirtualSourceField<typename TypeEquation::Complexe, Dimension3>* HarmonicMaxwell_3D_Eq<TypeEquation>
  ::GetNewVolumeSourceFunction(const IVect& ref, const VectString& param,
                               int& nb, Vector<Complexe>& polar,
                               VolumetricSource_Base<Complexe, Dimension3>& var) const
  {
    for (int k = nb; k < param.GetM(); k++)
      if (param(k) == "Polarization")
        {
          for (int p = 0; p < polar.GetM(); p++)
            polar(p) = to_num<Complexe>(param(k+p+1));
          
          break;
        }
    
    int pos = nb;
    if ((param(nb) == "GAUSSIAN") || (param(nb) == "CURL_GAUSSIAN"))
      {
        nb++;
        MaxwellGaussianSource<Complexe>* fgauss = new MaxwellGaussianSource<Complexe>();
        R3 origin; bool presence_curl = false;
        Real_wp RadiusSource(0), RadiusSourceCutOff(0);
        for (int k = 0; k < 3; k++)
          origin(k) = to_num<Real_wp>(param(nb++));
	
        RadiusSource = to_num<Real_wp>(param(nb++));
        if (param.GetM() == nb)
          RadiusSourceCutOff = RadiusSource;
        else
          RadiusSourceCutOff = to_num<Real_wp>(param(nb++));
        
        if (param(pos) == "CURL_GAUSSIAN")
          presence_curl = true;
        else
          presence_curl = false;

        fgauss->Init(origin, RadiusSource, RadiusSourceCutOff, presence_curl, polar);
        return fgauss;
      }
    
    return VarHarmonic<TypeEquation>::GetNewVolumeSourceFunction(ref, param, nb, polar, var);
  }


  template<class TypeEquation>
  VirtualSourceField<typename TypeEquation::Complexe, Dimension3>* HarmonicMaxwell_3D_Eq<TypeEquation>
  ::GetNewSurfaceSourceFunction(const IVect& ref_surf, const VectString& param,
				int& nb, Vector<Complexe>& polar,
				VolumetricSource_Base<Complexe, Dimension3>& var) const
  {
    if (param(nb) == "CIRCULAR_MODE")
      {
	nb++;
	CircularCoaxialModeField<Complexe>* fmode;
	fmode = new CircularCoaxialModeField<Complexe>(ref_surf, *this, *this, param, nb);
	return fmode;
      }

    return VarHarmonic<TypeEquation>::GetNewSurfaceSourceFunction(ref_surf, param, nb, polar, var);
  }


  //! allocation of a modal source
  template<class TypeEquation>
  ModalSourceBoundary_Dim<typename TypeEquation::Complexe, Dimension3>* 
  HarmonicMaxwell_3D_Eq<TypeEquation>::GetNewModalSourceEquation() const
  {
    return new ModalSourceMaxwell3D<Complexe>(this->GetLeafClass());
  }


  template<class TypeEquation>
  void HarmonicMaxwell_3D_Eq<TypeEquation>
  ::FindIntervalDofSignSymmetry(int& i0, int& i1, int& j0, int& j1) const
  {
    HarmonicMaxwell_3D_Base::FindIntervalDofSignSymmetry(i0, i1, j0, j1);
  }

  template<class TypeEquation>
  void HarmonicMaxwell_3D_Cplx<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension3, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension3, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    HarmonicMaxwell_3D<typename TypeEquation::Complexe>::
      GetVaryingIndices(rho_complex, num_ref, num_index,
			num_component, compute_grad, compute_hess);
  }
  

  /*************************************
   * CondensationBlockSolver_Maxwell3D *
   *************************************/

  
  template<class T, class TypeEquation>
  CondensationBlockSolver_Maxwell3D<T, TypeEquation>
  ::CondensationBlockSolver_Maxwell3D(EllipticProblem<TypeEquation>& var)
    : CondensationBlockSolver_Fem<T>(var), var_problem(var)
  {
    
  }

  template<>
  void CondensationBlockSolver_Maxwell3D<Real_wp, HarmonicMaxwellEquation_3D>
  ::Init(const GlobalGenericMatrix<Real_wp>& nat_mat, bool compress_matrix)
  {
    abort();
  }

  
  template<class T, class TypeEquation>
  void CondensationBlockSolver_Maxwell3D<T, TypeEquation>
  ::Init(const GlobalGenericMatrix<T>& nat_mat, bool compress_matrix)
  {
    CondensationBlockSolver_Fem<T>::Init(nat_mat, compress_matrix);
    if (!var_problem.FirstOrderFormulation())
      return;

    if (var_problem.GetNbGlobalEltPML() <= 0)
      return;

    if (var_problem.LightStaticCondensation())
      return;
    
    IVect& MatchingNumber_Pb = var_problem.GetProcMatchingNeighbor();
    Vector<IVect>& MatchingDofOrig_Pb = var_problem.GetOriginalMatchingDofNeighbor();
    
    int ns = MatchingNumber_Pb.GetM();
    IVect nb_dof(ns); nb_dof.Zero();
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    ns = 0;
    for (int i = 0; i < MatchingNumber_Pb.GetM(); i++)
      {
	for (int k = 0; k < MatchingDofOrig_Pb(i).GetM(); k++)
	  if (MatchingDofOrig_Pb(i)(k) >= Nvol)
	    nb_dof(i)++;

	if (nb_dof(i) > 0)
	  ns++;
      }

    IVect SharingProcPML; Vector<IVect> SharingRowsPML;
    SharingProcPML.Reallocate(ns);
    SharingRowsPML.Reallocate(ns);
    ns = 0;
    for (int i = 0; i < MatchingNumber_Pb.GetM(); i++)
      if (nb_dof(i) > 0)
	{
	  SharingRowsPML(ns).Reallocate(nb_dof(i));
	  nb_dof(i) = 0;
	  for (int k = 0; k < MatchingDofOrig_Pb(i).GetM(); k++)
	    if (MatchingDofOrig_Pb(i)(k) >= Nvol)
	      {
		SharingRowsPML(ns)(nb_dof(i)) = MatchingDofOrig_Pb(i)(k) - Nvol;
		nb_dof(i)++;
	      }

	  SharingProcPML(ns) = MatchingNumber_Pb(i);
	  ns++;
	}

    // retrieving tau on pml dofs
    num_u.Reallocate(mesh_num.GetNbDofPML()); num_u.Fill(-1);
    tau.Reallocate(mesh_num.GetNbDofPML());
    int offset = var_problem.mesh.GetNbElt() - var_problem.GetNbEltPML();
    coor_u.Reallocate(mesh_num.GetNbDofPML());
    coor_u.Fill(3);
    for (int i = offset; i < var_problem.mesh.GetNbElt(); i++)
      {
        const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(i);
        HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass = var_problem.Glob_matMass_elem(i);
        for (int k = 0; k < mesh_num.Element(i).GetNbDof(); k++)
          {
            int num_dof = mesh_num.Element(i).GetNumberDof(k);
            int npml = mesh_num.GetDofPML(num_dof);
            if (npml >= 0)
              {
                tau(npml) = mass.GetTauPML(k);
                coor_u(npml) = Fb.GetCoordinateDof(k);
              }
          }
      }
    
    for (int i = 0; i < mesh_num.GetNbDof(); i++)
      {
	int npml = mesh_num.GetDofPML(i);
	if (npml >= 0)
	  num_u(npml) = i;
      }

    //num_u.WriteText("num_pml.dat");
    
    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    ratio_mass.Reallocate(mesh_num.GetNbDofPML());
    mass_pml.Reallocate(mesh_num.GetNbDofPML());
    ratio_mass.Zero(); mass_pml.Zero();
    
    for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
      {
        int ref = var_problem.mesh.Element(i).GetReference();
        Complexe epsilon = var_problem.ref_epsilon(ref).GetConstant()(0, 0);
        const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(i);
        HarmonicMaxwell3D_PhysGeomInfo<Complexe>& mass = var_problem.Glob_matMass_elem(i);
        if (var_problem.InsidePML(i))
          {
            if (!mass.IsDiagonalMass())
              {
                cout << "Static condensation working only for diagonal mass with PML" << endl;
                abort();
              }
            
            for (int k = 0; k < mesh_num.Element(i).GetNbDof(); k++)
              {
                int num_dof = mesh_num.Element(i).GetNumberDof(k);
                int npml = mesh_num.GetDofPML(num_dof);
                if (npml >= 0)
                  mass_pml(npml) += mass.GetDh(k) / (m_iomega*epsilon);
              }
          }
        
        if (mass.IsDiagonalMass())
          {
            for (int k = 0; k < mesh_num.Element(i).GetNbDof(); k++)
              {
                int num_dof = mesh_num.Element(i).GetNumberDof(k);
                int npml = mesh_num.GetDofPML(num_dof);
                if (npml >= 0)
                  ratio_mass(npml) += mass.GetDh(k) / m_iomega;
              }
          }
        else
          {
            int N = Fb.GetNbDof();
            /* int Nquad = Fb.GetNbPointsQuadratureInside();
               Vector<TinyMatrix<Complexe, General, 3, 3> > Bmass(Nquad);
                   for (int k = 0; k < Nquad; k++)
                   Bmass(k) = mass.GetBh(k);
                   
                   Matrix<Complexe, Symmetric, RowSymPacked> mat_mass;
                   mat_mass.Reallocate(N, N); mat_mass.Zero();
                   Fb.AddVariableMassMatrix(0, 0, Bmass, mat_mass);
            */
            
            for (int k = 0; k < N; k++)
              {
                int num_dof = mesh_num.Element(i).GetNumberDof(k);
                int npml = mesh_num.GetDofPML(num_dof);
                if (npml >= 0)
                  {
                    cout << "non diagonal mass matrix not handled" << endl;
                    abort();
                    /*
                      for (int j = 0; j < N; j++)
                      if (j != k)
                      if (abs(mat_mass(j, k)) > 1e-12)
                      {
                      DISP(j); DISP(k); DISP(mat_mass(j, k));
                      cout << "non diagonal mass matrix not handled" << endl;
                      abort();
                      }
                      
                      ratio_mass(npml) += mat_mass(k, k) / m_iomega; */
                  }
              }
          }
      }

    // assembling between processors
#ifdef SELDON_WITH_MPI
    AssembleVector(mass_pml, MPI_SUM, SharingProcPML, SharingRowsPML,
                   var_problem.comm_group_mode, mesh_num.GetNbDofPML(), 1, 21);
    
    AssembleVector(ratio_mass, MPI_SUM, SharingProcPML, SharingRowsPML,
                   var_problem.comm_group_mode, mesh_num.GetNbDofPML(), 1, 22);
#endif
    
    for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
      ratio_mass(i) /= mass_pml(i);
    
#ifdef SELDON_WITH_MPI
    if (var_problem.GetNbProcPerMode() > 1)
      {
        VectReal_wp tau_scal(tau.GetM());
        for (int k = 0; k < 3; k++)
          {
            for (int i = 0; i < tau.GetM(); i++)
              tau_scal(i) = tau(i)(k);
    
            AssembleVector(tau_scal, MPI_MAX, SharingProcPML, SharingRowsPML,
                           var_problem.comm_group_mode, mesh_num.GetNbDofPML(),
                           1, 21);
            
            for (int i = 0; i < tau.GetM(); i++)
              tau(i)(k) = tau_scal(i);
          }

	Vector<int> coor_scal(coor_u.GetM());
	for (int i = 0; i < coor_u.GetM(); i++)
	  coor_scal(i) = coor_u(i);
	
	AssembleVector(coor_scal, MPI_MIN, SharingProcPML, SharingRowsPML,
		       var_problem.comm_group_mode, mesh_num.GetNbDofPML(),
		       1, 22);

	for (int i = 0; i < coor_u.GetM(); i++)
	  coor_u(i) = coor_scal(i);

      }
#endif

  }
  

  template<>
  void CondensationBlockSolver_Maxwell3D<Real_wp, HarmonicMaxwellEquation_3D>
  ::RecomposeSolution(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }


  template<class T, class TypeEquation>
  void CondensationBlockSolver_Maxwell3D<T, TypeEquation>
  ::RecomposeSolution(const SeldonTranspose& trans, Vector<T>& x_sol, const GlobalGenericMatrix<T>& nat_mat) const
  {
    IVect num_ddl;
    Matrix<T> mat_elt_unsym;

    int nb_elt = var_problem.GetNbElt();
    Matrix<T, Symmetric, RowSymPacked> a22_sym;
    Matrix<T> a12, a21, a22_unsym;

    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    bool sym = this->UseSymmetrization();

    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    // on forme E* (for PML layers)
    Vector<T> Estar;
    if ((var_problem.GetNbGlobalEltPML() > 0) && var_problem.LightStaticCondensation()
        && var_problem.FirstOrderFormulation())
      {
        int Nvol = mesh_num.GetNbDof();
        Estar.Reallocate(Nvol);
        for (int i = 0; i < Nvol; i++)
          {
            Estar(i) = x_sol(i);
            int npml = mesh_num.GetDofPML(i);
            if (npml >= 0)
              Estar(i) = x_sol(Nvol + npml);
          }
      }
    
    var_problem.SetLeafStaticCondensation(true);
    T sig = nat_mat.GetCoefDamping();
    T m = nat_mat.GetCoefMass()*m_iomega;
    T s = nat_mat.GetCoefStiffness();
    T one; SetComplexOne(one);
    
    for (int i = 0; i < nb_elt; i++)
      {	
        const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(i);
        int nb_dof_loc = Fb.GetNbDof();
        
	if (!this->store_block)
	  var_problem.ComputeElementaryMatrix(i, num_ddl, mat_elt_unsym,
					      const_cast<CondensationBlockSolver_Maxwell3D<T, TypeEquation>& >(*this),
					      nat_mat);

	// solution on internal nodes is recovered
	if (this->symmetric_elem_matrix)
	  {
	    if (this->store_block)
	      this->RecomposeSolGen(trans, x_sol, i, this->block_sym(i).num_ddl, mat_elt_unsym,
				     this->block_sym(i).a12, this->block_sym(i).a21, this->block_sym(i).inv_a22);
	    else
	      this->RecomposeSolGen(trans, x_sol, i, num_ddl, mat_elt_unsym, a12, a21, a22_sym);
	  }
	else
	  {
	    if (this->store_block)
	      this->RecomposeSolGen(trans, x_sol, i, this->block_unsym(i).num_ddl, mat_elt_unsym,
				    this->block_unsym(i).a12, this->block_unsym(i).a21, this->block_unsym(i).inv_a22);
	    else
	      this->RecomposeSolGen(trans, x_sol, i, num_ddl, mat_elt_unsym, a12, a21, a22_unsym);
	  }

	// solution on vectorial dofs is recovered (H and H*)
        if (var_problem.FirstOrderFormulation())
	  {
	    bool pml_elt = var_problem.InsidePML(i);
	    
	    // we extract E*
	    TinyVector<Vector<T>, 1> x;
	    x(0).Reallocate(nb_dof_loc); x(0).Zero();
            if (var_problem.LightStaticCondensation())
              var_problem.GetLocalUnknownVector(Estar, i, x);
            else
              var_problem.GetLocalUnknownVector(x_sol, i, x);
            
            if (trans.Trans())
              for (int j = 0; j < nb_dof_loc; j++)
                {
                  int num_dof = mesh_num.Element(i).GetNumberDof(j);
                  if (num_dof >= 0)
                    if (var_problem.IsDofDirichlet(num_dof))
                      SetComplexZero(x(0)(j));
                }
            
	    int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
	    Vector<T> y(3*nb_pts_quad);

	    // curl of E*
	    Fb.ApplyRhTranspose(x(0), y);
	    
	    bool affine = var_problem.mesh.IsElementAffine(i); 
	    int ref_domain = var_problem.mesh.Element(i).GetReference();
            
            bool mu_drude = false; T coef_drude; SetComplexOne(coef_drude);
            if (var_problem.ref_drude(ref_domain).IsEnabled())
              if (var_problem.ref_drude(ref_domain).IsModeTE())
                {
                  mu_drude = true;
                  int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                  Real_wp gamma(0), omega0_2(0); T invDenom(0);
                  SetComplexZero(coef_drude);
                  for (int kp = 0; kp < nPole; kp++)
                    {
                      gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                      omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                      Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                      Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp);
                      invDenom = one / (m*(m+gamma*s) + omega0_2*s*s);
                      coef_drude += (s*s*coef_eps_inf + s*m*coef_sig) * invDenom;
                    }
                  
                  coef_drude = one / ( var_problem.ref_drude(ref_domain).eps_inf
                                       + coef_drude );
                }
            
	    TinyVector<T, 3> vec_u, vec_v, vec_w;
	    int offset = var_problem.offset_dof_unknown(1) + var_problem.GetOffsetDofV(i);
	    int offset_Hs = offset + 3*nb_pts_quad;
	    T coef = nat_mat.GetCoefStiffness()/(m_iomega*nat_mat.GetCoefMass());
	    T coef_mass = Real_wp(1) / (m_iomega*nat_mat.GetCoefMass());
	    T coef1, coef2;
	    int i1 = i - var_problem.mesh.GetNbElt() + var_problem.GetNbEltPML();	    
	    for (int k = 0; k < nb_pts_quad; k++)
	      {
		CopyVector(y, k, vec_u);
		if (affine)
		  {
		    Mlt(var_problem.Glob_DFj(i)(0), vec_u, vec_v);
		    coef1 = coef / var_problem.Glob_jacobian(i)(0);
		    coef2 = coef_mass / (var_problem.Glob_jacobian(i)(0)*Fb.WeightsND(k));
		  }
		else
		  {
		    Mlt(var_problem.Glob_DFj(i)(k), vec_u, vec_v);
		    coef1 = coef*Fb.WeightsND(k) / var_problem.Glob_jacobian(i)(k);
		    coef2 = coef_mass / var_problem.Glob_jacobian(i)(k);
		  }

		if (sym)
		  coef2 = -coef2;
                else if (trans.Trans())
                  coef1 = -coef1;
                
                if (trans.Trans() && pml_elt)
                  ExtractVector(x_sol, offset_Hs + 3*k, vec_w);
                else
                  ExtractVector(x_sol, offset + 3*k, vec_w);
		
		vec_u = -coef1*vec_v + coef2*vec_w;
                
                if (trans.Trans() && pml_elt)
                  vec_v = vec_u;
                else
                  {
                    if (mu_drude)
                      vec_v = coef_drude*vec_u;
                    else
                      var_problem.ref_invMu(ref_domain).MltMatrix(var_problem, i, k, vec_u, vec_v);
                  }
                
		// computing H* (and modifying H for PML layers)
		if (pml_elt)
		  {
		    Real_wp tau0 = realpart(var_problem.GetTauPML(i1, k, 0));
		    Real_wp tau1 = realpart(var_problem.GetTauPML(i1, k, 1));
		    Real_wp tau2 = realpart(var_problem.GetTauPML(i1, k, 2));
		    Real_wp poids = var_problem.Glob_jacobian(i)(0)*Fb.WeightsND(k);
                    if (!affine)
                      poids = var_problem.Glob_jacobian(i)(k);
                    
                    if (trans.Trans())
                      {
                        vec_v(0) *= m / (m + sig*tau2);
                        vec_v(1) *= m / (m + sig*tau0);
                        vec_v(2) *= m / (m + sig*tau1);

                        ExtractVector(vec_v, offset_Hs + 3*k, x_sol);

                        ExtractVector(x_sol, offset + 3*k, vec_u);
                        
                        vec_u(0) = (vec_u(0)/poids + (m + sig*tau0)*vec_v(0)) / (m + sig*tau1);
                        vec_u(1) = (vec_u(1)/poids + (m + sig*tau1)*vec_v(1)) / (m + sig*tau2);
                        vec_u(2) = (vec_u(2)/poids + (m + sig*tau2)*vec_v(2)) / (m + sig*tau0);
                        
                        var_problem.ref_invMu(ref_domain).MltMatrix(var_problem, i, k, vec_u, vec_v);
                      }
                    else
                      {
                        vec_v(0) *= m / (m + sig*tau1);
                        vec_v(1) *= m / (m + sig*tau2);
                        vec_v(2) *= m / (m + sig*tau0);
                        
                        ExtractVector(x_sol, offset_Hs + 3*k, vec_u);
                        
                        vec_u(0) /= (m + sig*tau2) * poids;
                        vec_u(1) /= (m + sig*tau0) * poids;
                        vec_u(2) /= (m + sig*tau1) * poids;
                        
                        vec_u(0) += (m + sig*tau0) / (m + sig*tau2) * vec_v(0);
                        vec_u(1) += (m + sig*tau1) / (m + sig*tau0) * vec_v(1);
                        vec_u(2) += (m + sig*tau2) / (m + sig*tau1) * vec_v(2);
                        
                        ExtractVector(vec_u, offset_Hs + 3*k, x_sol);
                      }
		  }

		ExtractVector(vec_v, offset + 3*k, x_sol);
	      }

            // P and Q are recovered (Drude's model)
            if (var_problem.ref_drude(ref_domain).IsEnabled())
	      {
                int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                Vector<T> coef_E(nPole), coef_p(nPole);
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                    Real_wp omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                    T invDenom = Real_wp(1) / (m*(m+gamma*s) + omega0_2*s*s);
                    Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                    Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp);
                    
                    T coefE = s*m*invDenom;                
                    T coefp = s/m;
                    if (var_problem.use_symm_drude)
                      {
                        coefE *= coef_eps_inf;
                        if (trans.Trans() && (!sym))
                          {
                            coefE = -coefE;                
                            coefp = -coefp;
                          }
                      }
                    else
                      {
                        if (trans.NoTrans())
                          {
                            coefE *= coef_eps_inf;
                            coefE += m*m*invDenom*coef_sig;
                          }
                        else
                          {
                            if (!sym)
                              {
                                coefE = -coefE;                          
                                coefp *= -omega0_2;
                              }
                            else
                              coefp *= omega0_2;
                          }
                      }

                    coef_E(kp) = coefE; coef_p(kp) = coefp;
                  }
                
                TinyVector<T, 3> vec_P, vec_Q;
                
                int offset_Q = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                int offset_P = offset_Q - var_problem.GetNbVectorialDofDrude();
                offset_P += var_problem.OffsetDofDrudeV(i);
                offset_Q += var_problem.OffsetDofDrudeV(i);

                if (var_problem.ref_drude(ref_domain).IsModeTE())
                  {
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        for (int k = 0; k < nb_pts_quad; k++)
                          {
                            ExtractVector(x_sol, offset + 3*k, vec_v);
                            ExtractVector(x_sol, offset_P + 3*k, vec_P);
                            ExtractVector(x_sol, offset_Q + 3*k, vec_Q);
                            
                            Add(coef_E(kp), vec_v, vec_Q);
                            Add(coef_p(kp), vec_Q, vec_P);
                            
                            ExtractVector(vec_P, offset_P+3*k, x_sol);
                            ExtractVector(vec_Q, offset_Q+3*k, x_sol);                            
                          }

                        offset_P += 3*nb_pts_quad;
                        offset_Q += 3*nb_pts_quad;
                      }
                  }
                else
                  {
                    Fb.ApplyChTranspose(x(0), y);
                    
                    for (int k = 0; k < nb_pts_quad; k++)
                      {
                        CopyVector(y, k, vec_u);
                        if (affine)
                          {
                            MltTrans(var_problem.Glob_DFjm1(i)(0), vec_u, vec_v);
                            vec_v *= Real_wp(1) / var_problem.Glob_jacobian(i)(0);
                          }
                        else
                          {
                            Real_wp jacob = var_problem.Glob_jacobian(i)(k)/Fb.WeightsND(k);
                            MltTrans(var_problem.Glob_DFjm1(i)(k), vec_u, vec_v);
                            vec_v *= Real_wp(1) / jacob;
                          }

                        for (int kp = 0; kp < nPole; kp++)
                          {
                            ExtractVector(x_sol, offset_P + 3*(kp*nb_pts_quad+k), vec_P);
                            ExtractVector(x_sol, offset_Q + 3*(kp*nb_pts_quad+k), vec_Q);
                            
                            Add(coef_E(kp), vec_v, vec_Q);
                            Add(coef_p(kp), vec_Q, vec_P);
                            
                            ExtractVector(vec_P, offset_P+3*(kp*nb_pts_quad+k), x_sol);
                            ExtractVector(vec_Q, offset_Q+3*(kp*nb_pts_quad+k), x_sol);
                          }
                      }
                  }
              }
          }
        else
          {
            int ref_domain = var_problem.mesh.Element(i).GetReference();
            if (var_problem.ref_drude(ref_domain).IsEnabled() && var_problem.ref_drude(ref_domain).IsModeTM() && var_problem.linearize_drude)
              {
                TinyVector<Vector<T>, 1> x;
                x(0).Reallocate(nb_dof_loc); x(0).Zero();                
                var_problem.GetLocalUnknownVector(x_sol, i, x);
                int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                Vector<T> coef_E(nPole), coefDiag(nPole);
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                    Real_wp omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                    Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                    Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp); 
                    if (trans.NoTrans())
                      coef_E(kp) = coef_eps_inf*s + m_iomega*sig*coef_sig;
                    else
                      coef_E(kp) = -m_iomega*m;

                    coefDiag(kp) = m_iomega*m + m_iomega*sig*gamma + s*omega0_2;
                    coef_E(kp) /= coefDiag(kp);
                  }
                
                int offset_P = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                offset_P += var_problem.OffsetDofDrudeV(i);

                bool affine = var_problem.mesh.IsElementAffine(i); 
                
                int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
                Vector<T> y(3*nb_pts_quad);
                Fb.ApplyChTranspose(x(0), y);
                
                TinyVector<T, 3> vec_u, vec_v, vec_P;
                for (int k = 0; k < nb_pts_quad; k++)
                  {
                    CopyVector(y, k, vec_u);
                    if (affine)
                      {
                        MltTrans(var_problem.Glob_DFjm1(i)(0), vec_u, vec_v);
                        vec_v *= Real_wp(1) / var_problem.Glob_jacobian(i)(0);
                      }
                    else
                      {
                        Real_wp jacob = var_problem.Glob_jacobian(i)(k)/Fb.WeightsND(k);
                        MltTrans(var_problem.Glob_DFjm1(i)(k), vec_u, vec_v);
                        vec_v *= Real_wp(1) / jacob;
                      }
                    
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        ExtractVector(x_sol, offset_P + 3*(kp*nb_pts_quad+k), vec_P);
                        Add(coef_E(kp), vec_v, vec_P);
                        ExtractVector(vec_P, offset_P + 3*(kp*nb_pts_quad+k), x_sol);
                      }
                  }
              }                          
          }
      }

    var_problem.SetLeafStaticCondensation(false);

    if (!var_problem.FirstOrderFormulation())
      return;

    if (var_problem.LightStaticCondensation())
      return;
    
    // recovering E and E* (for PML)
    int Nvol = mesh_num.GetNbDof();
    TinyVector<int, 3> permut201(2, 0, 1);
    TinyVector<int, 3> permut120(1, 2, 0);
    for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
      {
        int iglob = num_u(i);
        if (trans.Trans())
          {
            T E = x_sol(iglob);
            T fE = x_sol(Nvol + i);
            T Estar = ((m + sig*tau(i)(permut120(coor_u(i)))) * ratio_mass(i) * E - fE / mass_pml(i))
              / (m + sig*tau(i)(coor_u(i)));
            
            x_sol(Nvol +i) = Estar;
          }
        else
          {
            T Estar = x_sol(iglob);
            T fE_star = x_sol(Nvol + i);
            T E = (m + sig*tau(i)(permut201(coor_u(i)))) / (m + sig*tau(i)(coor_u(i))) * Estar;
            E -= fE_star / (mass_pml(i) * (m + sig*tau(i)(coor_u(i)))); 
            x_sol(iglob) = E;
            x_sol(Nvol + i) = Estar;
          }
      }
  }


  template<>
  void CondensationBlockSolver_Maxwell3D<Real_wp, HarmonicMaxwellEquation_3D>
  ::ModifyRhsStaticCondensation(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }


  template<class T, class TypeEquation>
  void CondensationBlockSolver_Maxwell3D<T, TypeEquation>
  ::ModifyRhsStaticCondensation(const SeldonTranspose& trans, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const
  {
    var_problem.SetLeafStaticCondensation(true);

    Vector<T> b_dir;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    if (trans.NoTrans())
      {
        b_dir.Reallocate(var_problem.GetNbDirichletDof());
        if ((var_problem.GetNbGlobalEltPML() <= 0) || !var_problem.FirstOrderFormulation()
            || var_problem.LightStaticCondensation())
          {
            for (int i = 0; i < b_dir.GetM(); i++)
              {
                int num_dof = var_problem.GetDirichletDofNumber(i);	  
                b_dir(i) = b_source(num_dof);
              }
          }
        else
          for (int i = 0; i < b_dir.GetM(); i++)
            {
              int num_dof = var_problem.GetDirichletDofNumber(i);
              if (num_dof < mesh_num.GetNbDof())
                {
                  int npml = mesh_num.GetDofPML(num_dof);
                  if (npml < 0)
                    b_dir(i) = b_source(num_dof);
                  else
                    {
                      T bstar = b_source(Nvol+npml);
                      b_source(Nvol+npml) = b_source(num_dof);
                      b_source(num_dof) = bstar;
                      b_dir(i) = bstar;
                    }
                }
            }
      }
    
    IVect num_ddl;
    Matrix<T> mat_elt_unsym;
    Matrix<T, Symmetric, RowSymPacked> a22_sym;
    Matrix<T> a12, a21, a22_unsym;

    //b_source.Write("source_before.dat");
    
    Complexe m_iomega; var_problem.GetMiomega(m_iomega);
    T sig = nat_mat.GetCoefDamping();
    T m = nat_mat.GetCoefMass()*m_iomega;
    T s = nat_mat.GetCoefStiffness();

    // part due to E*
    TinyVector<int, 3> permut120(1, 2, 0);
    TinyVector<int, 3> permut201(2, 0, 1);
    if (var_problem.FirstOrderFormulation() && !var_problem.LightStaticCondensation())
      for (int i = 0; i < mesh_num.GetNbDofPML(); i++)
	{
	  int iglob = num_u(i);
          if (trans.Trans())
            {
              T fE = b_source(iglob);
              b_source(iglob) = b_source(Nvol+i) + fE * (m + sig*tau(i)(permut201(coor_u(i)))) / (m + sig*tau(i)(coor_u(i)));
              b_source(Nvol+i) = fE;
            }
          else
            {
              T f_Estar = b_source(Nvol+i);
              b_source(iglob) += ratio_mass(i) * f_Estar * (m + sig*tau(i)(permut120(coor_u(i)))) / (m + sig*tau(i)(coor_u(i)));
            }
	}    
    
    bool sym = this->UseSymmetrization();
    int nb_elt = var_problem.mesh.GetNbElt();
    Vector<T> xdir;
    T one; SetComplexOne(one);
    for (int i = 0; i < nb_elt; i++)
      {		
        const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(i);
        int nb_dof_loc = Fb.GetNbDof();
        
	// part due to elimination of H
	if (var_problem.FirstOrderFormulation())
	  {	    	    
	    int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
	    Vector<T> y(3*nb_pts_quad);
	    
	    TinyVector<T, 3> vec_u, vec_v, vec_w, vec_H;
	    bool affine = var_problem.mesh.IsElementAffine(i);
	    bool pml_elt = var_problem.InsidePML(i);
	    int ref_domain = var_problem.mesh.Element(i).GetReference();
	    int offset = var_problem.offset_dof_unknown(1) + var_problem.GetOffsetDofV(i);
	    T coef = nat_mat.GetCoefStiffness() / m;
	    if (sym)
	      coef = -coef;
	    else if (trans.Trans())
	      coef = -coef;
            
	    int i1 = i - var_problem.mesh.GetNbElt() + var_problem.GetNbEltPML();
            bool mu_drude = false; T coef_drude; SetComplexOne(coef_drude);

            // part due to elimination of P and Q (Drude's model)
            if (var_problem.ref_drude(ref_domain).IsEnabled())
              if (var_problem.ref_drude(ref_domain).IsModeTE())
                {
                  mu_drude = true;
                  Real_wp gamma(0), omega0_2(0); T invDenom(0);
                  int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                  Vector<T> coef_p(nPole), coef_q(nPole);
                  VectReal_wp coef_p2(nPole), coef_q2(nPole);
                  SetComplexZero(coef_drude);
                  for (int kp = 0; kp < nPole; kp++)
                    {
                      gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                      omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                      Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                      Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp);
                      invDenom = one / (m*(m+gamma*s) + omega0_2*s*s);
                      coef_drude += (s*s*coef_eps_inf + s*m*coef_sig) * invDenom;

                      coef_q2(kp) = Real_wp(1);
                      coef_p2(kp) = Real_wp(1);
                      if (var_problem.use_symm_drude)
                        {
                          coef_p(kp) = s*s*coef_eps_inf * invDenom;
                          coef_q(kp) = -m*s * invDenom * coef_eps_inf;
                          coef_q2(kp) = -Real_wp(1);
                          coef_p2(kp) = coef_eps_inf / omega0_2;

                          if (trans.Trans() || sym)
                            {
                              coef_q(kp) = -coef_q(kp);
                              coef_q2(kp) = Real_wp(1);
                            }
                          
                          if (sym)
                            coef_p2(kp) = -coef_p2(kp);
                        }
                      else
                        {
                          if (trans.Trans())
                            {
                              coef_p(kp) = s*s * coef_eps_inf * invDenom;
                              coef_q(kp) = m*s* coef_eps_inf * invDenom;
                              coef_q2(kp) = Real_wp(1) / coef_eps_inf;
                              if (sym)
                                coef_p2(kp) = Real_wp(-1);
                            }
                          else
                            {
                              coef_p(kp) = omega0_2 * s*s * invDenom;
                              coef_q(kp) = -m*s * invDenom;
                              coef_q2(kp) = Real_wp(-1);
                              if (sym)
                                {
                                  coef_q(kp) = -coef_q(kp);
                                  coef_q2(kp) = Real_wp(1);
                                  coef_p2(kp) = Real_wp(-1);
                                }
                            }
                        }                      
                    }
                  
                  coef_drude = one / ( var_problem.ref_drude(ref_domain).eps_inf
                                       + coef_drude );                  
                  
                  TinyVector<T, 3> vec_P, vec_Q; Real_wp poids(1);
                  int offset_Q = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                  int offset_P = offset_Q - var_problem.GetNbVectorialDofDrude();
                  offset_P += var_problem.OffsetDofDrudeV(i);
                  offset_Q += var_problem.OffsetDofDrudeV(i);
                  for (int k = 0; k < nb_pts_quad; k++)
                    {
                      if (affine)
                        poids = var_problem.Glob_jacobian(i)(0)*Fb.WeightsND(k);
                      else
                        poids = var_problem.Glob_jacobian(i)(k);

                      ExtractVector(b_source, offset + 3*k, vec_H);
                      vec_v = vec_H;
                      for (int kp = 0; kp < nPole; kp++)
                        {
                          ExtractVector(b_source, offset_P + 3*(kp*nb_pts_quad+k), vec_P);
                          ExtractVector(b_source, offset_Q + 3*(kp*nb_pts_quad+k), vec_Q);
                          
                          vec_u = coef_p(kp)*vec_P + coef_q(kp)*vec_Q;
                          vec_v += vec_u;                          
                          
                          vec_u *= coef_q2(kp) / (s*poids);
                          vec_P *= coef_p2(kp) / (m*poids);
                      
                          ExtractVector(vec_P, offset_P+3*(kp*nb_pts_quad+k), b_source);
                          ExtractVector(vec_u, offset_Q+3*(kp*nb_pts_quad+k), b_source);
                        }

                      ExtractVector(vec_v, offset + 3*k, b_source);
                    }
                }
            
	    for (int k = 0; k < nb_pts_quad; k++)
	      {
		ExtractVector(b_source, offset + 3*k, vec_u);

		if (mu_drude)
                  vec_v = coef_drude*vec_u;
                else
                  var_problem.ref_invMu(ref_domain).MltMatrix(var_problem, i, k, vec_u, vec_v);
                
                vec_u = vec_v;
                
		if (pml_elt)
		  {
		    // we extract fH*
		    ExtractVector(b_source, offset + 3*nb_pts_quad + 3*k, vec_w);
                    
		    Real_wp tau0 = realpart(var_problem.GetTauPML(i1, k, 0));
		    Real_wp tau1 = realpart(var_problem.GetTauPML(i1, k, 1));
		    Real_wp tau2 = realpart(var_problem.GetTauPML(i1, k, 2));

		    T dx = m + sig*tau0, dy = m + sig*tau1, dz = m + sig*tau2;

                    if (trans.Trans())
                      {
                        vec_u(0) *= m/dy;
                        vec_u(1) *= m/dz;
                        vec_u(2) *= m/dx;
                        
                        var_problem.ref_invMu(ref_domain).MltMatrix(var_problem, i, k, vec_w, vec_v);
                        vec_u(0) += m*vec_v(0) * dx / (dy*dz);
                        vec_u(1) += m*vec_v(1) * dy / (dx*dz);
                        vec_u(2) += m*vec_v(2) * dz / (dx*dy);
                      }
                    else
                      {
                        vec_u(0) *= (m*dx) / (dy*dz);
                        vec_u(1) *= (m*dy) / (dx*dz);
                        vec_u(2) *= (m*dz) / (dx*dy);
                        
                        vec_u(0) += m*vec_w(0) / dz;
                        vec_u(1) += m*vec_w(1) / dx;
                        vec_u(2) += m*vec_w(2) / dy;
                      }
		  }

		if (affine)
		  {
		    MltTrans(var_problem.Glob_DFj(i)(0), vec_u, vec_v);
		    vec_v *= coef / var_problem.Glob_jacobian(i)(0);
		  }
		else
		  {
		    MltTrans(var_problem.Glob_DFj(i)(k), vec_u, vec_v);
		    vec_v *= coef*Fb.WeightsND(k) / var_problem.Glob_jacobian(i)(k);
		  }

		CopyVector(vec_v, k, y);
	      }	    
            
	    TinyVector<Vector<T>, 1> x;
	    x(0).Reallocate(nb_dof_loc); x(0).Zero();
	    Fb.ApplyRh(y, x(0));
            
            // part due to elimination of P and Q (Drude's model)
            coef = nat_mat.GetCoefStiffness() / m;
            if (sym)
	      coef = -coef;
            else if (trans.Trans())
              coef = -coef;
            
            if (var_problem.ref_drude(ref_domain).IsModeTM())
              {
                int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                Vector<T> coef_p(nPole), coef_q(nPole);
                VectReal_wp coef_p2(nPole), coef_q2(nPole);
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp gamma(0), omega0_2(0); T invDenom(0);
                    coef_p(kp) = T(0); coef_q(kp) = T(0);
                    gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                    omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                    Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                    
                    invDenom = one / (m*(m+gamma*s) + omega0_2*s*s);
                    
                    coef_q2(kp) = Real_wp(1); coef_p2(kp) = Real_wp(1);
                    if (var_problem.use_symm_drude)
                      {
                        coef_p(kp) = s*s*coef_eps_inf * invDenom;
                        coef_q(kp) = -m*s * invDenom * coef_eps_inf;
                        coef_q2(kp) = -Real_wp(1);
                        coef_p2(kp) = coef_eps_inf / omega0_2;
                        if (sym)
                          coef_q(kp) = -coef_q(kp);
                        else if (trans.Trans())
                          {
                            coef_q(kp) = -coef_q(kp);                    
                            coef_q2(kp) = Real_wp(1);
                          }
                      }
                    else
                      {
                        if (trans.Trans())
                          {
                            coef_p(kp) = s*s * coef_eps_inf * invDenom;
                            coef_q(kp) = m*s* coef_eps_inf * invDenom;
                            coef_q2(kp) = Real_wp(1) / coef_eps_inf;
                            if (sym)
                              coef_q2(kp) = -coef_q2(kp);
                          }
                        else
                          {
                            coef_p(kp) = omega0_2 * s*s * invDenom;
                            coef_q(kp) = -m*s * invDenom;
                            coef_q2(kp) = Real_wp(-1);
                            if (sym)
                              coef_q(kp) = -coef_q(kp);
                          }
                      }
                  }
                
                TinyVector<T, 3> vec_P, vec_Q; Real_wp poids;
                int offset_Q = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                int offset_P = offset_Q - var_problem.GetNbVectorialDofDrude();
                offset_P += var_problem.OffsetDofDrudeV(i);
                offset_Q += var_problem.OffsetDofDrudeV(i);
                for (int k = 0; k < nb_pts_quad; k++)
                  {
                    vec_H.Zero();
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        ExtractVector(b_source, offset_P + 3*(kp*nb_pts_quad+k), vec_P);
                        ExtractVector(b_source, offset_Q + 3*(kp*nb_pts_quad+k), vec_Q);
                        
                        vec_u = coef_p(kp)*vec_P + coef_q(kp)*vec_Q;
                        if (affine)
                          {
                            poids = var_problem.Glob_jacobian(i)(0)*Fb.WeightsND(k);
                            Mlt(var_problem.Glob_DFjm1(i)(0), vec_u, vec_v);
                            vec_v *= Real_wp(1) / var_problem.Glob_jacobian(i)(0);
                          }
                        else
                          {
                            poids = var_problem.Glob_jacobian(i)(k);
                            Real_wp jacob = var_problem.Glob_jacobian(i)(k)/Fb.WeightsND(k);
                            Mlt(var_problem.Glob_DFjm1(i)(k), vec_u, vec_v);
                            vec_v *= Real_wp(1) / jacob;
                          }
                        
                        vec_H += vec_v;
                        
                        vec_u *= coef_q2(kp) / (s*poids);
                        vec_P *= coef_p2(kp) / (m*poids);
                        
                        ExtractVector(vec_P, offset_P+3*(kp*nb_pts_quad+k), b_source);
                        ExtractVector(vec_u, offset_Q+3*(kp*nb_pts_quad+k), b_source);
                      }

                    CopyVector(vec_H, k, y);
                  }
                
                Vector<T> ChY(Fb.GetNbDof()); 
                Fb.ApplyCh(y, ChY);

                x(0) += ChY;
              }
            
	    var_problem.AddLocalUnknownVector(one, x, i, b_source);
	  }
        else
          {
            int ref_domain = var_problem.mesh.Element(i).GetReference();
            if (var_problem.ref_drude(ref_domain).IsEnabled() && var_problem.ref_drude(ref_domain).IsModeTM() && var_problem.linearize_drude)
              {
                int nPole = var_problem.ref_drude(ref_domain).gamma.GetM();
                Vector<T> coef_E(nPole), coefDiag(nPole);
                for (int kp = 0; kp < nPole; kp++)
                  {
                    Real_wp gamma = var_problem.ref_drude(ref_domain).gamma(kp);
                    Real_wp omega0_2 = var_problem.ref_drude(ref_domain).omega_02(kp);
                    Real_wp coef_eps_inf = var_problem.ref_drude(ref_domain).eps_omega_p2(kp);
                    Real_wp coef_sig = var_problem.ref_drude(ref_domain).eps_sigma(kp); 
                    if (trans.NoTrans())
                      coef_E(kp) = -m_iomega*m;
                    else
                      coef_E(kp) = coef_eps_inf*s + m_iomega*sig*coef_sig;

                    coefDiag(kp) = m_iomega*m + m_iomega*sig*gamma + s*omega0_2;
                    coef_E(kp) /= coefDiag(kp);
                  }

                int offset_P = var_problem.GetNbDof() - var_problem.GetNbVectorialDofDrude();
                offset_P += var_problem.OffsetDofDrudeV(i);

                int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
                Vector<T> y(3*nb_pts_quad);
                
                bool affine = var_problem.mesh.IsElementAffine(i); 
                TinyVector<T, 3> vec_u, vec_v, vec_P, vec_H; Real_wp poids;
                for (int k = 0; k < nb_pts_quad; k++)
                  {
                    vec_H.Zero();
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        ExtractVector(b_source, offset_P + 3*(kp*nb_pts_quad+k), vec_P);
                        
                        vec_u = coef_E(kp)*vec_P;
                        if (affine)
                          {
                            poids = var_problem.Glob_jacobian(i)(0)*Fb.WeightsND(k);
                            Mlt(var_problem.Glob_DFjm1(i)(0), vec_u, vec_v);
                            vec_v *= Real_wp(1) / var_problem.Glob_jacobian(i)(0);
                          }
                        else
                          {
                            poids = var_problem.Glob_jacobian(i)(k);
                            Real_wp jacob = var_problem.Glob_jacobian(i)(k)/Fb.WeightsND(k);
                            Mlt(var_problem.Glob_DFjm1(i)(k), vec_u, vec_v);
                            vec_v *= Real_wp(1) / jacob;
                          }
                        
                        vec_H += vec_v;
                        
                        vec_P *= Real_wp(1) / (coefDiag(kp)*poids);
                        
                        ExtractVector(vec_P, offset_P+3*(kp*nb_pts_quad+k), b_source);
                      }

                    CopyVector(vec_H, k, y);
                  }
                
                TinyVector<Vector<T>, 1> x;
                x(0).Reallocate(nb_dof_loc); x(0).Zero();
                
                Fb.ApplyCh(y, x(0));
                var_problem.AddLocalUnknownVector(one, x, i, b_source);
              }
          }
	
        if (trans.Trans())
          {            
            xdir.Reallocate(nb_dof_loc);
            if (var_problem.IsHomogeneousDirichlet())
              xdir.Zero();
            else
              for (int j = 0; j < nb_dof_loc; j++)
                {
                  int num_dof = mesh_num.Element(i).GetNumberDof(j);
                  if (num_dof >= 0)
                    if (var_problem.IsDofDirichlet(num_dof))
                      xdir(j) = b_source(num_dof);
                }
          }
        
	// part due to elimination of internal nodes of u
	if (this->symmetric_elem_matrix)
	  {
	    if (this->store_block)
	      this->ModifyRhsGen(trans, b_source, i, this->block_sym(i).num_ddl, mat_elt_unsym,
				 this->block_sym(i).a12, this->block_sym(i).a21, this->block_sym(i).inv_a22);
	    else
	      this->ModifyRhsGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_sym);
	  }
	else
	  {
	    if (this->store_block)
	      this->ModifyRhsGen(trans, b_source, i, this->block_unsym(i).num_ddl, mat_elt_unsym,
				 this->block_unsym(i).a12, this->block_unsym(i).a21, this->block_unsym(i).inv_a22);
	    else
	      this->ModifyRhsGen(trans, b_source, i, num_ddl, mat_elt_unsym, a12, a21, a22_unsym);
          }

        if (trans.Trans())
          for (int j = 0; j < nb_dof_loc; j++)
            {
              int num_dof = mesh_num.Element(i).GetNumberDof(j);
              if (num_dof >= 0)
                if (var_problem.IsDofDirichlet(num_dof))
                  b_source(num_dof) = xdir(j);
            }
        
      }
    
    //b_source.Write("source_after.dat");
    
    if (trans.NoTrans())
      {
        if ((var_problem.GetNbGlobalEltPML() <= 0) || !var_problem.FirstOrderFormulation())
          {
            for (int i = 0; i < b_dir.GetM(); i++)
              b_source(var_problem.GetDirichletDofNumber(i)) = b_dir(i);
          }
        else
          for (int i = 0; i < b_dir.GetM(); i++)
            {
              int num_dof = var_problem.GetDirichletDofNumber(i);
              if (num_dof < mesh_num.GetNbDof())
                b_source(num_dof) = b_dir(i);
            }
      }
    
    var_problem.SetLeafStaticCondensation(false);

  }


  //! returns the memory used by the object in bytes
  template<class T, class TypeEquation>
  size_t CondensationBlockSolver_Maxwell3D<T, TypeEquation>::GetMemorySize() const
  {
    size_t taille = CondensationBlockSolver_Fem<T>::GetMemorySize();
    return taille;
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    HarmonicMaxwell_3D<Complex_wp>::
      ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, this->GetReferenceElementHcurl(i));
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    HarmonicMaxwell_3D<Complex_wp>::
      AddElementaryFluxesSipg(mat_sp, nat_mat, offset_row, offset_col);
  }
  

  void FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquation_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddHcurl3D(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		   X, Complex_wp(1, 0), Y, false);
  }


  /*************
   * Impedance *
   *************/

  
  //! impedance coefficient
  template<class T, class Complexe>
  void ImpedanceABC_Maxwell3D<T, Complexe>::
  EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                       const GlobalGenericMatrix<T>& nat_mat, int ref,
                       const SetPoints<Dimension3>& Pts, const SetMatrices<Dimension3>& Mat)
  {
    // variable impedance => allocating arrays
    if (k == 0)
      {
	int nb_points = Mat.GetNbPointsQuadratureBoundary();
	this->stored_coef_phi(0).Reallocate(nb_points);
	this->stored_coef_phi(0).Fill(0);

        this->vec_normale.Reallocate(Mat.GetNbPointsQuadratureBoundary());
        for (int kp = 0; kp < Mat.GetNbPointsQuadratureBoundary(); kp++)
          this->vec_normale(kp) = Mat.GetNormaleQuadratureBoundary(kp);
      }    
    
    T coef;
    SetComplexZero(coef);
    if (var_problem.FirstOrderFormulation())
      coef = var_maxwell.coefficient_impedance_absorbing(ref)*nat_mat.GetCoefDamping();
    else
      {
        Complexe m_iomega; 
        var_problem.GetMiomega(m_iomega);
        coef = m_iomega*var_maxwell.coefficient_impedance_absorbing(ref)*nat_mat.GetCoefDamping();
        if (var_boundary.take_into_account_curvature_for_abc)
          {
            Real_wp k1 = Mat.GetK1QuadratureBoundary(k);
            Real_wp k2 = Mat.GetK2QuadratureBoundary(k);
            Real_wp h = 0.5*(k1+k2);
            coef += var_maxwell.ref_invMu(ref).GetConstant()(0,0)*h*nat_mat.GetCoefStiffness();
          }
      }
    
    this->stored_coef_phi(0)(k) = coef;
  }
  
  
  /*********************************
   * HarmonicMaxwellEquation_3D_DG *
   *********************************/
  
  bool HarmonicMaxwellEquation_3D_DG::store_dfjm1(true);
  
  bool HarmonicMaxwellEquation_3D_DG::SymmetricGlobalMatrix()
  {
    return true;
  }

  bool HarmonicMaxwellEquation_3D_DG::SymmetricElementaryMatrix()
  {
    return true;
  }

  //! empty method
  template<class TypeEquation>
  void HarmonicMaxwellEquation_3D_DG::
  ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		    int i, const ElementReference_Dim<Dimension3>&)
  {
  }
   

  //! returns the matrix M, in the integral \f$ \int_K M u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] mass matrix M
  */  
  template<class TypeEquation, class T0, class MatMass>
  void HarmonicMaxwellEquation_3D_DG::
  GetTensorMass(const EllipticProblem<TypeEquation>& var,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		int ref, MatMass& mass)
  {
    FillZero(mass);
    Complex_wp coef = var.GetOmega()*nat_mat.GetCoefMass();
    mass(0, 0) = -var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 0)*coef;
    mass(1, 1) = -var.ref_epsilon(ref).GetCoefficient(var, i, j, 1, 1)*coef;
    mass(2, 2) = -var.ref_epsilon(ref).GetCoefficient(var, i, j, 2, 2)*coef;
    
    mass(3, 3) = -var.ref_mu(ref).GetCoefficient(var, i, j, 0, 0)*coef;
    mass(4, 4) = -var.ref_mu(ref).GetCoefficient(var, i, j, 1, 1)*coef;
    mass(5, 5) = -var.ref_mu(ref).GetCoefficient(var, i, j, 2, 2)*coef;
    
  }
  

  //! Applies matrix M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U vector to be multiplied by M
    \param[out] V result vector V = M U
  */  
  template<class TypeEquation, class T0, class Vector1>
  void HarmonicMaxwellEquation_3D_DG::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref,
		  Vector1& U, Vector1& V)
  {
    V = U;
    Complex_wp coef = var.GetOmega()*nat_mat.GetCoefMass();
    V(0) *= -var.ref_epsilon(ref).GetCoefficient(var, i, j, 0, 0)*coef;
    V(1) *= -var.ref_epsilon(ref).GetCoefficient(var, i, j, 1, 1)*coef;
    V(2) *= -var.ref_epsilon(ref).GetCoefficient(var, i, j, 2, 2)*coef;
    
    V(3) *= -var.ref_mu(ref).GetCoefficient(var, i, j, 0, 0)*coef;
    V(4) *= -var.ref_mu(ref).GetCoefficient(var, i, j, 1, 1)*coef;
    V(5) *= -var.ref_mu(ref).GetCoefficient(var, i, j, 2, 2)*coef;
  }
  

  //! which derivatives to evaluate during matrix-vector product ?  
  template<class TypeEquation, class T0, class Vector1>
  void HarmonicMaxwellEquation_3D_DG::
  GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                      const GlobalGenericMatrix<T0>& nat_mat,
                      Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    // on a besoin de ne deriver que la premiere inconnue, le champ electrique E
    unknown_to_derive(0) = true;
    unknown_to_derive(1) = true;
    unknown_to_derive(2) = true;
    
    unknown_to_derive(3) = false;
    unknown_to_derive(4) = false;
    unknown_to_derive(5) = false;
    

    fct_test_to_derive(0) = true;
    fct_test_to_derive(1) = true;
    fct_test_to_derive(2) = true;
        
    fct_test_to_derive(3) = false;
    fct_test_to_derive(4) = false;
    fct_test_to_derive(5) = false;
  }
  

  //! fills tensors D and E appearing in the variational formulation
  /*!
    \param[in] vars considered problem
    \param[in] num_elem element where D and E must be computed
    \param[in] jloc quadrature point where D and E must be computed
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the element
    \param[out] Ephi_grad tensor E
    \param[out] Dgrad_phi tensor D    
    The tensors D and E are appearing in the terms
    \int_K D \nabla u v + E u \nabla v dx
    of the variational formulation
   */  
  template<class TypeEquation, class T0, class MatStiff>
  void HarmonicMaxwellEquation_3D_DG::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
                   int num_elem, int jloc,
		   const GlobalGenericMatrix<T0>& nat_mat, int ref,
                   MatStiff& Ephi_grad, MatStiff& Dgrad_phi)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    Complex_wp c = nat_mat.GetCoefStiffness();
    Dgrad_phi(4,2)(0) = -c; Dgrad_phi(5,1)(0) = c;
    Dgrad_phi(3,2)(1) = c; Dgrad_phi(5,0)(1) = -c;
    Dgrad_phi(3,1)(2) = -c; Dgrad_phi(4,0)(2) = c;    

    Ephi_grad(2,4)(0) = -c; Ephi_grad(1,5)(0) = c;
    Ephi_grad(2,3)(1) = c; Ephi_grad(0,5)(1) = -c;
    Ephi_grad(1,3)(2) = -c; Ephi_grad(0,4)(2) = c;

  }
  

  //! Applying the tensor D to grad(v)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Vn gradient of the unknown vector V
    \param[out] Un result to D grad(v)
  */  
  template<class TypeEquation,
           class T0, class Vector1, class Vector2>
  void HarmonicMaxwellEquation_3D_DG::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                       int ref, Vector1& Vn, Vector2& Un)
  {
    Un.Fill(0);
    Un(3) = Vn(2)(1) - Vn(1)(2);
    Un(4) = -Vn(2)(0) + Vn(0)(2);
    Un(5) = Vn(1)(0) - Vn(0)(1);
    Un *= nat_mat.GetCoefStiffness();
  }
  

  //! Applying the tensor E to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un unknown vector U
    \param[out] Vn result E Un
  */  
  template<class TypeEquation,
           class T0, class Vector1, class Vector2>
  void HarmonicMaxwellEquation_3D_DG::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
                       int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
                       int ref, Vector1& Un, Vector2& Vn)
  {
    Vn.Fill(Complex_wp(0));
    Vn(2)(1) = -Un(3); Vn(1)(2) = Un(3);
    Vn(2)(0) = Un(4); Vn(0)(2) = -Un(4);
    Vn(1)(0) = -Un(5); Vn(0)(1) = Un(5);
    Vn *= -nat_mat.GetCoefStiffness();
  }
  
  
  //! computation of matrix N associated to the boundary condition
  /*!
    \param[out] Nabc matrix N
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref_d reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Matrix1, class GenericPb, class T0>
  void HarmonicMaxwellEquation_3D_DG::
  GetNabc(Matrix1& Nabc, R3& normale, int ref, int iquad, int k,
	  const GlobalGenericMatrix<T0>& nat_mat, int ref2, const GenericPb& vars,
	  const ElementReference<Dimension3, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Nabc.Fill(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	Nabc(0,4) = -normale(2); Nabc(0,5) = normale(1);
	Nabc(1,5) = -normale(0); Nabc(1,3) = normale(2);
	Nabc(2,3) = -normale(1); Nabc(2,4) = normale(0);
	Nabc(3,1) = normale(2); Nabc(3,2) = -normale(1);
	Nabc(4,2) = normale(0); Nabc(4,0) = -normale(2);
	Nabc(5,0) = normale(1); Nabc(5,1) = -normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	Nabc(0,4) = normale(2); Nabc(0,5) = -normale(1);
	Nabc(1,5) = normale(0); Nabc(1,3) = -normale(2);
	Nabc(2,3) = normale(1); Nabc(2,4) = -normale(0);
	Nabc(3,1) = -normale(2); Nabc(3,2) = normale(1);
	Nabc(4,2) = -normale(0); Nabc(4,0) = normale(2);
	Nabc(5,0) = -normale(1); Nabc(5,1) = normale(0);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	Nabc(0,0) = -1.0; Nabc(0,2) = normale(0)*normale(2); Nabc(0,1) = normale(0)*normale(1);
	Nabc(1,1) = -1.0; Nabc(1,0) = Nabc(0,1); Nabc(1,2) = normale(2)*normale(1);
	Nabc(2,2) = -1.0; Nabc(2,0) = Nabc(0,2); Nabc(2,1) = Nabc(1,2);
	
	Nabc(3,3) = -1.0; Nabc(3,4) = Nabc(0,1); Nabc(3,5) = Nabc(0,2);
	Nabc(4,4) = -1.0; Nabc(4,3) = Nabc(0,1); Nabc(4,5) = Nabc(1,2);
	Nabc(5,5) = -1.0; Nabc(5,3) = Nabc(0,2); Nabc(5,4) = Nabc(1,2);
	Mlt(Iwp, Nabc);
      }

    Nabc *= nat_mat.GetCoefStiffness();
  }
  
  
  //! Computation of "penalization" matrices C
  /*!
    \param[out] Nabc penalization matrix C
    \param[in] normale outward normale
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Matrix1, class GenericPb, class T0>
  void HarmonicMaxwellEquation_3D_DG::
  GetPenalDG(Matrix1& Nabc, R3& normale, int iquad, int k, int nf,
             const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
	     const GenericPb& vars, const ElementReference<Dimension3, 1>& Fb)
  {
    Nabc.Fill(0);
    Nabc(0,0) = -1.0; Nabc(0,2) = normale(0)*normale(2); Nabc(0,1) = normale(0)*normale(1);
    Nabc(1,1) = -1.0; Nabc(1,0) = Nabc(0,1); Nabc(1,2) = normale(2)*normale(1);
    Nabc(2,2) = -1.0; Nabc(2,0) = Nabc(0,2); Nabc(2,1) = Nabc(1,2);
    
    Nabc(3,3) = -1.0; Nabc(3,4) = Nabc(0,1); Nabc(3,5) = Nabc(0,2);
    Nabc(4,4) = -1.0; Nabc(4,3) = Nabc(0,1); Nabc(4,5) = Nabc(1,2);
    Nabc(5,5) = -1.0; Nabc(5,3) = Nabc(0,2); Nabc(5,4) = Nabc(1,2);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	{
	  Nabc(i,j) *= Iwp*vars.alpha_penalization;
	  Nabc(3+i,3+j) *= Iwp*vars.delta_penalization;
	}

    Nabc *= nat_mat.GetCoefStiffness();
  }
   
  
  //! Multiplication by penalization matrices
  /*!
    \param[in] normale outward normale
    \param[in] Vn vector to multiply
    \param[out] Un result vector Un = C*Vn
    \param[in] iquad element number
    \param[in] k local quadrature point number    
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Vector1, class Vector2, class GenericPb, class T0>
  void HarmonicMaxwellEquation_3D_DG::
  MltPenalDG(const R3& normale, const Vector1& Vn, Vector2& Un,
             int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
	     int ref, int ref2, const GenericPb& vars,
	     const ElementReference<Dimension3, 1>& Fb)
  {
    // int ref_elem = vars->mesh.elements(i).GetReference();
    // condition absorbante sur les sauts
    // Un.Fill(0);
    Un(0) = -(Vn(0)-normale(0)*(normale(2)*Vn(2)+normale(1)*Vn(1)))*vars.alpha_penalization;
    Un(1) = -(Vn(1)-normale(1)*(normale(0)*Vn(0)+normale(2)*Vn(2)))*vars.alpha_penalization;
    Un(2) = -(Vn(2)-normale(2)*(normale(1)*Vn(1)+normale(0)*Vn(0)))*vars.alpha_penalization;
    
    Un(3) = (Vn(3)-normale(0)*(normale(2)*Vn(5)+normale(1)*Vn(4)))*vars.delta_penalization;
    Un(4) = (Vn(4)-normale(1)*(normale(0)*Vn(3)+normale(2)*Vn(5)))*vars.delta_penalization;
    Un(5) = (Vn(5)-normale(2)*(normale(1)*Vn(4)+normale(0)*Vn(3)))*vars.delta_penalization;
    
    Mlt(Iwp*nat_mat.GetCoefStiffness(), Un);
  }
  
  
  //! multiplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] num_elem1 element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Vector1, class TypeEquation, class T0>
  void HarmonicMaxwellEquation_3D_DG::
  MltNabc(R3& normale, int ref, const Vector1& Vn, Vector1& Un,
          int num_elem1, int k, const GlobalGenericMatrix<T0>& nat_mat,
	  int ref2, const EllipticProblem<TypeEquation>& vars,
	  const ElementReference<Dimension3, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Un.Fill(0);
    if (cond == vars.LINE_DIRICHLET)
      {
	Un(0) = -Vn(4)*normale(2) + Vn(5)*normale(1);
	Un(1) = -Vn(5)*normale(0) + Vn(3)*normale(2);
	Un(2) = -Vn(3)*normale(1) + Vn(4)*normale(0);
	Un(3) = Vn(1)*normale(2) - Vn(2)*normale(1);
	Un(4) = Vn(2)*normale(0) - Vn(0)*normale(2);
	Un(5) = Vn(0)*normale(1) - Vn(1)*normale(0);
      }
    else if (cond == vars.LINE_NEUMANN)
      {
	Un(0) = Vn(4)*normale(2) - Vn(5)*normale(1);
	Un(1) = Vn(5)*normale(0) - Vn(3)*normale(2);
	Un(2) = Vn(3)*normale(1) - Vn(4)*normale(0);
	Un(3) = -Vn(1)*normale(2) + Vn(2)*normale(1);
	Un(4) = -Vn(2)*normale(0) + Vn(0)*normale(2);
	Un(5) = -Vn(0)*normale(1) + Vn(1)*normale(0);
      }
    else if (cond == vars.LINE_ABSORBING)
      {
	// int ref_elem = vars->mesh.elements(num_elem1).GetReference();
	// condition absorbante
	Un(0) = -(Vn(0)-normale(0)*(normale(2)*Vn(2)+normale(1)*Vn(1)));
	Un(1) = -(Vn(1)-normale(1)*(normale(0)*Vn(0)+normale(2)*Vn(2)));
	Un(2) = -(Vn(2)-normale(2)*(normale(1)*Vn(1)+normale(0)*Vn(0)));
	
	Un(3) = (Vn(3)-normale(0)*(normale(2)*Vn(5)+normale(1)*Vn(4)));
	Un(4) = (Vn(4)-normale(1)*(normale(0)*Vn(3)+normale(2)*Vn(5)));
	Un(5) = (Vn(5)-normale(2)*(normale(1)*Vn(4)+normale(0)*Vn(3)));
	Mlt(Iwp, Un);
      }
    
    Un *= nat_mat.GetCoefStiffness();
    // DISP(Un); DISP(Vn);
  }    
  
  
  /*****************************
   * TransparencySolverMaxwell *
   *****************************/
  

  TransparencySolver<HarmonicMaxwellEquation_3D>
  ::TransparencySolver(EllipticProblem<HarmonicMaxwellEquation_3D>& var,
		       All_LinearSolver& solver) :
    TransparencySolver_Maxwell3D(var, solver) 
  {
  }


  TransparencySolver<HarmonicMaxwellEquation_3D_DG>
  ::TransparencySolver(EllipticProblem<HarmonicMaxwellEquation_3D_DG>& var,
		       All_LinearSolver& solver) :
    TransparencySolver_Maxwell3D(var, solver)
  {
  }


  TransparencySolver<HarmonicMaxwellEquationHdg_3D>
  ::TransparencySolver(EllipticProblem<HarmonicMaxwellEquationHdg_3D>& var,
		       All_LinearSolver& solver) :
    TransparencySolver_Maxwell3D(var, solver)
  {
  }


  //! computation of E and H from values defined on a closed surface, for Maxwell equations
  /*!
    \param[in] trace_En values of E \times n on quadrature points of the surface \Gamma
    \param[in] trace_Hn values of H \times n on quadrature points of the surface \Gamma
    \param[in] mesh_ surface mesh representing \Gamma
    \param[in] pointX point where E and H are computed
    \param[in] normaleX not used
    \param[out] E_pot value of E on pointX
    \param[out] H_pot value of H on pointX
   */
  void TransparencySolver_Maxwell3D::
  ComputeIntegralRepresentation(const Vector<Complex_wp>& trace_En,
				const Vector<Complex_wp>& trace_Hn,
				const MeshInterpolationFEM<Dimension3>& mesh_,
				const R3& pointX, const R3& normaleX,
				R3_Complex_wp& E_pot, R3_Complex_wp& H_pot) const
  {   
    R3 pointY;
    R3_Complex_wp En, Hn, contrib1, contrib2;
    
    // initialization to zero
    E_pot.Fill(0); H_pot.Fill(0);
    
    Complex_wp phi, poids; R3_Complex_wp grad_phi;
    Matrix3_3sym_Complex_wp hessian_phi,dyadic_G;

    Real_wp omega = var_problem.GetOmega();
    Real_wp omega2 = omega*omega;      
    // loop on quadrature points
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // normaleY = mesh_.GetQuadratureNormale(k); Mlt(-1, normaleY);
        pointY =  mesh_.GetQuadraturePoint(k);
        poids = mesh_.GetQuadratureWeight(k);
        
	CopyVector(trace_En, k, En);
	CopyVector(trace_Hn, k, Hn);
        
        this->ComputeGreenKernel(pointX, pointY, omega, phi, grad_phi, hessian_phi);
        // G = \phi(x,y) I + 1/k^2 \nabla_y \nabla_y \phi(x,y)
        Copy(hessian_phi, dyadic_G); Mlt(Complex_wp(1.0/omega2), dyadic_G);
        
        dyadic_G(0, 0) += phi;
        dyadic_G(1, 1) += phi;
        dyadic_G(2, 2) += phi;
        
        // contrib1 = ik G n \times H
        Seldon::Mlt(dyadic_G, Hn, contrib1); contrib1 *= Iwp*omega;
	
        // contrib2 = (n \times E) \times (\nabla_y \phi)
        Seldon::TimesProd(En, grad_phi, contrib2);
	
        // DISP(phi); DISP(grad_phi); DISP(hessian_phi); DISP(contrib1); DISP(contrib2);
        contrib2 += contrib1; E_pot += poids*contrib2;
	
        // contrib1 = ik G n \times E
        Seldon::Mlt(dyadic_G, En, contrib1); contrib1 *= -Iwp*omega;
	
        // contrib2 = (n \times H) \times (\nabla_y \phi)
        Seldon::TimesProd(Hn, grad_phi, contrib2);
	
        // DISP(phi); DISP(grad_phi); DISP(hessian_phi); DISP(contrib1); DISP(contrib2);
        contrib2 += contrib1; Seldon::Add(poids, contrib2, H_pot); 
	
      }
    
  }
  
  
  //! computation of transparent source for Maxwell equations
  void TransparencySolver_Maxwell3D
  ::GetSource(const VectComplex_wp& traceEn, const VectComplex_wp& traceHn,
	      int n, const Real_wp& k_inf, const R3& point, const R3& normale,
	      Vector<VectComplex_wp>& scal_g, int j) const
  {
    if (var_problem.FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
	R3_Complex_wp tmp, En, Hn_times_n, Epot, Hpot;
	CopyVector(traceEn, n, Epot);
	CopyVector(traceHn, n, Hpot);
	
	TimesProd(Epot, normale, En);
	TimesProd(Hpot, normale, tmp); TimesProd(tmp, normale, Hn_times_n);
	
	tmp = En - Hn_times_n;
	
	TimesProd(normale, tmp, Epot);
        
	Mlt(-Iwp*k_inf, Epot);
	CopyVector(Epot, j, scal_g(0));
      }
    else
      {
	// not implemented
        abort();
      }  
  }


  /*****************************
   * VarComputationRCS_Maxwell *
   *****************************/


  void VarComputationRCS_Maxwell3D::ComputeRCS(const VectReal_wp& U0)
  {
    cout << "not possible" <<endl;
    abort();
  }


  //! computation of radar cross section
  /*!
    \param[in] U0 solution vector
    \param[in] infoRCS parameters for radar cross section
    The radar cross section is computed and written on output files 
   */
  
  void VarComputationRCS_Maxwell3D::ComputeRCS(const VectComplex_wp& U0)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    int nbPointsRCS = this->nb_angles_RCS;
   
    this->InitComputationRCS(true);
    
    VectComplex_wp trace_En;
    VectComplex_wp trace_Hn;
    var_problem.ComputeEnHnOnBoundary(this->var_mesh, U0, trace_En, trace_Hn);
    
    R3 ur; R3_Complex_wp En, Hn, Hn_u, Einf; 
    Real_wp teta;
    Real_wp step_angle = (this->last_angle_RCS-this->first_angle_RCS) / (nbPointsRCS-1);
    R3_Complex_wp rcs_teta;
    
    int nb_points_effective = nbPointsRCS;
    int num_angle0 = 0;
    GetParallelDistributionPoints(nb_proc, rank_proc, nbPointsRCS, nb_points_effective, num_angle0);
    
    VectReal_wp RCSField(nb_points_effective);
    VectReal_wp TetaField(nb_points_effective);
    
    R3 Etilde, ktilde;
    var_problem.GetPolarization(Etilde);
    
    Real_wp omega = var_problem.GetOmega();
    Mlt(1.0/Norm2(Etilde), Etilde);
    ktilde = var_problem.GetWaveVector(); Mlt(1.0/Norm2(ktilde), ktilde);
    
    // loop on each angle
    for (int num_angle = num_angle0; num_angle < num_angle0+nb_points_effective; num_angle++)
      {
	rcs_teta.Fill(0); Einf.Fill(0);
	teta = this->first_angle_RCS + step_angle*num_angle;
	ur(0) = cos(teta)*ktilde(0) + sin(teta)*Etilde(0); 
	ur(1) = cos(teta)*ktilde(1) + sin(teta)*Etilde(1); 
	ur(2) = cos(teta)*ktilde(2) + sin(teta)*Etilde(2);
	
	// integration over the boundary
	for (int k = 0; k < this->var_mesh.GetNbAllQuadraturePoints(); k++)
          {
	    CopyVector(trace_En, k, En);
	    CopyVector(trace_Hn, k, Hn);
	    
	    Real_wp kr = DotProd(ur, this->var_mesh.GetQuadraturePoint(k))*omega;
            Complex_wp arg = Complex_wp(cos(kr), sin(kr));
	    
            TimesProd(ur, Hn, Hn_u);
	    Hn_u += En;
            Add(this->var_mesh.GetQuadratureWeight(k)*arg, Hn_u, Einf);
	  }
	
	Mlt(Iwp*omega, Einf);
	TimesProd(ur, Einf, rcs_teta);
	
	Real_wp module_h = Norm2(rcs_teta); module_h *= module_h;
	
        int m = num_angle - num_angle0;
	TetaField(m) = 180.0*teta/pi_wp;
        RCSField(m) = 10.0*log(1.0/(4.0*pi_wp)*module_h)/log(10.0);
        
      }
    
    ofstream file_out;
    if (rank_proc == 0)
      {
        file_out.open(this->file_RCS.data());
        file_out.precision(cout.precision());
      }
    
    VectReal_wp AllRCS, AllTeta;

    if (nb_proc > 1)
      {

#ifdef SELDON_WITH_MPI
	Vector<int64_t> xtmp;
	IVect NbPointsPerProc, OffsetPointsProc;
	// nmax is the maximum number of points (among all proc)
	int nmax = nb_points_effective;
	NbPointsPerProc.Reallocate(nb_proc);
	OffsetPointsProc.Reallocate(nb_proc);
	
        MPI_Gather(&nb_points_effective, 1, MPI_INTEGER,
                   NbPointsPerProc.GetData(), 1, MPI_INTEGER, 0, var_problem.comm_group_mode);
	
        MPI_Gather(&num_angle0, 1, MPI_INTEGER,
                   OffsetPointsProc.GetData(), 1, MPI_INTEGER, 0, var_problem.comm_group_mode);
	
	nmax = NbPointsPerProc.GetNormInf();
        
        MPI_Bcast(&nmax, 1, MPI_INTEGER, 0, var_problem.comm_group_mode);
	
        if (nmax != nb_points_effective)
          {
            RCSField.Resize(nmax);
            TetaField.Resize(nmax);
          }
	
	int Nall = nmax*nb_proc;
        AllRCS.Reallocate(Nall);
        AllTeta.Reallocate(Nall);
        
	MpiGather(var_problem.comm_group_mode, RCSField, xtmp, AllRCS, nmax, 0);
	MpiGather(var_problem.comm_group_mode, TetaField, xtmp, AllTeta, nmax, 0);
	
        if (rank_proc == 0)
	  for (int p = 0; p < nb_proc; p++)
	    for (int k = 0; k < NbPointsPerProc(p); k++)
	      {
		int num_angle = OffsetPointsProc(p) + k;
		file_out<<AllTeta(num_angle)<<"  "<<AllRCS(num_angle)<<'\n';
	      }
#endif
        
      }
    else
      {
        for (int num_angle = 0; num_angle < nbPointsRCS; num_angle++)
          file_out<<TetaField(num_angle)<<"  "<<RCSField(num_angle)<<'\n';
      }
    
    if (rank_proc == 0)
      file_out.close();

    if ((rank_proc == 0)||(this->print_level >= 10))
      if (this->print_level >= 0)
	cout<<"Radar cross section OK"<<endl;
  }
  

  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
			const GridInterpolation<Dimension3>& var_interp,
			int iquad, bool compute_grad) const
  {
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
			const GridInterpolation<Dimension3>& var_interp,
			int iquad, bool compute_grad) const
  {
    if (compute_grad)
      {
	Complex_wp coef = -Iwp/this->omega;
	grad_u(0) *= coef;
        grad_u(1) *= coef;
        grad_u(2) *= coef;
      }
  }
  

  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			int i, bool u_component, bool grad_component) const
  {
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			int i, bool u_component, bool grad_component) const
  {
    if (grad_component)
      {
	Complex_wp coef = -Iwp/this->omega;
	for (int k = 0; k < GradNodal(0).GetM(); k++)
	  {
	    GradNodal(0)(k) *= coef;
	    GradNodal(1)(k) *= coef;
	    GradNodal(2)(k) *= coef;
	  }
      }
  }
  

  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& curl_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const
  {
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& curl_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const
  {
    HarmonicMaxwell_3D_Base::ComputeEnHnQuadratureHcurl(u_quadrature, curl_quadrature,
                                                        num_elem, pts, normale, compute_H, En_quad, Hn_quad);
  }
  

  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ComputeEnHnNodal(Vector<VectReal_wp >& u_nodal, Vector<VectReal_wp >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<VectReal_wp >& En_nodal, Vector<VectReal_wp >& Hn_nodal) const
  {
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D>
  ::ComputeEnHnNodal(Vector<VectComplex_wp >& u_nodal, Vector<VectComplex_wp >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<VectComplex_wp >& En_nodal, Vector<VectComplex_wp >& Hn_nodal) const
  {
    HarmonicMaxwell_3D_Base::ComputeEnHnNodal(u_nodal, grad_nodal, num_elem, pts, normale,
					      En_nodal, Hn_nodal);
  }


  IncidentWaveProjector<Complex_wp, Dimension3>*
  EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::GetNewIncidentProjector(int n, const Vector<VectString>& param,
			    IncidentWaveField<Complex_wp, Dimension3>& u_inc) const
  {
    return new IncidentWaveProjector_Maxwell3D_DG(*this, u_inc);
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
			  const GlobalGenericMatrix<Real_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
			  const GlobalGenericMatrix<Complex_wp>& nat_mat,
			  int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }

  
  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>& solver,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u, int i,
			const GridInterpolation<Dimension3>& var_interp,
			int iquad, bool compute_grad) const
  {
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u, int i,
			const GridInterpolation<Dimension3>& var_interp,
			int iquad, bool compute_grad) const
  {
    val_u(3) *= -this->GetOmega();
    val_u(4) *= -this->GetOmega();
    val_u(5) *= -this->GetOmega();    
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ModifyOutputUnknown(Vector<VectReal_wp>& Unodal, Vector<VectReal_wp>& GradNodal,
			int i, bool u_component, bool grad_component) const
  {
  }


  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>& Unodal, Vector<VectComplex_wp>& GradNodal,
			int i, bool u_component, bool grad_component) const
  {
    if (u_component)
      {
	for (int k = 0; k < Unodal(0).GetM(); k++)
	  {
	    Unodal(3)(k) *= -this->GetOmega();
	    Unodal(4)(k) *= -this->GetOmega();
	    Unodal(5)(k) *= -this->GetOmega();	    
	  }
      }
  }

  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& curl_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const
  {
  }

  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& curl_quadrature,
			  int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const
  {
    HarmonicMaxwell_3D_Base::ComputeEnHnQuadratureDG(u_quadrature, curl_quadrature,
                                                     num_elem, pts, normale, compute_H, En_quad, Hn_quad);
  }
  
  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ComputeEnHnNodal(Vector<VectReal_wp >& u_nodal, Vector<VectReal_wp >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<VectReal_wp >& En_nodal, Vector<VectReal_wp >& Hn_nodal) const
  {
  }
  
  void EllipticProblem<HarmonicMaxwellEquation_3D_DG>
  ::ComputeEnHnNodal(Vector<VectComplex_wp >& u_nodal, Vector<VectComplex_wp >& grad_nodal,
		     int num_elem, const Vector<R3>& pts, const Vector<R3>& normale,
		     Vector<VectComplex_wp >& En_nodal, Vector<VectComplex_wp >& Hn_nodal) const
  {
    HarmonicMaxwell_3D_Base::ComputeEnHnNodal(u_nodal, grad_nodal, num_elem, pts, normale,
					      En_nodal, Hn_nodal);
  }


  template<class T>
  bool MaxwellEquationHdg_3D<T>::store_dfjm1(true);


  template<class T> template<class TypeEquation>
  void MaxwellEquationHdg_3D<T>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var_problem,
                      int iquad, const ElementReference_Dim<Dimension3>& Fb0)
  {
    const ElementReference<Dimension3, 2>& Fb = dynamic_cast<const ElementReference<Dimension3, 2>& >(Fb0);

    if (!Fb.UsePiolaTransform())
      return;

    HarmonicMaxwell3D_PhysGeomInfo<T>& mass = var_problem.Glob_matMass_elem(iquad);
    int ref_domain = var_problem.mesh.Element(iquad).GetReference();
    bool variable = var_problem.UseNumericalIntegration(iquad);
    bool affine = var_problem.mesh.IsElementAffine(iquad);
	
    bool diag_mass = false;
    if ((var_problem.OrthogonalElement(iquad) == 0) && Fb.MassLumpingOrthogonalElement() && affine)
      diag_mass = true;
    
    if ((var_problem.ref_epsilon(ref_domain).GetAnisotropy() == TensorPhysicalIndice<Dimension3, 3, T>::ANISOTROPE) ||
	(var_problem.ref_sigma(ref_domain).GetAnisotropy() == TensorPhysicalIndice<Dimension3, 3, T>::ANISOTROPE) )
      diag_mass = false;

    if (!diag_mass)
      mass.ClearDiagonalDh();

    int N = Fb.GetNbPointsQuadratureInside();
    if (!variable)
      N = 1;

    T m_iomega;
    var_problem.GetMiomega(m_iomega);
    
    Matrix3_3 dfjm1; Real_wp jacob, jacob_weighted;
    TinyMatrix<T, Symmetric, 3, 3> epsilon, sigma, mu, Bmass, Cdamp, Astiff;
    TinyMatrix<T, General, 3, 3> A_tmp;
    int i1 = iquad - var_problem.mesh.GetNbElt() + var_problem.GetNbEltPML();
    
    if (diag_mass)
      {
	epsilon = var_problem.ref_epsilon(ref_domain).GetConstant();
	sigma = var_problem.ref_sigma(ref_domain).GetConstant();
	mu = var_problem.ref_mu(ref_domain).GetConstant()(0, 0);
	if (!affine)
	  { cout << "impossible " << endl; abort(); }

	jacob = var_problem.Glob_jacobian(iquad)(0);
	dfjm1 = var_problem.Glob_DFjm1(iquad)(0);                
	Mlt(1.0/jacob, dfjm1);

	int Ndof = Fb.GetNbDof();
	mass.ReallocateDhDiag(Ndof, Ndof);
	
	MltTrans(epsilon, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Bmass);
	Mlt(m_iomega*jacob, Bmass);

	MltTrans(sigma, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Cdamp); 
	Mlt(jacob, Cdamp);

	MltTrans(mu, dfjm1, A_tmp);
	Mlt(dfjm1, A_tmp, Astiff);
 	Mlt(m_iomega*jacob, Astiff);

	const VectReal_wp& weights_dof = Fb.WeightsDofND();

	for (int i = 0; i < Ndof; i++)
	  {
	    int n = Fb.GetCoordinateDof(i);
	    
	    mass.SetDh(i, Bmass(n, n)*weights_dof(i));
	    mass.SetDhSigma(i, Cdamp(n, n)*weights_dof(i));
	    mass.SetDhStiff(i, Astiff(n, n)*weights_dof(i));
	  }
      }
    else
      {
	mass.ReallocateAhBh(N, N, N);
	for (int i = 0; i < N; i++)
	  {
	    epsilon = var_problem.ref_epsilon(ref_domain).GetCoefficient(var_problem, iquad, i);
	    sigma = var_problem.ref_sigma(ref_domain).GetCoefficient(var_problem, iquad, i);
	    mu = var_problem.ref_mu(ref_domain).GetCoefficient(var_problem, iquad, i);
	    // indices are modified inside PML
	    if ((var_problem.InsidePML(iquad)) && (!var_problem.FirstOrderFormulation()))
	      var_problem.ModifyPhysicalCoefPML(epsilon, mu, sigma, i1, i);
	    
	    var_problem.GetInverseJacobianMatrix(iquad, i, affine, Fb.GetGeometricElement(),
						 dfjm1, jacob, jacob_weighted);
	    
	    if (!variable)
	      jacob_weighted = jacob;
	    
	    MltTrans(epsilon, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Bmass); 
	    Mlt(m_iomega*jacob_weighted, Bmass);

	    mass.SetBh(i, Bmass);
	    
	    MltTrans(sigma, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Cdamp); 
	    Mlt(jacob_weighted, Cdamp);
	    
	    mass.SetBhSigma(i, Cdamp);
	    
	    MltTrans(mu, dfjm1, A_tmp);
	    Mlt(dfjm1, A_tmp, Astiff);
	    Mlt(m_iomega*jacob_weighted, Astiff);
	    mass.SetAh(i, Astiff);
	  }
      }
    
    N = 0;
    for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
      N += Fb.GetNbQuadBoundary(num_loc);

    mass.ReallocateBhSurf(N);
    N = 0;
    Matrix3_3 dfj_trans;
    T coef_tau = -var_problem.alpha_penalization;
    if (var_problem.upwind_fluxes)
      coef_tau = -var_problem.alpha_penalization*var_problem.coefficient_impedance_absorbing(ref_domain);
    
    TinyMatrix<T, Symmetric, 2, 2> mat_DF_DFt;
    TinyMatrix<Real_wp, Symmetric, 3, 3> DF_DFt;
    for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
      {
	int num_face, ref_boundary, rf, rot; bool new_face; T phase, phase_conj;
	const Matrix<int>& FacesQuadRotation =
	  var_problem.GetGeometryPhaseData(iquad, num_loc,
					   num_face, ref_boundary, rf, new_face, rot,
					   phase, phase_conj);

	const ElementReference<Dimension2, 2>& Fb_s = dynamic_cast<const ElementReference<Dimension2, 2>& >(var_problem.GetSurfaceFiniteElement(num_face));
	
	// for hexas only
	int nx = 0, ny = 1;
	if ((num_loc == 0) || (num_loc == 5))
	  { nx = 1; ny = 2; }
	else if ((num_loc == 1) || (num_loc == 4))
	  { nx = 0; ny = 2; }	

	Real_wp dsj;
	for (int k = 0; k < Fb.GetNbQuadBoundary(num_loc); k++)
	  {
	    int krot = FacesQuadRotation(rot, k);
	    
	    if (new_face)
	      dsj = var_problem.Glob_dsj(num_face)(k);
	    else
	      dsj = var_problem.Glob_dsj(num_face)(krot);
	    
	    int kvol = Fb.GetQuadNumber(num_loc, k);
	    var_problem.GetInverseJacobianMatrix(iquad, kvol, affine,
						 Fb.GetGeometricElement(),
						 dfjm1, jacob, jacob_weighted);
	    
	    GetInverse(dfjm1, dfj_trans); Transpose(dfj_trans);
	    MltTrans(dfj_trans, dfj_trans,  DF_DFt);
	    Mlt(Fb_s.WeightsND(k)/dsj, DF_DFt);
	    
	    mat_DF_DFt(0, 0) = coef_tau*DF_DFt(ny, ny);
	    mat_DF_DFt(0, 1) = -coef_tau*DF_DFt(nx, ny);
	    mat_DF_DFt(1, 1) = coef_tau*DF_DFt(nx, nx);

	    mass.SetBhSurf(N + k, mat_DF_DFt);
	  }
	
	N += Fb.GetNbQuadBoundary(num_loc);
      }
  }
  
	      
  template<class T>
  size_t CondensedBlockSolver_MaxwellHdg3D<T>::GetMemorySize() const
  {
    size_t taille = invBtilde.GetMemorySize() + invBh.GetMemorySize() + SnD.GetMemorySize() + Mint.GetMemorySize();
    return taille;
  }
  

  template<class T, class Complexe>
  void CondensationBlockSolver_MaxwellHdg3D<T, Complexe>
  ::SetOptimizedCondensation(int i, bool optim)
  {
    if (optim_elt.GetM() <= 0)
      optim_elt.Reallocate(this->GetNbCondensedElt());
    
    optim_elt(i).flag = optim;
  }


  template<>
  void CondensationBlockSolver_MaxwellHdg3D<Real_wp, Complex_wp>
  ::ModifyRhsStaticCondensation(const SeldonTranspose& trans, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }

  
  template<class T, class Complexe>
  void CondensationBlockSolver_MaxwellHdg3D<T, Complexe>
  ::ModifyRhsStaticCondensation(const SeldonTranspose& trans, Vector<T>& b_source, const GlobalGenericMatrix<T>& nat_mat) const
  {
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nl = mesh_num.GetNbDof();
    Matrix<T> mat_elt_unsym;
    int Nvol = var_problem.GetNbMainUnknownDof();

    Vector<T> Fu, Fv, Ft;
    TinyVector<Vector<T>, 1> Fl_;
    Vector<T>& Fl = Fl_(0);

    this->level_volume.SetLevel();
    int nb_elt_lvl = this->level_volume.GetNbElt();
    T coef_s = Real_wp(1) / nat_mat.GetCoefStiffness();
    for (int i0 = 0; i0 < nb_elt_lvl; i0++)
      {
	int iquad = this->level_volume.GetElementNumber(i0);
	if (!optim_elt(i0).flag)
	  {
	    if (!this->store_block)
	      {
		cout << "Not implemented" << endl;
		abort();
	      }

	    if (this->symmetric_elem_matrix)
	      CondensationBlockSolver_Fem<T>::
		ModifyRhsGen(trans, b_source, iquad, this->block_sym(i0).num_ddl, mat_elt_unsym,
			     this->block_sym(i0).a12, this->block_sym(i0).a21, this->block_sym(i0).inv_a22);
	    else
	      CondensationBlockSolver_Fem<T>::
		ModifyRhsGen(trans, b_source, iquad, this->block_unsym(i0).num_ddl, mat_elt_unsym,
			     this->block_unsym(i0).a12, this->block_unsym(i0).a21, this->block_unsym(i0).inv_a22);
	    
	    continue;
	  }

	const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(iquad);
	const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	
	int nb_dof_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	int offset = mesh_num.Element(iquad).GetNbDof();
	int offset_u = Nl + var_problem.GetOffsetDofV(iquad);
	    
	// etapes de modification du second membre
	Fu.Reallocate(nb_dof_elt); Fl.Reallocate(offset); Fl.Zero();
	Fv.Reallocate(nb_dof_elt);
	Ft.Reallocate(nb_dof_elt);
	TinyVector<T, 3> vec_u, vec_v;

	const Vector<TinyMatrix<T, Symmetric, 3, 3> > & invBh = this->GetInverseBh(i0);

	// on recupere Fu et Fv	
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    Fu(i) = b_source(offset_u + i);
	    Fv(i) = b_source(offset_u + Nvol + i);
	  }

	// on multiplie Fv par B_h^{-1}
	for (int i = 0; i < nb_points_quad; i++)
	  {
	    int ix, iy, iz;
	    Fb_hex.GetDofNumber_FromPointNode(i, ix, iy, iz);

	    vec_u.Init(Fv(ix), Fv(iy), Fv(iz));
	    Mlt(invBh(i), vec_u, vec_v);
	    Fv(ix) = coef_s*vec_v(0); Fv(iy) = coef_s*vec_v(1); Fv(iz) = coef_s*vec_v(2);
	  }

	// on calcule Ft = Fu - (R^T - Sn) B_h^{-1} Fv
	const Matrix<Real_wp, General, ArrayRowSparse>& R = Fb_hex.GetPermutedRh();
	Mlt(SeldonTrans, R, Fv, Ft);
	Ft += Fu;

	// on multiplie Ft par inv(Btilde)
	Mlt(this->GetBtildeMatrix(i0), Ft, Fu);

	// on calcule Fh = Fl + S_d Bh^{-1} Fv +/- Mint Btilde^{-1} Ft 
	if (this->UseSymmetrization())
	  {
	    MltAdd(Real_wp(-1), this->GetMint(i0), Fu, Real_wp(1), Fl);
	    MltAdd(Real_wp(-1), this->GetSnD(i0), Fv, Real_wp(1), Fl);
	  }
	else
	  {
	    MltAdd(Real_wp(1), this->GetMint(i0), Fu, Real_wp(1), Fl);
	    MltAdd(Real_wp(1), this->GetSnD(i0), Fv, Real_wp(1), Fl);
	  }
	
	// on met a jour Fl, Fu et Fv
	var_problem.AddLocalUnknownVector(T(1), Fl_, iquad, b_source);
	
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    b_source(offset_u + i) = Ft(i);
	    b_source(offset_u + Nvol + i) = coef_s*Fv(i);
	  }
      }
  }


  template<>
  void CondensationBlockSolver_MaxwellHdg3D<Real_wp, Complex_wp>
  ::RecomposeSolution(const SeldonTranspose&, Vector<Real_wp>& b_source, const GlobalGenericMatrix<Real_wp>& nat_mat) const
  {
    abort();
  }

  
  template<class T, class Complexe>
  void CondensationBlockSolver_MaxwellHdg3D<T, Complexe>
  ::RecomposeSolution(const SeldonTranspose& trans, Vector<T>& x_sol, const GlobalGenericMatrix<T>& nat_mat) const
  {
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);
    int Nl = mesh_num.GetNbDof();
    Matrix<T> mat_elt_unsym;
    int Nvol = var_problem.GetNbMainUnknownDof();

    TinyVector<Vector<T>, 1> Lambda_;
    Vector<T> Fu, Fv, Uloc;
    Vector<T>& Lambda = Lambda_(0);
    
    this->level_volume.SetLevel();
    int nb_elt_lvl = this->level_volume.GetNbElt();
    T coef_s = Real_wp(1) / nat_mat.GetCoefStiffness();

    for (int i0 = 0; i0 < nb_elt_lvl; i0++)
      {
	int iquad = this->level_volume.GetElementNumber(i0);
	if (!optim_elt(i0).flag)
	  {
	    if (!this->store_block)
	      {
		cout << "Not implemented" << endl;
		abort();
	      }
	    
	    if (this->symmetric_elem_matrix)
	      CondensationBlockSolver_Fem<T>::
		RecomposeSolGen(trans, x_sol, iquad, this->block_sym(i0).num_ddl, mat_elt_unsym,
				this->block_sym(i0).a12, this->block_sym(i0).a21, this->block_sym(i0).inv_a22);
	    else
	      CondensationBlockSolver_Fem<T>::
		RecomposeSolGen(trans, x_sol, iquad, this->block_unsym(i0).num_ddl, mat_elt_unsym,
				this->block_unsym(i0).a12, this->block_unsym(i0).a21, this->block_unsym(i0).inv_a22);
	    
	    continue;
	  }
	
	const ElementReference<Dimension3, 2>& Fb = var_problem.GetReferenceElementHcurl(iquad);
	const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	
	int nb_dof_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	int offset = mesh_num.Element(iquad).GetNbDof();
	int offset_u = Nl + var_problem.GetOffsetDofV(iquad);
	    
	// etapes de reconstruction de U et V
	Fu.Reallocate(nb_dof_elt); Lambda.Reallocate(offset);
	Fv.Reallocate(nb_dof_elt); Uloc.Reallocate(nb_dof_elt);
	TinyVector<T, 3> vec_u, vec_v;
	
	const Vector<TinyMatrix<T, Symmetric, 3, 3> > & invBh = this->GetInverseBh(i0);

	// on recupere Lambda, Fu et Fv
	var_problem.GetLocalUnknownVector(x_sol, iquad, Lambda_);
	
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    Fu(i) = x_sol(offset_u + i);
	    Fv(i) = x_sol(offset_u + Nvol + i);
	  }

	// on calcule U = Btilde^{-1} ( F_tilde - Mint^T Lambda)
	MltAdd(Real_wp(-1), SeldonTrans, this->GetMint(i0), Lambda, Real_wp(1), Fu);
	Mlt(this->GetBtildeMatrix(i0), Fu, Uloc);

	// on calcule F = Fv_tilde + Bh^{-1} (R U + Sd^T Lambda)
	const Matrix<Real_wp, General, ArrayRowSparse>& R = Fb_hex.GetPermutedRh();
	Mlt(R, Uloc, Fu);
	MltAdd(Real_wp(-1), SeldonTrans, this->GetSnD(i0), Lambda, Real_wp(1), Fu);
	for (int i = 0; i < nb_points_quad; i++)
	  {
	    int ix, iy, iz;
	    Fb_hex.GetDofNumber_FromPointNode(i, ix, iy, iz);

	    vec_u.Init(Fu(ix), Fu(iy), Fu(iz));
	    Mlt(invBh(i), vec_u, vec_v);
	    if (!this->UseSymmetrization())
	      vec_v = -vec_v;
	    
	    Fv(ix) += coef_s*vec_v(0); Fv(iy) += coef_s*vec_v(1); Fv(iz) += coef_s*vec_v(2);
	  }
	
	for (int i = 0; i < nb_dof_elt; i++)
	  {
	    x_sol(offset_u + i) = Uloc(i);
	    x_sol(offset_u + Nvol + i) = Fv(i);
	  }
      }
  }
  
  
  template<class T, class Complexe>
  size_t CondensationBlockSolver_MaxwellHdg3D<T, Complexe>::GetMemorySize() const
  {
    size_t taille = CondensationBlockSolver_Fem<T>::GetMemorySize();
    taille += optim_elt.GetM()*sizeof(void*);
    for (int i = 0; i < optim_elt.GetM(); i++)
      taille += optim_elt(i).GetMemorySize();
    
    return taille;
  }
  
  
  EllipticProblem<HarmonicMaxwellEquationHdg_3D>::EllipticProblem()
  {
    this->dg_formulation = ElementReference_Base::HDG;
  }

  
  void EllipticProblem<HarmonicMaxwellEquationHdg_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicMaxwellEquationHdg_3D>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>& cond_solver,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    this->ComputeElementaryMatrixHdg(i, num_dof, mat_elem,
				     cond_solver, nat_mat, this->GetReferenceElementHcurl(i));
  }


  void FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquationHdg_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Not possible" << endl;
    abort();
  }

  
  void FemMatrixFreeClass<Complex_wp, HarmonicMaxwellEquationHdg_3D>
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
               const SeldonTranspose& trans, int lvl, 
               const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddHcurlHdg3D(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		     X, Complex_wp(1, 0), Y, false);
  }

} // end namespace

#define MONTJOIE_FILE_HARMONIC_MAXWELL_3D_CXX
#endif
