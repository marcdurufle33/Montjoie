#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HELMHOLTZ_CXX

namespace Montjoie
{

  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, HelmholtzEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "not possible" << endl;
    abort();
  }


#ifdef MONTJOIE_WITH_TWO_DIM
  template<class T>
  void ApplyRhSplitElement(const ElementReference<Dimension2, 1>& Fb,
			   const Vector<T>& gradU_quad, Vector<T>& ProdUquad,
			   Vector<T>& ProdUtilde_quad,Vector<T>& ProdUterce_quad,
			   Vector<T>& ProdUloc, Vector<T>& ProdUtilde, Vector<T>& ProdUterce)
  {
    if (Fb.UseQuadraturePointsForRh())
      {
	Fb.ApplyRhQuadratureSplit(gradU_quad, ProdUquad, ProdUtilde_quad);
	T vx, vy;
	for (int i = 0; i < ProdUquad.GetM(); i++)
	  {
	    vx = ProdUquad(i); vy = ProdUtilde_quad(i);
	    ProdUquad(i) = vx + vy;
	    ProdUtilde_quad(i) = vx - vy;
	  }
      }
    else
      {
	Fb.ApplyRhSplit(gradU_quad, ProdUloc, ProdUtilde);
	T vx, vy;
	for (int i = 0; i < ProdUloc.GetM(); i++)
	  {
	    vx = ProdUloc(i); vy = ProdUtilde(i);
	    ProdUloc(i) = vx + vy;
	    ProdUtilde(i) = vx - vy;
	  }
      }
  }

  template<class T>
  void ApplyRhTransSplitElement(const ElementReference<Dimension2, 1>& Fb,
				int nb_points_quad, const T& s, const Matrix2_2& dfjm1, 
				const Vector<T>& gradUtilde_quad,
				const Vector<T>& gradUterce_quad,
				TinyVector<Vector<T>, 2>& ProdVquad)
  {
    TinyVector<T, 2> vec_u, grad_U;
    for (int k = 0; k < nb_points_quad; k++)
      {
	T poids = Fb.WeightsND(k)*s;
	CopyVector(gradUtilde_quad, k, vec_u);
	MltTrans(dfjm1, vec_u, grad_U);
	grad_U *= poids;	

	ProdVquad(0)(k) += grad_U(0);
	ProdVquad(1)(k) -= grad_U(1);
      }
  }
  
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  template<class T>
  void ApplyRhSplitElement(const ElementReference<Dimension3, 1>& Fb,
			   const Vector<T>& gradU_quad, Vector<T>& ProdUquad,
			   Vector<T>& ProdUtilde_quad,Vector<T>& ProdUterce_quad,
			   Vector<T>& ProdUloc, Vector<T>& ProdUtilde, Vector<T>& ProdUterce)
  {
    if (Fb.UseQuadraturePointsForRh())
      {
	Fb.ApplyRhQuadratureSplit(gradU_quad, ProdUquad, ProdUtilde_quad, ProdUterce_quad);
	T vx, vy, vz;
	for (int i = 0; i < ProdUquad.GetM(); i++)
	  {
	    vx = ProdUquad(i); vy = ProdUtilde_quad(i); vz = ProdUterce_quad(i);
	    ProdUquad(i) = vx + vy + vz;
	    ProdUtilde_quad(i) = vx - vy;
	    ProdUterce_quad(i) = vx - vz;
	  }
      }
    else
      {
	Fb.ApplyRhSplit(gradU_quad, ProdUloc, ProdUtilde, ProdUterce);
	T vx, vy, vz;
	for (int i = 0; i < ProdUloc.GetM(); i++)
	  {
	    vx = ProdUloc(i); vy = ProdUtilde(i); vz = ProdUterce(i);
	    ProdUloc(i) = vx + vy + vz;
	    ProdUtilde(i) = vx - vy;
	    ProdUterce(i) = vx - vz;
	  }
      }
  }

  template<class T>
  void ApplyRhTransSplitElement(const ElementReference<Dimension3, 1>& Fb,
				int nb_points_quad, const T& s, const Matrix3_3& dfjm1, 
				const Vector<T>& gradUtilde_quad,
				const Vector<T>& gradUterce_quad,
				TinyVector<Vector<T>, 3>& ProdVquad)
  {
    TinyVector<T, 3> vec_u, vec_v, grad_U, grad_V;
    for (int k = 0; k < nb_points_quad; k++)
      {
	T poids = Fb.WeightsND(k)*s;
	CopyVector(gradUtilde_quad, k, vec_u);
	CopyVector(gradUterce_quad, k, vec_v);
	MltTrans(dfjm1, vec_u, grad_U);
	MltTrans(dfjm1, vec_v, grad_V);

	ProdVquad(0)(k) += poids*(grad_U(0) + grad_V(0));
	ProdVquad(1)(k) -= poids*grad_U(1);
	ProdVquad(2)(k) -= poids*grad_V(2);
      }
  }
#endif

  
  //! generic for matrix vector product for continuous elements and Helmholtz equation (mixed formulation)
  template<class T, class T0, class TypeEquation>
  void MltAddScalarHelmholtzH1(const T& alpha, const GlobalGenericMatrix<T0>& nat_mat,
			       const SeldonTranspose& trans, int lvl,
			       const FemMatrixFreeClass<T0, TypeEquation>& A, 
			       const Vector<T>& X, const T& beta, Vector<T>& Y, bool assemble)
  {
    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (beta == zero)
      Y.Zero();
    else if (beta != one)
      Mlt(beta, Y);

    const EllipticProblem<TypeEquation>& var = A.var;
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);

    Complexe m_iomega; var.GetMiomega(m_iomega);
    int Nvol = mesh_num.GetNbDof();
    int Nscal = var.offset_dof_unknown(1);
    
    TinyVector<Vector<T>, 1> Uloc, ProdUloc;
    TinyVector<Vector<T>, Dimension::dim_N> Vloc, ProdVquad;
    Vector<T> Uquad, gradU_quad, ProdUquad, gradUtilde_quad, gradUterce_quad;
    Vector<T> Utilde, Uterce, Utilde_quad, Uterce_quad;
    Vector<T> ProdUtilde, ProdUterce, ProdUtilde_quad, ProdUterce_quad;
    
    T m = nat_mat.GetCoefMass()*alpha;
    T s = nat_mat.GetCoefStiffness()*alpha;
    T sig = nat_mat.GetCoefDamping()*alpha;
    TinyVector<T, Dimension::dim_N> vec_v, grad_U, vec_u, vec_w;
    bool sym = A.GetSymmetrizationUse();
    bool stiff = false, mass = true;
    if (s != zero)
      stiff = true;
    
    if (m != zero)
      mass = true;

    if (sig != zero)
      mass = true;

    T coef_v = m_iomega;
    bool signU = true, signV = true;
    if (sym)
      coef_v = -coef_v;

    if (!sym)
      {
	if (trans.Trans())
	  signU = false;
	else
	  signV = false;
      }

    Real_wp coef_sym(1);
    if (sym)
      coef_sym = Real_wp(-1);
    
    int d = Dimension::dim_N;
    int Npml = mesh_num.GetNbDofPML();

    // loop over elements
    for (int i0 = 0; i0 < var.mesh.GetNbElt(); i0++)
      {
	//int i = list_level.GetElementNumber(i0);
	int i = i0;

	const ElementReference<Dimension, 1>& Fb = var.GetReferenceElementH1(i);

	int nb_points_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	bool affine = mesh.IsElementAffine(i);
	int ref = mesh.Element(i).GetReference();
	bool pml_elt = var.InsidePML(i);
	bool ortho = (var.OrthogonalElement(i) == 0);
        bool drude_model = var.ref_drude(ref).IsEnabled();
	//bool ortho_z = (var.OrthogonalElement(i) == 1);
	int i1 = i - mesh.GetNbElt() + var.GetNbEltPML();
	if (!pml_elt)
	  ortho = true;
	
	if (!ortho)
	  {
	    cout << "Non-orthogonal PML not implemented for the matrix-vector product" << endl;
	    abort();
	  }
	
	ProdUloc(0).Reallocate(nb_points_elt); ProdUloc(0).Zero();
	ProdUquad.Reallocate(nb_points_quad);
	ProdUquad.Zero();
	for (int n = 0; n < Dimension::dim_N; n++)
	  {
	    ProdVquad(n).Reallocate(nb_points_quad);
	    ProdVquad(n).Zero();
	  }

	if (pml_elt)
	  {
	    Utilde.Reallocate(nb_points_elt); Utilde.Zero();
	    ProdUtilde.Reallocate(nb_points_elt);
	    ProdUtilde_quad.Reallocate(nb_points_quad);
	    Utilde_quad.Reallocate(nb_points_quad);
	    ProdUtilde.Zero(); 
	    ProdUtilde_quad.Zero(); 

	    if (d == 3)
	      {
		Uterce.Reallocate(nb_points_elt); Uterce.Zero();
		Uterce_quad.Reallocate(nb_points_quad);
		ProdUterce.Reallocate(nb_points_elt);
		ProdUterce_quad.Reallocate(nb_points_quad);
		ProdUterce.Zero();
		ProdUterce_quad.Zero();
	      }
	  }
	
	// we retrieve values of u on the element
	int nb_dof_elt = mesh_num.Element(i).GetNbDof();
	int offset = Nscal + var.GetOffsetDofV(i);
	Uloc(0).Reallocate(nb_points_elt);
	var.GetLocalUnknownVector(X, i, Uloc);
	if (pml_elt)
	  {
	    for (int j = 0; j < nb_dof_elt; j++)
	      {
		int num_dof = mesh_num.Element(i).GetNumberDof(j);
		if (num_dof >= 0)
		  {
		    int num_pml = Nvol + mesh_num.GetDofPML(num_dof);
		    Utilde(j) = X(num_pml);
		    if (d == 3)
		      Uterce(j) = X(Npml+num_pml);
		  }
	      }		  
	  }
	
 	// and values of v
	for (int k = 0; k < Dimension::dim_N; k++)
	  Vloc(k).Reallocate(nb_points_quad);
	
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    ExtractVector(X, offset + Dimension::dim_N*j, vec_v);
	    CopyVector(vec_v, j, Vloc);
	  }
	
	// we compute u on quadrature points
	Uquad.Reallocate(nb_points_quad);
	Fb.ApplyChTranspose(Uloc(0), Uquad);
	if (pml_elt)
	  {
	    Fb.ApplyChTranspose(Utilde, Utilde_quad);
	    if (d == 3)
	      Fb.ApplyChTranspose(Uterce, Uterce_quad);
	  }

	// we add stiffness terms if present
	if (stiff)
	  {
	    gradU_quad.Reallocate(nb_points_quad*Dimension::dim_N);
	    
	    if (Fb.UseQuadraturePointsForRh())
	      Fb.ApplyRhQuadratureTranspose(Uquad, gradU_quad);
	    else
	      Fb.ApplyRhTranspose(Uloc(0), gradU_quad);

	    if (trans.Trans() && pml_elt)
	      {
		gradUtilde_quad.Reallocate(nb_points_quad*Dimension::dim_N);
		if (Fb.UseQuadraturePointsForRh())
		  Fb.ApplyRhQuadratureTranspose(Utilde_quad, gradUtilde_quad);
		else
		  Fb.ApplyRhTranspose(Utilde, gradUtilde_quad);

		if (Dimension::dim_N == 3)
		  {
		    gradUterce_quad.Reallocate(nb_points_quad*Dimension::dim_N);
		    if (Fb.UseQuadraturePointsForRh())
		      Fb.ApplyRhQuadratureTranspose(Uterce_quad, gradUterce_quad);
		    else
		      Fb.ApplyRhTranspose(Uterce, gradUterce_quad);
		  }
	      }
	    
	    if (affine)
	      {
		typename Dimension::MatrixN_N dfjm1 = var.Glob_DFjm1(i)(0);
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    T poids = Fb.WeightsND(k)*s;
		    
		    // part \nabla u \cdot psi
		    CopyVector(gradU_quad, k, vec_u);
		    MltTrans(dfjm1, vec_u, grad_U);
		    grad_U *= poids;
		    if (!signV)
		      grad_U = -grad_U;
		    
		    CopyVector(grad_U, k, ProdVquad);

		    // part v \cdot \nabla \phi
		    CopyVector(Vloc, k, vec_u);
		    Mlt(dfjm1, vec_u, vec_v);
		    vec_v *= poids;
		    if (!signU)
		      vec_v = -vec_v;

		    CopyVector(vec_v, k, gradU_quad);
		  }

		if (pml_elt && trans.Trans())
		  ApplyRhTransSplitElement(Fb, nb_points_quad, s, dfjm1, gradUtilde_quad,
					   gradUterce_quad, ProdVquad);
	      }
	    else
	      for (int k = 0; k < nb_points_quad; k++)
		{
		  
		  T poids = Fb.WeightsND(k)*s;
		  
		  // part \nabla u \cdot psi
		  CopyVector(gradU_quad, k, vec_u);
		  MltTrans(var.Glob_DFjm1(i)(k), vec_u, grad_U);
		  grad_U *= poids;
		  if (!signV)
		    grad_U = -grad_U;
		  
		  CopyVector(grad_U, k, ProdVquad);
		  
		  // part v \cdot \nabla \phi
		  CopyVector(Vloc, k, vec_u);
		  Mlt(var.Glob_DFjm1(i)(k), vec_u, vec_v);
		  vec_v *= poids;
		  if (!signU)
		    vec_v = -vec_v;
		  
		  CopyVector(vec_v, k, gradU_quad);
		}

	    if (pml_elt && (trans.NoTrans()))
	      {
		ApplyRhSplitElement(Fb, gradU_quad, ProdUquad, ProdUtilde_quad, ProdUterce_quad,
				    ProdUloc(0), ProdUtilde, ProdUterce);
	      }
	    else
	      {
		if (Fb.UseQuadraturePointsForRh())
		  Fb.ApplyRhQuadrature(gradU_quad, ProdUquad);
		else
		  Fb.ApplyRh(gradU_quad, ProdUloc(0));
	      }
            
            if (drude_model)
              {
                if (Fb.UseQuadraturePointsForRh())
                  {
                    cout << "Case not implemented" << endl;
                    abort();
                  }

		int offset_P, offset_Q;
                offset_Q = var.GetNbDof() - var.GetNbVectorialDofDrude();
                offset_P = offset_Q - var.GetNbVectorialDofDrude();
                offset_P += var.OffsetDofDrudeV(i);
                offset_Q += var.OffsetDofDrudeV(i);
		
                VectReal_wp gamma = var.ref_drude(ref).gamma;
                VectReal_wp coef_omega2 = var.ref_drude(ref).omega_02;
                VectReal_wp coef_eps_inf = var.ref_drude(ref).eps_omega_p2;

                int nPole = gamma.GetM();
                VectReal_wp coef_Q(nPole);
                VectReal_wp coef_P(nPole);
                if (!var.use_symm_drude)
                  {
                    // pas de symmetrisation pour traiter omega_0 = 0
                    coef_Q.Fill(coef_sym); coef_P.Fill(Real_wp(1));
                  }
                else
                  {
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        coef_Q(kp) = coef_sym / coef_eps_inf(kp);
                        coef_P(kp) = coef_omega2(kp) / coef_eps_inf(kp);
                      }
                  }
                
                if (var.ref_drude(ref).IsModeTE())
                  {
                    coef_P *= coef_sym;
                    coef_Q *= coef_sym;
                  }

                VectReal_wp coef_u_q(nPole), coef_q_u;
                coef_u_q.Fill(Real_wp(1)); coef_q_u = -coef_Q*coef_eps_inf;
                if (var.ref_drude(ref).IsModeTE())
                  coef_u_q.Fill(coef_sym);
                
                VectReal_wp coef_q_q; coef_q_q = coef_Q*gamma;
                VectReal_wp coef_q_p; coef_q_p = coef_Q*coef_omega2;
                VectReal_wp coef_p_q; coef_p_q = -coef_P;
                if (!trans.NoTrans())
                  {
                    SwapPointer(coef_u_q, coef_q_u);
                    SwapPointer(coef_p_q, coef_q_p);
                  }

                for (int j = 0; j < nb_dof_elt; j++)
                  {
                    Real_wp poids;
                    if (affine)
                      poids = Fb.WeightsND(j) * var.Glob_jacobian(i)(0);
                    else
                      poids = var.Glob_jacobian(i)(j);
                    
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        if (var.ref_drude(ref).IsModeTM())
                          {              
                            int num_dof = mesh_num.Element(i).GetNumberDof(j);
                            int num_dof_P = offset_P + j*nPole + kp;
                            int num_dof_Q = offset_Q + j*nPole + kp;
                            
                            T q_val = X(num_dof_Q), p_val = X(num_dof_P), u_val = X(num_dof);
                            Y(num_dof) += s*poids*coef_u_q(kp)*q_val;
                            Y(num_dof_P) += s*poids*coef_p_q(kp)*q_val;
                            Y(num_dof_Q) += s*poids*(coef_q_p(kp)*p_val + coef_q_q(kp)*q_val + coef_q_u(kp)* u_val);
                          }
                        else
                          {
                            int num_dofHx = offset + 2*j, num_dofHy = num_dofHx + 1;
                            int num_dofPx = offset_P + 2*j*nPole+kp, num_dofPy = num_dofPx + nPole;
                            int num_dofQx = offset_Q + 2*j*nPole+kp, num_dofQy = num_dofQx + nPole;
                            T qx_val = X(num_dofQx), qy_val = X(num_dofQy);
                            T px_val = X(num_dofPx), py_val = X(num_dofPy);
                            T hx_val = X(num_dofHx), hy_val = X(num_dofHy);
                            
                            Y(num_dofHx) += s*poids*coef_u_q(kp)*qx_val;
                            Y(num_dofHy) += s*poids*coef_u_q(kp)*qy_val;
                            
                            Y(num_dofPx) += s*poids*coef_p_q(kp)*qx_val;
                            Y(num_dofPy) += s*poids*coef_p_q(kp)*qy_val;
                            
                            Y(num_dofQx) += s*poids*(coef_q_p(kp)*px_val + coef_q_q(kp)*qx_val + coef_q_u(kp)*hx_val);
                            Y(num_dofQy) += s*poids*(coef_q_p(kp)*py_val + coef_q_q(kp)*qy_val + coef_q_u(kp)*hy_val);			
                          }
                      }
                  }            
              }
	  }
        
	// mass terms
	if (mass)
	  {
	    T poids(zero), poids1(zero), poids2(zero), poids3(zero); Real_wp jacobian;
	    bool variable = var.UseNumericalIntegration(i);
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		if (variable)
                  poids = var.Glob_matMass_Dh(i)(k)*m + var.Glob_matMass_DhSigma(i)(k)*sig;
                else
                  poids = Fb.WeightsND(k)*(var.Glob_matMass_Dh(i)(0)*m + var.Glob_matMass_DhSigma(i)(0)*sig);

		ProdUquad(k) += poids*Uquad(k);

		if (affine)
		  jacobian = var.Glob_jacobian(i)(0)*Fb.WeightsND(k);
		else
		  jacobian = var.Glob_jacobian(i)(k);

		CopyVector(Vloc, k, vec_w);
		var.ref_invMu(ref).MltMatrix(var, i, k, vec_w, vec_v);
		vec_v *= coef_v*m*jacobian;

		CopyVector(ProdVquad, k, vec_u);
		vec_v += vec_u;

		if (pml_elt)
		  {
		    vec_u = var.GetTauPML(i1, k);
		    vec_u *= jacobian*sig;
		    var.ref_invMu(ref).MltMatrix(var, i, k, vec_u, grad_U);
		    vec_v += grad_U*vec_w;
		    
		    if (d == 2)
		      {
			ProdUtilde_quad(k) += poids*Utilde_quad(k);
			
			poids = var.Glob_matMass_DhSigmaDiff(i1)(k)*sig;		    
			ProdUquad(k) += poids*Utilde_quad(k);
			ProdUtilde_quad(k) += poids*Uquad(k);
		      }
		    else
		      {
			ProdUtilde_quad(k) += poids*Utilde_quad(k);
			ProdUterce_quad(k) += poids*Uterce_quad(k);

			poids1 = var.Glob_matMass_DhSigmaDiff(i1)(k)*sig;
			poids2 = var.Glob_matMass_DhSigmaDiff2(i1)(k)*sig;
			poids3 = var.Glob_matMass_DhSigmaDiff3(i1)(k)*sig;

			ProdUtilde_quad(k) += poids3*Utilde_quad(k);
			ProdUterce_quad(k) -= poids3*Uterce_quad(k);
			
			ProdUquad(k) += poids1*Utilde_quad(k);			
			ProdUtilde_quad(k) += poids1*Uquad(k);

			ProdUquad(k) += poids2*Uterce_quad(k);
			ProdUterce_quad(k) += poids2*Uquad(k);
			
			if (trans.Trans())
			  {
			    ProdUtilde_quad(k) -= poids3*Uquad(k);
			    ProdUterce_quad(k) += poids3*Uquad(k);
			    ProdUterce_quad(k) += poids1*Utilde_quad(k);
			    ProdUtilde_quad(k) += poids2*Uterce_quad(k);
			  }
			else
			  {
			    ProdUquad(k) -= poids3*Utilde_quad(k);
			    ProdUquad(k) += poids3*Uterce_quad(k);
			    ProdUtilde_quad(k) += poids1*Uterce_quad(k);
			    ProdUterce_quad(k) += poids2*Utilde_quad(k);
			  }
		      }
		  }

		CopyVector(vec_v, k, ProdVquad);
	      }

            if (drude_model)
              {
		int offset_P, offset_Q;
                offset_Q = var.GetNbDof() - var.GetNbVectorialDofDrude();
                offset_P = offset_Q - var.GetNbVectorialDofDrude();
                offset_P += var.OffsetDofDrudeV(i);
                offset_Q += var.OffsetDofDrudeV(i);
                
                int nPole = var.ref_drude(ref).gamma.GetM();
                VectReal_wp coef_omega2 = var.ref_drude(ref).omega_02;
                VectReal_wp coef_eps_inf = var.ref_drude(ref).eps_omega_p2;
                VectReal_wp coef_u_q(nPole), coef_q_u(nPole); 
                coef_u_q.Zero();
                
                VectReal_wp coef_Q(nPole), coef_P(nPole);
                if (!var.use_symm_drude)
                  {
                    // pas de symmetrisation pour traiter omega_0 = 0
                    coef_Q.Fill(coef_sym); coef_P.Fill(Real_wp(1));
                    coef_q_u = -coef_sym*var.ref_drude(ref).eps_sigma;
                  }
                else
                  {
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        coef_Q(kp) = coef_sym / coef_eps_inf(kp);
                        coef_P(kp) = coef_omega2(kp) / coef_eps_inf(kp);
                        coef_q_u(kp) = -coef_sym*var.ref_drude(ref).eps_sigma(kp) / coef_eps_inf(kp);
                      }
                  }
                
                if (var.ref_drude(ref).IsModeTE())
                  {
                    coef_P *= coef_sym;
                    coef_Q *= coef_sym;
                    coef_q_u *= coef_sym;
                  }

                if (!trans.NoTrans())
                  SwapPointer(coef_u_q, coef_q_u);
                
                for (int j = 0; j < nb_dof_elt; j++)
                  {
                    Real_wp poids;
                    if (affine)
                      poids = Fb.WeightsND(j) * var.Glob_jacobian(i)(0);
                    else
                      poids = var.Glob_jacobian(i)(j);

                    T coef_m = m*m_iomega*poids;
                    for (int kp = 0; kp < nPole; kp++)
                      {
                        if (var.ref_drude(ref).IsModeTM())
                          {
                            int num_dof = mesh_num.Element(i).GetNumberDof(j);
                            int num_dof_P = offset_P + j*nPole + kp;
                            int num_dof_Q = offset_Q + j*nPole + kp;
                            T q_val = X(num_dof_Q), p_val = X(num_dof_P), u_val = X(num_dof);
                            Y(num_dof) += coef_m*coef_u_q(kp)*q_val;
                            Y(num_dof_P) += coef_m*coef_P(kp)*p_val;
                            Y(num_dof_Q) += coef_m*(coef_Q(kp)*q_val + coef_q_u(kp)*u_val);
                          }
                        else
                          {
                            int num_dofHx = offset+2*j, num_dofHy = num_dofHx+1;
                            int num_dofPx = offset_P + 2*j*nPole + kp, num_dofPy = num_dofPx + nPole;
                            int num_dofQx = offset_Q + 2*j*nPole + kp, num_dofQy = num_dofQx + nPole;
                            T qx_val = X(num_dofQx), qy_val = X(num_dofQy);
                            T px_val = X(num_dofPx), py_val = X(num_dofPy);
                            T hx_val = X(num_dofHx), hy_val = X(num_dofHy);
                            
                            Y(num_dofHx) += coef_m*coef_u_q(kp)*qx_val;
                            Y(num_dofHy) += coef_m*coef_u_q(kp)*qy_val;
                            
                            Y(num_dofPx) += coef_m*coef_P(kp)*px_val;
                            Y(num_dofPy) += coef_m*coef_P(kp)*py_val;
                            
                            Y(num_dofQx) += coef_m*(coef_Q(kp)*qx_val + coef_q_u(kp)*hx_val);
                            Y(num_dofQy) += coef_m*(coef_Q(kp)*qy_val + coef_q_u(kp)*hy_val);
                          }
                      }
                  }
              }
	  }
        
	if (Fb.UseQuadraturePointsForRh())
	  Fb.ApplyCh(ProdUquad, ProdUloc(0));
	else
	  {
	    if (mass)
	      {
		Fb.ApplyCh(ProdUquad, Uloc(0));
		ProdUloc(0) += Uloc(0);
	      }
	  }
	
	if (pml_elt)
	  {
	    if (Fb.UseQuadraturePointsForRh())
	      {
		Fb.ApplyCh(ProdUtilde_quad, ProdUtilde);
		if (d == 3)
		  Fb.ApplyCh(ProdUterce_quad, ProdUterce);
	      }
	    else
	      {
		if (mass)
		  {
		    Fb.ApplyCh(ProdUtilde_quad, Uloc(0));
		    ProdUtilde += Uloc(0);
		    if (d == 3)
		      {
			Fb.ApplyCh(ProdUterce_quad, Uloc(0));
			ProdUterce += Uloc(0);
		      }
		  }		
	      }
	  }
	
	// contributions are added to Y
	var.AddLocalUnknownVector(one, ProdUloc, i, Y);
	if (pml_elt)
	  for (int j = 0; j < nb_dof_elt; j++)
	    {
	      int num_dof = mesh_num.Element(i).GetNumberDof(j);
	      if (num_dof >= 0)
		{
		  int num_pml = Nvol + mesh_num.GetDofPML(num_dof);
		  Y(num_pml) += ProdUtilde(j);
		  if (d == 3)
		    Y(Npml + num_pml) += ProdUterce(j);
		}
	    }
	
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    CopyVector(ProdVquad, j, vec_v);
	    AddVector(vec_v, offset + Dimension::dim_N*j, Y);
	  }
      }
    
    A.AddExtraBoundaryTerms(alpha, X, Y);

    if (!A.DirichletDofIgnored())
      var.ImposeNullDirichletCondition(Y);    
  }

  
  template<class T, class TypeEquation>
  void MltAddDrudeHelmholtzH1(const GlobalGenericMatrix<T>& nat_mat, const SeldonTranspose& trans,
                              int lvl, const FemMatrixFreeClass<T, TypeEquation>& A,
                              const Vector<T>& X, Vector<T>& Y)
  {
    const EllipticProblem<TypeEquation>& var = A.var;
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);

    Complexe m_iomega; var.GetMiomega(m_iomega);
    Complexe m_omega2 = m_iomega*m_iomega;
    
    T m = nat_mat.GetCoefMass();
    T s = nat_mat.GetCoefStiffness();
    T sig = nat_mat.GetCoefDamping();

    if (var.linearize_drude)
      for (int i0 = 0; i0 < var.mesh.GetNbElt(); i0++)
        {
          //int i = list_level.GetElementNumber(i0);
          int i = i0;
          int ref = mesh.Element(i).GetReference();
          bool drude_model = var.ref_drude(ref).IsEnabled();
          Complexe invRho0 = Real_wp(1) / var.ref_rho(ref).GetConstant();
          if (drude_model && var.ref_drude(ref).IsModeTM())
            {
              const ElementReference<Dimension, 1>& Fb = var.GetReferenceElementH1(i);
              
              int nb_points_quad = Fb.GetNbPointsQuadratureInside();
              bool affine = mesh.IsElementAffine(i);
              
              int nPole = var.ref_drude(ref).gamma.GetM(); T val;
              Vector<T> coef_E(nPole), coef_P(nPole), coef_diag(nPole);
              for (int k = 0; k < nPole; k++)
                {
                  val = -s*var.ref_drude(ref).eps_omega_p2(k) - sig*m_iomega*var.ref_drude(ref).eps_sigma(k);
                  if (trans.NoTrans())
                    {
                      coef_E(k) = nat_mat.GetCoefMass()*m_omega2;
                      coef_P(k) = val;
                    }
                  else
                    {
                      coef_P(k) = nat_mat.GetCoefMass()*m_omega2;
                      coef_E(k) = val;
                    }
                  
                  coef_diag(k) = m*m_omega2 + sig*m_iomega*var.ref_drude(ref).gamma(k) + s*var.ref_drude(ref).omega_02(k);
                }
              
              int offset_P = var.GetNbDof() - var.nb_dof_drude_vec;
              offset_P += var.OffsetDofDrudeV(i);
              Real_wp poids(0); T coefE, coefP;
              for (int j = 0; j < nb_points_quad; j++)
                {
                  int num_dof = mesh_num.Element(i).GetNumberDof(j);
                  if (affine)
                    poids = Fb.WeightsND(j) * var.Glob_jacobian(i)(0);
                  else
                    poids = var.Glob_jacobian(i)(j);                
                  
                  for (int k = 0; k < nPole; k++)
                    {
                      coefE = poids*coef_E(k); coefP = poids*coef_P(k);
                      if (var.InsidePML(i))
                        {
                          if (trans.NoTrans())
                            coefE = var.Glob_matMass_Dh(i)(j)*m*invRho0;
                          else
                            coefP = var.Glob_matMass_Dh(i)(j)*m*invRho0;
                        }
                      
                      int dofP = offset_P + j*nPole + k;
                      Y(num_dof) += coefE*X(dofP);
                      Y(dofP) += poids*coef_diag(k)*X(dofP) + coefP*X(num_dof);
                    }
                }
            }
        }
  }
  

  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, HelmholtzEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    if (this->var.FirstOrderFormulation())
      MltAddScalarHelmholtzH1(Complex_wp(1, 0), nat_mat, trans, lvl, *this,
			      X, Complex_wp(1, 0), Y, false);
    else
      {
        MltAddScalarH1(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
                       X, Complex_wp(1, 0), Y, false);
        
        MltAddDrudeHelmholtzH1(nat_mat, trans, lvl, *this, X, Y);
      }
  }
  
  
  //! generic for matrix vector product for discontinuous elements (HDG) and Helmholtz equation
  template<class T, class T0, class TypeEquation>
  void MltAddScalarHelmholtz(const T& alpha, const GlobalGenericMatrix<T0>& nat_mat,
			     const SeldonTranspose& trans, int lvl,
			     const FemMatrixFreeClass<T0, TypeEquation>& A, 
			     const Vector<T>& X, const T& beta, Vector<T>& Y, bool assemble)
  {
    T zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (beta == zero)
      Y.Zero();
    else if (beta != one)
      Mlt(beta, Y);

    const EllipticProblem<TypeEquation>& var = A.var;
    if (var.FormulationDG() != ElementReference_Base::HDG)
      {
	cout << "Only implemented for HDG formulation" << endl;
	abort();
      }

    if (trans.ConjTrans())
      {
	cout << "Impossible" << endl;
	abort();
      }
    
    typedef typename TypeEquation::Complexe Complexe;
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    typedef typename Dimension::DimensionBoundary DimensionB;

    int nb_dof_L = var.GetOffsetDofUnknown(1);
    int Nvol = var.GetNbMainUnknownDof();

    TinyVector<Vector<T>, 1> Lloc;
    Vector<T> ProdUloc, ProdUquad, feval, Uquad;
    TinyVector<Vector<T>, Dimension::dim_N> Vloc, Vquad, ProdVloc, ProdVquad, face_V, fevalV;
    Vector<T> Lambda_quad, face_U, Uloc, contrib, dU, fevalL;

    // penalization and coefficients (depending on signs and formulation used)
    Complexe m_iomega; Real_wp coef_ipp;
    T coef_u, coef_v, coef_u_damped;    

    T s = nat_mat.GetCoefStiffness();
    TinyVector<T, Dimension::dim_N> vec_v, grad_U, vec_u;
    bool sym = A.GetSymmetrizationUse();
    
    // loop over elements
    for (int i0 = 0; i0 < var.mesh.GetNbElt(); i0++)
      {
	//int i = list_level.GetElementNumber(i0);
	int i = i0;
	int offset_u = nb_dof_L + var.GetOffsetDofV(i);
	
	const ElementReference<Dimension, 1>& Fb = var.GetReferenceElementH1(i);

        int ref_domain = mesh.Element(i).GetReference();
	int nb_points_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	ProdUloc.Reallocate(nb_points_elt); ProdUloc.Zero();
	ProdUquad.Reallocate(nb_points_quad);
	ProdUquad.Zero();
	for (int n = 0; n < Dimension::dim_N; n++)
	  {
	    ProdVloc(n).Reallocate(nb_points_elt);
	    ProdVquad(n).Reallocate(nb_points_quad);
	    ProdVloc(n).Zero(); ProdVquad(n).Zero();
	  }

        // ipp coefficients
        var.GetCoefficientIPP(var, m_iomega, coef_ipp, coef_u, coef_v,
                              coef_u_damped);
        
        if (sym)
          {
            coef_ipp = 1.0;
            coef_v = -coef_v;
          }
        
        T coef_v_orig(coef_v);        
        coef_u *= nat_mat.GetCoefMass();
        coef_v *= nat_mat.GetCoefMass();
        coef_u_damped *= nat_mat.GetCoefDamping();
        
        Real_wp coef_trans = 1.0, coef_ipp_orig(coef_ipp);
        if (trans.Trans())
          {
            coef_trans = coef_ipp;
            coef_ipp = 1.0;
          }
        
        T c_ipp = coef_ipp*nat_mat.GetCoefStiffness();
        T c_ipp_orig = coef_ipp_orig*nat_mat.GetCoefStiffness();
        T c_trans = coef_trans*nat_mat.GetCoefStiffness();
        
	// we retrieve values of lambda on the element
	mesh_num.number_map.GetLocalUnknownVector(mesh_num, X, i, Lloc);

	// u and v as well
	Uloc.Reallocate(nb_points_elt);
	Uquad.Reallocate(nb_points_quad);
	for (int n = 0; n < Dimension::dim_N; n++)
	  {
	    Vloc(n).Reallocate(nb_points_elt);
	    Vquad(n).Reallocate(nb_points_quad);
	  }
	
	for (int j = 0; j < nb_points_elt; j++)
	  {
	    Uloc(j) = X(offset_u + j);
	    for (int n = 0; n < Dimension::dim_N; n++)
	      Vloc(n)(j) = X(offset_u + Nvol*(n+1) + j);
	  }

	bool affine = mesh.IsElementAffine(i);
	bool compute_u_quad = false;
	bool compute_v_quad = false;
	if (Fb.UseQuadraturePointsForRh() || Fb.UseQuadraturePointsForSh() || (!affine) )
	  {
	    compute_u_quad = true;
	    compute_v_quad = true;
	  }

	if (var.ref_mu(ref_domain).IsVarying())
	  compute_v_quad = true;

	if (var.ref_rho(ref_domain).IsVarying() || var.ref_sigma(ref_domain).IsVarying())
	  compute_u_quad = true;
	
	if (compute_u_quad)
	  Fb.ApplyChTranspose(Uloc, Uquad);

	if (compute_v_quad)
	  for (int n = 0; n < Dimension::dim_N; n++)
	    Fb.ApplyChTranspose(Vloc(n), Vquad(n));
	
	// loop over boundaries of the element
	int offset = 0;
	//int offset_face = mesh_num.OffsetQuadElementNumber(i);
	for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
	    int num_face = mesh.Element(i).numBoundary(num_loc);
            
            T coef_tau;
            var.GetPenalizationTauIPP(var, num_face, coef_tau);
            coef_tau *= nat_mat.GetCoefStiffness();
            
	    int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    const ElementReference<DimensionB, 1>& Fb_s = var.GetSurfaceFiniteElementH1(num_face);
	    
	    contrib.Reallocate(Fb_s.GetNbDof());
	    fevalL.Reallocate(nb_points_face);
	    Lambda_quad.Reallocate(nb_points_face);
	    
	    // values of Lambda on quadrature points are computed
	    for (int j = 0; j < contrib.GetM(); j++)
	      contrib(j) = Lloc(0)(offset + j);
	    
	    Fb_s.ApplyChTranspose(contrib, Lambda_quad);

	    Real_wp dsj; typename Dimension::R_N normale;
	    
	    // values of u and v on quadrature points are retrieved
	    face_U.Reallocate(nb_points_face);
	    for (int k = 0; k < Dimension::dim_N; k++)
	      face_V(k).Reallocate(nb_points_face);

	    if (Fb.UseQuadraturePointsForSh())
	      {
		Fb.ApplyShQuadratureTranspose(num_loc, Uquad, face_U, rf);
		for (int k = 0; k < Dimension::dim_N; k++)
		  Fb.ApplyShQuadratureTranspose(num_loc, Vquad(k), face_V(k), rf);
	      }
	    else
	      {
		Fb.ApplyShTranspose(num_loc, Uloc, face_U, rf);
		for (int k = 0; k < Dimension::dim_N; k++)
		  Fb.ApplyShTranspose(num_loc, Vloc(k), face_V(k), rf);
	      }
	    
	    bool new_face = var.IsNewFace(i)(num_loc);
	    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
	    
	    int num_elem2 = var.mesh.Boundary(num_face).numElement(0);
	    if ((num_elem2 == i)&&(var.mesh.Boundary(num_face).GetNbElements()==2))
	      num_elem2 = var.mesh.Boundary(num_face).numElement(1);
	    
	    int rot1 = mesh.Element(i).GetOrientationBoundary(num_loc), rot2 = 0;
	    int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0;
	    if (num_elem2 != i)
	      {
		int num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(num_face);
		if (num_pos2_face < 0)
		  {
		    int nf2 = mesh_num.GetPeriodicBoundary(num_face);
		    num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(nf2);
		  }
                
		rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
		rot = mesh_num.GetRotationFace(rot1, rot2, nv);
	      }
	    
	    int ref_boundary = mesh.Boundary(num_face).GetReference();
	    bool dirichlet_boundary = false;
	    if ((ref_boundary != 0) &&
		(mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_DIRICHLET))
	      dirichlet_boundary = true;
	    
	    feval.Reallocate(nb_points_face);
	    for (int d = 0; d < Dimension::dim_N; d++)
	      fevalV(d).Reallocate(nb_points_face);

	    Real_wp poids; T diff_u, v_dot_n;
	    for (int k = 0; k < nb_points_face; k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  {
		    normale = var.Glob_normale(num_face)(k);
		    dsj = var.Glob_dsj(num_face)(k);
		  }
		else
		  {
		    normale = var.Glob_normale(num_face)(krot);
		    dsj = var.Glob_dsj(num_face)(krot);
		    normale = -normale;
		  }

		poids = dsj*Fb_s.WeightsND(k);
		SetComplexZero(diff_u);

		CopyVector(face_V, k, vec_v);
		v_dot_n = DotProd(vec_v, normale);
		  
		feval(k) = poids* (coef_tau * (face_U(k)-coef_trans*Lambda_quad(k)) - c_trans * v_dot_n);
		diff_u = c_ipp*face_U(k);
		
		if (!dirichlet_boundary)
		  {
		    diff_u -= c_ipp_orig*Lambda_quad(k);
		    fevalL(k) = c_ipp_orig*poids*v_dot_n;
		  }
		else
		  SetComplexZero(fevalL(k));
		
		vec_v = -poids*normale*diff_u;
		CopyVector(vec_v, k, fevalV);

		T coef_abc(coef_tau*coef_ipp_orig);
		var.ModifyCoefficientBC_HDG(var, var, mesh, ref_boundary,
					    0, s, coef_ipp_orig, coef_v_orig,
					    ref_domain, i, num_loc, coef_abc);
		
		fevalL(k) += poids*(coef_abc*Lambda_quad(k) - coef_tau*coef_ipp*face_U(k));
	      }
	    
	    // contributions for Lambda is put in Lloc
	    Fb_s.ApplyCh(fevalL, contrib);
	    for (int j = 0; j < contrib.GetM(); j++)
	      Lloc(0)(offset + j) = contrib(j);
	    
	    // then contributions for u and v are treated
	    if (Fb.UseQuadraturePointsForSh())
	      {
		Fb.ApplyShQuadrature(Real_wp(1), num_loc, feval, ProdUquad, rf);
		for (int d = 0; d < Dimension::dim_N; d++)
		  Fb.ApplyShQuadrature(Real_wp(1), num_loc, fevalV(d), ProdVquad(d), rf);
	      }
	    else
	      {
		Fb.ApplySh(Real_wp(1), num_loc, feval, ProdUloc, rf);
		for (int d = 0; d < Dimension::dim_N; d++)
		  Fb.ApplySh(Real_wp(1), num_loc, fevalV(d), ProdVloc(d), rf);
	      }
	    
	    offset += Fb_s.GetNbDof();
	    //offset_face += nb_points_face;
	  }

	// values for Lambda are added to Y
	mesh_num.number_map.AddLocalUnknownVector(mesh_num, alpha, Lloc, i, Y);

	// now we treat volume integrals
	// part for ProdVloc
	if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
	  nb_points_quad = nb_points_elt;
	
	dU.Reallocate(nb_points_quad*Dimension::dim_N);
	contrib.Reallocate(nb_points_elt);
	feval.Reallocate(nb_points_quad);
	for (int d = 0; d < Dimension::dim_N; d++)
	  fevalV(d).Reallocate(nb_points_quad);
	
	// gradient of u on reference element
	if (Fb.UseQuadraturePointsForSh())
	  Fb.ApplyRhQuadratureTranspose(Uquad, dU);
	else
	  {
	    if (Fb.UseQuadraturePointsForRh())
	      Fb.ApplyRhQuadratureTranspose(Uquad, dU);
	    else
	      {
		if (affine)
		  Fb.ApplyConstantRhTranspose(Uloc, dU);
		else
		  Fb.ApplyRhTranspose(Uloc, dU);
	      }
	  }

	// then application of DF
	if (affine)
	  {
	    const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N>&
	      dfjm1 = var.Glob_DFjm1(i)(0);
	    
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		CopyVector(dU, k, vec_v); Mlt(c_ipp, vec_v);
		MltTrans(dfjm1, vec_v, grad_U);
		CopyVector(grad_U, k, fevalV);
	      }
	    
	    if (Fb.UseQuadraturePointsForRh())
	      for (int k = 0; k < nb_points_quad; k++)
		MltVector(Fb.WeightsND(k), k, fevalV);
	  }
	else
	  for (int k = 0; k < nb_points_quad; k++)
	    {
	      CopyVector(dU, k, vec_v);
	      MltTrans(var.Glob_DFjm1(i)(k), vec_v, grad_U);	  
	      grad_U *= c_ipp*Fb.WeightsND(k);
	      CopyVector(grad_U, k, fevalV);
	    }

	// we add mass part
	if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
	  {
	    if (var.ref_mu(ref_domain).IsVarying())
	      {
		for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
		  {
		    CopyVector(Vquad, k, vec_u);
		    var.ref_invMu(ref_domain).MltMatrix(var, i, k, vec_u, vec_v);
		    vec_v *= coef_v*var.Glob_jacobian(i)(0)*Fb.WeightsND(k);
		    CopyVector(vec_v, k, ProdVquad);
		  }

		for (int d = 0; d < Dimension::dim_N; d++)
		  {
		    Fb.ApplyCh(ProdVquad(d), contrib);
		    fevalV(d) += contrib;
		  }
	      }
	    else
	      {
		for (int d = 0; d < Dimension::dim_N; d++)
		  face_V(d).Reallocate(nb_points_elt);
		
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    CopyVector(Vloc, k, vec_u);
		    var.ref_invMu(ref_domain).MltMatrix(var, i, 0, vec_u, vec_v);
		    vec_v *= coef_v*var.Glob_jacobian(i)(0);
		    CopyVector(vec_v, k, face_V);
		  }

		for (int d = 0; d < Dimension::dim_N; d++)
		  {
		    Fb.MltMassMatrix(face_V(d));
		    fevalV(d) += face_V(d);
		  }
	      }
	  }
	else
	  {
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		CopyVector(Vquad, k, vec_u);
		var.ref_invMu(ref_domain).MltMatrix(var, i, k, vec_u, vec_v);
		vec_v *= coef_v*var.Glob_jacobian(i)(k);
		CopyVector(fevalV, k, vec_u);
		vec_v += vec_u;
		CopyVector(vec_v, k, fevalV);
	      }
	  }
		
	// integration against basis functions
	if (Fb.UseQuadraturePointsForSh())
	  for (int d = 0; d < Dimension::dim_N; d++)
	    {
	      ProdVquad(d) += fevalV(d);
	      Fb.ApplyCh(ProdVquad(d), ProdVloc(d));
	    }
	else
	  for (int d = 0; d < Dimension::dim_N; d++)
	    {
	      if ( affine && !Fb.UseQuadraturePointsForRh())
		ProdVloc(d) += fevalV(d);
	      else
		{
		  Fb.ApplyCh(fevalV(d), contrib);
		  ProdVloc(d) += contrib;
		}
	    }

	// part for ProdUloc
	
	// application of transformation DF_i^{-1}	
	if (affine)
	  {
	    const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N>&
	      dfjm1 = var.Glob_DFjm1(i)(0);
	    
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		if (Fb.UseQuadraturePointsForRh())
		  {
		    CopyVector(Vquad, k, vec_u);		
		    Mlt(dfjm1, vec_u, vec_v);
		    vec_v *= c_trans*Fb.WeightsND(k);
		  }
		else
		  {
		    CopyVector(Vloc, k, vec_u);
		    Mlt(dfjm1, vec_u, vec_v);
		    vec_v *= c_trans;   
		  }
		
		CopyVector(vec_v, k, dU);
	      }
	  }
	else
	  {
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		CopyVector(Vquad, k, vec_u);
		
		Mlt(var.Glob_DFjm1(i)(k), vec_u, vec_v);
		vec_v *= c_trans*Fb.WeightsND(k);
		CopyVector(vec_v, k, dU);
	      }
	  }

	// we add mass part
	Complexe rho, sigma; T vloc;
	if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
	  {
	    if (var.ref_rho(ref_domain).IsVarying()
		|| var.ref_sigma(ref_domain).IsVarying())
	      {
		for (int k = 0; k < Fb.GetNbPointsQuadratureInside(); k++)
		  {
		    rho = var.ref_rho(ref_domain).GetCoefficient(var, i, k);
		    sigma = var.ref_sigma(ref_domain).GetCoefficient(var, i, k);
		    vloc = Uquad(k) * (rho*coef_u + sigma*coef_u_damped)
		      * var.Glob_jacobian(i)(0)*Fb.WeightsND(k);
		    
		    feval(k) = vloc;
		  }

		Fb.ApplyCh(feval, contrib);
		ProdUloc += contrib;
	      }
	    else
	      {
		rho = var.ref_rho(ref_domain).GetConstant();
		sigma = var.ref_sigma(ref_domain).GetConstant();
		vloc = (rho*coef_u + sigma*coef_u_damped) * var.Glob_jacobian(i)(0);
		contrib = vloc * Uloc; 
		
		Fb.MltMassMatrix(contrib);
		ProdUloc += contrib;
	      }
	  }
	else
	  {
	    for (int k = 0; k < nb_points_quad; k++)
	      {
		rho = var.ref_rho(ref_domain).GetCoefficient(var, i, k);
		sigma = var.ref_sigma(ref_domain).GetCoefficient(var, i, k);
		vloc = Uquad(k) * (rho*coef_u + sigma*coef_u_damped) * var.Glob_jacobian(i)(k);
		feval(k) = vloc;
	      }
	  }

	// integration against grad phi
	if (Fb.UseQuadraturePointsForRh())
	  {
	    Fb.ApplyRhQuadrature(dU, Uquad);	    
	    if (Fb.UseQuadraturePointsForSh())
	      {
		ProdUquad += Uquad + feval;
		Fb.ApplyCh(ProdUquad, ProdUloc);
	      }
	    else
	      {
		Uquad += feval;
		Fb.ApplyCh(Uquad, Uloc);
		ProdUloc += Uloc;
	      }
	  }
	else
	  {	
	    if (affine)
	      Fb.ApplyConstantRh(dU, Uloc);
	    else
	      {
		Fb.ApplyRh(dU, Uloc);
		Fb.ApplyCh(feval, contrib);
		Uloc += contrib;
	      }
	    
	    ProdUloc += Uloc;
	  }

	// contributions are added to ProdU and ProdV
	for (int j = 0; j < ProdUloc.GetM(); j++)
	  Y(offset_u + j) += alpha*ProdUloc(j);

	offset_u += Nvol;
	for (int d = 0; d < Dimension::dim_N; d++)
	  {
	    for (int j = 0; j < ProdVloc(d).GetM(); j++)
	      Y(offset_u + j) += alpha*ProdVloc(d)(j);

	    offset_u += Nvol;
	  }
      }
  }
  
  
  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, HelmholtzEquationDG<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "not possible" << endl;
    abort();
  }


  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, HelmholtzEquationDG<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddScalarHelmholtz(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
			  X, Complex_wp(1, 0), Y, false);
  }


  /*******************
   * LaplaceEquation *
   *******************/

  
  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    if (this->var.FirstOrderFormulation())
      MltAddScalarHelmholtzH1(Real_wp(1), nat_mat, trans, lvl, *this,
			      X, Real_wp(1), Y, false);
    else
      MltAddScalarH1(Real_wp(1), nat_mat, trans, lvl, *this, 
		     X, Real_wp(1), Y, false);
  }
  

  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    if (this->var.FirstOrderFormulation())
      MltAddScalarHelmholtzH1(Complex_wp(1), nat_mat, trans, lvl, *this,
			      X, Complex_wp(1), Y, false);
    else
      MltAddScalarH1(Complex_wp(1), nat_mat, trans, lvl, *this, 
		     X, Complex_wp(1), Y, false);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, LaplaceEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Types not compatible" << endl;
    abort();
  }
  

  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, LaplaceEquation<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    if (this->var.FirstOrderFormulation())
      MltAddScalarHelmholtzH1(Complex_wp(1, 0), nat_mat, trans, lvl, *this,
			      X, Complex_wp(1, 0), Y, false);
    else
      MltAddScalarH1(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
		     X, Complex_wp(1, 0), Y, false);
  }
  

  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, LaplaceEquationDG<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddScalarHelmholtz(Real_wp(1), nat_mat, trans, lvl, *this, 
			  X, Real_wp(1), Y, false);
  }


  template<class Dimension>
  void FemMatrixFreeClass<Real_wp, LaplaceEquationDG<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, LaplaceEquationDG<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    cout << "Types not compatible " << endl;
    abort();
  }


  template<class Dimension>
  void FemMatrixFreeClass<Complex_wp, LaplaceEquationDG<Dimension> >
  ::MltAddFree(const GlobalGenericMatrix<Complex_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    MltAddScalarHelmholtz(Complex_wp(1, 0), nat_mat, trans, lvl, *this, 
			  X, Complex_wp(1, 0), Y, false);
  }

}

#define MONTJOIE_FILE_PROD_MAT_VECT_HELMHOLTZ_CXX
#endif
