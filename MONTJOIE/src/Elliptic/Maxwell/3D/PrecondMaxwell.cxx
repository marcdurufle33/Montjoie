#ifndef MONTJOIE_FILE_PRECOND_MAXWELL_CXX

namespace Montjoie
{

  /**********************************
   * PreconditionerHelmholtzMaxwell *
   **********************************/
  
  
  //! Default constructor
  template<class T>
  PreconditionerHelmholtzMaxwell<T>::PreconditionerHelmholtzMaxwell()
  {      
    type_smoother = JACOBI;
    omega_ssor = 0.5;
    mat_ssor_sym = NULL;
    mat_ssor_unsym = NULL;
    matCSR_ssor_sym = NULL;
    matCSR_ssor_unsym = NULL;
  }
  
  
  //! returns size of memory used by the object
  template<class T>
  size_t PreconditionerHelmholtzMaxwell<T>::GetMemorySize() const
  {
    size_t taille = mat_stored_sym.GetMemorySize() + mat_stored_unsym.GetMemorySize();
    taille += Projection_NodalToEdge.GetMemorySize() + mat_galerkin_H1.GetMemorySize();
    taille += sizeof(T)*diagonal_galerkin_H1.GetM();
    taille += sizeof(T)*diagonal_jacobi.GetM() + sizeof(*this);
    taille += sizeof(int)*(DirichletHelmholtz_dof.GetM() + DirichletMaxwell_dof.GetM());
    return taille;
  }

  
  template<class T>
  void PreconditionerHelmholtzMaxwell<T>::SetParameterRelaxation(const Real_wp& w)
  {
    omega_ssor = w;
  }  
  
  
  //! constructing of preconditioning for Helmholtz preconditioner
  template<class T>
  void PreconditionerHelmholtzMaxwell<T>
  ::ConstructPreconditioner(int type, const VarComputationProblem& var_,
			    const FemMatrixFreeClass_Base<T>& A,
			    const GlobalGenericMatrix<T>& nature, bool same_nature)
  {
    type_smoother = type;
    FemMatrixFreeClass_Base<T>& Avar = const_cast<FemMatrixFreeClass_Base<T>& >(A);
    VarComputationProblem& var = const_cast<VarComputationProblem& >(var_);
    
    // Helmholtz problem is computed and its preconditioning
    this->ComputeHelmPrecond(var, this->mat_galerkin_H1,
			     this->Projection_NodalToEdge, 
			     this->diagonal_galerkin_H1, nature);
    
    // computes the iterative matrix for relaxation algorithm
    if (type_smoother == SSOR)
      {
	// computing matrix for SSOR iteration
        this->mat_stored_sym.Clear();
        this->mat_stored_unsym.Clear();
        if (IsSymmetricMatrix(A))
          {
            bool compute_mat = true;
            if (same_nature)
              if (Avar.SucceedInAffectingPointer(this->mat_ssor_sym, this->matCSR_ssor_sym))
                compute_mat = false;
            
            if (compute_mat)
              {
                var.AddMatrixWithBC(this->mat_stored_sym, nature);
                this->mat_ssor_sym = &this->mat_stored_sym;
              }
	    
	    //this->mat_ssor_sym->WriteText("mat_ssor.dat");
          }
        else
          {
            bool compute_mat = true;
            if (same_nature)
              if (Avar.SucceedInAffectingPointer(this->mat_ssor_unsym, this->matCSR_ssor_unsym))
                compute_mat = false;
            
            if (compute_mat)
              {
                var.AddMatrixWithBC(this->mat_stored_unsym, nature);
                this->mat_ssor_unsym = &this->mat_stored_unsym;
              }
          }
      }
    else
      {
	// we compute inverse of diagonals of mat_galerkin_H1 and A
	for (int i = 0; i < this->diagonal_galerkin_H1.GetM(); i++)
	  this->diagonal_galerkin_H1(i) = this->omega_ssor/this->diagonal_galerkin_H1(i);
	
	if (same_nature)
	  var.ComputeDiagonalMatrix(this->diagonal_jacobi, A, nature);
	else
	  var.ComputeDiagonalMatrix(this->diagonal_jacobi, nature);

	for (int i = 0; i < this->diagonal_jacobi.GetM(); i++)
	  this->diagonal_jacobi(i) = this->omega_ssor/this->diagonal_jacobi(i);
      }
  }
  

  //! construction of preconditioner based on Helmholtz decomposition for Maxwell equations
  /*!
    \param[in] var finite element object
    \param[in] mat_iterative matrix that has to be preconditioned
    \param[out] mat_gal Galerkin matrix (with nodal elements)
    \param[out] mat_csr_proj projector from nodal finite element to edge finite element
    \param[out] diag_iter diagonal of iterative matrix
    \param[out] diag_gal diagonal of Galkerin matrix
   */
  template<class T>
  void PreconditionerHelmholtzMaxwell<T>
  ::ComputeHelmPrecond(const VarComputationProblem& var_,
		       typename GenericStorage<T>::SparseSymMatrix& mat_gal,
		       Matrix<Real_wp, General, RowSparse>& mat_csr_proj,
                       Vector<T>& diag_gal, const GlobalGenericMatrix<T>& nature)
  {
    Matrix<Real_wp, General, ArrayRowSparse> mat_proj;
    const EllipticProblem<HarmonicMaxwellEquation_3D>* var_ptr = NULL;
    
    try 
      {
        var_ptr = dynamic_cast<const EllipticProblem<HarmonicMaxwellEquation_3D>* >(&var_);
      }
    catch (const std::bad_cast&)
      {
	cout << "Preconditioning only available for H(curl) elements" << endl;
	abort();
      }

    const EllipticProblem<HarmonicMaxwellEquation_3D>& var = *var_ptr;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    
    int rank_proc = var.GetRankProcMode();
    
    // we construct P
    if (mesh_num.GetOrder() == 0)
      {
        int nb_edges = var.mesh.GetNbEdges();
        int nb_vertices = var.mesh.GetNbVertices();
	mat_proj.Reallocate(nb_edges, nb_vertices);
	for (int i = 0; i < nb_edges; i++)
	  {
	    int n1 = var.mesh.GetEdge(i).numVertex(0);
	    int n2 = var.mesh.GetEdge(i).numVertex(1);
	    mat_proj.AddInteraction(i, n1, -1.0);
	    mat_proj.AddInteraction(i, n2, 1.0);
	  }
      }
    else
      {
	if ((rank_proc == 0)||(var.print_level >= 10))
	  cout<<"computation helmholtz problem "<<endl;
	
	// constructing nodal finite element    
	EllipticProblem<HelmholtzEquation<Dimension3> > var_helm; 
	MeshNumbering<Dimension3>& mesh_num_helm = var_helm.GetMeshNumbering(0);
        
        var_helm.mesh = var.mesh;
	mesh_num_helm.SetOrder(mesh_num.GetOrder());
	var_helm.SetFrequency(0.0);
	var_helm.InitIndices(var.ref_epsilon.GetM()-1);
	var_helm.ref_mu = var.ref_epsilon;
	
        var_helm.grid_to_be_computed = false;
        var_helm.SetTypeEquation("HELMHOLTZ");
	var_helm.ComputeMeshAndFiniteElement(string("TETRAHEDRON_LOBATTO"));
	//var_helm.output_grid_param = var.output_grid_param;
	//var_helm.var_grid = var.var_grid;
	//var_helm.type_output_file = var.type_output_file;

        GlobalGenericMatrix<T> nat_mat;
        T zero, one; SetComplexZero(zero);
	SetComplexOne(one);
        nat_mat.SetCoefStiffness(-var.GetSquareOmega()*nature.GetCoefMass());
	nat_mat.SetCoefMass(zero);
	var_helm.ComputeMassMatrix();
        
	// first-order abc
	ImpedanceFunction_Base<T, Dimension3> fimped(var_helm);
	T coef_imped;
	to_complex(-Iwp*var.GetOmega()*sqrt(nature.GetCoefMass()), coef_imped);
	fimped.SetCoefficient(zero, coef_imped);
        
	mat_gal.Reallocate(var_helm.GetNbDof(), var_helm.GetNbDof());
	var_helm.AddMatrixImpedanceBoundary(one, var_helm.mesh.GetBoundaryCondition(),
					    BoundaryConditionEnum::LINE_ABSORBING,
					    nat_mat, mat_gal, 0, 0, fimped);
        
        var_helm.AddMatrixWithBC(mat_gal, nat_mat);
	
        //mat_gal.WriteText("mat_helm.dat");
        var_helm.ComputeDiagonalMatrix(diag_gal, mat_gal, nat_mat);
	
	DirichletHelmholtz_dof = var_helm.GetDirichletDofNumber();
        DirichletMaxwell_dof = var.GetDirichletDofNumber();
        
	mat_proj.Clear();
	mat_proj.Reallocate(var.GetNbDof(), var_helm.GetNbDof());
        
	// now asssembling the matrix P
        VectReal_wp InvDh(var.GetNbDof());
        InvDh.Fill(0);
		
	for (int i = 0; i < var.mesh.GetNbElt(); i++)
	  {
	    const ElementReference<Dimension3, 2>& Fb_ref = var.GetReferenceElementHcurl(i);
	    const ElementReference<Dimension3, 1>& Fb_h1 = var_helm.GetReferenceElementH1(i);
            
            const HexahedronHcurlFirstFamily* Fb_ptr = NULL;
            const HexahedronLobatto* hex_lob_ptr = NULL;
	    try
	      {
		Fb_ptr = &dynamic_cast<const HexahedronHcurlFirstFamily& >(Fb_ref);
		hex_lob_ptr = &dynamic_cast<const HexahedronLobatto& >(Fb_h1);
	      }
	    catch (const std::bad_cast&)
	      {
		cout << "Projection defined only for HexahedronHcurlFirstFamily" << endl;
		abort();
	      }
            	    
	    const HexahedronHcurlFirstFamily& Fb = *Fb_ptr;	    	    
	    const HexahedronLobatto& hex_lob = *hex_lob_ptr;
            
            const Array3D<int>& NumDofs_X = Fb.GetNumDofsX();
            const Array3D<int>& NumDofs_Y = Fb.GetNumDofsY();
            const Array3D<int>& NumDofs_Z = Fb.GetNumDofsZ();
            const Matrix<Real_wp>& dGL_G = Fb.GetGradGL_G();
            Vector<bool> NegativeDof(Fb.GetNbDof());
            NegativeDof.Fill(false);
            const Vector<int>& neg_dof = mesh_num.Element(i).GetNegativeDofNumber();
            for (int j = 0; j < neg_dof.GetM(); j++)
              NegativeDof(neg_dof(j)) = true;
            
            for(int j = 0; j < Fb.GetNbDof(); j++)
              InvDh(mesh_num.Element(i).GetNumberDof(j)) += 1.0;
	    
	    int r = Fb.GetOrder();
	    for (int i1 = 0; i1 <= r; i1++)
	      for (int i2 = 0; i2 <= r; i2++)
		for (int i3 = 0; i3 <= r; i3++)
		  {
		    int num_loc = hex_lob.GetNumNodes3D(i1, i2, i3);
		    int num_dof = mesh_num_helm.Element(i).GetNumberDof(num_loc);
                    
		    for (int j = 0; j < r; j++)
		      {
			// for dofs oriented along x
			int num_loc_vec = NumDofs_X(j, i2, i3);
			int num_dof_vec = mesh_num.Element(i).GetNumberDof(num_loc_vec);
			Real_wp vloc = dGL_G(i1, j);
                        if (!var.IsDofDirichlet(num_dof_vec) &&
                            !var_helm.IsDofDirichlet(num_dof))
                          {
                            if (NegativeDof(num_loc_vec))
                              mat_proj.AddInteraction(num_dof_vec, num_dof, -vloc);
                            else
                              mat_proj.AddInteraction(num_dof_vec, num_dof, vloc);
                          }
                        
			// for dofs oriented along y
			num_loc_vec = NumDofs_Y(i1, j, i3);
			num_dof_vec = mesh_num.Element(i).GetNumberDof(num_loc_vec);
			vloc = dGL_G(i2, j);
                        if (!var.IsDofDirichlet(num_dof_vec) &&
                            !var_helm.IsDofDirichlet(num_dof))
                          {
                            if (NegativeDof(num_loc_vec))
                              mat_proj.AddInteraction(num_dof_vec, num_dof, -vloc);
                            else
                              mat_proj.AddInteraction(num_dof_vec, num_dof, vloc);
                          }
                        
			// for dofs oriented along z
			num_loc_vec = NumDofs_Z(i1, i2, j);
			num_dof_vec = mesh_num.Element(i).GetNumberDof(num_loc_vec);
			vloc = dGL_G(i3, j);
                        if (!var.IsDofDirichlet(num_dof_vec) &&
                            !var_helm.IsDofDirichlet(num_dof))
                          {
                            if (NegativeDof(num_loc_vec))
                              mat_proj.AddInteraction(num_dof_vec, num_dof, -vloc);
                            else
                              mat_proj.AddInteraction(num_dof_vec, num_dof, vloc);
                          }
                      }
                  }
	  }
	
        for (int i = 0; i < InvDh.GetM(); i++)
          InvDh(i) = 1.0/InvDh(i);
	
	// ponderation
	for (int i = 0; i < var.GetNbDof(); i++)
	  for (int j = 0; j < mat_proj.GetRowSize(i); j++)
	    mat_proj.Value(i,j) *= InvDh(i);
      }
    
    //mat_proj.WriteText("ProjNodal.dat");
    Seldon::Copy(mat_proj, mat_csr_proj);
  }
  


  //! applies preconditioning
  template<class T> template<class T0>
  void PreconditionerHelmholtzMaxwell<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<T0>& r, Vector<T0>& z)
  {
    if (trans.Trans())
      {
	cout << "Not implemented" << endl;
	abort();
      }
    
    z.Fill(0);
    if (type_smoother == SSOR)
      {
        if (A.IsSymmetric())
          {
	    if (mat_ssor_sym != NULL)
	      this->SolveHelmholtzDecomposition(*mat_ssor_sym, r, z, true);
	    else
	      this->SolveHelmholtzDecomposition(*matCSR_ssor_sym, r, z, true);
	  }
        else
          {
	    if (mat_ssor_unsym != NULL)
	      this->SolveHelmholtzDecomposition(*mat_ssor_unsym, r, z, true);
	    else
	      this->SolveHelmholtzDecomposition(*matCSR_ssor_unsym, r, z, true);
	  }
      }
    else
      this->SolveHelmholtzDecomposition(A, r, z, true);
  }
  
  
  //! applies preconditioning based on Helmholtz decomposition
  /*!
    \param[in] A matrix that is preconditioned
    \param[in] r vector to multiply
    \param[in,out] z result vector after mutiplication by the preconditioning operator
    \param[in] init_guess_null true if z is set to 0 initially
    The preconditioning may be unsymmetric for Dirichlet conditions
   */
  template<class T> template<class T0>
  void PreconditionerHelmholtzMaxwell<T>
  ::SolveHelmholtzDecomposition(const VirtualMatrix<T>& A, const Vector<T0>& r,
                                Vector<T0>& z, bool init_guess_null)
  {
    T0 one, zero;
    SetComplexOne(one); SetComplexZero(zero);
    int nphi(0), n = A.GetM();
    if (this->type_smoother == JACOBI)
      nphi = this->diagonal_galerkin_H1.GetM();
    else
      nphi = this->mat_galerkin_H1.GetM();
    
    Vector<T0> xphi(nphi), rphi(nphi), w(n);
    xphi.Zero(); rphi.Zero(); w.Zero();
    
    //VectComplexe& r_modif = const_cast<VectComplexe&>(r);
    //for (int i = 0; i < DirichletMaxwell_dof.GetM(); i++)
    //r_modif(DirichletMaxwell_dof(i)) = 0;
    
    // projection on nodal elements
    // r_\phi \; = \; P^t (r - A z)
    if (init_guess_null)
      Mlt(SeldonTrans, this->Projection_NodalToEdge, r, rphi);
    else
      {
        w.Copy(r);
        A.MltAddVector(-one, z, one, w);
        Mlt(SeldonTrans, this->Projection_NodalToEdge, w, rphi);
      }

    for (int i = 0; i < DirichletHelmholtz_dof.GetM(); i++)
      rphi(DirichletHelmholtz_dof(i)) = 0;
    
    if (this->type_smoother == JACOBI)
      {
        for (int i = 0; i < nphi; i++)
          xphi(i) = this->diagonal_galerkin_H1(i)*rphi(i);
      }
    else
      {
        // descente Gauss-Seidel sur A_\phi x_\phi \; = \; r_\phi
        this->mat_galerkin_H1.ApplySor(SeldonNoTrans, xphi, rphi, this->omega_ssor, 1, 2);
        // xphi = rphi; SolveLU(mat_glob_lu, xphi);
      }
    
    // rajout de la relaxation : z = z + P x_\phi
    MltAdd(one, this->Projection_NodalToEdge, xphi, one, z);

    //for (int i = 0; i < DirichletMaxwell_dof.GetM(); i++)
    //z(DirichletMaxwell_dof(i)) = 0;
    
    if (this->type_smoother == JACOBI)
      {
        w.Copy(r);
        A.MltAddVector(-one, z, one, w);
	
        for (int i = 0; i < n; i++)
          z(i) += this->diagonal_jacobi(i)*w(i);
      }
    else
      {
        // descente et remontee de Gauss-Seidel sur A z = r
        if (A.IsSymmetric())
          {
	    if (this->mat_ssor_sym != NULL)
	      SOR(*this->mat_ssor_sym, z, r, this->omega_ssor, 1, 6);
	    else
	      SOR(*this->matCSR_ssor_sym, z, r, this->omega_ssor, 1, 6);
	  }
        else
          {
	    if (this->mat_ssor_unsym != NULL)
	      SOR(*this->mat_ssor_unsym, z, r, this->omega_ssor, 1, 6);
	    else
	      SOR(*this->matCSR_ssor_unsym, z, r, this->omega_ssor, 1, 6);
	  }
      }
    
    // calcul de r_\phi \; = \; P^t (r - A z)
    w.Copy(r);
    A.MltAddVector(-one, z, one, w);
    MltAdd(one, SeldonTrans, this->Projection_NodalToEdge, w, zero, rphi);
    
    for (int i = 0; i < DirichletHelmholtz_dof.GetM(); i++)
      rphi(DirichletHelmholtz_dof(i)) = 0;
    
    if (this->type_smoother == JACOBI)
      {
        for (int i = 0; i < nphi; i++)
          xphi(i) = this->diagonal_galerkin_H1(i)*rphi(i);
      }
    else
      {
        xphi.Fill(0); 
        // Remontee de Gauss-Seidel sur A_\phi x_\phi = r_\phi
        this->mat_galerkin_H1.ApplySor(SeldonNoTrans, xphi, rphi, this->omega_ssor, 1, 3);
        // xphi = rphi; SolveLU(var_helmholtz.mat_glob_lu, xphi);
      }
    
    // rajout de la contribution de x_\phi
    MltAdd(one, this->Projection_NodalToEdge, xphi, one, z);

    //for (int i = 0; i < DirichletMaxwell_dof.GetM(); i++)
    //z(DirichletMaxwell_dof(i)) = 0;
  }
  
  
  //! applying pre-smoothing step
  /*!
    \param[in] A matrix that is preconditioned
    \param[in] r vector to multiply
    \param[out] z result vector after mutiplication by the preconditioning operator
    \param[in] init_guess_null 
   */
  template<class T> template<class T0>
  void PreconditionerHelmholtzMaxwell<T>
  ::ApplyPreSmoothing(const VirtualMatrix<T>& A, const Vector<T0>& r,
                      Vector<T0>& z, bool init_guess_null)
  {
    T0 one, zero;
    SetComplexOne(one); SetComplexZero(zero);

    int nphi(0), n = A.GetM();
    if (this->type_smoother == JACOBI)
      nphi = this->diagonal_galerkin_H1.GetM();
    else
      nphi = this->mat_galerkin_H1.GetM();
	
    Vector<T0> xphi(nphi), rphi(nphi), w(n);
    xphi.Zero(); rphi.Zero();
    // projection on nodal elements
    // r_\phi \; = \; P^t (r - A z)
    if (init_guess_null)
      Mlt(SeldonTrans, this->Projection_NodalToEdge, r, rphi);
    else
      {
	w.Copy(r);
	A.MltAddVector(-one, z, one, w);
        Mlt(SeldonTrans, this->Projection_NodalToEdge, w, rphi);
      }
    
    for (int i = 0; i < DirichletHelmholtz_dof.GetM(); i++)
      rphi(DirichletHelmholtz_dof(i)) = zero;
        
    if (this->type_smoother == JACOBI)
      {
	for (int i = 0; i < nphi; i++)
	  xphi(i) = this->diagonal_galerkin_H1(i)*rphi(i);
      }
    else
      {
	// descente Gauss-Seidel sur A_\phi x_\phi \; = \; r_\phi
	this->mat_galerkin_H1.ApplySor(SeldonNoTrans, xphi, rphi, this->omega_ssor, 1, 2);
      }
    
    // rajout de la relaxation : z = z + P x_\phi
    MltAdd(one, this->Projection_NodalToEdge, xphi, one, z);
	
    if (this->type_smoother == JACOBI)
      {
	w.Copy(r);
	A.MltAddVector(-one, z, one, w);
	
	for (int i = 0; i < n; i++)
	  z(i) += this->diagonal_jacobi(i)*w(i);
      }
    else
      {
	// descente de Gauss-Seidel sur A z = r
	if (A.IsSymmetric())
          {
	    if (this->mat_ssor_sym != NULL)
	      SOR(*this->mat_ssor_sym, z, r, this->omega_ssor, 1, 2);
	    else
	      SOR(*this->matCSR_ssor_sym, z, r, this->omega_ssor, 1, 2);
	  }
        else
          {
	    if (this->mat_ssor_unsym != NULL)
	      SOR(*this->mat_ssor_unsym, z, r, this->omega_ssor, 1, 2);
	    else
	      SOR(*this->matCSR_ssor_unsym, z, r, this->omega_ssor, 1, 2);
	  }
      }
  }
  
  
  //! applying post-smoothing step
  /*!
    \param[in] A matrix that is preconditioned
    \param[in] r vector to multiply
    \param[out] z result vector after mutiplication by the preconditioning operator
    \param[in] init_guess_null 
   */
  template<class T> template<class T0>
  void PreconditionerHelmholtzMaxwell<T>
  ::ApplyPostSmoothing(const VirtualMatrix<T>& A, const Vector<T0>& r,
                       Vector<T0>& z, bool init_guess_null)
  { 
    T0 one, zero;
    SetComplexOne(one); SetComplexZero(zero);

    int nphi(0), n = A.GetM();
    if (this->type_smoother == JACOBI)
      nphi = this->diagonal_galerkin_H1.GetM();
    else
      nphi = this->mat_galerkin_H1.GetM();
    
    Vector<T0> xphi(nphi), rphi(nphi), w(n);
    xphi.Zero(); rphi.Zero();
        
    if (this->type_smoother == JACOBI)
      {
	w.Copy(r);
	A.MltAddVector(-one, z, one, w);
	
	for (int i = 0; i < n; i++)
	  z(i) += this->diagonal_jacobi(i)*w(i);
      }
    else
      {
	// remontee de Gauss-Seidel sur A z = r
	if (A.IsSymmetric())
          {
	    if (this->mat_ssor_sym != NULL)
	      SOR(*this->mat_ssor_sym, z, r, this->omega_ssor, 1, 3);
	    else
	      SOR(*this->matCSR_ssor_sym, z, r, this->omega_ssor, 1, 3);
	  }
        else
          {
	    if (this->mat_ssor_unsym != NULL)
	      SOR(*this->mat_ssor_unsym, z, r, this->omega_ssor, 1, 3);
	    else
	      SOR(*this->matCSR_ssor_unsym, z, r, this->omega_ssor, 1, 3);
	  }
      }
    
    // calcul de r_\phi \; = \; P^t (r - A z)
    w.Copy(r);
    A.MltAddVector(-one, z, one, w);
    Mlt(SeldonTrans, this->Projection_NodalToEdge, w, rphi);
    
    for (int i = 0; i < DirichletHelmholtz_dof.GetM(); i++)
      rphi(DirichletHelmholtz_dof(i)) = zero;
    
    if (this->type_smoother == JACOBI)
      {
	for (int i = 0; i < nphi; i++)
	  xphi(i) = this->diagonal_galerkin_H1(i)*rphi(i);
      }
    else
      {
	// Remontee de Gauss-Seidel sur A_\phi x_\phi = r_\phi
	xphi.Fill(0);
	Seldon::SOR(this->mat_galerkin_H1, xphi, rphi, this->omega_ssor, 1, 3);
      }
    
    // rajout de la contribution de x_\phi
    MltAdd(one, this->Projection_NodalToEdge, xphi, one, z);
  }


  /***********************************
   * LowOrderPreconditioning_Maxwell *
   ***********************************/
  


  //! correspondance locale
  template<class T>
  void LowOrderPreconditioning_Maxwell<T>
  ::GetLocalPermutation_LowOrder(const ElementReference<Dimension3, 2>& FaceFine_,
				 const ElementReference<Dimension3, 2>& FaceCoarse_,
                                 const VectReal_wp& Points,
                                 IVect& num_element, IVect& num_dof_local,
                                 VectReal_wp& LocalScaling)
  {
    const HexahedronHcurlFirstFamily* FaceCoarse_ptr = NULL;
    const HexahedronHcurlFirstFamily* FaceFine_ptr = NULL;
    try
      {
	FaceCoarse_ptr
	  = &dynamic_cast<const HexahedronHcurlFirstFamily& >(FaceCoarse_);
	
	FaceFine_ptr
	  = &dynamic_cast<const HexahedronHcurlFirstFamily& >(FaceFine_);
      }
    catch (const std::bad_cast&)
      {
	cout << "Projection defined only for HexahedronHcurlFirstFamily" << endl;
	abort();
      }

    const HexahedronHcurlFirstFamily& FaceCoarse = *FaceCoarse_ptr;
    const HexahedronHcurlFirstFamily& FaceFine = *FaceFine_ptr;
    
    int nb_dof_loc = FaceFine.GetNbDof();
    int r = FaceFine.GetOrder();
    LocalScaling.Reallocate(nb_dof_loc); LocalScaling.Fill(0);
    num_element.Reallocate(nb_dof_loc); num_dof_local.Reallocate(nb_dof_loc);
    num_element.Fill(-1); num_dof_local.Fill(-1);
    int num_elem = 0; IVect num(4); num.Zero();
    R3 ext1, ext2;
    const Array3D<int>& NumDofs_X = FaceFine.GetNumDofsX();
    const Array3D<int>& NumDofs_Y = FaceFine.GetNumDofsY();
    const Array3D<int>& NumDofs_Z = FaceFine.GetNumDofsZ();
    const Matrix<int>& CoordinateDofs = FaceCoarse.GetCoordinateDofs();    
    Real_wp eps = epsilon_machine;
    // boucle sur tous les petits cubes du gros cube divise en r intervalles
    for (int i = 0; i < r; i++)
      for (int j = 0; j < r; j++)
	for (int k = 0; k < r; k++)
	  {
            // ext1 et ext 2 sont les points opposés du petit cube
	    ext1(0) = Points(i)-eps; ext1(1) = Points(j)-eps; ext1(2) = Points(k)-eps;
	    ext2(0) = Points(i+1)+eps; ext2(1) = Points(j+1)+eps; ext2(2) = Points(k+1)+eps;

	    // boucle sur tous les ddl du petit cube (on utilise Q1 sur chaque petit cube)
	    for (int p_c = 0; p_c < 12; p_c++)
	      {
		int num_dof_coarse = -1;
		for (int l = 0; l < 4; l++)
		  num(l) = CoordinateDofs(p_c, l);
                
                if (num(3) == 0)
                  num_dof_coarse = NumDofs_X(i, j+num(1), k+num(2));
                else if (num(3) == 1)
                  num_dof_coarse = NumDofs_Y(i+num(0), j, k+num(2));
		else
                  num_dof_coarse = NumDofs_Z(i+num(0), j+num(1), k);
		
                num_element(num_dof_coarse) = num_elem;
		num_dof_local(num_dof_coarse) = p_c;
                
                LocalScaling(num_dof_coarse) = 1.0/(ext2(num(3)) - ext1(num(3)));
              }
            
            num_elem++;
          }
  }
  

  //! correspondance maillage ordre eleve et maillage d'ordre bas
  template<class T>
  void LowOrderPreconditioning_Maxwell<T>
  ::ComputeLocalProlongation(const VectReal_wp& Points,
                             DistributedProblem<Dimension3>& var_fine_,
                             DistributedProblem<Dimension3>& var_low_)
  {
    IVect num_element, num_dof_local;
    VectReal_wp LocalScaling;
        
    VarProblem<Dimension3>* var_fine_ptr = NULL;
    VarProblem<Dimension3>* var_low_ptr = NULL;

    try
      {
	var_fine_ptr
	  = &dynamic_cast<VarProblem<Dimension3>& >(var_fine_);
	
        var_low_ptr
	  = &dynamic_cast<VarProblem<Dimension3>& >(var_low_);
      }
    catch (const std::bad_cast&)
      {
	cout << "Preconditioning only available for H(curl) elements" << endl;
	abort();
      }
    
    VarProblem<Dimension3>& var_fine = *var_fine_ptr;
    VarProblem<Dimension3>& var_low = *var_low_ptr;
    
    const ElementReference<Dimension3, 2>& Fb_coarse = var_low.GetReferenceElementHcurl(0);
    const ElementReference<Dimension3, 2>& Fb_fine = var_fine.GetReferenceElementHcurl(0);

    this->GetLocalPermutation_LowOrder(Fb_fine, Fb_coarse, Points,
                                       num_element, num_dof_local, LocalScaling);

    const MeshNumbering<Dimension3>& mesh_num_fine = var_fine.GetMeshNumbering(0);
    const MeshNumbering<Dimension3>& mesh_num_low = var_low.GetMeshNumbering(0);
        
    int Nvol = mesh_num_fine.GetNbDof();
    this->NumberDof_LowOrder.Reallocate(Nvol);
    this->ScalingDof_LowOrder.Reallocate(Nvol);
    this->ScalingDof_LowOrder.Fill(1.0);
    int offset_elt = 0;
    for (int i = 0; i < var_fine.mesh.GetNbElt(); i++)
      {
        // retrieving sign of dofs
        int nb_dof_loc = mesh_num_fine.Element(i).GetNbDof();
        Vector<bool> NegativeFine(nb_dof_loc);
        NegativeFine.Fill(false);
        const Vector<int>& neg_dof = mesh_num_fine.Element(i).GetNegativeDofNumber();
        for (int j = 0; j < neg_dof.GetM(); j++)
          NegativeFine(neg_dof(j)) = true;
        
	int r = Fb_fine.GetOrder();
        Matrix<bool> NegativeLow(r*r*r, 12);
        NegativeLow.Fill(false);
        for (int k = 0; k < r*r*r; k++)
          {
            const Vector<int>& neg_low = mesh_num_low.Element(offset_elt+k).GetNegativeDofNumber();
            for (int j = 0; j < neg_low.GetM(); j++)
              NegativeLow(k, neg_low(j)) = true;
          }
        
        // then filling ScalingDof_LowOrder and NumberDof_LowOrder
        for (int j = 0; j < nb_dof_loc; j++)
          {
            int num_dof_fine = mesh_num_fine.Element(i).GetNumberDof(j);
            int num_low = mesh_num_low.Element(offset_elt+num_element(j))
              .GetNumberDof(num_dof_local(j));
            
            this->NumberDof_LowOrder(num_dof_fine) = num_low;
            if (NegativeFine(j) == NegativeLow(num_element(j), num_dof_local(j)))
              this->ScalingDof_LowOrder(num_low) = LocalScaling(j);
            else
              this->ScalingDof_LowOrder(num_low) = -LocalScaling(j);
          }
        
        offset_elt += r*r*r;
      }
    
    //this->NumberDof_LowOrder.WriteText("num_low.dat");
    //this->ScalingDof_LowOrder.WriteText("scaling_low.dat");
    //abort();
  }
  
  
  //! applies preconditioning
  template<class T> template<class T0>
  void LowOrderPreconditioning_Maxwell<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	     const Vector<T0>& b, Vector<T0>& x)
  {
    Vector<T0> z(b);
    for (int i = 0; i < z.GetM(); i++)
      {
        int j = this->NumberDof_LowOrder(i);
        z(j) = this->ScalingDof_LowOrder(j)*b(i);
      }
    
    this->mat_lu.Solve(trans, z);
    
    for (int i = 0; i < z.GetM(); i++)
      {
        int j = this->NumberDof_LowOrder(i);
        x(i) = this->ScalingDof_LowOrder(j)*z(j);
      }    
  }


  //! applies preconditioning
  template<class T>
  void LowOrderPreconditioning_Maxwell<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonNoTrans, A, b, x);
  }

  
  //! applies transpose preconditioning
  template<class T>
  void LowOrderPreconditioning_Maxwell<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& b, Vector<T>& x)
  {
    SolveGen(SeldonTrans, A, b, x);
  }


  template<>
  void LowOrderPreconditioning_Maxwell<Complex_wp>
  ::Solve(const SeldonTranspose&, const VirtualMatrix<Complex_wp>& A, const Vector<Real_wp>& b, Vector<Real_wp>& x)
  {
    cout << "Incompatibles types" << endl;
    abort();
  }

  
  template<class T>
  void LowOrderPreconditioning_Maxwell<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x)
  {
    SolveGen(trans, A, b, x);
  }

  
  template<class T>
  void LowOrderPreconditioning_Maxwell<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x)
  {
    SolveGen(trans, A, b, x);
  }

  
  /************************************
   * MultigridPreconditioning_Maxwell *
   ************************************/


  //! returns size of memory used in bytes
  template<class T>
  size_t MultigridPreconditioning_Maxwell<T>::GetMemorySize() const
  {
    size_t taille = MultigridPreconditioning_Dim<T, Dimension3>::GetMemorySize();
    taille += smoother_helm.GetMemorySize();
    return taille;
  }
  
  
  //! additional parameters
  template<class T>
  void MultigridPreconditioning_Maxwell<T>
  ::SetInputData(const string& keyword, const VectString& param)
  {
    if (!keyword.compare("Smoother"))
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of MultigridPreconditioning" << endl;
	    cout << "Smoother needs more parameters, for instance :" << endl;
	    cout << "Smoother = JACOBI" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

	VectString other_param(param.GetM()-1);
	for (int i = 1; i < param.GetM(); i++)
	  other_param(i-1) = param(i);
	
	if (!param(0).compare("JACOBI"))
	  {
	    this->type_smoother = JACOBI;
	    this->smoother_jacobi.SetInputPreconditioning(param(0), other_param);
	  }
	else if (!param(0).compare("BLOCK_JACOBI"))
	  {
	    this->type_smoother = BLOCK_JACOBI;
	    this->smoother_jacobi.SetInputPreconditioning(param(0), other_param);
	  }
	else if (!param(0).compare("SSOR"))
	  {
	    if (param.GetM() <= 1)
	      {
		cout << "In SetInputData of MultigridPreconditioning" << endl;
		cout << "Smoother needs more parameters, for instance :" << endl;
		cout << "Smoother = SSOR omega" << endl;
		cout << "Current parameters are : " << endl << param << endl;
		abort();
	      }

	    this->type_smoother = SSOR;
	    this->smoother_ssor.SetParameterRelaxation(to_num<Real_wp>(param(1)));
	    if (param.GetM() > 2)
	      this->smoother_ssor.SetNumberIterations(to_num<int>(param(2)));
	  }
	else if ((!param(0).compare("HELMHOLTZ_SSOR")) || (!param(0).compare("HELMHOLTZ")))
	  {
	    if (param.GetM() <= 1)
	      {
		cout << "In SetInputData of MultigridPreconditioning" << endl;
		cout << "Smoother needs more parameters, for instance :" << endl;
		cout << "Smoother = HELMHOLTZ_SSOR omega" << endl;
		cout << "Current parameters are : " << endl << param << endl;
		abort();
	      }

	    this->type_smoother = HELMHOLTZ_SSOR;
	    smoother_helm.SetParameterRelaxation(to_num<Real_wp>(param(1)));
	  }
	else if (!param(0).compare("HELMHOLTZ_JACOBI"))
	  {
	    if (param.GetM() <= 1)
	      {
		cout << "In SetInputData of MultigridPreconditioning" << endl;
		cout << "Smoother needs more parameters, for instance :" << endl;
		cout << "Smoother = HELMHOLTZ_JACOBI omega" << endl;
		cout << "Current parameters are : " << endl << param << endl;
		abort();
	      }

	    this->type_smoother = HELMHOLTZ_JACOBI;
	    smoother_helm.SetParameterRelaxation(to_num<Real_wp>(param(1)));
	  }        
      }
    else 
      MultigridPreconditioning_Dim<T, Dimension3>::SetInputData(keyword, param);
  }

  
  //! computes smoother
  template<class T>
  void MultigridPreconditioning_Maxwell<T>
  ::ConstructSmoother(const FemMatrixFreeClass_Base<T>& A,
		      const GlobalGenericMatrix<T>& nature)
  {
    if (this->type_smoother == HELMHOLTZ_SSOR)
      {
        this->smoother_helm.ConstructPreconditioner(this->smoother_helm.SSOR,
                                                    this->var_fine, A, nature, true);
      }
    else if (this->type_smoother == HELMHOLTZ_JACOBI)
      {
        this->smoother_helm.ConstructPreconditioner(this->smoother_helm.JACOBI,
                                                    this->var_fine, A, nature, true);
      }
    else
      MultigridPreconditioning_Dim<T, Dimension3>::ConstructSmoother(A, nature);
  }
  
  
  //! applying smoothing iteration for multigrid algorithm
  /*!
    \param[in] A preconditioned matrix
    \param[out] x iterate after smoothing
    \param[in] b right-hand-side (Ax = b)
    \param[in,out] r intermediate vector
    \param[in] init_guess_null if true x is assumed equal to 0
   */
  template<class T> template<class T0>
  void MultigridPreconditioning_Maxwell<T>
  ::ApplyPreSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			Vector<T0>& r, bool init_guess_null)
  {
    if ((this->type_smoother == HELMHOLTZ_JACOBI) || (this->type_smoother == HELMHOLTZ_SSOR))
      smoother_helm.ApplyPreSmoothing(A, b, x, init_guess_null);
    else
      MultigridPreconditioning_Dim<T, Dimension3>::
        ApplyPreSmoother(A, x, b, r, init_guess_null);
  }


  //! applying smoothing iteration for multigrid algorithm
  /*!
    \param[in] A preconditioned matrix
    \param[out] x iterate after smoothing
    \param[in] b right-hand-side (Ax = b)
    \param[in,out] r intermediate vector
    \param[in] init_guess_null if true x is assumed equal to 0
   */
  template<class T> template<class T0>
  void MultigridPreconditioning_Maxwell<T>
  ::ApplyPostSmootherGen(const VirtualMatrix<T>& A, Vector<T0>& x, const Vector<T0>& b,
			 Vector<T0>& r, bool init_guess_null)
  {
    if ((this->type_smoother == HELMHOLTZ_JACOBI) || (this->type_smoother == HELMHOLTZ_SSOR))
      smoother_helm.ApplyPostSmoothing(A, b, x, init_guess_null);
    else
      MultigridPreconditioning_Dim<T, Dimension3>::
        ApplyPostSmoother(A, x, b, r, init_guess_null);
  }


  template<>
  void MultigridPreconditioning_Maxwell<Complex_wp>
  ::ApplyPreSmoother(const VirtualMatrix<Complex_wp>& A, Vector<Real_wp>& x, const Vector<Real_wp>& b,
		     Vector<Real_wp>& r, bool init_guess_null)
  {
    cout << "Incompatibles types" << endl;
    abort();
  }
  

  template<class T>
  void MultigridPreconditioning_Maxwell<T>
  ::ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
		     Vector<Treal>& r, bool init_guess_null)
  {
    ApplyPreSmootherGen(A, x, b, r, init_guess_null);
  }

  
  template<class T>
  void MultigridPreconditioning_Maxwell<T>
  ::ApplyPreSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
                           Vector<Tcplx>& r, bool init_guess_null)
  {
    ApplyPreSmootherGen(A, x, b, r, init_guess_null);
  }
  

  template<>
  void MultigridPreconditioning_Maxwell<Complex_wp>
  ::ApplyPostSmoother(const VirtualMatrix<Complex_wp>& A, Vector<Real_wp>& x, const Vector<Real_wp>& b,
		      Vector<Real_wp>& r, bool init_guess_null)
  {
    cout << "Incompatibles types" << endl;
    abort();
  }

  
  template<class T>
  void MultigridPreconditioning_Maxwell<T>
  ::ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Treal>& x, const Vector<Treal>& b,
                          Vector<Treal>& r, bool init_guess_null)
  {
    ApplyPostSmootherGen(A, x, b, r, init_guess_null);
  }
  
  
  template<class T>
  void MultigridPreconditioning_Maxwell<T>
  ::ApplyPostSmoother(const VirtualMatrix<T>& A, Vector<Tcplx>& x, const Vector<Tcplx>& b,
                           Vector<Tcplx>& r, bool init_guess_null)
  {
    ApplyPostSmootherGen(A, x, b, r, init_guess_null);
  }
  

  /******************************
   * All_Preconditioner_Maxwell *
   ******************************/

  
  //! returns size of memory used by the object in bytes
  template<class T>
  size_t All_Preconditioner_Maxwell<T>
  ::GetMemorySize() const
  {
    size_t taille = All_Preconditioner_Base<T>::GetMemorySize();
    taille += prec_helm.GetMemorySize();
    return taille;
  }
  

  //! sets parameters in the line TypeResolution = GMRES HELMHOLTZ
  template<class T>
  void All_Preconditioner_Maxwell<T>
  ::SetInputPreconditioning(const string& keyword, const Vector<string>& parameters)
  {
    if (keyword == "HELMHOLTZ")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputPreconditioning of All_Preconditioner_Maxwell" << endl;
	    cout << "HELMHOLTZ needs more parameters, for instance :" << endl;
	    cout << "HELMHOLTZ omega" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	this->type_precond = HELMHOLTZ_SSOR;
        Real_wp omega = to_num<Real_wp>(parameters(0));
        prec_helm.SetParameterRelaxation(omega);
      }
    else if (keyword == "HELMHOLTZ_JACOBI")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputPreconditioning of All_Preconditioner_Maxwell" << endl;
	    cout << "HELMHOLTZ_JACOBI needs more parameters, for instance :" << endl;
	    cout << "HELMHOLTZ_JACOBI omega" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	this->type_precond = HELMHOLTZ_JACOBI;
        Real_wp omega = to_num<Real_wp>(parameters(0));
        prec_helm.SetParameterRelaxation(omega);
      }
    else if (keyword == "HELMHOLTZ_SSOR")
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputPreconditioning of All_Preconditioner_Maxwell" << endl;
	    cout << "HELMHOLTZ_SSOR needs more parameters, for instance :" << endl;
	    cout << "HELMHOLTZ_SSOR omega" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	this->type_precond = HELMHOLTZ_SSOR;
        Real_wp omega = to_num<Real_wp>(parameters(0));
        prec_helm.SetParameterRelaxation(omega);
      }
    else
      All_Preconditioner_Base<T>::SetInputPreconditioning(keyword, parameters);
  }
  
  
  //! construction of preconditioning operator
  template<class T>
  void All_Preconditioner_Maxwell<T>
  ::ConstructPreconditioner(VarComputationProblem& var, All_LinearSolver& solver,
			    const FemMatrixFreeClass_Base<T>& A,
			    const GlobalGenericMatrix<T>& nature,
			    CondensationBlockSolver_Fem<T>& cond_solver)
  {    
    GlobalGenericMatrix<T> nat_mat;
    nat_mat.SetCoefMass(nature.GetCoefMass()*this->alpha_laplace);
    nat_mat.SetCoefStiffness(nature.GetCoefStiffness());
    nat_mat.SetCoefDamping(nature.GetCoefDamping());

    T one; SetComplexOne(one);
    if ((this->type_precond == HELMHOLTZ_JACOBI) || (this->type_precond == HELMHOLTZ_SSOR))
      {
        int type_smoother = prec_helm.JACOBI;
        if (this->type_precond == HELMHOLTZ_SSOR)
          type_smoother = prec_helm.SSOR;
        
        bool same_nature = false;
        if (this->alpha_laplace == one)
          same_nature = true;

        prec_helm.ConstructPreconditioner(type_smoother, var, A, nat_mat, same_nature);
      }
    else
      All_Preconditioner_Base<T>::ConstructPreconditioner(var, solver, A, nature, cond_solver);
  }
  
  
  //! applies preconditioning 
  /*!
    \param[in] A matrix that is preconditioned
    \param[in] r vector to multiply
    \param[out] z result vector after mutiplication by the preconditioning operator
   */
  template<class T>
  void All_Preconditioner_Maxwell<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
  {
    SolveGen(SeldonNoTrans, A, r, z);
  }

  
  template<class T>
  void All_Preconditioner_Maxwell<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
  {
    SolveGen(SeldonTrans, A, r, z);
  }

  
  template<class T> template<class T0>
  void All_Preconditioner_Maxwell<T>
  ::SolveGen(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	     const Vector<T0>& r, Vector<T0>& z)
  {
    if ((this->type_precond == HELMHOLTZ_JACOBI)
        ||(this->type_precond == HELMHOLTZ_SSOR))
      prec_helm.Solve(trans, A, r, z);
    else
      All_Preconditioner_Base<T>::Solve(trans, A, r, z);

  }

  
  template<>
  void All_Preconditioner_Maxwell<Complex_wp>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<Complex_wp>& A,
	  const Vector<Real_wp>& b, Vector<Real_wp>& x)
  {
    cout << "Incompatible types" << endl;
    abort();
  }
  
    
  template<class T>
  void All_Preconditioner_Maxwell<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Treal>& b, Vector<Treal>& x)
  {
    SolveGen(trans, A, b, x);
  }


  template<class T>
  void All_Preconditioner_Maxwell<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A, const Vector<Tcplx>& b, Vector<Tcplx>& x)
  {
    SolveGen(trans, A, b, x);
  }

}

#define MONTJOIE_FILE_PRECOND_MAXWELL_CXX
#endif

