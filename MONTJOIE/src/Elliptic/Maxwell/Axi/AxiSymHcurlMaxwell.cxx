#ifndef MONTJOIE_FILE_AXISYM_HCURL_MAXWELL_CXX

namespace Montjoie
{
  
  //! default constructor
  template<class TypeEquation>
  HarmonicMaxwellAxi_Hcurl<TypeEquation>::HarmonicMaxwellAxi_Hcurl()
    : HarmonicMaxwellAxi_Eq<TypeEquation>()
  {
    modified_formulation = false;
    this->mixed_formulation = false;
  }
  
  
  //! reading a line of the data file
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::SetInputData(const string& keyword, const Vector<string>& param)
  {
    HarmonicMaxwellAxi_Eq<TypeEquation>::SetInputData(keyword, param);
    
    if (keyword == "ModifiedFormulation")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of HarmonicMaxwellAxi_Hcurl" << endl;
	    cout << "ModifiedFormulation needs at least one parameter, for instance :" << endl;
	    cout << "ModifiedFormulation = YES" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "YES")
          modified_formulation = true;
        else
          modified_formulation = false;
      }
  }


  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::SetTypeEquation(const string&)
  {
    this->other_dg_formulation.Reallocate(1);
    this->other_dg_formulation(0) = ElementReference_Base::CONTINUOUS;
    
    if (this->other_mesh_num.GetM() > 0)
      for (int i = 0; i < this->other_mesh_num.GetM(); i++)
        delete this->other_mesh_num(i);
    
    this->other_mesh_num.Reallocate(1);
    this->other_mesh_num(0) = new MeshNumbering<Dimension2>(this->mesh);

    this->all_mesh_num.Reallocate(2);
    this->all_mesh_num(0) = &this->mesh_num;
    this->all_mesh_num(1) = this->other_mesh_num(0);
    
    this->mesh_num_unknown.Reallocate(2);
    this->mesh_num_unknown(0) = 0;
    this->mesh_num_unknown(1) = 1;
    this->InitPolarization();
  }


  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::SetTypeElement(string& name_elt)
  {
    string name_element(name_elt);
    if (name_element == "QUADRANGLE_HCURL_AXI")
      {
        name_elt = "QUADRANGLE_GAUSS_FIRST_FAMILY";
        this->name_other_elements.Reallocate(1);
        this->name_other_elements(0) = "TRIANGLE_CLASSICAL";
      }
    else
      {
        cout << "Unknown finite element" << endl;
        abort();
      }
  }
  
  
  //! computes the number of degrees of freedom for the problem
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::ComputeNumberOfDofs()
  {
    HarmonicMaxwellAxi_Eq<TypeEquation>::ComputeNumberOfDofs();

    this->OffsetDofH.Reallocate(this->mesh.GetNbElt()+1);
    this->OffsetDofH(0) = 0;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        const ElementReference<Dimension2, 2>& Fb_hcurl = this->GetReferenceElementHcurl(i);
        const ElementReference<Dimension2, 1>& Fb_h1 = this->GetReferenceElementH1(i, 1);
	int nb_dof = Fb_hcurl.GetNbDofH() + 2*Fb_h1.GetNbDof();
        OffsetDofH(i+1) = OffsetDofH(i) + nb_dof;
      }
    
    this->nodl = this->offset_dof_unknown(2);
    if (this->FirstOrderFormulation())
      {
	this->OffsetDofV = OffsetDofH;
        this->nodl += this->OffsetDofV(this->mesh.GetNbElt());
      }        
  }
    

  //! adding dofs for H and Hteta for first order formulation
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::PerformOtherInitializations()
  {            
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    this->coefficient_volumic_source = Real_wp(1) / this->GetOmega();
    
    // dofs located on the axis are retrieved
    // if IsDofOnAxe(i) is true, dof i is located on the z-axis
    Vector<bool> IsDofOnAxe(this->offset_dof_unknown(2));
    IsDofOnAxe.Fill(false);
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
	int n1 = this->mesh.BoundaryRef(i).numVertex(0),
	  n2 = this->mesh.BoundaryRef(i).numVertex(1);
	
	if (this->Vertex_On_Axe(n1)&&this->Vertex_On_Axe(n2))
	  {
            // scalar dofs on this edge
            int ne = i;
            int num_elem = this->mesh.Boundary(ne).numElement(0);
            int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(ne);
            int r = this->mesh_num.GetOrderElement(num_elem);
            
            const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem, 1);
	    for (int j = 0; j <= r; j++)
	      {
		int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
		int num_dof = this->other_mesh_num(0)->Element(num_elem).GetNumberDof(num_dof_loc);
		IsDofOnAxe(this->offset_dof_unknown(1) + num_dof) = true;
	      }
	  }
      }
    
    int nb_dof_axe = 0;
    for (int i = 0; i < IsDofOnAxe.GetM(); i++)
      if ((IsDofOnAxe(i)) && (!this->is_dof_dirichlet(i)))
        nb_dof_axe++;
    
    this->Dof_On_Axe.Reallocate(nb_dof_axe);
    nb_dof_axe = 0;
    for (int i = 0; i < IsDofOnAxe.GetM(); i++)
      if ((IsDofOnAxe(i)) && (!this->is_dof_dirichlet(i)))
        this->Dof_On_Axe(nb_dof_axe++) = i;    

#ifdef SELDON_WITH_MPI
    // arrays needed to assemble Dirichlet dofs are constructed    
    this->ConstructDirichletComm();
#endif

    if ((this->print_level >= 0) && (rank_proc == 0))
      cout << "Number of degrees of freedom with H, Hteta = " << this->nodl << endl;      
  }
  

  //! affectation of Dirichlet dofs
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::UpdateDirichlet(int n)
  {
    if (!modified_formulation)
      return;
    
    if (this->number_mode == 0)
      {
        // for m = 0, E_theta = 0 on the axis
        for (int i = 0; i < this->Dof_On_Axe.GetM(); i++)
          this->is_dof_dirichlet(this->Dof_On_Axe(i)) = true;        
      }
    else
      for (int i = 0; i < this->Dof_On_Axe.GetM(); i++)
        this->is_dof_dirichlet(this->Dof_On_Axe(i)) = false;
    
    this->nb_dof_dirichlet = 0;
    for (int i = 0; i < this->is_dof_dirichlet.GetM(); i++)
      if (this->is_dof_dirichlet(i))
        this->nb_dof_dirichlet++;
      
    this->Dirichlet_dof.Reallocate(this->nb_dof_dirichlet);
    this->nb_dof_dirichlet = 0;
    for (int i = 0; i < this->is_dof_dirichlet.GetM(); i++)
      if (this->is_dof_dirichlet(i))
        this->Dirichlet_dof(this->nb_dof_dirichlet++) = i;    

    this->nb_dof_dirichlet_all = this->nb_dof_dirichlet;
    
#ifdef SELDON_WITH_MPI
    MPI_Allreduce(&this->nb_dof_dirichlet, &this->nb_dof_dirichlet_all, 1, 
                  MPI_INTEGER, MPI_SUM, this->comm_group_mode);
#endif
  }
  
  
  //! allocation arrays needed to compte elementary matrices
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::AllocateMassMatrices()
  {
    HarmonicMaxwellAxi_Eq<TypeEquation>::AllocateMassMatrices();
    
    matmass_Hrz.Reallocate(this->mesh.GetNbElt());
    matmass_Erz.Reallocate(this->mesh.GetNbElt());
    matmass_Hteta.Reallocate(this->mesh.GetNbElt());
    matmass_Eteta.Reallocate(this->mesh.GetNbElt());
    if (this->modified_formulation)
      {
	matmass_Eteta_m0.Reallocate(this->mesh.GetNbElt());
	matmass_Eteta_rr.Reallocate(this->mesh.GetNbElt());
	matmass_Erz_r3.Reallocate(this->mesh.GetNbElt());
	matmass_Erz_m2.Reallocate(this->mesh.GetNbElt());
	matmass_Erz_zz.Reallocate(this->mesh.GetNbElt());
	matmass_Hteta_r3.Reallocate(this->mesh.GetNbElt());
	matmass_Hrz_r3.Reallocate(this->mesh.GetNbElt());
        vecstiff_Erz.Reallocate(this->mesh.GetNbElt());
        vecstiff_Eteta.Reallocate(this->mesh.GetNbElt());
	vecstiff_Eteta_r2.Reallocate(this->mesh.GetNbElt());
        //vec_coupling_E.Reallocate(this->mesh.GetNbElt());
        //mat_coupling_E.Reallocate(this->mesh.GetNbElt());
      }
    else
      {
        matstiff_radius.Reallocate(this->mesh.GetNbElt());
        matstiff_a10.Reallocate(this->mesh.GetNbElt());
        matstiff_a11.Reallocate(this->mesh.GetNbElt());
      }
  }
  
  
  //! computation of needed coefficients for an element of the mesh,
  //! in order to compute elementary matrix
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeHcurlMassMatrix(int i, const ElementReference_Dim<Dimension2>& Fb)
  {
    const EllipticProblem<TypeEquation> & var_leaf = this->GetLeafClass();

    int N = Fb.GetNbPointsQuadratureInside();

    matmass_Hrz(i).Reallocate(N);
    matmass_Erz(i).Reallocate(N);
    matmass_Hteta(i).Reallocate(N);
    matmass_Eteta(i).Reallocate(N);
    if (this->modified_formulation)
      {
	matmass_Eteta_m0(i).Reallocate(N);
	matmass_Eteta_rr(i).Reallocate(N);
	matmass_Erz_r3(i).Reallocate(N);
	matmass_Erz_m2(i).Reallocate(N);
	matmass_Erz_zz(i).Reallocate(N);
	matmass_Hteta_r3(i).Reallocate(N);
	matmass_Hrz_r3(i).Reallocate(N);
        vecstiff_Erz(i).Reallocate(N);
        vecstiff_Eteta(i).Reallocate(N);
        vecstiff_Eteta_r2(i).Reallocate(N);
        //vec_coupling_E(i).Reallocate(N);
        //mat_coupling_E(i).Reallocate(N);
      }
    else
      {
        matstiff_radius(i).Reallocate(N);
        matstiff_a10(i).Reallocate(N);
        matstiff_a11(i).Reallocate(N);
      }
        
    
    TinyMatrix<Complexe, Symmetric, 2, 2> mu_rz, eps_rz, invMu_rz;
    Complexe eps_teta, mu_teta, invMu_teta, miomega;
    
    Matrix2_2 dfjm1;
	
    Real_wp coef_mass = TypeEquation::GetCoefficientMassMatrix(var_leaf);
    Real_wp coef_stiff = TypeEquation::GetCoefficientStiffnessMatrix(var_leaf);

    int ref = this->mesh.Element(i).GetReference();    
    this->GetMiomega(miomega);
        
    R2 point; Real_wp jacobian, poids; Complexe r2, r3, rtilde;
    TinyMatrix<Complexe, Symmetric, 3, 3> epsilon_3D, mu_3D;
    TinyMatrix<Complexe, General, 2, 2> A_tmp, eps_dfjm1, invMu_dfjm1;
    TinyMatrix<Complexe, Symmetric, 2, 2> Eps_Tmp, InvMu_tmp, invMu_tilde, InvMu_rev;
    bool affine = this->mesh.IsElementAffine(i);
    int i1 = i - this->mesh.GetNbElt() + this->GetNbEltPML();
    for (int j = 0; j < N; j++)
      {
	epsilon_3D = this->ref_epsilon(ref).GetCoefficient(var_leaf, i, j);
	mu_3D = this->ref_mu(ref).GetCoefficient(var_leaf, i, j);
	Complexe isigma_omega = this->ref_sigma(ref).GetCoefficient(var_leaf, i, j)/miomega;
	
	eps_rz(0,0) = epsilon_3D(0, 0) + isigma_omega; eps_rz(0, 1) = epsilon_3D(0, 2);
	eps_rz(1,0) = epsilon_3D(2, 0); eps_rz(1, 1) = epsilon_3D(2, 2) + isigma_omega;
	eps_teta = epsilon_3D(1, 1) + isigma_omega;
	
	mu_rz(1,0) = mu_3D(2, 0); mu_rz(1, 1) = mu_3D(2, 2);
	mu_rz(0,0) = mu_3D(0, 0); mu_rz(0, 1) = mu_3D(0, 2);
	mu_teta = mu_3D(1, 1);     
 
	point = this->Glob_PointsQuadrature(i)(j);
        if (affine)
	  {
	    dfjm1 = this->Glob_DFjm1(i)(0);
	    jacobian = this->Glob_jacobian(i)(0);
	    Mlt(1.0/jacobian, dfjm1);
	  }
	else
	  {
	    dfjm1 = this->Glob_DFjm1(i)(j);
	    jacobian = this->Glob_jacobian(i)(j) / Fb.WeightsND(j);
	    Mlt(1.0/jacobian, dfjm1);
	  }
	
        //radius = point(0);
	poids = Fb.WeightsND(j);
	rtilde = this->Glob_rtilde(i)(j);	
        Complexe dr; SetComplexOne(dr);
        Complexe dz; SetComplexOne(dz);
        
        // damping factors for PML layers
        if (this->InsidePML(i))
	  this->ModifyPMLCoefficient(eps_rz, eps_teta, mu_rz, mu_teta, dr, dz, i1, j);
	
	invMu_teta = 1.0/mu_teta;
	GetInverse(mu_rz, invMu_rz);
	
        if (this->modified_formulation)
          {
            Complexe coef;
	    
	    InvMu_rev(0, 0) = mu_rz(1, 1); InvMu_rev(0, 1) = -mu_rz(0, 1); InvMu_rev(1, 1) = mu_rz(0, 0);
	    GetInverse(InvMu_rev);	    

            MltTrans(eps_rz, dfjm1, eps_dfjm1); Mlt(dfjm1, eps_dfjm1, Eps_Tmp);
            MltTrans(invMu_rz, dfjm1, invMu_dfjm1); Mlt(dfjm1, invMu_dfjm1, InvMu_tmp);
            MltTrans(InvMu_rev, dfjm1, invMu_dfjm1); Mlt(dfjm1, invMu_dfjm1, InvMu_rev);
            r2 = rtilde*rtilde; r3 = r2*rtilde;
	    poids = Fb.WeightsND(j);
            
	    /* Terms needed for mode=0 */
	    
	    // term - omega^2 r epsilon E \cdot \phi
	    coef = -this->GetSquareOmega()*rtilde*jacobian*poids;
	    matmass_Erz(i)(j) = Eps_Tmp;
	    Mlt(coef, matmass_Erz(i)(j));
	    
	    // term r / mu rot E  rot phi
	    matmass_Hteta(i)(j) = rtilde*invMu_teta*poids/jacobian;
                
	    // term r / mu  grad(E_theta) \cdot grad(psi)
	    coef = jacobian*poids*rtilde;
	    matmass_Hrz(i)(j) = InvMu_rev;
	    matmass_Hrz(i)(j) *= coef;
	    
	    // term 1/mu (e_r E_\theta \cdot grad(psi) + e_r psi \cdot grad(E_\theta) )
	    vecstiff_Eteta(i)(j)(0) = dfjm1(0, 0)*invMu_rz(1, 1) - dfjm1(0, 1)*invMu_rz(0, 1);
	    vecstiff_Eteta(i)(j)(1) = dfjm1(1, 0)*invMu_rz(1, 1) - dfjm1(1, 1)*invMu_rz(0, 1);
	    vecstiff_Eteta(i)(j) *= jacobian*poids*dr;
	    
	    // term -omega^2 r \epsilon E_\theta psi
	    matmass_Eteta(i)(j) = -this->GetSquareOmega()*jacobian*eps_teta*rtilde*poids;
	    
	    // term 1/(r mu) E_\theta psi
	    matmass_Eteta_m0(i)(j) = jacobian*poids*invMu_rz(1, 1)/rtilde*(dr*dr);
	    
	    /* Terms needed for modes different from 0 */
	    
	    // term - omega^2 r^3 epsilon U \cdot V
	    coef = -this->GetSquareOmega()*r3*jacobian*poids;
	    matmass_Erz_r3(i)(j) = Eps_Tmp; Mlt(coef, matmass_Erz_r3(i)(j));
	    
	    // term m^2 r / mu U \cdot V
	    coef = jacobian*poids*rtilde;
	    matmass_Erz_m2(i)(j) = InvMu_rev; matmass_Erz_m2(i)(j) *= coef;
	    
	    // term - m^2  \omega^2 r \epsilon E_theta psi
	    matmass_Eteta(i)(j) = -this->GetSquareOmega()*jacobian*eps_teta*rtilde*poids;	    
	    
	    // term r^3 / mu rot U rot V
	    matmass_Hteta_r3(i)(j) = r3*poids*invMu_teta/jacobian;
	    
	    // term dr r^2 / mu (U \cdot e_z rot(V) + V \cdot e_z rot(U) )
	    vecstiff_Erz(i)(j)(0) = dr*r2*poids*invMu_teta*dfjm1(0, 1);
	    vecstiff_Erz(i)(j)(1) = dr*r2*poids*invMu_teta*dfjm1(1, 1);
            
	    // term  + dr^2 r/mu U \cdot e_z  V \cdot e_z
	    coef = dr*dr*rtilde*jacobian*poids*invMu_teta;
	    matmass_Erz_zz(i)(j)(0, 0) = dfjm1(0, 1)*dfjm1(0, 1)*coef;
	    matmass_Erz_zz(i)(j)(1, 1) = dfjm1(1, 1)*dfjm1(1, 1)*coef;
	    matmass_Erz_zz(i)(j)(1, 0) = dfjm1(0, 1)*dfjm1(1, 1)*coef;
	    
	    // term - omega^2 r^3 epsilon grad(E_\theta) grad(psi)
	    coef = -this->GetSquareOmega()*r3*jacobian*poids;
	    matmass_Hrz_r3(i)(j) = Eps_Tmp; Mlt(coef, matmass_Hrz_r3(i)(j));
	    
	    // term - omega^2 dr r^2 epsilon (E_\theta e_r \cdot grad(psi) + psi e_r \cdot grad(E_theta))
	    vecstiff_Eteta_r2(i)(j)(0) = -this->GetSquareOmega()*r2*dr*jacobian*poids
	      *(eps_rz(0, 0)*dfjm1(0, 0) + eps_rz(0, 1)*dfjm1(0, 1));
	    
	    vecstiff_Eteta_r2(i)(j)(1) = -this->GetSquareOmega()*r2*dr*jacobian*poids
	      *(eps_rz(0, 0)*dfjm1(1, 0) + eps_rz(0, 1)*dfjm1(1, 1));
	    
	    // term -omega^2 r dr^2 epsilon E_\theta psi
	    matmass_Eteta_rr(i)(j) = -this->GetSquareOmega()*dr*dr*eps_rz(0, 0)
	      *rtilde*jacobian*poids;
	    
	    // term -\omega^2 r^3 epsilon (U \cdot grad(psi) + V \cdot grad (E_\theta))
	    // already present in matmass_Hrz_r3
	    //coef = -this->GetSquareOmega()*r3*jacobian*poids;
	    //mat_coupling_E(i)(j) = Eps_Tmp; Mlt(coef, mat_coupling_E(i)(j));
	    
	    // term -\omega^2 r^2 dr epsilon (U \cdot e_r psi + V \cdot e_r E_theta)
	    // already present in vecstiff_Eteta_r2
	    //vec_coupling_E(i)(j)(0) = -this->GetSquareOmega()*r2*dr*jacobian
	    //  *poids*(eps_tilde(0, 0)*dfjm1(0, 0) + eps_tilde(0, 1)*dfjm1(0, 1));
	    
	    //vec_coupling_E(i)(j)(1) = -this->GetSquareOmega()*r2*dr*jacobian
	    //  *poids*(eps_tilde(0, 0)*dfjm1(1, 0) + eps_tilde(0, 1)*dfjm1(1, 1));
	    
          }
        else
          {
            // mass matrix related to E
            // matmass_Erz = -\omega \tilde{r} J_i DF_i^{-1} \epsilon_{r,z} (1/d 0; 0 d) DF_i^{-t} 
            MltTrans(eps_rz, dfjm1, A_tmp); Mlt(dfjm1, A_tmp, matmass_Erz(i)(j));
            Mlt(coef_mass*rtilde*Fb.WeightsND(j)*jacobian, matmass_Erz(i)(j));
            
            // mass matrix related to H
            // matmass_Hrz = -\omega \tilde{r} J_i DF_i^{-1} \mu_{r,z} (1/d 0; 0 d) DF_i^{-t} 
            MltTrans(mu_rz, dfjm1, A_tmp); Mlt(dfjm1, A_tmp, matmass_Hrz(i)(j));
            Mlt(-rtilde/coef_stiff*Fb.WeightsND(j)*jacobian, matmass_Hrz(i)(j));
            
            // mass matrix related to Eteta and H teta
            // matmass_Eteta = - \omega \epsilon r J_i
            // matmass_Hteta = - \omega \mu r J_i
            matmass_Eteta(i)(j) = coef_mass*eps_teta*rtilde*jacobian*Fb.WeightsND(j);
            matmass_Hteta(i)(j) = -1.0/coef_stiff*mu_teta*rtilde*jacobian*Fb.WeightsND(j);
            
            // rtilde
            matstiff_radius(i)(j) = rtilde*Fb.WeightsND(j);
            
            // term containing derivative of damping
            // a = \tilde{r}/d dd/dr J_i DF_i^{-1} e_z
            // d = d J_i DF_i^{-1} e_z
            matstiff_a10(i)(j) = dfjm1(0, 1)*jacobian*dr*Fb.WeightsND(j);
            matstiff_a11(i)(j) = dfjm1(1, 1)*jacobian*dr*Fb.WeightsND(j);
          }
      }
  }


  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::GetLocalUnknownVectorE(const Vector<Complexe>& U0, int iquad, const ElementReference<Dimension2, 2>& Fb_hcurl,
                           const ElementReference<Dimension2, 1>& Fb_h1, Vector<Complexe>& E, Vector<Complexe>& Eteta) const
  {
    int nb_dof_E = Fb_hcurl.GetNbDof();
    int nb_dof_Eteta = Fb_h1.GetNbDof();

    E.Zero(); Eteta.Zero();
    for (int j = 0; j < nb_dof_E; j++)
      {
        int num_dof = this->mesh_num.Element(iquad).GetNumberDof(j);
        if (num_dof >= 0)
          E(j) = U0(num_dof);
      }
    
    this->mesh_num.number_map.ModifyLocalComponentVector(this->mesh_num, E, iquad);
    
    for (int j = 0; j < nb_dof_Eteta; j++)
      {
        int num_dof = this->other_mesh_num(0)->Element(iquad).GetNumberDof(j);
        if (num_dof >= 0)
          Eteta(j) = U0(num_dof + this->offset_dof_unknown(1));
      }
    
    this->other_mesh_num(0)->number_map.ModifyLocalComponentVector(*this->other_mesh_num(0), Eteta, iquad);
  }


  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::GetLocalUnknownVectorE(const Vector<Vector<Complexe> >& U0, int iquad, const ElementReference<Dimension2, 2>& Fb_hcurl,
                           const ElementReference<Dimension2, 1>& Fb_h1, Vector<Complexe>& E, Vector<Complexe>& Eteta) const
  {
    int nb_dof_E = Fb_hcurl.GetNbDof();
    int nb_dof_Eteta = Fb_h1.GetNbDof();

    E.Zero(); Eteta.Zero();
    for (int j = 0; j < nb_dof_E; j++)
      {
        int num_dof = this->mesh_num.Element(iquad).GetNumberDof(j);
        if (num_dof >= 0)
          E(j) = U0(0)(num_dof);
      }
    
    this->mesh_num.number_map.ModifyLocalComponentVector(this->mesh_num, E, iquad);
    
    for (int j = 0; j < nb_dof_Eteta; j++)
      {
        int num_dof = this->other_mesh_num(0)->Element(iquad).GetNumberDof(j);
        if (num_dof >= 0)
          Eteta(j) = U0(1)(num_dof);
      }
    
    this->other_mesh_num(0)->number_map.ModifyLocalComponentVector(*this->other_mesh_num(0), Eteta, iquad);
  }


  //! computation of electromagnetical energy \int \epsilon E^2 + mu H^2
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::ComputeEnergy(const Vector<Complexe>& U0, VectReal_wp& normE, VectReal_wp& normH)
  {
    int nb_ref = this->ref_epsilon.GetM();
    normE.Reallocate(nb_ref);
    normH.Reallocate(nb_ref);
    normE.Fill(0); normH.Fill(0);
    for (int iquad = 0; iquad < this->mesh.GetNbElt(); iquad++)
      {
 	const ElementReference<Dimension2, 2>& Fb_hcurl = this->GetReferenceElementHcurl(iquad);
        const ElementReference<Dimension2, 1>& Fb_h1 = this->GetReferenceElementH1(iquad, 1);
        
        int nb_dof_E = Fb_hcurl.GetNbDof();
	int nb_dof_Eteta = Fb_h1.GetNbDof();
	int nb_dof_H = 2*Fb_h1.GetNbDof();
	int nb_dof_Hteta = Fb_hcurl.GetNbDofH();
	
	VectComplexe E(nb_dof_E), Eteta(nb_dof_Eteta), Hx(nb_dof_Eteta), Hy(nb_dof_Eteta), Hteta(nb_dof_Hteta);
	int N = Fb_h1.GetNbPointsQuadratureInside();
	VectComplexe Equad(2*N), Eteta_quad(N), Hx_quad(N), Hy_quad(N), Hteta_quad(N);
	int offset = this->offset_dof_unknown(2) + this->OffsetDofH(iquad);
	int ref = this->mesh.Element(iquad).GetReference();

        this->GetLocalUnknownVectorE(U0, iquad, Fb_hcurl, Fb_h1, E, Eteta);
        
	if (this->FirstOrderFormulation())
	  {
	    for (int i = 0; i < nb_dof_Eteta; i++)
	      {
                Hx(i) = U0(offset + 2*i);
                Hy(i) = U0(offset + 2*i+1);
              }
            
	    for (int i = 0; i < nb_dof_Hteta; i++)
	      Hteta(i) = U0(offset + nb_dof_H + i);
	    
	    Fb_h1.ApplyChTranspose(Hx, Hx_quad);
            Fb_h1.ApplyChTranspose(Hy, Hy_quad);
            Fb_hcurl.ApplyChTransposeH(Hteta, Hteta_quad);
	  }
	
	Fb_hcurl.ApplyChTranspose(E, Equad);
	Fb_h1.ApplyChTranspose(Eteta, Eteta_quad);
        
	TinyVector<Complexe, 2> vecE, vecH, vecHt, vecEt;
	Matrix2_2 dfjm1;
	for (int i = 0; i < N; i++)
	  {
	    if (this->Glob_DFjm1(iquad).GetM() == 1)
	      dfjm1 = this->Glob_DFjm1(iquad)(0);
	    else
	      dfjm1 = this->Glob_DFjm1(iquad)(i);
	    
	    vecE(0) = Equad(2*i); vecE(1) = Equad(2*i+1);
	    if (this->FirstOrderFormulation())
	      {
		vecH(0) = Hx_quad(i); vecH(1) = Hy_quad(i);
		MltTrans(dfjm1, vecH, vecHt);
		Real_wp Hsquare = absSquare(vecHt(0)) + absSquare(vecHt(1))
                  + absSquare(Hteta_quad(i));
		normH(ref) += Fb_h1.WeightsND(i)*abs(this->Glob_rtilde(iquad)(i))*Hsquare;
	      }
            
	    MltTrans(dfjm1, vecE, vecEt);
	    
	    Real_wp Esquare = absSquare(vecEt(0)) + absSquare(vecEt(1))
              + absSquare(Eteta_quad(i));
	    
	    normE(ref) += Fb_h1.WeightsND(i)*abs(this->Glob_rtilde(iquad)(i))*Esquare;	  
	  }
      }
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(this->comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(this->comm_group_mode, &nb_proc);
    if (nb_proc > 1)
      {
        VectReal_wp normEloc(normE), normHloc(normH);
        Vector<int64_t> xtmp;
        MpiReduce(this->comm_group_mode, normEloc, xtmp, normE, normE.GetM(),
                  MPI_SUM, 0);
	
        MpiReduce(this->comm_group_mode, normHloc, xtmp, normH, normH.GetM(),
                  MPI_SUM, 0);
      }
#else
    int rank_proc(0);
#endif    
    
    if (rank_proc == 0)
      for (int ref = 0; ref < nb_ref; ref++)
        {
          normE(ref) = sqrt(normE(ref));
          normH(ref) = sqrt(normH(ref));
        }    
  }
  
    
  //! Computation of E(point_glob) from values on degrees of freedom
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeEpolar(const Vector<VectComplexe>& U0, int num_elem, const R2& point_loc, const R2& pt_glob,
		const Matrix2_2& dfjm1, TinyVector<Complexe, 3>& val_E,
                const ElementReference<Dimension2, 2>& Fb_hcurl,
                const ElementReference<Dimension2, 1>& Fb_h1) const
  {
    TinyVector<Complexe, 2> Erz; Complexe Eteta; R2 tmp;
    VectR2 val_phi_hcurl(Fb_hcurl.GetNbDof());
    Vector<TinyVector<Real_wp, 1> > val_psi_h1(Fb_h1.GetNbDof());
    
    Erz.Fill(0); Eteta = 0.0;
    // we compute phi
    Fb_hcurl.ComputeValuesPhi(point_loc, val_phi_hcurl, dfjm1, this->mesh_num, num_elem);
    Fb_h1.ComputeValuesPhi(point_loc, val_psi_h1, dfjm1, *this->other_mesh_num(0), num_elem);
    
    if ((this->modified_formulation) && (this->number_mode != 0))
      {
	// and curl phi
        VectR2 grad_psi(Fb_h1.GetNbDof());
        Fb_h1.ComputeValuesGradientPhi(point_loc, grad_psi, dfjm1, *this->other_mesh_num(0), num_elem);
        
        Real_wp radius = pt_glob(0);
        for (int j = 0; j < Fb_hcurl.GetNbDof(); j++)
          {
	    int num_dof = this->mesh_num.Element(num_elem).GetNumberDof(j);
            Erz(0) += radius*val_phi_hcurl(j)(0)*U0(0)(num_dof);
            Erz(1) += radius*val_phi_hcurl(j)(1)*U0(0)(num_dof);
          }

        for (int j = 0; j < Fb_h1.GetNbDof(); j++)
          {
            int num_dof = this->other_mesh_num(0)->Element(num_elem).GetNumberDof(j);
            Complexe uloc = U0(1)(num_dof);
            Eteta += val_psi_h1(j)(0)*uloc;
            Erz(0) += (radius*grad_psi(j)(0) + val_psi_h1(j)(0))*uloc;
            Erz(1) += radius*grad_psi(j)(1)*uloc;
          }
        
        Mlt(Real_wp(1.0/abs(this->number_mode)), Erz);
      }
    else
      {
        for (int j = 0; j < Fb_hcurl.GetNbDof(); j++)
          {
	    int num_dof = this->mesh_num.Element(num_elem).GetNumberDof(j);
            Erz(0) += val_phi_hcurl(j)(0)*U0(0)(num_dof);
            Erz(1) += val_phi_hcurl(j)(1)*U0(0)(num_dof);
          }

        for (int j = 0; j < Fb_h1.GetNbDof(); j++)
          {
            int num_dof = this->other_mesh_num(0)->Element(num_elem).GetNumberDof(j);
            Complexe uloc = U0(1)(num_dof);
            Eteta += val_psi_h1(j)(0)*uloc;
          }
      }
    
    Complexe Icplx;
    SetImaginaryI(Icplx);
    
    if (this->number_mode >= 0)
      Eteta *= -Icplx;
    else
      Eteta *= Icplx;
    
    val_E(0) = Erz(0);
    val_E(1) = Eteta;
    val_E(2) = Erz(1);
  }
  
  
  //! Computation of H(point_glob) from values on degrees of freedom
  template<class TypeEquation> 
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeHpolar(const Vector<VectComplexe>& U0, int num_elem, const R2& point_loc,
		const R2& point_glob, const Matrix2_2& dfjm1,
                TinyVector<Complexe, 3>& val_H,
                const ElementReference<Dimension2, 2>& Fb_hcurl,
                const ElementReference<Dimension2, 1>& Fb_h1) const
  {
    TinyVector<Complexe, 2> Hrz; Complexe Hteta; R2 tmp;
    int num_dof;
    Hrz.Fill(0); Hteta = 0.0;

    Complexe Icplx;
    SetImaginaryI(Icplx);

    if (this->modified_formulation)
      {
	int ref = this->mesh.Element(num_elem).GetReference();
	VectR2 val_phi_hcurl(Fb_hcurl.GetNbDof()), grad_psi(Fb_h1.GetNbDof());
        Vector<TinyVector<Real_wp, 1> > val_psi_h1(Fb_h1.GetNbDof()), curl_phi(Fb_hcurl.GetNbDof());
        
	Fb_h1.ComputeValuesPhi(point_loc, val_psi_h1, dfjm1, *this->other_mesh_num(0), num_elem);
        Fb_h1.ComputeValuesGradientPhi(point_loc, grad_psi, dfjm1, *this->other_mesh_num(0), num_elem);
        Fb_hcurl.ComputeValuesPhi(point_loc, val_phi_hcurl, dfjm1, this->mesh_num, num_elem);
	Fb_hcurl.ComputeValuesGradientPhi(point_loc, curl_phi, dfjm1, this->mesh_num, num_elem);
	
	Real_wp r = point_glob(0);
	TinyVector<Complexe, 3> rotE, vecH;
	if (this->number_mode == 0)
	  {
            for (int j = 0; j < Fb_h1.GetNbDof(); j++)
              {
                num_dof = this->other_mesh_num(0)->Element(num_elem).GetNumberDof(j);
                rotE(0) += grad_psi(j)(1)*U0(1)(num_dof);
                rotE(2) -= (val_psi_h1(j)(0)/r + grad_psi(j)(0))*U0(1)(num_dof);
              }

            for (int j = 0; j < Fb_hcurl.GetNbDof(); j++)
              {
                num_dof = this->mesh_num.Element(num_elem).GetNumberDof(j);
                rotE(1) += curl_phi(j)(0)*U0(0)(num_dof);
              }                
          }
	else
	  {
            for (int j = 0; j < Fb_hcurl.GetNbDof(); j++)
              {
                num_dof = this->mesh_num.Element(num_elem).GetNumberDof(j);
                rotE(0) -= val_phi_hcurl(j)(1)*U0(0)(num_dof);
                rotE(1) += (val_phi_hcurl(j)(1) + r*curl_phi(j)(0))*U0(0)(num_dof);
                rotE(2) += val_phi_hcurl(j)(0)*U0(0)(num_dof);
              }                
            
	    rotE(1) /= Real_wp(abs(this->number_mode));
	  }
        
	if (this->ref_mu(ref).IsVarying())
	  {
	    TinyMatrix<Complexe, Symmetric, 3, 3> mu;
	    TinyMatrix<Complexe, Symmetric, 2, 2> invMu_rz;
	    
	    VectReal_wp phi;
	    this->mesh.ComputeValuesPhiNodalRef(num_elem, point_loc, phi);
	    mu = this->ref_mu(ref).GetMatrixPoint(this->ElementRho(num_elem), phi, point_glob);
	    
	    invMu_rz(0, 0) = mu(0, 0); invMu_rz(0, 1) = mu(0, 2); invMu_rz(1, 1) = mu(2, 2);
	    GetInverse(invMu_rz);
	    
	    vecH(0) = invMu_rz(0, 0)*rotE(0) + invMu_rz(0, 1)*rotE(2);
	    vecH(2) = invMu_rz(0, 1)*rotE(0) + invMu_rz(1, 1)*rotE(2);
	    vecH(1) = rotE(1)/mu(1, 1);
	  }
	else
	  this->ref_invMu(ref).MltMatrix(this->GetLeafClass(), num_elem, 0, rotE, vecH);
	
	Complexe coef = 1.0/this->GetOmega();
	Hteta = Icplx*coef*vecH(1);
	Hrz(0) = coef*vecH(0);
	Hrz(1) = coef*vecH(2);
      }
    else
      {
	VectReal_wp val_phi(Fb_h1.GetNbDof());
	VectReal_wp scal_phi(Fb_hcurl.GetNbDofH());
	
	// we compute the vectorial values of basis fonctions at point_loc
	Fb_h1.ComputeValuesPhiRef(point_loc, val_phi);
	
	for (int j = 0; j < Fb_h1.GetNbDof(); j++)
	  {
	    // we add the contribution DFi*-1 phi_j
	    num_dof = OffsetDofH(num_elem) + 2*j;
            Hrz(0) += val_phi(j)*(U0(0)(num_dof)*dfjm1(0, 0) + U0(0)(num_dof+1)*dfjm1(1, 0));
            Hrz(1) += val_phi(j)*(U0(0)(num_dof)*dfjm1(0, 1) + U0(0)(num_dof+1)*dfjm1(1, 1));
	  }
        
	// we compute the scalar values of basis fonctions at point_loc
	Fb_hcurl.ComputeValuesPhiHRef(point_loc, scal_phi);
	
	for (int j = 0; j < Fb_hcurl.GetNbDofH(); j++)
	  {
	    num_dof = OffsetDofH(num_elem) + 2*Fb_h1.GetNbDof() + j;	    
	    Hteta += scal_phi(j)*U0(0)(num_dof);
	  }
	
        Hteta *= -Icplx;
      }
    
    if (this->number_mode < 0)
      Hrz *= -1.0;    
    
    val_H(0) = Hrz(0);
    val_H(1) = Hteta;
    val_H(2) = Hrz(1);
  }
  

  template<class TypeEquation> template<class T>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeInterpolationUlocGen(const Vector<Vector<T> > & U0,
                              const GridInterpolation<Dimension2> & var_interp,
                              Vector<T>& trace_vec, Vector<T>& trace_grad_vec,
                              const IVect& list_points, int nnz, bool compute_grad) const
  {
    cout << "Not implemented" << endl;
    abort();
  }


  //! computation of E on a predefined grid
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeInterpolationUlocGen(const Vector<VectComplexe> & U0,
                              const GridInterpolation<Dimension2> & var_interp,
                              Vector<Complexe>& trace_vec, Vector<Complexe>& trace_grad_vec,
                              const IVect& list_points, int nnz, bool compute_grad) const
  {   
    VectR2 s; R2 point_loc;
    TinyVector<Complexe, 3> val_E, val_H;
    int iquad; Real_wp teta;
    
    Matrix2_2 dfjm1, dfj;
    
    trace_vec.Reallocate(nnz*3);
    trace_vec.Fill(0);
    Vector<VectComplexe> EvalH_Nodal;
    if (compute_grad)
      {
	trace_grad_vec.Reallocate(nnz*3);
	trace_grad_vec.Fill(0);

	if (!modified_formulation)
	  this->EvaluateH_MixedFormulation(U0, EvalH_Nodal);
	else
	  EvalH_Nodal.SetData(U0.GetM(), U0.GetData());	
      }
    
    int nb = 0;
    // Loop on all the points of the grid
    for (int i1 = 0; i1 < list_points.GetM(); i1++)
      {
	int i = list_points(i1);
	val_E.Zero(); val_H.Zero();
	if (i >= 0)
	  {
	    // iquad is the element where the point i is
	    iquad = var_interp.GetElementNumber(i);
	    // local coordinates in the element
	    point_loc = var_interp.GetLocalCoordinate(i);
            R2 point_glob = var_interp.GetGlobalCoordinate(i);
	    // if the point has been found in the grid
	    if ((iquad >= 0)&&(iquad < this->mesh.GetNbElt()))
	      {
		// we get the vertices of the element
		dfjm1 = var_interp.GetDFjm1(i);

		this->ComputeEpolar(U0, iquad, point_loc, point_glob,
				    dfjm1, val_E, this->GetReferenceElementHcurl(iquad),
                                    this->GetReferenceElementH1(iquad, 1));
		
		if (compute_grad)
		  this->ComputeHpolar(EvalH_Nodal, iquad, point_loc, point_glob,
				      dfjm1, val_H, this->GetReferenceElementHcurl(iquad),
                                      this->GetReferenceElementH1(iquad, 1));
                
		teta = var_interp.GetTheta(i);
		Complexe Ex = val_E(0);
		val_E(0) = val_E(0)*cos(teta) - val_E(1)*sin(teta);
		val_E(1) = Ex*sin(teta) + val_E(1)*cos(teta);

		Complexe exp_mteta;
 		this->GetFourierMode(teta, exp_mteta);
		Mlt(exp_mteta, val_E); 
		
		if (compute_grad)
		  {
		    Complexe Hx = val_H(0);
		    val_H(0) = val_H(0)*cos(teta) - val_H(1)*sin(teta);
		    val_H(1) = Hx*sin(teta) + val_H(1)*cos(teta);
		    Mlt(exp_mteta, val_H); 
		  }
              }
	    
	    trace_vec(nb) = val_E(0);
	    trace_vec(nb+nnz) = val_E(1);
	    trace_vec(nb+2*nnz) = val_E(2);
	    if (compute_grad)
	      {		
		trace_grad_vec(nb) = val_H(0);
		trace_grad_vec(nb+nnz) = val_H(1);
		trace_grad_vec(nb+2*nnz) = val_H(2);
	      }
	    
	    nb++;
	  }	
      }

    if (modified_formulation)
      EvalH_Nodal.Nullify();
  }
  


  template<class TypeEquation> 
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::ComputeInterpolationUloc(const Vector<VectReal_wp> & U0,
                             const GridInterpolation<Dimension2> & var_interp,
                             Vector<Real_wp>& trace_vec, Vector<Real_wp>& trace_grad_vec,
                             const IVect& list_points, int nnz, bool compute_grad) const
  {
    ComputeInterpolationUlocGen(U0, var_interp, trace_vec, trace_grad_vec, list_points, nnz, compute_grad);
  }
  
    
  template<class TypeEquation> 
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::ComputeInterpolationUloc(const Vector<VectComplex_wp> & U0,
                             const GridInterpolation<Dimension2> & var_interp,
                             Vector<Complex_wp>& trace_vec, Vector<Complex_wp>& trace_grad_vec,
                             const IVect& list_points, int nnz, bool compute_grad) const
  {
    ComputeInterpolationUlocGen(U0, var_interp, trace_vec, trace_grad_vec, list_points, nnz, compute_grad);
  }

  template<class TypeEquation> template<class T>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  AddVolumeSourceGen(const T& alpha, Vector<Vector<T> > & b_source,
                     Vector<VirtualSourceFEM<T, Dimension2>* >& f) const
  {
    cout << "Not implemented for this type of T" << endl;
    abort();
  }

  
  //! adding the source term \int f \varphi to a vector
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  AddVolumeSourceGen(const Complexe& alpha, Vector<VectComplexe> & b_source,
                     Vector<VirtualSourceFEM<Complexe, Dimension2>* >& f) const
  {
    VirtualMaxwellAxiSourceFEM<Complexe>& f_axi = dynamic_cast<VirtualMaxwellAxiSourceFEM<Complexe>& >(*f(0));    
    VectR2 s;

    SetPoints<Dimension2> PointsElem;
    SetMatrices<Dimension2> MatricesElem;
    Vector<Complexe> f_eval(3), g_eval(3);
    Vector<Vector<Complexe> > gradf_eval(1), curlf_eval(1);
    Vector<Vector<Complexe> > feval_h1(1), feval_hcurl(1);
    Vector<Vector<Complexe> > contrib_h1(1), contrib_hcurl(1);
    
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	// we get the vertices of the face
	this->mesh.GetVerticesElement(i, s);

        const ElementReference<Dimension2, 2>& Fb_hcurl = this->GetReferenceElementHcurl(i);
	const ElementReference<Dimension2, 1>& Fb_h1 = this->GetReferenceElementH1(i, 1);
	
	// number of dofs in the face
	int nb_dof_scal = Fb_h1.GetNbDof();
        int nb_dof_vec = Fb_hcurl.GetNbDof();
	int nb_points_quadrature = Fb_h1.GetNbPointsQuadratureInside();	
	
	f_axi.InitElement(i, s);
	
	// we compute the transformation Fi on quadrature points
	Fb_h1.FjElemQuadrature(s, PointsElem, this->mesh, i);
	Fb_h1.DFjElemQuadrature(s, PointsElem, MatricesElem, this->mesh, i);
	
	f_eval.Fill(0); g_eval.Fill(0);
	feval_hcurl(0).Reallocate(2*nb_points_quadrature);
        feval_h1(0).Reallocate(nb_points_quadrature);
	if (f_axi.IsNonNullVolumetricSource(s))
	  {
	    // \int_{K_i} f \phi_i = \int_\hat{K} J_i f \phi_i
	    for (int j = 0; j < nb_points_quadrature; j++)
	      {
		f_axi.EvaluateVolumetricSource(i, j, PointsElem.
					       GetPointQuadrature(j), f_eval);
		
		feval_hcurl(0)(2*j) = f_eval(0);
		feval_hcurl(0)(2*j+1) = f_eval(1);
		feval_h1(0)(j) = f_eval(2);
	      }
            
	    // we compute the integrals int( f phi_j dx ) j=0..nb_dof_elt-1 
	    Fb_h1.ComputeIntegral(MatricesElem, feval_h1, contrib_h1,
                                  *this->other_mesh_num(0), i);

            Fb_hcurl.ComputeIntegral(MatricesElem, feval_hcurl, contrib_hcurl,
                                     this->mesh_num, i);

	    // we add contributions to b_source
	    for (int j = 0; j < nb_dof_vec; j++)
	      {
		int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
		b_source(0)(num_dof) += alpha*contrib_hcurl(0)(j);
	      }

            for (int j = 0; j < nb_dof_scal; j++)
	      {
		int num_dof = this->other_mesh_num(0)->Element(i).GetNumberDof(j);
                num_dof += this->offset_dof_unknown(1);
		b_source(0)(num_dof) += alpha*contrib_h1(0)(j);
	      }
	  }
        
	gradf_eval(0).Reallocate(2*nb_points_quadrature);
        curlf_eval(0).Reallocate(nb_points_quadrature);
	if ((f_axi.IsNonNullGradientSource(s)) && (this->modified_formulation))
	  {
	    // \int_{K_i} f \phi_i = \int_\hat{K} J_i f \phi_i
	    for (int j = 0; j < nb_points_quadrature; j++)
	      {
		f_axi.EvaluateGradientSource(i, j, PointsElem.
					     GetPointQuadrature(j), g_eval);
		
		gradf_eval(0)(2*j) = -g_eval(1);
                gradf_eval(0)(2*j+1) = g_eval(0);
		curlf_eval(0)(j) = g_eval(2);
	      }
            
	    // we compute the integrals int( f phi_j dx ) j=0..nb_dof_elt-1 
	    Fb_h1.ComputeIntegralGradient(MatricesElem, gradf_eval, contrib_h1,
                                          *this->other_mesh_num(0), i);

            Fb_hcurl.ComputeIntegralGradient(MatricesElem, curlf_eval, contrib_hcurl,
                                             this->mesh_num, i);
            
	    // we add contributions to b_source
            for (int j = 0; j < nb_dof_vec; j++)
              {
                int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
		b_source(0)(num_dof) += alpha*contrib_hcurl(0)(j);
	      }
            
            for (int j = 0; j < nb_dof_scal; j++)
	      {
		int num_dof = this->other_mesh_num(0)->Element(i).GetNumberDof(j);
                num_dof += this->offset_dof_unknown(1);
		b_source(0)(num_dof) += alpha*contrib_h1(0)(j);
	      }
	  }
        
	if (f_axi.IsNonNull_SourceH(s))
	  {
            Vector<Complexe> feval_hx(nb_points_quadrature);
            Vector<Complexe> feval_hz(nb_points_quadrature);
            Vector<Complexe> feval_hteta(nb_points_quadrature);
	    Vector<Complexe> contrib_Hx(Fb_h1.GetNbDof()), contrib_Hz(Fb_h1.GetNbDof());
            Vector<Complexe> contrib_Hteta(Fb_hcurl.GetNbDofH());
	    
	    // source in H \int_{K_i} r \, f \cdot \varphi 
	    //    \; = \; \int_\hat{K} r \, J_i DF_i^{-1} f \cdot \hat{\varphi} 
	    // source in H_\theta \int_{K_i} r \, f \psi \; = \; \int_\hat{K} r \, J_i
	    bool zero_source = true;
            TinyVector<Complexe, 2> phi, fvec; Matrix2_2 dfjm1;
	    for (int j = 0; j < nb_points_quadrature; j++)
	      {
		f_axi.EvaluateSource_H(i, j, PointsElem.GetPointQuadrature(j), f_eval);

                phi(0) = f_eval(0); phi(1) = f_eval(1);
                GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);
                Mlt(dfjm1, phi, fvec);
                Real_wp radius = PointsElem.GetPointQuadrature(j)(0);
                Real_wp poids = radius*Fb_h1.WeightsND(j)*Det(MatricesElem.GetPointQuadrature(j));
                feval_hx(j) = fvec(0)*poids;
                feval_hz(j) = fvec(1)*poids;

                feval_hteta(j) = f_eval(2)*poids;
		if (Norm2(f_eval) > epsilon_machine)
		  zero_source = false;
	      }

	    if (!zero_source)
	      {
                Fb_h1.ApplyCh(feval_hx, contrib_Hx);
                Fb_h1.ApplyCh(feval_hz, contrib_Hz);
                Fb_hcurl.ApplyChH(feval_hteta, contrib_Hteta);

		this->ComputeSourceH_ComplementSchur(contrib_Hx, contrib_Hz,
                                                     contrib_Hteta, contrib_hcurl(0),
                                                     contrib_h1(0), Fb_hcurl, Fb_h1, i);

                for (int j = 0; j < nb_dof_vec; j++)
                  {
                    int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
                    b_source(0)(num_dof) += alpha*contrib_hcurl(0)(j);
                  }
                
                for (int j = 0; j < nb_dof_scal; j++)
                  {
                    int num_dof = this->other_mesh_num(0)->Element(i).GetNumberDof(j);
                    num_dof += this->offset_dof_unknown(1);
                    b_source(0)(num_dof) += alpha*contrib_h1(0)(j);
                  }
              }
          }
      }
    //b_source(0).Write("source.dat"); exit(0);
  }
  

  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::AddVolumeSource(const Real_wp& alpha, Vector<VectReal_wp> & b_source,
                    Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const
  {
    AddVolumeSourceGen(alpha, b_source, f);
  }
  
  
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  AddVolumeSource(const Complex_wp& alpha, Vector<VectComplex_wp> & b_source,
                  Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const
  {
    AddVolumeSourceGen(alpha, b_source, f);
  }


  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::AddDiracSource(const Real_wp& alpha, Vector<Vector<Real_wp> > & b_source,
		   Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const
  {
    abort();
  }

  
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>
  ::AddDiracSource(const Complex_wp& alpha, Vector<Vector<Complex_wp> > & b_source,
		   Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const
  {
    VectR3 pt_source3D(1);
    pt_source3D(0) = this->origine_phase3D;

    if (f.GetM() != 1)
      {
        cout << "Only one source for Dirac" << endl;
        abort();
      }

    int m = this->number_mode;
    if (f(0)->IsDiracSource())
      {	        
        R3 pt = pt_source3D(0);

        // conversion in cylindrical coordinates of the origin
        Real_wp r0, z0(pt(2)), theta0;  //cos_theta0, sin_theta0;
        CartesianToPolar(pt(0), pt(1), r0, theta0);
        //cos_theta0 = cos(theta0); sin_theta0 = sin(theta0);

        // checking if the point belongs to the axis Oz
        bool point_on_axis = false;
        if (r0 <= epsilon_machine)
          point_on_axis = true;

        // beta = 1/(2 pi) if the point is on the axis, exp(im theta0) / (2pi) otherwise
        Complex_wp beta = 1.0/(2.0*pi_wp);
        if (!point_on_axis)
          to_complex(beta*exp(Iwp*Real_wp(m)*theta0), beta);

        Vector<Complex_wp> polar = f(0)->GetPolarization();
        Mlt(beta, polar);
        
        if (this->modified_formulation)
          {
            cout << "Not implemented" << endl;
            abort();
          }

        R2 pt2D(r0, z0); VectR2 val_phi; Vector<TinyVector<Real_wp, 1> > val_psi; int nb_eval = 0;
        Vector<Complex_wp> b_add(this->mesh_num.GetNbDof());
        TinyVector<R2, 2> enveloppe;
        b_add.Zero();
        for (int num_elem = 0; num_elem < this->mesh.GetNbElt(); num_elem++)
          {
            const ElementReference<Dimension2, 1>& Fb_h1 = this->GetReferenceElementH1(num_elem, 1);
            const ElementReference<Dimension2, 2>& Fb_hcurl = this->GetReferenceElementHcurl(num_elem);
            
            // initialization of F_i^-1
            FjInverseProblem<Dimension2> inverseFj(this->mesh, num_elem);
            
            SetPoints<Dimension2>& PointsElem = inverseFj.GetSetPoints();
            VectR2& s = inverseFj.GetVertices();
            
            this->mesh.GetBoundingBox(num_elem, s, PointsElem, enveloppe);
            IVect Nodle = this->GetDofNumberOnElement(num_elem);
            
            // localizing the point
            if (PointInsideBoundingBox(pt2D, enveloppe))
              {
                R2 pt_loc;
                bool pt_inside = inverseFj.Solve(pt2D, pt_loc);
                Real_wp distance_boundary_elt
                  = this->mesh.GetDistanceToBoundary(pt_loc, num_elem);
                
                if (distance_boundary_elt < 0)
                  this->mesh.ProjectPointOnBoundary(pt_loc, num_elem);
                
                Matrix2_2 dfj, dfjm1;
                if (pt_inside)
                  {
                    Fb_h1.DFj(s, PointsElem, pt_loc, dfj, this->mesh, num_elem);
		  
                    GetInverse(dfj, dfjm1);
                    // Fb_h1.ComputeValuesPhi(pt_loc, val_psi, dfjm1, *this->other_mesh_num(0), num_elem);
                    Fb_hcurl.ComputeValuesPhi(pt_loc, val_phi, dfjm1, this->mesh_num, num_elem);
		  
                    for (int i = 0; i < Fb_hcurl.GetNbDof(); i++)
                      {
                        int num_dof = Nodle(i);
                        if (num_dof >= 0)
                          {
                            if (m == 0)
                              b_add(num_dof) += polar(2)*val_phi(i)(1);
                          }
                      }
		  
                    nb_eval++;
                  }
              }
          }
    
        int nb_eval_all(nb_eval);
#ifdef SELDON_WITH_MPI
        MPI_Allreduce(&nb_eval, &nb_eval_all,
                      1, MPI_INTEGER, MPI_SUM, this->comm_group_mode);
#endif    
        
        if (nb_eval_all > 1)
          {
            Real_wp coef = Real_wp(1)/nb_eval_all;
            Mlt(coef, b_add);
          }
        
        Complex_wp one; SetComplexOne(one);
        Add(one, b_add, b_source(0));
      }    
  }

  //! computation of sources if using the mixed formulation
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeSourceH_ComplementSchur(Vector<Complexe>& contrib_Hx, Vector<Complexe>& contrib_Hz,
                                 Vector<Complexe>& contrib_Hteta, Vector<Complexe>& contrib_hcurl,
				 Vector<Complexe>& contrib_h1, const ElementReference<Dimension2, 2>& Fb_hcurl,
                                 const ElementReference<Dimension2, 1>& Fb_h1, int num_elem) const
  {
    //Complexe mode_m = abs(TypeEquation::number_mode)*TypeEquation::GetCoefMode(var_leaf);
    
    contrib_hcurl.Reallocate(Fb_hcurl.GetNbDof());
    contrib_h1.Reallocate(Fb_h1.GetNbDof());

    // complement schur to have sources in E and E_\theta
    // in E, we have -C_h \, (D_h^1)^{-1} F_h - R_h^2 (D_h^2)^{-1} F_h_\theta
    // in E_\theta, we have - R_h^2 (D_h^1)^{-1} F_h
    Matrix<Complexe> Ch, Rh1, Rh2, Dh1, Dh2, Bh1, Bh2;
    this->ComputeLocalMatrices_MixedFormulation(num_elem, Ch, Rh1, Rh2, Dh1, Dh2, Bh1, Bh2, Fb_hcurl, Fb_h1);
    
    IVect pivot1(Dh1.GetM()), pivot2(Dh2.GetM());
    Vector<Complexe> contrib_H(Dh1.GetM());
    for (int i = 0; i < contrib_Hx.GetM(); i++)
      {
        contrib_H(2*i) = contrib_Hx(i);
        contrib_H(2*i+1) = contrib_Hz(i);
      }
    
    Seldon::GetLU(Dh1, pivot1); Seldon::GetLU(Dh2, pivot2);
    Seldon::SolveLU(Dh1, pivot1, contrib_H);
    Seldon::SolveLU(Dh2, pivot2, contrib_Hteta);
    
    contrib_hcurl.Fill(0);
    MltAdd(Complexe(-1), Ch, contrib_H, Complexe(1), contrib_hcurl);
    MltAdd(Complexe(-1), Rh1, contrib_Hteta, Complexe(1), contrib_hcurl);
    contrib_h1.Fill(0);
    MltAdd(Complexe(-1), Rh2, contrib_H, Complexe(1), contrib_h1);
    
    this->mesh_num.number_map.
      ModifyLocalUnknownVector(this->mesh_num, contrib_hcurl, num_elem);

    this->other_mesh_num(0)->number_map.
      ModifyLocalUnknownVector(*this->other_mesh_num(0), contrib_hcurl, num_elem);
  }
  

  //! computation of matrices of mixed formulation, necessary for elements near the axis
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  ComputeLocalMatrices_MixedFormulation(int iquad, Matrix<Complexe>& Ch, Matrix<Complexe>& Rh1,
                                        Matrix<Complexe>& Rh2, Matrix<Complexe>& Dh1,
                                        Matrix<Complexe>& Dh2, Matrix<Complexe>& Bh1,
                                        Matrix<Complexe>& Bh2,
					const ElementReference<Dimension2, 2>& Fb_hcurl,
                                        const ElementReference<Dimension2, 1>& Fb_h1) const
  {
    //Complexe val, vloc1, vloc2, zero(0);
    // we store local matrices of mixed formulation, ie C_h, R_h^1, R_h^2, D_h^1 and D_h^2
    int nb_dof_E = Fb_hcurl.GetNbDof();
    int nb_dof_Eteta = Fb_h1.GetNbDof();
    int nb_dof_H = 2*Fb_h1.GetNbDof();
    int nb_dof_Hteta = Fb_hcurl.GetNbDofH();
    int nb_pts_quad = Fb_h1.GetNbPointsQuadratureInside();
    
    Ch.Reallocate(nb_dof_E, nb_dof_H);
    Rh1.Reallocate(nb_dof_E, nb_dof_Hteta);
    Rh2.Reallocate(nb_dof_Eteta, nb_dof_H);
    Dh1.Reallocate(nb_dof_H, nb_dof_H);
    Dh2.Reallocate(nb_dof_Hteta, nb_dof_Hteta);
    Bh1.Reallocate(nb_dof_E, nb_dof_E);
    Bh2.Reallocate(nb_dof_Eteta, nb_dof_Eteta);
    
    Ch.Fill(0); Rh1.Fill(0); Rh2.Fill(0);
    Dh1.Fill(0); Dh2.Fill(0);
    Bh1.Fill(0); Bh2.Fill(0);
    
    Complexe mode_m = abs(this->number_mode)*TypeEquation::GetCoefMode(this->GetLeafClass());
    TinyVector<Complexe, 2> vec_u, vec_v;
    
    VectReal_wp Ones(nb_dof_H/2), val_Hx(nb_pts_quad), val_Hy(nb_pts_quad);
    VectComplexe fevalCh(2*nb_pts_quad), fevalDh_x(nb_pts_quad), fevalDh_y(nb_pts_quad);
    VectComplexe fevalRh2_psi(nb_pts_quad), fevalRh2_gradPsi(2*nb_pts_quad);
    VectComplexe contribCh(nb_dof_E), contribDh_x(nb_dof_H/2), contribDh_y(nb_dof_H/2);
    VectComplexe contribRh2(nb_dof_Eteta), contribRh2_bis(nb_dof_Eteta);
    contribCh.Fill(0); contribDh_x.Fill(0); contribRh2.Fill(0); contribRh2_bis.Fill(0);
    fevalRh2_psi.Fill(0); fevalRh2_gradPsi.Fill(0);
    
    // loop over dofs of H
    for (int i = 0; i < nb_dof_H; i++)
      {
        Ones.Fill(0); Ones(i/2) = 1.0;
        if (i%2 == 0)
          {
            Fb_h1.ApplyChTranspose(Ones, val_Hx);
            val_Hy.Zero();
          }
        else
          {
            Fb_h1.ApplyChTranspose(Ones, val_Hy);
            val_Hx.Zero();
          }
        
        for (int j = 0; j < nb_pts_quad; j++)
          {
            // Ch part
            fevalCh(2*j) = val_Hy(j)*mode_m*Fb_h1.WeightsND(j);
            fevalCh(2*j+1) = -val_Hx(j)*mode_m*Fb_h1.WeightsND(j);
            
            // Rh^2 part
            fevalRh2_psi(j) = -matstiff_a10(iquad)(j)*val_Hx(j) 
              - matstiff_a11(iquad)(j)*val_Hy(j);
            
            fevalRh2_gradPsi(2*j) = -matstiff_radius(iquad)(j)*val_Hy(j);
            fevalRh2_gradPsi(2*j+1) = matstiff_radius(iquad)(j)*val_Hx(j);
            
            // Dh^1 part
            vec_u(0) = val_Hx(j); vec_u(1) = val_Hy(j);
            Mlt(matmass_Hrz(iquad)(j), vec_u, vec_v);
            fevalDh_x(j) = vec_v(0);
            fevalDh_y(j) = vec_v(1);
          }
        
        Fb_hcurl.ApplyCh(fevalCh, contribCh);
        Fb_h1.ApplyCh(fevalDh_x, contribDh_x);
        Fb_h1.ApplyCh(fevalDh_y, contribDh_y);
        
        Fb_h1.ApplyCh(fevalRh2_psi, contribRh2);
        Fb_h1.ApplyRh(fevalRh2_gradPsi, contribRh2_bis);
        
        for (int j = 0; j < nb_dof_E; j++)
          Ch(j, i) += contribCh(j);
        
        for (int j = 0; j < nb_dof_H/2; j++)
          {
            Dh1(2*j, i) += contribDh_x(j);
            Dh1(2*j+1, i) += contribDh_y(j);
          }
        
        for (int j = 0; j < nb_dof_Eteta; j++)
          Rh2(j, i) += contribRh2(j) + contribRh2_bis(j);        
      }

    VectComplexe fevalDh(2*nb_pts_quad);
    Ones.Reallocate(nb_dof_E);
    VectReal_wp val_H(2*nb_pts_quad);
    for (int i = 0; i < nb_dof_E; i++)
      {
        Ones.Fill(0); Ones(i) = 1.0;
        Fb_hcurl.ApplyChTranspose(Ones, val_H);
        for (int j = 0; j < nb_pts_quad; j++)
          {
            // Bh^1 part
            vec_u(0) = val_H(2*j); vec_u(1) = val_H(2*j+1);
            Mlt(matmass_Erz(iquad)(j), vec_u, vec_v);
            fevalDh(2*j) = vec_v(0);
            fevalDh(2*j+1) = vec_v(1);
          }
        
        Fb_hcurl.ApplyCh(fevalDh, contribCh);
        for (int j = 0; j < nb_dof_E; j++)
          Bh1(j, i) += contribCh(j);        
      }
    
    // loop over dofs of Hteta
    Ones.Reallocate(nb_dof_Hteta); val_H.Reallocate(nb_pts_quad);
    VectComplexe fevalRh1(nb_pts_quad), fevalDh2(nb_pts_quad);
    contribRh2.Reallocate(nb_dof_Hteta);
    for (int i = 0; i < nb_dof_Hteta; i++)
      {
        Ones.Fill(0); Ones(i) = 1.0;
        Fb_hcurl.ApplyChTransposeH(Ones, val_H);
        
        for (int j = 0; j < nb_pts_quad; j++)
          {
            fevalRh1(j) = -matstiff_radius(iquad)(j)*val_H(j);
            fevalDh2(j) = matmass_Hteta(iquad)(j)*val_H(j);
          }
        
        Fb_hcurl.ApplyRh(fevalRh1, contribCh);
        Fb_hcurl.ApplyChH(fevalDh2, contribRh2);
        
        for (int j = 0; j < nb_dof_E; j++)
          Rh1(j, i) += contribCh(j);
        
        for (int j = 0; j < nb_dof_Hteta; j++)
          Dh2(j, i) += contribRh2(j);
      }
    
    contribRh2.Reallocate(nb_dof_Eteta);
    Ones.Reallocate(nb_dof_Eteta);
    for (int i = 0; i < nb_dof_Eteta; i++)
      {
        Ones.Fill(0); Ones(i) = 1.0;
        Fb_h1.ApplyChTranspose(Ones, val_H);
        
        for (int j = 0; j < nb_pts_quad; j++)
          fevalDh2(j) = matmass_Eteta(iquad)(j)*val_H(j);
          
        Fb_h1.ApplyCh(fevalDh2, contribRh2);
        for (int j = 0; j < nb_dof_Eteta; j++)
          Bh2(j, i) += contribRh2(j);
      }
    
  }
  
  
  //! recomposition of H using mixed formulation
  template<class TypeEquation>
  void HarmonicMaxwellAxi_Hcurl<TypeEquation>::
  EvaluateH_MixedFormulation(const Vector<Vector<Complexe> >& U0, Vector<Vector<Complexe> >& EvalH_Nodal) const
  {
    if (modified_formulation)
      {
	EvalH_Nodal = U0;
	return;
      }

    EvalH_Nodal.Reallocate(1);
    EvalH_Nodal(0).Reallocate(OffsetDofH(this->mesh.GetNbElt()));
    for (int num_elem = 0; num_elem < this->mesh.GetNbElt(); num_elem++)
      {
	const ElementReference<Dimension2, 2>& Fb_hcurl = this->GetReferenceElementHcurl(num_elem);
        const ElementReference<Dimension2, 1>& Fb_h1 = this->GetReferenceElementH1(num_elem, 1);
        int nb_dof_H = 2*Fb_h1.GetNbDof(), nb_dof_Hteta = Fb_hcurl.GetNbDofH();
        int nb_dof_E = Fb_hcurl.GetNbDof(), nb_dof_Eteta = Fb_h1.GetNbDof();
	if (this->FirstOrderFormulation())
	  {
	    int offset = OffsetDofH(num_elem);
	    for (int i = 0; i < nb_dof_H; i++)
	      EvalH_Nodal(0)(OffsetDofH(num_elem) + i) = U0(2)(offset + i);
	    
	    for (int i = 0; i < nb_dof_Hteta; i++)
              EvalH_Nodal(0)(OffsetDofH(num_elem) + nb_dof_H + i) = U0(2)(offset+nb_dof_H+i);
	  }
	else
	  {
	    Vector<Complexe> E(nb_dof_E), Eteta(nb_dof_Eteta);
	    Vector<Complexe> H(nb_dof_H), Hteta(nb_dof_Hteta);
	    
	    Matrix<Complexe> Ch, Rh1, Rh2, Dh1, Dh2, Bh1, Bh2;
	    
	    IVect pivot1(nb_dof_H), pivot2(nb_dof_Hteta);
	    
	    // we retrieve E and Eteta
	    this->GetLocalUnknownVectorE(U0, num_elem, Fb_hcurl, Fb_h1, E, Eteta);
            
	    // we compute Dh^1, Dh^2, Ch, Rh1, Rh2
	    this->ComputeLocalMatrices_MixedFormulation(num_elem, Ch, Rh1, Rh2,
							Dh1, Dh2, Bh1, Bh2, Fb_hcurl, Fb_h1);
	    
	    Seldon::GetLU(Dh1, pivot1);
	    Seldon::GetLU(Dh2, pivot2);
	    
	    MltAdd(Complexe(-1), SeldonTrans, Rh2, Eteta, Complexe(0), H);
	    MltAdd(Complexe(-1), SeldonTrans, Ch, E, Complexe(1), H);
	    Seldon::SolveLU(Dh1, pivot1, H);
	    
	    MltAdd(Complexe(-1), SeldonTrans, Rh1, E, Complexe(0), Hteta);
	    Seldon::SolveLU(Dh2, pivot2, Hteta);
	    
	    // we fill EvalH_Nodal
	    for (int i = 0; i < nb_dof_H; i++)
	      EvalH_Nodal(0)(OffsetDofH(num_elem) + i) = H(i);
	    
	    for (int i = 0; i < nb_dof_Hteta; i++)
	      EvalH_Nodal(0)(OffsetDofH(num_elem) + nb_dof_H + i) = Hteta(i);    
	    
	  }
      }
  }

  
  IncidentWaveProjector<Complex_wp, Dimension2>*
  EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>
  ::GetNewIncidentProjector(int n, const Vector<VectString>&,
			    IncidentWaveField<Complex_wp, Dimension2>& u_inc) const
  {
    return new IncidentWaveProjector_MaxwellAxi(*this, u_inc);
  }

    //! Computation of elementary matrix 
  /*! \param[in] iquad element number
    \param[out] num_dof dof numbers
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat coefficients alpha and beta
    \param[in] vars given problem
  */
  template<class Complexe, class Prop, class Storage, class TypeEquation>
  void ComputeElementaryMatrix(int iquad, IVect& num_dof, Matrix<Complexe, Prop, Storage>& mat,
			       const GlobalGenericMatrix<Complexe>& nat_mat,
			       const EllipticProblem<TypeEquation>& vars,
                               const ElementReference<Dimension2, 2>& Fb_hcurl,
                               const ElementReference<Dimension2, 1>& Fb_h1)
  {
    int nb_dof_face = Fb_hcurl.GetNbDof() + Fb_h1.GetNbDof();
    int nb_dof_H = 2*Fb_h1.GetNbDof();
    int nb_dof_Hteta = Fb_hcurl.GetNbDofH();
    int nb_dof_E = Fb_hcurl.GetNbDof();
    int nb_dof_Eteta = Fb_h1.GetNbDof();
    int nb_dof_all = nb_dof_face;
    if (vars.FirstOrderFormulation())
      nb_dof_all += nb_dof_H + nb_dof_Hteta;
    
    mat.Reallocate(nb_dof_all, nb_dof_all);
    mat.Fill(0);

    typedef typename TypeEquation::Complexe ComplexeEq;
    //const Mesh<Dimension2>& mesh = vars.mesh;
    const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
    const MeshNumbering<Dimension2>& mesh_num_h1 = vars.GetMeshNumbering(1);
    
    num_dof.Reallocate(nb_dof_all);
    for (int j = 0; j < nb_dof_E; j++)
      num_dof(j) = mesh_num.Element(iquad).GetNumberDof(j);

    for (int j = 0; j < nb_dof_Eteta; j++)
      num_dof(nb_dof_E + j) = vars.offset_dof_unknown(1) + mesh_num_h1.Element(iquad).GetNumberDof(j);
    
    int Nvol = vars.offset_dof_unknown(2);
    if (vars.FirstOrderFormulation())
      for (int j = 0; j < nb_dof_H + nb_dof_Hteta; j++)
        num_dof(nb_dof_face + j) = Nvol + vars.OffsetDofH(iquad) + j;
    
    if (vars.modified_formulation)
      {
	int N = Fb_h1.GetNbPointsQuadratureInside();
	// alternative computation of finite element matrix
        /*
	if (false)
	{
	  Vector<R2> ValueU(Fb.GetNbDof());
	  Vector<Real_wp> ValueEteta(Fb.GetNbDof());
	  
	  Vector<Real_wp> Curl_rU(Fb.GetNbDof()), Curl_U(Fb.GetNbDof());
	  Vector<R2> Grad_rEteta(Fb.GetNbDof());
	  
	  VectR3 val_phi(Fb.GetNbDof()), curl_phi(Fb.GetNbDof()); 
	  VectR2 s;
	  mesh.GetVerticesElement(iquad, s);
	  SetPoints<Dimension2> PtsElem; SetMatrices<Dimension2> MatElem;
	  
	  Fb.FjElemQuadrature(s, PtsElem, mesh, iquad);
	  Fb.DFjElemQuadrature(s, PtsElem, MatElem, mesh, iquad);
	  
	  Matrix2_2 mat_dfj, dfjm1;
	  int ref = mesh.Element(iquad).GetReference();
	  Complex_wp invMu = vars.ref_invMu(ref)(0, 0);
	  Complex_wp eps_om2 = (vars.ref_epsilon(ref)(0, 0)+Iwp*vars.ref_sigma(ref).GetConstant()/vars.GetOmega())*vars.GetSquareOmega();
	  for (int k = 0; k < N; k++)
	    {
	      R2 point_loc = Fb.PointsND(k);
	      R2 point_glob = PtsElem.GetPointQuadrature(k);
	      mat_dfj = MatElem.GetPointQuadrature(k);
	      GetInverse(mat_dfj, dfjm1);
	      Fb.ComputeValuesPhi(point_loc, val_phi, dfjm1, mesh_num, iquad);

	      Fb.ComputeValuesGradientPhi(point_loc, curl_phi, dfjm1, mesh_num, iquad);
	      
	      Real_wp poids = Det(mat_dfj)*Fb.WeightsND(k);
	      Real_wp r = point_glob(0);
	      Real_wp m2 = square(vars.GetCurrentModeNumber());
	      for (int i = 0; i < Fb.GetNbDof(); i++)
		{
		  ValueU(i).Init(val_phi(i)(0), val_phi(i)(2));
		  ValueEteta(i) = val_phi(i)(1);
		  
		  // curl (rU) = d/dr(r U_z) - d/dz(r U_r) = U_z + r curl U
		  Curl_rU(i) = val_phi(i)(2) + r*curl_phi(i)(1);
		  Curl_U(i) = curl_phi(i)(1);
		  
		  // grad(r Eteta) = [Eteta + r dEteta/dr, r dEteta/dz]
		  Grad_rEteta(i)(0) = val_phi(i)(1) - r*curl_phi(i)(2);
		  Grad_rEteta(i)(1) = r*curl_phi(i)(0);
		}
	      
	      Complex_wp vloc(0, 0);
	      for (int i = 0; i < Fb.GetNbDof(); i++)
		{
		  int j0 = 0;
		  if (IsSymmetricMatrix(mat))
		    j0 = i;
		  
		  for (int j = j0; j < Fb.GetNbDof(); j++)
		    {
		      if (vars.GetCurrentModeNumber() == 0)
			{
			  vloc = invMu/r*DotProd(Grad_rEteta(i), Grad_rEteta(j));
			  vloc += invMu*r*Curl_U(i)*Curl_U(j);		
			  vloc -= eps_om2*r*DotProd(ValueU(i), ValueU(j));
			  vloc -= eps_om2*r*ValueEteta(i)*ValueEteta(j);
			}
		      else
			{
			  vloc = invMu*r*m2*DotProd(ValueU(i), ValueU(j));
			  vloc += invMu*r*Curl_rU(i)*Curl_rU(j);
			  
			  vloc -= eps_om2*r*DotProd(r*ValueU(i) + Grad_rEteta(i),
						    r*ValueU(j) + Grad_rEteta(j));
			  
			  vloc -= eps_om2*m2*r*ValueEteta(i)*ValueEteta(j);
			}
		      
		      mat(i, j) += poids*vloc;
		    }
		  
		}
	    }
	}	  
        */
	//Matrix<Complexe, Prop, Storage> mat_ref = mat;
	//mat.Fill(0);
	
        Vector<TinyMatrix<ComplexeEq, Symmetric, 2, 2> > mat_mass(N), mat_stiff(N), mat_coup(N);
        Vector<TinyVector<ComplexeEq, 2> > vecE_stiff(N), vec_stiff(N), vec_coup(N);
	Vector<ComplexeEq> coef_stiff(N), coef_mass(N);
	
	Real_wp m2 = square(vars.GetCurrentModeNumber());
	for (int i = 0; i < N; i++)
	  {
	    if (vars.GetCurrentModeNumber() == 0)
	      {
		mat_mass(i) = vars.matmass_Erz(iquad)(i);
		mat_stiff(i) = vars.matmass_Hrz(iquad)(i);
		vec_stiff(i) = vars.vecstiff_Eteta(iquad)(i);
		coef_mass(i) = vars.matmass_Eteta(iquad)(i) + vars.matmass_Eteta_m0(iquad)(i);
		coef_stiff(i) = vars.matmass_Hteta(iquad)(i);
	      }
	    else
	      {
		mat_mass(i) = vars.matmass_Erz_r3(iquad)(i)
		  + m2*vars.matmass_Erz_m2(iquad)(i) + vars.matmass_Erz_zz(iquad)(i);
		
		mat_stiff(i) = vars.matmass_Hrz_r3(iquad)(i);
		mat_coup(i) = vars.matmass_Hrz_r3(iquad)(i);        
		
		vecE_stiff(i) = vars.vecstiff_Erz(iquad)(i);
		vec_stiff(i) = vars.vecstiff_Eteta_r2(iquad)(i);
		vec_coup(i) = vars.vecstiff_Eteta_r2(iquad)(i);
		coef_stiff(i) = vars.matmass_Hteta_r3(iquad)(i);
		coef_mass(i) = vars.matmass_Eteta(iquad)(i)*m2 + vars.matmass_Eteta_rr(iquad)(i);
	      }
	  }
        
        Vector<Complexe> contribE(nb_dof_E);
        Vector<Complexe> contribE2(nb_dof_E);
        Vector<Complexe> contribEteta(nb_dof_Eteta);
        Vector<Complexe> contribEteta2(nb_dof_Eteta);
        contribE.Fill(0); contribE2.Fill(0);
        contribEteta.Fill(0); contribEteta2.Fill(0);
        Vector<Complexe> feval_phi(2*N), feval_curlPhi(N), feval_psi(N), feval_gradPsi(2*N);
        feval_phi.Fill(0); feval_curlPhi.Fill(0);
        feval_psi.Fill(0); feval_gradPsi.Fill(0);
        VectReal_wp Ones(nb_dof_E);
        VectReal_wp valE(2*N), curlE(N);
        valE.Fill(0); curlE.Fill(0);
        TinyVector<Real_wp, 2> vecE, gradE;
        TinyVector<Complexe, 2> vec_u; Complexe vloc;
        // loop on columns of E
        for (int j = 0; j < nb_dof_E; j++)
          {
            Ones.Fill(0); Ones(j) = 1.0;
            
            // computation of E and rot E on quadrature points
            Fb_hcurl.ApplyChTranspose(Ones, valE);
            Fb_hcurl.ApplyRhTranspose(Ones, curlE);
            
            // forming mat_mass E + vecE_stiff curl E, and coef_stiff curl E + vecE_stiff E
            // and vec_coup E,  mat_coup E
            for (int k = 0; k < N; k++)
              {
                vecE(0) = valE(2*k); vecE(1) = valE(2*k+1);
                Mlt(mat_mass(k), vecE, vec_u);
                Add(curlE(k), vecE_stiff(k), vec_u);
                feval_phi(2*k) = vec_u(0); feval_phi(2*k+1) = vec_u(1);
                
                vloc = DotProd(vecE_stiff(k), vecE);
                vloc += coef_stiff(k)*curlE(k);
                feval_curlPhi(k) = vloc;
                
                Mlt(mat_coup(k), vecE, vec_u);
                feval_gradPsi(2*k) = vec_u(0); feval_gradPsi(2*k+1) = vec_u(1);
                
                vloc = DotProd(vec_coup(k), vecE);
                feval_psi(k) = vloc;
              }
            
            Fb_hcurl.ApplyCh(feval_phi, contribE);
            Fb_hcurl.ApplyRh(feval_curlPhi, contribE2);
            Fb_h1.ApplyCh(feval_psi, contribEteta);
            Fb_h1.ApplyRh(feval_gradPsi, contribEteta2);
            
            for(int i = 0; i < nb_dof_E; i++)
              mat(i, j) = contribE(i) + contribE2(i);
            
            for(int i = 0; i < nb_dof_Eteta; i++)
              mat(i + nb_dof_E, j) = contribEteta(i) + contribEteta2(i);
          }
        
        // loop on columns of Eteta
        Ones.Reallocate(nb_dof_Eteta);
        curlE.Clear(); valE.Clear();
        VectReal_wp valEteta(N), gradEteta(2*N);
        valEteta.Fill(0); gradEteta.Fill(0);
        for (int j = 0; j < nb_dof_Eteta; j++)
          {
            Ones.Fill(0); Ones(j) = 1.0;
            
            // computation of E and rot E on quadrature points
            Fb_h1.ApplyChTranspose(Ones, valEteta);
            Fb_h1.ApplyRhTranspose(Ones, gradEteta);
            
            // forming vec_coup Eteta + mat_coup grad Eteta 
            // and coef_mass Eteta + vec_stiff grad Eteta, mat_stiff grad Eteta + vec_stiff Eteta
            for (int k = 0; k < N; k++)
              {
                gradE(0) = gradEteta(2*k); gradE(1) = gradEteta(2*k+1);
                Mlt(mat_coup(k), gradE, vec_u);
                Add(valEteta(k), vec_coup(k), vec_u);
                feval_phi(2*k) = vec_u(0); feval_phi(2*k+1) = vec_u(1);
                
                Mlt(mat_stiff(k), gradE, vec_u);
                Add(valEteta(k), vec_stiff(k), vec_u);
                feval_gradPsi(2*k) = vec_u(0); feval_gradPsi(2*k+1) = vec_u(1);
                
                vloc = DotProd(vec_stiff(k), gradE) + coef_mass(k)*valEteta(k);
                feval_psi(k) = vloc;
              }
	    
            Fb_hcurl.ApplyCh(feval_phi, contribE);
            Fb_h1.ApplyCh(feval_psi, contribEteta);
            Fb_h1.ApplyRh(feval_gradPsi, contribEteta2);
            
            for(int i = 0; i < nb_dof_E; i++)
              mat(i, nb_dof_E + j) = contribE(i);
            
            for(int i = 0; i < nb_dof_Eteta; i++)
              mat(nb_dof_E + i, nb_dof_E + j) = contribEteta(i) + contribEteta2(i);
          }

        mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat, iquad, 1, 0);
        mesh_num_h1.number_map.ModifyLocalRowMatrix(mesh_num_h1, mat, iquad, 1, nb_dof_E);

        mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat, iquad, 1, 0);
        mesh_num_h1.number_map.ModifyLocalColumnMatrix(mesh_num_h1, mat, iquad, 1, nb_dof_E);
	
        return;
      }
    
    Matrix<Complexe> Ch, Rh1, Rh2, Dh1, Dh2, Bh1, Bh2;
    vars.ComputeLocalMatrices_MixedFormulation(iquad, Ch, Rh1, Rh2, Dh1, Dh2, Bh1, Bh2, Fb_hcurl, Fb_h1);
    
    Complexe coef_mass = nat_mat.GetCoefMass();
    Complexe coef_stiff = nat_mat.GetCoefStiffness();
    
    for (int i = 0; i < nb_dof_E; i++)
      for (int j = 0; j < nb_dof_E; j++)
        mat(i, j) = coef_mass*Bh1(i, j);

    for (int i = 0; i < nb_dof_Eteta; i++)
      for (int j = 0; j < nb_dof_Eteta; j++)
        mat(nb_dof_E + i, nb_dof_E + j) = coef_mass*Bh2(i, j);
    
    if (vars.FirstOrderFormulation())
      {
        for (int i = 0; i < nb_dof_E; i++)
          for (int j = 0; j < nb_dof_H; j++)
            {
              int i2 = i;
              int j2 = nb_dof_face + j;
              mat(i2, j2) = coef_stiff*Ch(i, j);
              mat(j2, i2) = coef_stiff*Ch(i, j);
            }
        
        for (int i = 0; i < nb_dof_E; i++)
          for (int j = 0; j < nb_dof_Hteta; j++)
            {
              int i2 = i;
              int j2 = nb_dof_face + nb_dof_H + j;
              mat(i2, j2) = coef_stiff*Rh1(i, j);
              mat(j2, i2) = coef_stiff*Rh1(i, j);
            }

        for (int i = 0; i < nb_dof_Eteta; i++)
          for (int j = 0; j < nb_dof_H; j++)
            {
              int i2 = nb_dof_E + i;
              int j2 = nb_dof_face + j;
              mat(i2, j2) = coef_stiff*Rh2(i, j);
              mat(j2, i2) = coef_stiff*Rh2(i, j);
            }
        
        for (int i = 0; i < nb_dof_H; i++)
          for (int j = 0; j < nb_dof_H; j++)
            {
              int i2 = nb_dof_face + i;
              int j2 = nb_dof_face + j;
              mat(i2, j2) = coef_mass*Dh1(i, j);
            }

        for (int i = 0; i < nb_dof_Hteta; i++)
          for (int j = 0; j < nb_dof_Hteta; j++)
            {
              int i2 = nb_dof_face + nb_dof_H + i;
              int j2 = nb_dof_face + nb_dof_H + j;
              mat(i2, j2) = coef_mass*Dh2(i, j);
            }

        mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat, iquad, 1, 0);
        mesh_num_h1.number_map.ModifyLocalRowMatrix(mesh_num_h1, mat, iquad, 1, nb_dof_E);

        mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat, iquad, 1, 0);
        mesh_num_h1.number_map.ModifyLocalColumnMatrix(mesh_num_h1, mat, iquad, 1, nb_dof_E);

        return;
      }
    
    // now we compute schur complement
    // For E :
    // (B_h^1 - C_h (D_h^1)^{-1} C_h^t) E - C_h (D_h^1)^{-1} (R_h^2)^t Eteta 
    //    - R_h^1 (D_h^2)^{-1} (R_h^1)^t E
    // B_h^1 is already put on the matrix mat, we compute the other sub-matrices
    
    // For Eteta :
    // B_h^2 Eteta - R_h^2 (D_h^1)^{-1} C_h^t E - R_h^2 (D_h^1)^{-1} (R_h^2)^t Eteta
    
    // we take in accounty symmetry of mat, and compute only one time the coupling term
    // C_h (D_h^1)^{-1} (R_h^2)^t
    
    Matrix<Complexe> matE_E(nb_dof_E, nb_dof_E),
      matE_Eteta(nb_dof_E, nb_dof_Eteta), matEteta_Eteta(nb_dof_Eteta, nb_dof_Eteta);
    
    Matrix<Complexe> A_tmp(nb_dof_H, nb_dof_E),
      B_tmp(nb_dof_Hteta, nb_dof_E), C_tmp(nb_dof_H, nb_dof_Eteta);
    
    // we compute (D_h^1)^{-1} and (D_h^2)^{-1}
    GetInverse(Dh1); GetInverse(Dh2); //DISP(Dh1); DISP(Dh2);
    // A_tmp = (D_h^1)^{-1} C_h^t
    MltAdd(Complexe(1), SeldonNoTrans, Dh1, SeldonTrans, Ch, Complexe(0), A_tmp);
    // matE_E = - C_h (D_h^1)^{-1} C_h^t
    MltAdd(Complexe(-1), Ch, A_tmp, Complexe(0), matE_E);
    
    // B_tmp = (D_h^2)^{-1} (R_h^1)^t
    MltAdd(Complexe(1), SeldonNoTrans, Dh2, SeldonTrans, Rh1, Complexe(0), B_tmp);
    // matE_E = - C_h (D_h^1)^{-1} C_h^t - R_h^1 (D_h^2)^{-1} (R_h^1)^t
    MltAdd(Complexe(-1), Rh1, B_tmp, Complexe(1), matE_E);
    
    // C_tmp = (D_h^1)^{-1} (R_h^2)^t
    MltAdd(Complexe(1), SeldonNoTrans, Dh1, SeldonTrans, Rh2, Complexe(0), C_tmp);
    // matEteta_Eteta = - R_h^2 (D_h^1)^{-1} (R_h^2)^t
    MltAdd(Complexe(-1), Rh2, C_tmp, Complexe(0), matEteta_Eteta);
    
    // C_tmp = (D_h^1)^{-1} (R_h^2)^t
    // matE_Eteta = - C_h (D_h^1)^{-1} (R_h^2)^t
    MltAdd(Complexe(-1), Ch, C_tmp, Complexe(0), matE_Eteta);
    
    if (IsSymmetricMatrix(mat))
      {
        for (int i = 0; i < nb_dof_E; i++)
          for (int j = i; j < nb_dof_E; j++)
            mat(i, j) += coef_stiff*matE_E(i,j);
        
        for (int i = 0; i < nb_dof_Eteta; i++)
          for (int j = i; j < nb_dof_Eteta; j++)
            mat(nb_dof_E +i, nb_dof_E+j) += coef_stiff*matEteta_Eteta(i,j);
	
        for (int i = 0; i < nb_dof_E; i++)
          for (int j = 0; j < nb_dof_Eteta; j++)
            mat(i, nb_dof_E+j) += coef_stiff*matE_Eteta(i,j);
    
      }
    else
      {
        for (int i = 0; i < nb_dof_E; i++)
          for (int j = 0; j < nb_dof_E; j++)
            mat(i, j) += coef_stiff*matE_E(i,j);
        
        for (int i = 0; i < nb_dof_Eteta; i++)
          for (int j = 0; j < nb_dof_Eteta; j++)
            mat(nb_dof_E +i, nb_dof_E+j) += coef_stiff*matEteta_Eteta(i,j);
	
        for (int i = 0; i < nb_dof_E; i++)
          for (int j = 0; j < nb_dof_Eteta; j++)
            {
              mat(i, nb_dof_E+j) += coef_stiff*matE_Eteta(i,j);
              mat(nb_dof_E+j, i) += coef_stiff*matE_Eteta(i,j);
            }
      }

    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat, iquad, 1, 0);
    mesh_num_h1.number_map.ModifyLocalRowMatrix(mesh_num_h1, mat, iquad, 1, nb_dof_E);
    
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat, iquad, 1, 0);
    mesh_num_h1.number_map.ModifyLocalColumnMatrix(mesh_num_h1, mat, iquad, 1, nb_dof_E);

  }


  void EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicMaxwellEquation_HcurlAxi>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, dynamic_cast<Matrix<Complex_wp>& >(mat_elem), nat_mat, *this,
				      this->GetReferenceElementHcurl(i),
                                      this->GetReferenceElementH1(i, 1));
  }

}

#define MONTJOIE_FILE_AXISYM_HCURL_MAXWELL_CXX
#endif
