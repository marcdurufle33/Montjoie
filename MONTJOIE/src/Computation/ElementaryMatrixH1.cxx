#ifndef MONTJOIE_FILE_ELEMENTARY_MATRIX_H1_CXX

namespace Montjoie
{  
  
  //! elementary matrix for any H1 finite element (standard quadrature)
  /*!
    \param[in] iquad element number
    \param[out] num_dof dof numbers
    \param[out] mat_interac elementary matrix
    \param[in] nat_mat coefficients mass and stiff
    \param[in] vars given problem
    \param[in] Fb finite element considered
    elementary matrix \f$  \int Amass * \varphi_j \varphi_i 
    + \int Cgrad_grad \nabla \varphi_j \nabla \varphi_i 
    + \int Dgrad_phi \varphi_j \nabla \varphi_i  + \int Ephi_grad \nabla \varphi_j \varphi_i \f$
    for any element, the operators Amass, Cgrad_grad, Dgrad_phi and Ephi_grad
    can vary inside the element and depend on the solved equation
  */
  template<class Complexe, class TypeEquation,
	   class Vector1, class Vector2, class Vector3>
  void ComputeVariableElementaryMatrix(int iquad, VirtualMatrix<Complexe> & mat,
				       Vector1& Amass,
				       Vector2& Cgrad_grad, Vector3& Dgrad_phi, Vector3& Ephi_grad,
				       const GlobalGenericMatrix<Complexe>& nat_mat,
				       const EllipticProblem<TypeEquation>& vars,
				       const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Dimension Dimension;
    typename Dimension::R_N grad_phi1, grad_phi2;

    int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
    int nb_dof_elt = Fb.GetNbDof();
    
    VectReal_wp val_phi(nb_points_quadrature);
    Vector<Complexe> feval(nb_points_quadrature), feval_bis(nb_points_quadrature);
    Vector<Complexe> feval_grad(Dimension::dim_N*nb_points_quadrature);
    typename Dimension::VectR_N grad_phi(nb_points_quadrature);
    Vector<Complexe> contrib(nb_dof_elt), contrib_grad(nb_dof_elt);
    contrib_grad.Zero();

    Complexe one; SetComplexOne(one);
    TinyVector<Complexe, Dimension::dim_N> vec_phi;
    
    // loop on rows of the matrix
    for (int i = 0; i < nb_dof_elt; i++)
      {
        // we compute phi_i(\xi_k) and nabla phi_i for all k
	Fb.GetGradientSinglePhiQuadrature(i, val_phi, grad_phi);
        
        // then loops over unknowns
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
            {
              // computation of grad_feval = E phi_i + C \nabla phi_i
              for (int k = 0; k < nb_points_quadrature; k++)
                for (int q = 0; q < Dimension::dim_N; q++) 
                  {
                    feval_grad(Dimension::dim_N*k+q) = Ephi_grad(k)(m, n)(q)*val_phi(k);
                    for (int p = 0; p < Dimension::dim_N; p++)
                      feval_grad(Dimension::dim_N*k+q) += Cgrad_grad(k)(m, n)(p, q)
                        *grad_phi(k)(p);                    
                  }
              
              // then integration against \nabla phi_j
              if (Fb.UseQuadraturePointsForRh())
                Fb.ApplyRhQuadrature(feval_grad, feval_bis);
              else
                Fb.ApplyRh(feval_grad, contrib_grad);
              
              // computation of feval = Amass phi_i + D \nabla phi_i
              for (int k = 0; k < nb_points_quadrature; k++)
                {
                  feval(k) = Amass(k)(m, n)*val_phi(k);
                  for (int p = 0; p < Dimension::dim_N; p++)
                    feval(k) += Dgrad_phi(k)(m, n)(p)*grad_phi(k)(p);                  
                }
              
              // integration against varphi_j and summation with contribution of grad_feval
              if (Fb.UseQuadraturePointsForRh())
                {
                  Add(one, feval_bis, feval);
                  Fb.ApplyCh(feval, contrib);
                }
              else
                {
                  Fb.ApplyCh(feval, contrib);
                  Add(one, contrib_grad, contrib);
                }
              
              // we got a row of the matrix
              for (int j = 0; j < nb_dof_elt; j++)
                {
                  int num_col = j + n*nb_dof_elt;
                  int num_row = i + m*nb_dof_elt;
                  mat.SetEntry(num_row, num_col, contrib(j));
                }
            }
      }
  }
  

  template<class Complexe, class Dimension>
  void ComputeVariableElementaryMatrix(Array2D<Vector<Complexe> >& Amass, 
                                       Array2D<Vector<TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N> > >& Cgrad_grad,
                                       Array2D<Vector<TinyVector<Complexe, Dimension::dim_N> > >& Dgrad_phi,
                                       Array2D<Vector<TinyVector<Complexe, Dimension::dim_N> > >& Ephi_grad,
                                       const Vector<int>& offset_loc,  const Vector<int>& num_unknown, 
                                       Vector<const ElementReference<Dimension, 1>* >& Fb,
                                       VirtualMatrix<Complexe>& mat)
  {
    int nb_unknowns = Amass.GetM();
    int nb_mesh_num = Fb.GetM();
    int nb_points_quadrature = Amass(0, 0).GetM();
    VectReal_wp val_phi(nb_points_quadrature);
    typename Dimension::VectR_N grad_phi(nb_points_quadrature);

    Vector<Complexe> feval(nb_points_quadrature), feval_bis(nb_points_quadrature);
    Vector<Complexe> feval_grad(Dimension::dim_N*nb_points_quadrature);

    Vector<Vector<Complexe> > contrib(nb_mesh_num), contrib_grad(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      {
        int nb_dof_elt = Fb(i)->GetNbDof();
        contrib(i).Reallocate(nb_dof_elt);
        contrib_grad(i).Reallocate(nb_dof_elt);
        contrib(i).Zero(); contrib_grad(i).Zero();
      }

    Complexe one; SetComplexOne(one);
    
    // loop on rows of the matrix
    for (int m = 0; m < nb_unknowns; m++)
      {
        int nm = num_unknown(m);        
        for (int i = 0; i < Fb(nm)->GetNbDof(); i++)
          {          
            // we compute phi_i(\xi_k) and nabla phi_i for all k
            Fb(nm)->GetGradientSinglePhiQuadrature(i, val_phi, grad_phi);
            
            // then loops over unknowns (for columns)
            for (int n = 0; n < nb_unknowns; n++)
              {
                int nn = num_unknown(n);
                
                // computation of grad_feval = E phi_i + C \nabla phi_i
                for (int k = 0; k < nb_points_quadrature; k++)
                  for (int q = 0; q < Dimension::dim_N; q++) 
                    {
                      feval_grad(Dimension::dim_N*k+q) = Ephi_grad(m, n)(k)(q)*val_phi(k);
                      for (int p = 0; p < Dimension::dim_N; p++)
                        feval_grad(Dimension::dim_N*k+q) += Cgrad_grad(m, n)(k)(p, q)
                          *grad_phi(k)(p);                    
                    }
                
                // then integration against \nabla phi_j
                if (Fb(nn)->UseQuadraturePointsForRh())
                  Fb(nn)->ApplyRhQuadrature(feval_grad, feval_bis);
                else
                  Fb(nn)->ApplyRh(feval_grad, contrib_grad(nn));
                
                // computation of feval = Amass phi_i + D \nabla phi_i
                for (int k = 0; k < nb_points_quadrature; k++)
                  {
                    feval(k) = Amass(m, n)(k)*val_phi(k);
                    for (int p = 0; p < Dimension::dim_N; p++)
                      feval(k) += Dgrad_phi(m, n)(k)(p)*grad_phi(k)(p);                  
                  }
              
                // integration against varphi_j and summation with contribution of grad_feval
                if (Fb(nn)->UseQuadraturePointsForRh())
                  {
                    Add(one, feval_bis, feval);
                    Fb(nn)->ApplyCh(feval, contrib(nn));
                  }
                else
                  {
                    Fb(nn)->ApplyCh(feval, contrib(nn));
                    Add(one, contrib_grad(nn), contrib(nn));
                  }
                
                // we got a row of the matrix
                for (int j = 0; j < Fb(nn)->GetNbDof(); j++)
                  {
                    int num_col = j + offset_loc(n);
                    int num_row = i + offset_loc(m);
                    mat.SetEntry(num_row, num_col, contrib(nn)(j));
                  }
              }
          }
      }
  }
  
  
  //! computation of volumetric integrals on an element for variational formulation
  /*!
    \param[in] iquad element number
    \param[out] num_dof numbers of degrees of freedom
    \param[out] mat_interac elementary matrix
    \param[in] vars considered problem
    \param[in] Fb finite element class
    This function is used for H1 transformation (therefore H1 and/or DG formulation)
    It computes the elementary matrix \f$  M_{i, j} \ , = \, \int A * \varphi_j \varphi_i 
    + \int C \nabla \varphi_j \nabla \varphi_i 
    + \int D \varphi_j \nabla \varphi_i  + \int E \nabla \varphi_j \varphi_i \f$
    Operators A, C, D and E are provided in object TypeEquation respectively through 
    static methods GetTensorMass, GetGradGradTensor and GetGradPhiTensor (both D and E)
   */
  template<class Complexe, class TypeEquation>
  void ComputeElementaryMatrix(int iquad, IVect& num_dof,
			       VirtualMatrix<Complexe>& mat_interac,
			       const GlobalGenericMatrix<Complexe>& nat_mat,
			       const EllipticProblem<TypeEquation>& vars,
			       const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Dimension Dimension;
    
    int nb_dof_elt = Fb.GetNbDof();
    
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    int nb_mesh_num = vars.GetNbMeshNumberings();
    
    // number of integration points
    int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
    int ref_domain = mesh.Element(iquad).GetReference();
    
    // dof numbers
    Vector<int> offset_loc;
    Vector<const ElementReference<Dimension, 1>* > all_Fb;
    if (nb_mesh_num == 1)
      {
        int nb_dof_all = TypeEquation::nb_unknowns*nb_dof_elt;
        
        mat_interac.Clear();
        mat_interac.Reallocate(nb_dof_all, nb_dof_all);
        mat_interac.Zero();
        
        num_dof.Reallocate(nb_dof_all);
        IVect Nodle = vars.GetDofNumberOnElement(iquad);
        for (int j = 0; j < nb_dof_elt; j++)
          for (int m = 0; m < TypeEquation::nb_unknowns; m++)
            {
              num_dof(j + nb_dof_elt*m) = Nodle(j) + vars.offset_dof_unknown(m);
              if (Nodle(j) < 0)
                num_dof(j + nb_dof_elt*m) = -1;
            }
      }
    else
      {
        // case with several numberings
        Vector<IVect> Nodle(nb_mesh_num);
        all_Fb.Reallocate(nb_mesh_num);
        for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            Nodle(nm) = vars.GetDofNumberOnElement(iquad, nm);
            all_Fb(nm) = &vars.GetReferenceElementH1(iquad, nm);
          }
        
        offset_loc.Reallocate(TypeEquation::nb_unknowns+1);
        int nb_dof_all = 0; offset_loc(0) = 0;
        for (int n = 0; n < TypeEquation::nb_unknowns; n++)
          {
            int nm = vars.mesh_num_unknown(n);            
            //const_cast<ElementReference<Dimension, 1>* >
            nb_dof_all += all_Fb(nm)->GetNbDof();
            offset_loc(n+1) = nb_dof_all;            
          }
        
        num_dof.Reallocate(nb_dof_all);
        for (int n = 0; n < TypeEquation::nb_unknowns; n++)
          {
            int nm = vars.mesh_num_unknown(n);
            int nb_dof_loc = all_Fb(nm)->GetNbDof();
            for (int j = 0; j < nb_dof_loc; j++)
              {
                num_dof(offset_loc(n) + j) = Nodle(nm)(j) + vars.offset_dof_unknown(n);
                if (Nodle(nm)(j) < 0)
                  num_dof(offset_loc(n) + j) = -1;
              }
          }
        
        mat_interac.Clear();
        mat_interac.Reallocate(nb_dof_all, nb_dof_all);
        mat_interac.Zero();
      }
    
    bool variable = vars.UseNumericalIntegration(iquad);
    bool affine = vars.mesh.IsElementAffine(iquad);

    // types used to store geometric operators
    typedef TinyMatrix<Complexe, General,
      TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> Matrixn_n;
    
    typedef TinyMatrix<TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N>,
      General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> Matrixn2_n2;
    
    typedef TinyMatrix<TinyVector<Complexe, Dimension::dim_N>,
      General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> MatrixVectn_n;
    
    if (nb_mesh_num > 1)
      {
        int nb_u = TypeEquation::nb_unknowns;
        Array2D<Vector<TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N> > > Cgrad_grad(nb_u, nb_u);
        Array2D<Vector<TinyVector<Complexe, Dimension::dim_N> > > Dgrad_phi(nb_u, nb_u), Ephi_grad(nb_u, nb_u);
        Array2D<Vector<Complexe> > Amass(nb_u, nb_u);

        for (int m = 0; m < nb_u; m++)
          for (int n = 0; n < nb_u; n++)
            {
              Cgrad_grad(m, n).Reallocate(nb_points_quadrature);
              Dgrad_phi(m, n).Reallocate(nb_points_quadrature);
              Ephi_grad(m, n).Reallocate(nb_points_quadrature);
              Amass(m, n).Reallocate(nb_points_quadrature);
            }
        
        Matrixn_n A0; Matrixn2_n2 C0; MatrixVectn_n D0, E0;
        for (int j = 0; j < nb_points_quadrature; j++)
	  {
            TypeEquation::GetMassPhiDFiTensor(vars, iquad, j, nat_mat, ref_domain, A0,
                                              true, affine, Fb);
            
	    TypeEquation::GetGradGradDFiTensor(vars, iquad, j, nat_mat, ref_domain, C0,
                                               true, affine, Fb);
            
	    TypeEquation::GetGradPhiDFiTensor(vars, iquad, j, nat_mat, ref_domain,
                                              D0, E0,
                                              true, affine, Fb);            

            for (int m = 0; m < nb_u; m++)
              for (int n = 0; n < nb_u; n++)
                {
                  Amass(m, n)(j) = A0(m, n);
                  Cgrad_grad(m, n)(j) = C0(m, n);
                  Dgrad_phi(m, n)(j) = D0(m, n);
                  Ephi_grad(m, n)(j) = E0(m, n);
                }                       
	  }
	
        
        ComputeVariableElementaryMatrix(Amass, Cgrad_grad, Dgrad_phi, Ephi_grad,
                                        offset_loc, vars.mesh_num_unknown, all_Fb, mat_interac);
        
        return;
      }
    else if (variable)
      {
        // variable case occurs when numerical integration has to be performed
        // usually when the transformation Fi is non-affine
        //  or when physical coefficients are varying inside the element
	Vector<Matrixn_n> Amass(nb_points_quadrature);
	Vector<Matrixn2_n2> Cgrad_grad(nb_points_quadrature);
	Vector<MatrixVectn_n> Dgrad_phi(nb_points_quadrature), Ephi_grad(nb_points_quadrature);
	
	// we evaluate mass coefficient A
        // and stiffness coefficients C, D and E on each quadrature point
	for (int j = 0; j < nb_points_quadrature; j++)
	  {
            TypeEquation::GetMassPhiDFiTensor(vars, iquad, j, nat_mat, ref_domain, Amass(j),
                                              variable, affine, Fb);
            
	    TypeEquation::GetGradGradDFiTensor(vars, iquad, j, nat_mat, ref_domain, Cgrad_grad(j),
                                               variable, affine, Fb);
            
	    TypeEquation::GetGradPhiDFiTensor(vars, iquad, j, nat_mat, ref_domain,
                                              Dgrad_phi(j), Ephi_grad(j),
                                              variable, affine, Fb);
            
	  }
	
	//DISP(iquad); DISP(Amass); DISP(Cgrad_grad); DISP(Dgrad_phi); DISP(Ephi_grad);
        if (Fb.OptimizedComputationElementaryMatrix())
	  {
            // computing elementary matrix with methods of the finite element class
            Vector<TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N> >
              C(nb_points_quadrature);
            
            Vector<TinyVector<Complexe, Dimension::dim_N> >
              D(nb_points_quadrature), E(nb_points_quadrature);
            
            Vector<Complexe> A(nb_points_quadrature);
            for (int p = 0; p < TypeEquation::nb_unknowns; p++)
              for (int q = 0; q < TypeEquation::nb_unknowns; q++)
                {
                  bool no_grad = true;
		  TinyVector<bool, 4> null_term(false, true, true, true);
                  for (int i = 0; i < nb_points_quadrature; i++)
                    {
                      A(i) = Amass(i)(p, q);
                      C(i) = Cgrad_grad(i)(p, q);
                      D(i) = Dgrad_phi(i)(p, q);
                      E(i) = Ephi_grad(i)(p, q);
		      
                      if (!C(i).IsZero())
			{
			  no_grad = false;
			  null_term(1) = false;
			}

		      if (!D(i).IsZero())
			{
			  no_grad = false;
			  null_term(2) = false;
			}
		      
		      if (!E(i).IsZero())
                        {
			  no_grad = false;
			  null_term(3) = false;
			}
                    }
                  
                  if (!no_grad)
		    Fb.AddVariableElemMatrix(nb_dof_elt*p, nb_dof_elt*q, A, C, D, E, null_term, mat_interac);
                  else
		    Fb.AddVariableMassMatrix(nb_dof_elt*p, nb_dof_elt*q, A, mat_interac);
                }
          }
        else
          ComputeVariableElementaryMatrix(iquad, mat_interac, Amass, 
                                          Cgrad_grad, Dgrad_phi, Ephi_grad, nat_mat, vars, Fb);
      }
    else
      {
	// linear Fj, DFj is constant, physical coefficients are constants
        // use of precomputed matrices
	Matrixn2_n2 Cgrad_grad;
	MatrixVectn_n Dgrad_phi, Ephi_grad;
	Matrixn_n Amass;
	
        // we get A, C, D and E for all the element
	TypeEquation::GetMassPhiDFiTensor(vars, iquad, 0, nat_mat, ref_domain, Amass,
                                          variable, affine, Fb);
        
	TypeEquation::GetGradGradDFiTensor(vars, iquad, 0, nat_mat, ref_domain, Cgrad_grad,
                                           variable, affine, Fb);
        
	TypeEquation::GetGradPhiDFiTensor(vars, iquad, 0, nat_mat, ref_domain, Dgrad_phi, Ephi_grad,
                                          variable, affine, Fb);
        
	TinyVector<bool, 4> null_term(false, false, false, false);
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  for (int n = 0; n < TypeEquation::nb_unknowns; n++)
	    {
              int icol = n*nb_dof_elt; int irow = m*nb_dof_elt;
              
	      Fb.AddConstantElemMatrix(irow, icol, Amass(m, n), Cgrad_grad(m, n),
				       Dgrad_phi(m, n), Ephi_grad(m, n), null_term, mat_interac);
            }
      }
    
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad);    
  }
  
  
  //! Adding boundary integrals (usually for DG formulation)
  /*!
    adding surfacic term  \f$ + \displaystyle \sum_{\mbox{e edge or face of the mesh} }
    \int_{Gamma_e} \{ C \nabla u n \} [ \varphi ]  + \{ C \nabla \varphi n\} [ u ] 
    + \{ D u n \} [ \varphi ] + \{ E^* \varphi n \} [ u ] - 1/2 [ P u ] [\varphi]  \f$
    where e is an edge shared by elements K- and K+, n is the outward normale of K-,
    and \{ u \} = 1/2 (u+ + u- \} ,  [ u ] = (u+ - u-) 
    we add also boundary condition term \int_Gamma N u phi (for first order formulation only)
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
  template<class Complexe, class TypeEquation>
  void AddElementaryFluxDG(VirtualMatrix<Complexe>& mat_sp, const GlobalGenericMatrix<Complexe>& nat_mat,
			   const EllipticProblem<TypeEquation>& vars,
			   int iquad, int num_pos1_face, int num_elem2,
			   int num_face, int ref, bool new_face, int& nb_neighbor,
                           int offset_row, int offset_col,
			   const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    // iquad is considered as element K-, and num_elem2 element K+
    // since terms C+, D+, E+ and P+ may be unavailable (i.e. on another processor)
    // we keep only terms with C-, D-, E- and P- : 
    // C- (nabla u n)- [ phi ] +  C- (nabla phi n)- [ u ] + D- (u n)- [phi]
    // + (E^* phi n)- [u] + (P u)- [ phi ]
    // where [u] = 1/2 (u+ - u-)
    // therefore we have : ( (C nabla u n)- + (D u n)- + (P u)- ) [phi]
    //  and :    ( (C nabla phi n)- + (E^* phi n)- ) [u]
    // usually D = E^*, C and P are symmetric, leading to a symmetric matrix
    typedef typename TypeEquation::Dimension Dimension;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::R_N R_N;
    
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    
    int rf = mesh_num.GetOrderQuadrature(num_face);
    int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
    VectReal_wp PoidsFlux = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
    
    R_N normale;
    typedef TinyMatrix<Complexe, General,
      TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> Matrixn_n;
    
    typedef TinyMatrix<TinyMatrix<Complexe, General, Dimension::dim_N, Dimension::dim_N>,
      General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> Matrixn2_n2;
  
    typedef TinyMatrix<TinyVector<Complexe, Dimension::dim_N>,
      General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> MatrixVectn_n;
    
    MatrixVectn_n Ei, Di; Matrixn2_n2 Ci; Matrixn_n Ncl;
    Real_wp dsj, coef; MatrixN_N dfjm1;
    TinyVector<TinyVector<Complexe, Dimension::dim_N>, Dimension::dim_N> bvec;
    bool first_order = TypeEquation::FirstOrderFormulation;
    bool affine = vars.mesh.IsElementAffine(iquad);
    
    int ref_d = mesh.Element(iquad).GetReference();
    TinyVector<Complexe, Dimension::dim_N> vec_u;
    IVect Nodle = vars.GetDofNumberOnElement(iquad);
    IVect Nodle2 = vars.GetDofNumberOnElement(num_elem2);
    
    Vector<Matrixn_n> En_quad(nb_points_face), Dn_quad(nb_points_face);
    Vector<MatrixVectn_n> Cn_quad(nb_points_face);
    VectReal_wp CoefInternal(nb_points_face), CoefExternal(nb_points_face);
    CoefInternal.Fill(1); CoefExternal.Fill(1);
    
    Complexe phase(1), phase_conj(1), cone(1);
    int rot1 = mesh.Element(iquad).GetOrientationBoundary(num_pos1_face), rot2 = 0;
            
    int cond = vars.mesh.GetBoundaryCondition(ref);
    bool face_on_gamma = false;
    if (cond != BoundaryConditionEnum::LINE_INSIDE)
      {
        if (cond != BoundaryConditionEnum::LINE_NEIGHBOR)
          face_on_gamma = true;
      }
    
    int num_pos2_face = -1;
    int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0, ref2 = -1;
    if (num_elem2 != iquad)
      {
	num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(num_face);
	if (num_pos2_face < 0)
	  {
	    int nf2 = mesh_num.GetPeriodicBoundary(num_face);
            vars.GetPeriodicPhase(num_face, phase);
            vars.GetPeriodicPhase(nf2, phase_conj);
            num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(nf2);
	  }
 
        ref2 = mesh.Element(num_elem2).GetReference();
	rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
	rot = mesh_num.GetRotationFace(rot1, rot2, nv);
      }
    
    int krot, num_point = Fb.GetNbPointsQuadratureInside();
    for (int k = 0; k < num_pos1_face; k++)
      num_point += mesh_num.GetNbPointsQuadratureBoundary(mesh.Element(iquad).numBoundary(k));

    const ElementReference<Dimension, 1>* Fb2_ptr = &vars.GetReferenceElementH1(num_elem2);

    bool neighbor_face = (cond == BoundaryConditionEnum::LINE_NEIGHBOR);
    IVect NumRotQuad(nb_points_face);
#ifdef SELDON_WITH_MPI
    int proc2 = -1;
    if (neighbor_face)
      {
	int type_per = mesh_num.GetPeriodicityTypeForBoundary(num_face);
	if (type_per >= 0)
	  {
	    vars.GetPeriodicPhase(num_face, phase);
	    phase_conj = conjugate(phase);
	  }

        proc2 = vars.GetProcessorNeighboringFace(num_face);
        num_pos2_face = vars.GetLocalPositionNeighboringFace(num_face);
        int pos_loc_neighbor = mesh_num.GetLocalEdgeNumberNeighborElement(num_face);
        Fb2_ptr = dynamic_cast<const ElementReference<Dimension, 1>* >(&vars.GetNeighborReferenceElement(pos_loc_neighbor));
        rot = vars.GetRotationNeighboringFace(num_face);
        Nodle2 = vars.GetNodleNeighboringFace(num_face);
        ref2 = vars.GetRefDomainNeighboringFace(num_face);
        NumRotQuad.Reallocate(nb_points_face);
      }
#endif

    const ElementReference<Dimension, 1>& Fb2 = *Fb2_ptr;
    
    // loop over quadrature points of the face
    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
    
    for (int k = 0; k < nb_points_face; k++)
      {
        krot = FacesQuadRotation(rot, k);
        
        // surface element of integration ds and normale
        if (new_face)
          {
            normale = vars.Glob_normale(num_face)(k);
            dsj = vars.Glob_dsj(num_face)(k);
          }
        else
          {
            normale = vars.Glob_normale(num_face)(krot);
            dsj = vars.Glob_dsj(num_face)(krot);
            Mlt(Real_wp(-1), normale);
          }
        
        // computation of D n and E n
        TypeEquation::GetGradPhiTensor(vars, iquad, num_point, nat_mat, ref_d, Di, Ei);
        if (vars.InsidePML(iquad))
          TypeEquation::ApplyPmlGradPhi(vars, iquad, num_point, nat_mat, ref_d, Di, Ei);
        
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
            {
              Dn_quad(k)(m, n) = DotProd(Di(m, n), normale);
              En_quad(k)(m, n) = DotProd(Ei(m, n), normale);
            }
        
        // P and N are added to Dn
        if (!face_on_gamma)
          {
            // we provide the matrix D = \sum n_i (E_i - D_i)
            // which is useful (since for upwind fluxes the absolue value of this matrix is required)
            Ncl = En_quad(k) - Dn_quad(k);
            TypeEquation::GetPenalDG(Ncl, normale, iquad, num_point,
				     num_face, nat_mat, ref_d, ref2, vars, Fb);
            
            if ((vars.automatic_choice_penalization) && (!first_order))
              {
                // interior penalty formulation, the penalty coefficient is well known
                // usually alpha r(r+1)/h
                Ncl *= vars.Glob_CoefPenalDG(num_face);
              }
            
            Dn_quad(k) += Ncl;
          }
        
        if (first_order)
          {
            if (face_on_gamma)
              {
                // boundary condition is given by D u = Ncl u
                // where D = \sum n_i (E_i - D_i)
                // we provide the matrix D and we get the matrix Ncl in return
                Ncl = En_quad(k) - Dn_quad(k);
                TypeEquation::GetNabc(Ncl, normale, ref, iquad,
                                      num_point, nat_mat, ref_d, vars, Fb);
                
                Mlt(-Real_wp(1), Ncl);
                Dn_quad(k) += Ncl;
              }
          }
        else
          {
            // computing \sum_{p,q} C^{m,n}_{p, q} b_{j, q} n_p  (where b = DF_i*-1)
            TypeEquation::GetGradGradTensor(vars, iquad, num_point, nat_mat, ref_d, Ci);
            if (vars.InsidePML(iquad))
              TypeEquation::ApplyPmlGradGrad(vars, iquad, num_point, nat_mat, ref_d, Ci);
            
            if (affine)
              {
                dfjm1 = vars.Glob_DFjm1(iquad)(0);
                Mlt(1.0/vars.Glob_jacobian(iquad)(0), dfjm1);
              }
            else
              dfjm1 = vars.Glob_DFjm1(iquad)(num_point);
            
            for (int j = 0; j < Dimension::dim_N; j++)
              GetRow(dfjm1, j, bvec(j));
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                {
                  for (int j = 0; j < Dimension::dim_N; j++)
                    {
                      Mlt(Ci(m, n), bvec(j), vec_u);
                      Cn_quad(k)(m, n)(j) = DotProd(vec_u, normale);
                    }
                }
            
          }
        
        // multiplying by 1/2 omega_k ds
        coef = dsj*PoidsFlux(k);
        Cn_quad(k) *= coef; Dn_quad(k) *= coef; En_quad(k) *= coef;	
        
        if (ElementReference_Base::use_warburton_trick)
          {
            if (!affine)
              {
                int num_scal = k + Fb.GetOffsetSh(num_pos1_face);
                int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();
                int num_pt_quad = offset_war + num_scal;
                Real_wp invSqrtJacob = mesh.Glob_invSqrtJacobian(iquad)(num_pt_quad);
                CoefExternal(k) = invSqrtJacob;
                CoefInternal(k) = invSqrtJacob*invSqrtJacob;
              }
            
            if (!vars.mesh.IsElementAffine(num_elem2))
              {
                int num_scal2 = krot + Fb2.GetOffsetSh(num_pos2_face);
                int num_pt_quad = num_scal2;
                Real_wp invSqrtJacob = mesh.Glob_invSqrtJacobianNeighbor(nb_neighbor)(num_pt_quad);
                CoefExternal(k) *= invSqrtJacob;
              }
          }
        
        NumRotQuad(k) = krot;
        num_point++;
      }
    
    Vector<Complexe> contrib(Fb.GetNbDof()), contrib2(Fb2.GetNbDof());
    contrib.Fill(0); contrib2.Fill(0);
    VectReal_wp val_phi(nb_points_face), val_phi_quad(Fb.GetNbPointsQuadratureInside());
    Vector<Complexe> feval(nb_points_face), feval_quad(Fb.GetNbPointsQuadratureInside());
    Vector<Complexe> feval_quad2(Fb2.GetNbPointsQuadratureInside());
    VectReal_wp Ones(Fb.GetNbDof());
    feval_quad.Fill(0); feval_quad2.Fill(0);
    feval.Fill(0); val_phi.Fill(0); val_phi_quad.Fill(0);
    
    VectReal_wp grad_phi;
    if (!first_order)
      grad_phi.Reallocate(Dimension::dim_N*nb_points_face);
    
    // loop over dofs
    for (int i = 0; i < Fb.GetNbDof(); i++)
      {
        Ones.Fill(0); Ones(i) = 1.0;
        if (Fb.UseQuadraturePointsForSh())
          {
            Fb.ApplyChTranspose(Ones, val_phi_quad);
            Fb.ApplyShQuadratureTranspose(num_pos1_face, val_phi_quad, val_phi, rf);
            if (!first_order)
              Fb.ApplyNablaShQuadratureTranspose(num_pos1_face, val_phi_quad, grad_phi, rf);
          }
        else
          {
            Fb.ApplyShTranspose(num_pos1_face, Ones, val_phi, rf);
            if (!first_order)
              Fb.ApplyNablaShTranspose(num_pos1_face, Ones, grad_phi, rf);
          }
        
        if (Norm2(val_phi) > 10.0*epsilon_machine)
          {
            // double loop over unknowns
            for (int p = 0; p < TypeEquation::nb_unknowns; p++)
              for (int q = 0; q < TypeEquation::nb_unknowns; q++)
                {
                  /*********************
                   * First order terms *
                   *********************/
                  
                  // internal contribution
                  
                  // first treating - (Dn_quad u)- phi-  - (En_quad phi)- u- 
                  for (int k = 0; k < nb_points_face; k++)
                    feval(k) = val_phi(k)*(Dn_quad(k)(p, q) + En_quad(k)(p, q));
                                        
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
                  
                  for (int j = 0; j < Fb.GetNbDof(); j++)
                    if (abs(contrib(j)) > vars.GetThresholdMatrix())
                      {
                        int dofp = offset_row + Nodle(i) + vars.offset_dof_unknown(p);
                        int dofq = offset_col + Nodle(j) + vars.offset_dof_unknown(q);
                        mat_sp.AddInteraction(dofp, dofq, -contrib(j));
                      }
                  
                  if (!face_on_gamma)
                    {
                      // external contribution
                      
                      // first treating (Dn_quad u)- phi+
                      for (int k = 0; k < nb_points_face; k++)
                        feval(NumRotQuad(k)) = val_phi(k)*Dn_quad(k)(p, q);
                      
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
                            if (abs(contrib2(j)) > vars.GetThresholdMatrix())
                              {
                                int dofp = offset_row + Nodle2(j) + vars.GetOffsetGlobalUnknownDof(p);
                                int dofq = offset_col + Nodle(i) + vars.offset_dof_unknown(q);
                                mat_sp.AddRowDistantInteraction(dofp, dofq, proc2,
                                                                phase_conj*contrib2(j));
                              }
#endif
                        }
                      else
                        {
                          for (int j = 0; j < Fb2.GetNbDof(); j++)
                            if (abs(phase_conj*contrib2(j)) > vars.GetThresholdMatrix())
                              {
                                int dofp = offset_row + Nodle2(j) + vars.offset_dof_unknown(p);
                                int dofq = offset_col + Nodle(i) + vars.offset_dof_unknown(q);
				mat_sp.AddInteraction(dofp, dofq, phase_conj*contrib2(j));
                              }
                        }
                      
                      // then (En_quad phi)- u+ 
                      for (int k = 0; k < nb_points_face; k++)
                        feval(NumRotQuad(k)) = val_phi(k)*En_quad(k)(p, q);
                      
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
                            if (abs(contrib2(j)) > vars.GetThresholdMatrix())
                              {
                                int dofp = offset_row + Nodle(i) + vars.offset_dof_unknown(p);
                                int dofq = offset_col + Nodle2(j) + vars.GetOffsetGlobalUnknownDof(q);
                                mat_sp.AddDistantInteraction(dofp, dofq, proc2, phase*contrib2(j));
                              }
#endif
                        }
                      else
                        {
                          for (int j = 0; j < Fb2.GetNbDof(); j++)
                            if (abs(phase*contrib2(j)) > vars.GetThresholdMatrix())
                              {
                                int dofp = offset_row + Nodle(i) + vars.offset_dof_unknown(p);
                                int dofq = offset_col + Nodle2(j) + vars.offset_dof_unknown(q);
				mat_sp.AddInteraction(dofp, dofq, phase*contrib2(j));
                              }
                        }
                    }
                }
          }
        
        
        if (Norm2(grad_phi) > 10.0*epsilon_machine)
          {
            // double loop over unknowns
            for (int p = 0; p < TypeEquation::nb_unknowns; p++)
              for (int q = 0; q < TypeEquation::nb_unknowns; q++)
                {
                  if (!face_on_gamma)
                    {
                      
                      /**********************
                       * Second order terms *
                       **********************/
                      
                      if (!first_order)
                        {                          
                          // internal contribution -(C nabla u n)-  phi-  and -(C nabla phi n)- u-
                          for (int k = 0; k < nb_points_face; k++)
                            {
                              CopyVector(grad_phi, k, vec_u);
                              feval(k) = DotProd(Cn_quad(k)(p, q), vec_u);
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
                          
                          for (int j = 0; j < Fb.GetNbDof(); j++)
                            if (abs(contrib(j)) > vars.GetThresholdMatrix())
                              {
                                int dofp = Nodle(i) + vars.offset_dof_unknown(q);
                                int dofq = Nodle(j) + vars.offset_dof_unknown(p);
                                if (vars.sipg_formulation)
                                  {
                                    // SIPG formulation
                                    mat_sp.AddInteraction(offset_row + dofp, offset_col + dofq,
                                                          -contrib(j));
                                    mat_sp.AddInteraction(offset_row + dofq, offset_col + dofp,
                                                          -contrib(j));
                                  }
                                else
                                  {
                                    // NIPG formulation
                                    mat_sp.AddInteraction(offset_row + dofp, offset_col + dofq,
                                                          contrib(j));
                                    mat_sp.AddInteraction(offset_row + dofq, offset_col + dofp,
                                                          -contrib(j));
                                  }
                              }
                          
                          // external contribution (C nabla u n)-  phi+  and (C nabla phi n)- u+
                          for (int k = 0; k < nb_points_face; k++)
                            {
                              CopyVector(grad_phi, k, vec_u);
                              feval(NumRotQuad(k)) = DotProd(Cn_quad(k)(p, q), vec_u);
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
                                if (abs(contrib2(j)) > vars.GetThresholdMatrix())
                                  {
                                    int dofp = Nodle(i) + vars.offset_dof_unknown(q);
                                    int dofq = Nodle2(j) + vars.GetOffsetGlobalUnknownDof(p);
                                    if (vars.sipg_formulation)
                                      {
                                        mat_sp.AddDistantInteraction(offset_row + dofp,
                                                                     offset_col + dofq, proc2,
                                                                     phase*contrib2(j));
                                        
                                        mat_sp.AddRowDistantInteraction(offset_row + dofq,
                                                                        offset_col + dofp,
                                                                        proc2,
                                                                        phase_conj*contrib2(j));
                                      }
                                    else
                                      {
                                        mat_sp.AddDistantInteraction(offset_row + dofp,
                                                                     offset_col + dofq, proc2,
                                                                     -phase*contrib2(j));
                                        
                                        mat_sp.AddRowDistantInteraction(offset_row + dofq,
                                                                        offset_col + dofp, proc2,
                                                                        phase_conj*contrib2(j));
                                      }
                                  }
#endif
                            }
                          else
                            {
                              for (int j = 0; j < Fb2.GetNbDof(); j++)
                                if (abs(contrib2(j)) > vars.GetThresholdMatrix())
                                  {
                                    int dofp = Nodle(i) + vars.offset_dof_unknown(q);
                                    int dofq = Nodle2(j) + vars.offset_dof_unknown(p);
                                    if (vars.sipg_formulation)
                                      {
                                        // SIPG formulation
                                        mat_sp.AddInteraction(offset_row + dofp,
                                                              offset_col+ dofq, phase*contrib2(j));
                                        
                                        mat_sp.AddInteraction(offset_row + dofq,
                                                              offset_col + dofp,
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
          }
      } 
    
    if (neighbor_face)
      nb_neighbor++;
  }
    
  
  //! fluxes for DG methods on triangles/quadrilaterals/hexahedrals
  /*!
    surfacic terms  \f$ + \displaystyle \sum \int_{\partial K} A_i n_i [u] \varphi 
    + \sum \int_{\partial K} B_i n_i {u} \varphi + \sum \int_\Gamma N u \varphi  \f$ 
    of the discontinuous Galerkin formulation
    \param[in] mat_sp the terms are added to this sparse matrix
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] vars considered problem
   */
  template<class Complexe, class TypeEquation>
  void AddElementaryFluxesDG(VirtualMatrix<Complexe>& mat_sp,
			     const GlobalGenericMatrix<Complexe>& nat_mat,
			     const EllipticProblem<TypeEquation>& vars,
                             int offset_row, int offset_col)
  {
    int nb_neighbor = 0;
    // loop over all the faces
    // for sake of simplicity, we do a loop on elements and then on local faces
    // Hence, faces may be selected twice (because they can belong to two elements)
    // but only interactions 11 and 12 are added (and not the interactions 21 and 22)
    // if 1 and 2 refers to the element 1 and element 2 adjoining the face
    // so that each value is added only once
    for (int iquad = 0; iquad < vars.mesh.GetNbElt(); iquad++)
      for (int num_pos1_face = 0; num_pos1_face < vars.mesh.Element(iquad).GetNbBoundary();
           num_pos1_face++)
	{
	  int num_face = vars.mesh.Element(iquad).numBoundary(num_pos1_face);
          
          if (vars.FaceHasToBeConsideredForBoundaryIntegral(num_face))
            {
              int num_elem2 = vars.mesh.Boundary(num_face).numElement(0);
              int ref = vars.mesh.Boundary(num_face).GetReference();
              bool new_face = vars.IsNewFace(iquad)(num_pos1_face);
              if ((num_elem2 == iquad)&&(vars.mesh.Boundary(num_face).GetNbElements()==2))
                num_elem2 = vars.mesh.Boundary(num_face).numElement(1);

	      bool neighbor_face = (vars.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR);
	      
	      if (!neighbor_face)
		{
		  // for second order formulation, boundary condition is set with AddMatrixImpedanceBoundary
		  if ( (vars.mesh.GetBoundaryCondition(ref) != 0) && (!TypeEquation::FirstOrderFormulation))
		    continue;
		}
	      
              // AddElementaryFluxDG is called for the selected face
	      AddElementaryFluxDG(mat_sp, nat_mat, vars, iquad, num_pos1_face, num_elem2,
				  num_face, ref, new_face, nb_neighbor,
				  offset_row, offset_col, vars.GetReferenceElementH1(iquad));
            }
        }
  }

} // end namespace

#define MONTJOIE_FILE_ELEMENTARY_MATRIX_H1_CXX
#endif
