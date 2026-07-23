#ifndef MONTJOIE_FILE_ELEMENTARY_MATRIX_HCURL_2D_CXX

namespace Montjoie
{
  
  //! computation of elementary matrix for edge finite elements
  /*!
    \param[in] iquad element number
    \param[in] B mass coefficient
    \param[in] A stiffness coefficient
    \param[in,out] mat matrix to fill
    \param[in] Fb finite element object
    \param[in] nat_mat additional coefficients
    \param[in] vars given problem
    computation of \f$ -\omega^2 \int_K \varepsilon \varphi_i \varphi_j
    + \int_K \mu^{-1} \nabla \times \varphi_i \nabla \times \varphi_j \f$
    After change of variables from K to \f$ \hat{K}\f$, we get :
    \f$ -\omega^2 \int_{\hat{K}} J_i DF_i^{-1} \varepsilon DF_i^{*-1} \varphi_i \varphi_j 
    + \int_{\hat{K}} \frac{J_i}{\mu} \nabla \times \varphi_i \nabla \times \varphi_j \f$
    mass coefficient (2x2 matrix) \f$ B =  J_i DF_i^{-1} \varepsilon DF_i^{*-1} \f$
    stiffness coefficient (single coefficient) \f$ A =  \frac{J_i}{\mu} \f$
   */
  template<class Vect1, class Vect2, class Matrix1, class Complexe, class TypeEquation>
  void ComputeVariableElementaryMatrix(int iquad, const Vect1& B, const Vect2& A, Matrix1& mat,
				       const ElementReference<Dimension2, 2>& Fb,
                                       const GlobalGenericMatrix<Complexe>& nat_mat,
                                       const EllipticProblem<TypeEquation>& vars)
    
  {
    int nb_dof_elt = Fb.GetNbDof();
    int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
    
    VectReal_wp val_phi(2*nb_points_quadrature), Ones(nb_dof_elt);
    Vector<Complexe> feval(2*nb_points_quadrature);
    Vector<Complexe> feval_curl(nb_points_quadrature);
    VectReal_wp curl_phi(nb_points_quadrature);
    Vector<Complexe> contrib(nb_dof_elt), contrib_curl(nb_dof_elt);
    TinyVector<Complexe, 2> vec_u, vec_v;
    Complexe one; SetComplexOne(one);    
    
    // loop on dofs and quadrature points
    for (int i = 0; i < nb_dof_elt; i++)
      {
        // we compute phi_i(\xi_k) and curl phi_i for all k
        Ones.Fill(0); Ones(i) = 1.0;
        Fb.ApplyChTranspose(Ones, val_phi);
        Fb.ApplyRhTranspose(Ones, curl_phi);
        
        // application of A and B
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
            {
              for (int k = 0; k < nb_points_quadrature; k++)
                {
                  vec_u(0) = val_phi(2*k);
                  vec_u(1) = val_phi(2*k+1);
                  
                  Mlt(B(n, m)(k), vec_u, vec_v);
                  
                  feval(2*k) = vec_v(0);
                  feval(2*k+1) = vec_v(1);
                  
                  feval_curl(k) = A(n, m)(k)*curl_phi(k);
                }
              
              // then integration against \nabla phi_j
              Fb.ApplyRh(feval_curl, contrib_curl);
              Fb.ApplyCh(feval, contrib);
              Add(one, contrib_curl, contrib);
              
              // we got a column of the matrix
              for (int j = 0; j < nb_dof_elt; j++)
                mat.SetEntry(j + n*nb_dof_elt, i + m*nb_dof_elt, contrib(j));
            }
      }
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
  template<class Complexe, class TypeEquation>
  void ComputeElementaryMatrix(int iquad, IVect& num_dof,
                               VirtualMatrix<Complexe>& mat_interac,
			       const GlobalGenericMatrix<Complexe>& nat_mat,
			       const EllipticProblem<TypeEquation>& vars,
                               const ElementReference<Dimension2, 2>& Fb)
  {
    // number of integration points
    int nb_points_quadrature = Fb.GetNbPointsQuadratureInside();
    // number of dofs in the face
    int nb_dof_face = Fb.GetNbDof();
    const MeshNumbering<Dimension2>& mesh_num = vars.GetMeshNumbering(0);
    int Nvol = mesh_num.GetNbDof();
    num_dof.Reallocate(nb_dof_face*TypeEquation::nb_unknowns);
    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      for (int j = 0; j < nb_dof_face; j++)
        num_dof(j + m*nb_dof_face) = mesh_num.Element(iquad).GetNumberDof(j) + m*Nvol;
    
    //bool variable = vars.UseNumericalIntegration(iquad);
    //bool affine = vars.mesh.IsElementAffine(iquad);
    
    if (vars.FirstOrderFormulation())
      {
	// for mixed formulation, dof numbers for unknown H
	int nb_dof_all = (nb_dof_face + Fb.GetNbPointsNodalElt())*TypeEquation::nb_unknowns;
	num_dof.Resize(nb_dof_all);
        mat_interac.Reallocate(nb_dof_all, nb_dof_all);
	mat_interac.Zero();
        int offset_loc = nb_dof_face*TypeEquation::nb_unknowns;
        int offset_glob = Nvol*TypeEquation::nb_unknowns;
        int Nh = Fb.GetNbPointsNodalElt()*vars.mesh.GetNbElt();
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
            num_dof(j + m*Fb.GetNbPointsNodalElt() + offset_loc)
              = offset_glob + iquad*Fb.GetNbPointsNodalElt() + j + m*Nh;
      }
    else
      {
	mat_interac.Reallocate(num_dof.GetM(), num_dof.GetM());
	mat_interac.Zero();
      }
    
    TinyMatrix<TinyMatrix<Complexe, General, 2, 2>, General,
      TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> mat_mass;
    TinyMatrix<Complexe, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> mat_stiff;
    //DISP(vars.UseNumericalIntegration(iquad));
    if (vars.UseNumericalIntegration(iquad))
      {
	// variable coefficients, need of a numerical quadrature
	
	// mass and stiffness coefficients are computed on each quadrature point
	// mass coefficient (2x2 matrix) B =  DF_i^{-1} \epsilon DF_i^{*-1}
	// stiffness coefficient (single coefficient) A =  J_i/mu
      
	TinyMatrix<Vector<TinyMatrix<Complexe, General, 2, 2> >, General,
          TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> B;
	TinyMatrix<Vector<Complexe>, General,
          TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> A;
		
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
            {
              B(m, n).Reallocate(nb_points_quadrature);
              A(m, n).Reallocate(nb_points_quadrature);
              for (int k = 0; k < nb_points_quadrature; k++)
                {
                  vars.GetMassMatrixCoef(iquad, k, nat_mat, m, n, B(m, n)(k));
                  vars.GetStiffMatrixCoef(iquad, k, nat_mat, m, n, A(m, n)(k));
                }
            }
        
	// computation of matrix, by using the coefficients A and B
	if (Fb.OptimizedComputationElementaryMatrix())
          {
            for (int p = 0; p < TypeEquation::nb_unknowns; p++)
              for (int q = 0; q < TypeEquation::nb_unknowns; q++)
                {
                  Fb.AddVariableMassMatrix(nb_dof_face*p, nb_dof_face*q, B(p, q), mat_interac);
                  Fb.AddVariableStiffnessMatrix(nb_dof_face*p, nb_dof_face*q,
                                                A(p, q), mat_interac);
                }
          }
        else
          ComputeVariableElementaryMatrix(iquad, B, A, mat_interac, Fb, nat_mat, vars);
      }
    else
      {
	// constant coefficients in the integral (for triangles)
	// we can use pre-computed matrices
	// evaluation of A and B only in the first point 
        // (because they are constant on all the element)
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
            {
              vars.GetMassMatrixCoef(iquad, 0, nat_mat, m, n, mat_mass(m, n));
              vars.GetStiffMatrixCoef(iquad, 0, nat_mat, m, n, mat_stiff(m, n));
            }
	
	// adding interactions to the matrix
        if (vars.FirstOrderFormulation())
          {
            TinyMatrix<Vector<TinyMatrix<Complexe, General, 2, 2> >, General,
              TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> B;
            TinyMatrix<Vector<Complexe>, General, TypeEquation::nb_unknowns,
              TypeEquation::nb_unknowns> A;
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                {
                  B(m, n).Reallocate(nb_points_quadrature);
                  A(m, n).Reallocate(nb_points_quadrature);
                  for (int k = 0; k < nb_points_quadrature; k++)
                    {
                      B(m, n)(k) = mat_mass(m, n); Mlt(Fb.WeightsND(k), B(m, n)(k));
                      A(m, n)(k) = mat_stiff(m, n)*Fb.WeightsND(k);
                    }
                }
            
            if (Fb.OptimizedComputationElementaryMatrix())
              {
                for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                  for (int q = 0; q < TypeEquation::nb_unknowns; q++)
                    {
                      Fb.AddVariableMassMatrix(nb_dof_face*p, nb_dof_face*q,
                                               B(p, q), mat_interac);
                      Fb.AddVariableStiffnessMatrix(nb_dof_face*p, nb_dof_face*q,
                                                    A(p, q), mat_interac);
                    }
              }
            else
              ComputeVariableElementaryMatrix(iquad, B, A, mat_interac, Fb, nat_mat, vars);
          }
        else
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                {
                  Fb.AddConstantMassMatrix(m*nb_dof_face, n*nb_dof_face,
                                           mat_mass(m, n), mat_interac);
                  Fb.AddConstantStiffnessMatrix(m*nb_dof_face, n*nb_dof_face,
                                                mat_stiff(m, n), mat_interac);
                }
          }
      }
    
    // modification of the matrix if global dofs are obtained by linear combination of local dofs
    mesh_num.number_map.ModifyLocalRowMatrix(mesh_num, mat_interac, iquad);
    mesh_num.number_map.ModifyLocalColumnMatrix(mesh_num, mat_interac, iquad);    
    //DISP(mat_interac);
  }    

} // end namespace

#define MONTJOIE_FILE_ELEMENTARY_MATRIX_HCURL_2D_CXX
#endif
