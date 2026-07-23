#ifndef MONTJOIE_FILE_PROD_MAT_VECT_MAXWELL2D_CXX

namespace Montjoie
{


  //! matrix vector product by stiffness matrix for Nedelec's first family on quadrilaterals
  /*!
    \param[in] alpha coefficient
    \param[in] level matrix-vector product can be done for a subset of elements
    \param[in] var_hyper given problem
    \param[in] B vector to be multiplied
    \param[in] beta coefficient
    \param[in,out] C result vector
    alpha*A*B + beta*C -> C
  */
  template<class TypeEquation>
  void MltAddStiffnessScalarMaxwell2D(const Real_wp& alpha, int level,
				      const HyperbolicProblem<TypeEquation>& var_hyper,
				      const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C)
  {
    // B is H field, C is Prod_E
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var_hyper.var_harmonic;
    
    if (beta == Real_wp(0))
      C.Fill(0);
    else
      Mlt(beta, C);
    
    // local E_h
    const Mesh<Dimension2>& mesh = vars.mesh;
    const MeshNumbering<Dimension2>& mesh_num = vars.mesh_num;
    VectReal_wp H, Prod_Eh;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension2, 2>& Fb = vars.GetReferenceElement(i);
	int nb_dof_face = Fb.GetNbDof(), num_dof;
	//int nb_points_gauss = Fb.GetNbPointsQuadratureStiffness();
	int nb_points_gauss = 0; abort();
	H.Reallocate(nb_points_gauss);
	Prod_Eh.Reallocate(nb_dof_face);
	
	// we get H
	for (int j = 0; j < nb_points_gauss; j++)
	  { //H(j) = B(var.GetOffsetDofV(i)+j)*Fb.WeightsStiffness(j);	
	    abort(); }
	
	// step Prod_Eh = R_h H
	Fb.ApplyRh(H, Prod_Eh);
	
	// we fill vector C, we modify the sign if necessary
	for (int j = 0; j < nb_dof_face; j++)
	  {
	    num_dof = mesh_num.Element(i).GetNumberDof(j);
	    if (vars.IsNegativeSign(i,j))
	      C(num_dof) += alpha*Prod_Eh(j);
	    else
	      C(num_dof) -= alpha*Prod_Eh(j);
	  }
      }
    
    vars.ImposeNullDirichletCondition(C);    
  }
  
  
  //! matrix vector product by stiffness matrix for Nedelec's first family on quadrilaterals
  /*!
    \param[in] alpha coefficient
    \param[in] level matrix-vector product can be done for a subset of elements
    \param[in] var_hyper given problem
    \param[in] B vector to be multiplied
    \param[in] beta coefficient
    \param[in,out] C result vector
    alpha*A*B + beta*C -> C
  */
  template<class TypeEquation>
  void MltAddStiffnessVectorialMaxwell2D(const Real_wp& alpha, int level,
					 const HyperbolicProblem<TypeEquation>& var_hyper,
					 const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C)
  {
    // B is E field, C is Prod_H
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var_hyper.var_harmonic;
    
    if (beta == Real_wp(0))
      C.Fill(0);
    else
      Mlt(beta, C);
    
    // local E_h
    const Mesh<Dimension2>& mesh = vars.mesh;
    const MeshNumbering<Dimension2>& mesh_num = vars.mesh_num;
    VectReal_wp Eh_dof, H;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension2, 2>& Fb = vars.GetReferenceElement(i);
	int nb_dof_face = Fb.GetNbDof(), num_dof;
	//int nb_points_gauss = Fb.GetNbPointsQuadratureStiffness();
	int nb_points_gauss = 0;
	abort();
	
	// we fill vector E_h, we modify the sign if necessary
	for (int j = 0; j < nb_dof_face; j++)
	  {
	    num_dof = mesh_num.Element(i).GetNumberDof(j);
	    if (vars.IsNegativeSign(i,j))
	      Eh_dof(j) = -B(num_dof);
	    else
	      Eh_dof(j) = B(num_dof);
	  }
	
	// step H = R_h^t E
	Fb.ApplyRhTranspose(Eh_dof, H);
	
	// we fill C
	for (int j = 0; j < nb_points_gauss; j++)
	  {
	    abort();
	    //C(var.GetOffsetDofV(i)+j) -= alpha*Fb.WeightsStiffness(j)*H(j);
	  }
      }    
  }
  
} // end namespace

#define MONTJOIE_FILE_PROD_MAT_VECT_MAXWELL2D_CXX
#endif
