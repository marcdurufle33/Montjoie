#ifndef MONTJOIE_FILE_PROD_MAT_VECT_HCURL_2D_CXX

namespace Montjoie
{

  //! matrix vector product for an element and edge finite elements
  template<class T0, class TypeEquation, class Complexe>
  void MltAdd_ElementHcurl(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			   const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			   int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                           const ElementReference<Dimension2, 2>& Fb)
  {
    const MeshNumbering<Dimension2>& mesh_num = A.var.GetMeshNumbering(0);
    const EllipticProblem<TypeEquation>& vars = A.var;
    
    // local E_h
    int nb_dof_face = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    TinyVector<Vector<Complexe>, 1> Eh_dof;
    Eh_dof(0).Reallocate(nb_dof_face);
    Vector<Complexe> Eh_mass(nb_dof_face);
    Vector<Complexe> Eh_quad(2*nb_points_quad);
    Vector<Complexe> H_quad(nb_points_quad);
    Vector<Complexe> Prod_Eh(2*nb_points_quad);
    Vector<Complexe> Prod_H(nb_points_quad);
    TinyVector<Complexe, 2> Eloc, vloc;
    FillZero(Eh_dof); H_quad.Fill(0); Prod_Eh.Fill(0); Prod_H.Fill(0);
    int Nvol = mesh_num.GetNbDof();
    if (vars.FirstOrderFormulation())
      {
	// we retrieve local values of E
	vars.GetLocalUnknownVector(B, i, Eh_dof);
	
	// we get H on quadrature points
	for (int j = 0; j < nb_points_quad; j++)
	  H_quad(j) = B(Nvol + i*nb_points_quad + j);
	
	// interpolation of E on quadrature points
	Fb.ApplyChTranspose(Eh_dof(0), Eh_quad);
	
	// Prod_H = curl E
	Fb.ApplyRhTranspose(Eh_dof(0), Prod_H); 
	
	Complexe Astiff; TinyMatrix<Complexe, Symmetric, 2, 2> Bmass;
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    // second equation in H : -mass*D_h^{-1} H + stiff* curl E 
	    vars.GetStiffMatrixCoef(i, j, nat_mat, 0, 0, Astiff); 
	    Astiff = -1.0/Astiff;
	    Prod_H(j) = Prod_H(j)*nat_mat.GetCoefStiffness() + Astiff*H_quad(j);
	    
	    // mass matrix term B_h E = \omega_k J_e DF_e^{-1} \varepsilon DF_e^{*-1} 
	    vars.GetMassMatrixCoef(i, j, nat_mat, 0, 0, Bmass);
            
	    Eloc(0) = Eh_quad(2*j); Eloc(1) = Eh_quad(2*j+1);
	    Mlt(Bmass, Eloc, vloc);
	    Prod_Eh(2*j) = vloc(0); Prod_Eh(2*j+1) = vloc(1);
	  }
	
	for (int j = 0; j < nb_points_quad; j++)
	  C(Nvol + i*nb_points_quad + j) += alpha*Prod_H(j);
	
	// computation of curl(H)
	Fb.ApplyRh(H_quad, Eh_dof(0));
	
	// integration agains basis function for mass term
	Fb.ApplyCh(Prod_Eh, Eh_mass);
	for (int j = 0; j < nb_dof_face; j++)
	  Eh_dof(0)(j) += Eh_mass(j);
	
	// we fill vector C, we modify the sign if necessary
	vars.AddLocalUnknownVector(alpha, Eh_dof, i, C);
      }
    else
      {
	// we retrieve local values of E
	vars.GetLocalUnknownVector(B, i, Eh_dof);
	
	// E on quadrature points
	Fb.ApplyChTranspose(Eh_dof(0), Eh_quad);
	
	// curl(E) on quadrature points
	Fb.ApplyRhTranspose(Eh_dof(0), H_quad);
	
	// we take into account geometry
	Complexe Astiff; TinyMatrix<Complexe, Symmetric, 2, 2> Bmass;
	for (int j = 0; j < nb_points_quad; j++)
	  {
	    vars.GetStiffMatrixCoef(i, j, nat_mat, 0, 0, Astiff);
	    H_quad(j) *= Astiff;
	    
	    vars.GetMassMatrixCoef(i, j, nat_mat, 0, 0, Bmass);
            
	    Eloc(0) = Eh_quad(2*j); Eloc(1) = Eh_quad(2*j+1);
	    Mlt(Bmass, Eloc, vloc);
	    Prod_Eh(2*j) = vloc(0); Prod_Eh(2*j+1) = vloc(1);
	  }
	
	// integration against curl(phi_i)
	Fb.ApplyRh(H_quad, Eh_dof(0));
	
	// integration against phi_i
	Fb.ApplyCh(Prod_Eh, Eh_mass);
	for (int j = 0; j < nb_dof_face; j++)
	  Eh_dof(0)(j) += Eh_mass(j);
	
	// we add contribution to vector C
	vars.AddLocalUnknownVector(alpha, Eh_dof, i, C);
      }
  }
  
  
  //! generic matrix-vector product for edge finite elements in 2-D
  //! and Maxwell equations
  /*
    this algorithm separates the geometry \f$ DF_i^{*-1} \f$ from the 
    spatial derivatives. This allows low-storage and fast-algorithm
    in \f$ O(r^4) \f$ instead of \f$ O(r^6) \f$ if the full matrix is stored, where
    r is the order of approximation
  */
  template <class T0, class Complexe, class TypeEquation>
  void MltAddHcurl2D(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
		     const SeldonTranspose& trans, int level,
		     const FemMatrixFreeClass<Complexe, TypeEquation>& A,
		     const Vector<T0>& B2, const T0& beta, Vector<T0>& C2, bool assemble)
  {
#ifdef SELDON_WITH_MPI
    Vector<T0> B(B2), C(C2);
#else
    Vector<T0> B, C;
#endif
    
    A.ApplyRightScaling(B2, C2, B, C);    

    if (beta == T0(0))
      C2.Fill(0);
    else
      Mlt(beta, C2);
        
    // loop over all elements of the mesh
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      {
	MltAdd_ElementHcurl(alpha, nat_mat, A, i, B, C, A.var.GetReferenceElementHcurl(i));
      }
    
    A.AddExtraBoundaryTerms(alpha, B, C);
    
    A.ApplyLeftScaling(B2, C2, B, C);    
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C2);
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_PROD_MAT_VECT_HCURL_2D_CXX
#endif
