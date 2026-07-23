#ifndef MONTJOIE_FILE_OPT_PROD_MAT_VECT_NODAL_H1_CXX

namespace Montjoie
{
  
  //! optimized matrix vector product for Helmholtz (or Laplace equation)
  template<class T0, class TypeEquation, class Complexe>
  void MltAdd_Optimized(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			const FemMatrixFreeClass<Complexe,TypeEquation>& A,
			int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                        const ElementReference<Dimension, 1>& Fb, const ElementReference<Dimension, 1>& Fbr)
  {
    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension>& mesh = var.mesh;
    typedef TinyMatrix<Complexe, Symmetric,
      Dimension::dim_N, Dimension::dim_N> MatrixN_Nsym_Complexe;
    
    TinyVector<Complexe, Dimension::dim_N> tmp, vh_loc;
    //Complexe mass;
    bool variable = var.UseNumericalIntegration(i);
    //bool affine = mesh.IsElementAffine(i);
        
    bool presence_mass = (nat_mat.GetCoefMass()!= Complexe(0));
    bool presence_stiff = (nat_mat.GetCoefStiffness() != Complexe(0));
    int ref = mesh.Element(i).GetReference();
        
    TinyVector<Complexe, TypeElt::nb_dof> Uh_dof, Prod_Uh;
    TinyVector<Complexe, TypeElt::nb_quad> Uh_quad;
    TinyVector<Complexe, TypeElt::nb_quad*Dimension::dim_N> Vh_quad;
    TinyVector<int, TypeElt::nb_dof> num_ddl;

    // we get values of u on the element
    const IVect& Nodle = mesh.Element(i).GetNodle();
    for (int j = 0; j < TypeElt::nb_dof; j++ )
      {
	num_ddl(j) = Nodle(j);
	Uh_dof(j) = B(num_ddl(j));
      }
    
    // interpolation on quadrature points (if necessary)
    if (variable)
      {
	if (presence_mass)
	  {
	    Fbr.ApplyChTranspose(Uh_dof, Uh_quad);
            const Vector<Complexe>& mass = 
              TypeEquation::GetMassCoefficient(A.var, i, ref);
            
	    for (int j = 0; j < TypeElt::nb_quad; j++)
              Uh_quad(j) *= mass(j);
            
	    // projection on dofs (if necessary)
	    Fbr.ApplyCh(Uh_quad, Prod_Uh);
	  }
	else
	  Prod_Uh.Zero();
	
	if (presence_stiff)
	  {
	    // gradient of u
	    Fbr.ApplyRhTranspose(Uh_dof, Vh_quad);
	    
	    const Vector<MatrixN_Nsym_Complexe>& mat_dfj
	      = TypeEquation::GetStiffCoefficient(A.var, i, ref);
	    // we apply geometrical transformation J_i DF_i^{-1} \mu DF_i^{*-1}
	    for (int j = 0; j < TypeElt::nb_quad; j++)
	      { 
		// TypeEquation::GetStiffMatrix(A.var, i, j, 0, 0, ref, Bhm1);
		TinyVectorLoop<Dimension::dim_N>::CopyTinyVector(Vh_quad, j, tmp);
		
		//Mlt(Bhm1, tmp, vh_loc);
		Mlt(mat_dfj(j), tmp, vh_loc);
		
		TinyVectorLoop<Dimension::dim_N>::CopyTinyVector(vh_loc, Vh_quad, j);
	      }
	    
	    // integration againt grad(phi)
	    Fbr.ApplyRh(Vh_quad, Uh_dof);
	    
	    Add(nat_mat.GetCoefStiffness(), Uh_dof, Prod_Uh);
	  }
      }
    else
      {
	if (presence_mass)
	  {
            const Vector<Complexe>& mass
              = TypeEquation::GetMassCoefficient(A.var, i, ref);
            
            Fbr.MltMassMatrix(Uh_dof, Prod_Uh);
	    Mlt(mass(0), Prod_Uh);
	  }
	else
	  Prod_Uh.Zero();
	
	if (presence_stiff)
	  {
	    // gradient of u
	    Fbr.ApplyRhTranspose(Uh_dof, Vh_quad);
	    
	    const Vector<MatrixN_Nsym_Complexe>& mat_dfj
	      = TypeEquation::GetStiffCoefficient(A.var, i, ref);
            
	    const MatrixN_Nsym_Complexe& Bhm1 = mat_dfj(0);  
	    
	    // we apply geometrical transformation J_i DF_i^{-1} \mu DF_i^{*-1}
	    for (int j = 0; j < TypeElt::nb_quad; j++)
	      { 		
		TinyVectorLoop<Dimension::dim_N>::CopyTinyVector(Vh_quad, j, tmp);
		
		Mlt(Bhm1, tmp, vh_loc);
		Mlt(Fbr.WeightsND(j), vh_loc);
		
		TinyVectorLoop<Dimension::dim_N>::CopyTinyVector(vh_loc, Vh_quad, j);
	      }
	    
	    // integration againt grad(phi)
	    Fbr.ApplyRh(Vh_quad, Uh_dof);
	    
	    Add(nat_mat.GetCoefStiffness(), Uh_dof, Prod_Uh);
	  }
      }
    
    for (int j = 0; j < TypeElt::nb_dof; j++)
      C(num_ddl(j)) += alpha*Prod_Uh(j);
  }
  
}

#define MONTJOIE_FILE_OPT_PROD_MAT_VECT_NODAL_H1_CXX
#endif
