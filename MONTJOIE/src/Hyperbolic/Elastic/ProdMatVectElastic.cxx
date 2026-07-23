#ifndef MONTJOIE_FILE_PROD_MAT_VECT_ELASTIC_CXX

namespace Montjoie
{

  //! Performs matrix vector product for stiffness matrix related to scalar unknowns
  //! (an element only)
  /*!
    \param[in] alpha coefficient (C = C + alpha*A*B)
    \param[in] A finite element matrix
    \param[in] i element number
    \param[in] B vector to multiply
    \param[in] C result
    \param[in] Uh intermediary vector
    \param[in] Vh intermediary vector
    \param[in] Fb finite element
   */
  template<class Complexe, class TypeEquation>
  void MltAdd_ElementElasH1_Scalar(const Complexe& alpha,
                                   const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                   int i, const Vector<Complexe>& B, Vector<Complexe>& C,
				   TinyVector<Vector<Complexe>, 2>& Uh,
                                   TinyVector<Vector<Complexe>, 2>& Vh,
				   const ElementReference<Dimension2, 1>& Fb)
  {
    // 2-D case
    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension2>& mesh = var.mesh;
    const MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
    TinyVector<Complexe, 2> tmp, vh_loc;
    int dim_N = tmp.GetM();
    Complexe coef;
    bool affine = mesh.IsElementAffine(i);
    bool pml = var.InsidePML(i);
    bool ortho = (var.OrthogonalElement(i) == 0);    
    if (!var.InsidePML(i))
      ortho = true;
    
    if (!ortho)
      {
        cout << "Not implemented" << endl;
        abort();
      }
    
    // number of degrees of freedom and quadrature points
    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    
    int Nvol = mesh_num.GetNbDof();
    int Nscal = var.offset_dof_unknown(1);
    for (int m = 0; m < 2; m++)
      {
	if (Uh(m).GetM() != nb_dof_elt)
	  {
	    Uh(m).Reallocate(nb_dof_elt);
	    Uh(m).Fill(0);
	  }
	
	Vh(m).Reallocate(dim_N*nb_quad);
      }

    int offset = var.GetOffsetDofV(i);    

    if (!affine)
      for (int j = 0; j < nb_quad; j++)
	{ 	  
	  coef = alpha*Fb.WeightsND(j);
	  const TinyMatrix<Real_wp, General, 2, 2> & dfjm1 = var.Glob_DFjm1(i)(j);
	  
	  int offset_elt = offset + 3*j;
	  // retrieving sigma_xx, sigma_xy
	  tmp(0) = B(offset_elt);
	  tmp(1) = B(offset_elt+1);
          
	  Mlt(dfjm1, tmp, vh_loc);
	  Mlt(coef, vh_loc);
          
	  Vh(0)(2*j) = vh_loc(0);
	  Vh(0)(2*j+1) = vh_loc(1);
          
	  // then sigma_xy, sigma_yy
	  tmp(0) = B(offset_elt+1);
	  tmp(1) = B(offset_elt+2);
          
	  Mlt(dfjm1, tmp, vh_loc);
	  Mlt(coef, vh_loc);
          
	  Vh(1)(2*j) = vh_loc(0);
	  Vh(1)(2*j+1) = vh_loc(1);
          
	  if (pml)
	    {
	      offset_elt = offset + 3*(nb_quad+j);
              
	      // retrieving sigma_xx, sigma_xy
	      tmp(0) = B(offset_elt);
	      tmp(1) = B(offset_elt+1);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(0)(2*j) += vh_loc(0);
	      Vh(0)(2*j+1) += vh_loc(1);
              
	      // then sigma_xy, sigma_yy
	      tmp(0) = B(offset_elt+1);
	      tmp(1) = B(offset_elt+2);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(1)(2*j) += vh_loc(0);
	      Vh(1)(2*j+1) += vh_loc(1);
	    }
	}
    else
      {
	const TinyMatrix<Real_wp, General, 2, 2> &
	  dfjm1 = var.Glob_DFjm1(i)(0);
	
	for (int j = 0; j < nb_quad; j++)
	  { 	  
	    coef = alpha*Fb.WeightsND(j);
            
	    int offset_elt = offset + 3*j;
	    tmp(0) = B(offset_elt);
	    tmp(1) = B(offset_elt+1);
            
	    Mlt(dfjm1, tmp, vh_loc);
	    Mlt(coef, vh_loc);
            
	    Vh(0)(2*j) = vh_loc(0);
	    Vh(0)(2*j+1) = vh_loc(1);
            
	    tmp(0) = B(offset_elt+1);
	    tmp(1) = B(offset_elt+2);
            
	    Mlt(dfjm1, tmp, vh_loc);
	    Mlt(coef, vh_loc);
            
	    Vh(1)(2*j) = vh_loc(0);
	    Vh(1)(2*j+1) = vh_loc(1);
            
	    if (pml)
	      {
		offset_elt = offset + 3*(nb_quad+j);
                
		// retrieving sigma_xx, sigma_xy
		tmp(0) = B(offset_elt);
		tmp(1) = B(offset_elt+1);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(0)(2*j) += vh_loc(0);
		Vh(0)(2*j+1) += vh_loc(1);
                
		// then sigma_xy, sigma_yy
		tmp(0) = B(offset_elt+1);
		tmp(1) = B(offset_elt+2);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(1)(2*j) += vh_loc(0);
		Vh(1)(2*j+1) += vh_loc(1);
	      }
	  }            
      }
    
    // integration against \nabla \varphi
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    
    if (pml)
      {
	Fb.ApplyRhSplit(Vh(0), Uh(0), Uh(1));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		int num_pml = mesh_num.GetDofPML(num_dof);
		C(num_dof) -= Uh(0)(j) + Uh(1)(j);
		C(num_pml+Nvol) -= Uh(0)(j) - Uh(1)(j);
	      }
	  }
	
	Fb.ApplyRhSplit(Vh(1), Uh(0), Uh(1));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		int num_pml = mesh_num.GetDofPML(num_dof);
		C(num_dof+Nscal) -= Uh(0)(j) + Uh(1)(j);
		C(num_pml+Nvol+Nscal) -= Uh(0)(j) - Uh(1)(j);
	      }
	  }
      }
    else
      {
	Fb.ApplyRh(Vh(0), Uh(0));
	Fb.ApplyRh(Vh(1), Uh(1));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		C(num_dof) -= Uh(0)(j);
		C(num_dof+Nscal) -= Uh(1)(j);
	      }
	  }
	
      }        
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  template<class Complexe, class TypeEquation>
  void MltAdd_ElementElasH1_Scalar(const Complexe& alpha,
                                   const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                   int i, const Vector<Complexe>& B, Vector<Complexe>& C,
				   TinyVector<Vector<Complexe>, 3>& Uh,
                                   TinyVector<Vector<Complexe>, 3>& Vh,
				   const ElementReference<Dimension3, 1>& Fb)
  {
    // 3-D case
    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    TinyVector<Complexe, 3> tmp, vh_loc;
    int dim_N = tmp.GetM();
    Complexe coef;
    bool affine = mesh.IsElementAffine(i);
    bool pml = var.InsidePML(i);
    bool ortho = (var.OrthogonalElement(i) == 0);    
    if (!var.InsidePML(i))
      ortho = true;
    
    if (!ortho)
      {
        cout << "Not implemented" << endl;
        abort();
      }
    
    // number of degrees of freedom and quadrature points
    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    
    int Nvol = mesh_num.GetNbDof();
    int Nvol_pml = mesh_num.GetNbDofPML();
    int Nscal = var.offset_dof_unknown(1);
    for (int m = 0; m < 3; m++)
      {
	if (Uh(m).GetM() != nb_dof_elt)
	  {
	    Uh(m).Reallocate(nb_dof_elt);
	    Uh(m).Fill(0);
	  }
	
	Vh(m).Reallocate(dim_N*nb_quad);
      }

    int offset = var.GetOffsetDofV(i);    
    if (!affine)
      for (int j = 0; j < nb_quad; j++)
	{ 	  
	  coef = alpha*Fb.WeightsND(j);
	  const TinyMatrix<Real_wp, General, 3, 3>& dfjm1 = var.Glob_DFjm1(i)(j);
	  
	  int offset_elt = offset + 6*j;
	  // retrieving sigma_xx, sigma_xy, sigma_xz
	  tmp(0) = B(offset_elt);
	  tmp(1) = B(offset_elt+1);
	  tmp(2) = B(offset_elt+2);
          
	  Mlt(dfjm1, tmp, vh_loc);
	  Mlt(coef, vh_loc);
          
	  Vh(0)(3*j) = vh_loc(0);
	  Vh(0)(3*j+1) = vh_loc(1);
	  Vh(0)(3*j+2) = vh_loc(2);
          
	  // then sigma_xy, sigma_yy, sigma_yz
	  tmp(0) = B(offset_elt+1);
	  tmp(1) = B(offset_elt+3);
	  tmp(2) = B(offset_elt+4);
          
	  Mlt(dfjm1, tmp, vh_loc);
	  Mlt(coef, vh_loc);
          
	  Vh(1)(3*j) = vh_loc(0);
	  Vh(1)(3*j+1) = vh_loc(1);
	  Vh(1)(3*j+2) = vh_loc(2);
	  
	  // then sigma_xz, sigma_yz, sigma_zz
	  tmp(0) = B(offset_elt+2);
	  tmp(1) = B(offset_elt+4);
	  tmp(2) = B(offset_elt+5);
          
	  Mlt(dfjm1, tmp, vh_loc);
	  Mlt(coef, vh_loc);
          
	  Vh(2)(3*j) = vh_loc(0);
	  Vh(2)(3*j+1) = vh_loc(1);
	  Vh(2)(3*j+2) = vh_loc(2);
          
	  if (pml)
	    {
	      offset_elt = offset + 6*(nb_quad+j);
              
	      // retrieving sigma_xx, sigma_xy, sigma_xz
	      tmp(0) = B(offset_elt);
	      tmp(1) = B(offset_elt+1);
	      tmp(2) = B(offset_elt+2);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(0)(3*j) += vh_loc(0);
	      Vh(0)(3*j+1) += vh_loc(1);
	      Vh(0)(3*j+2) += vh_loc(2);
              
	      // then sigma_xy, sigma_yy, sigma_yz
	      tmp(0) = B(offset_elt+1);
	      tmp(1) = B(offset_elt+3);
	      tmp(2) = B(offset_elt+4);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(1)(3*j) += vh_loc(0);
	      Vh(1)(3*j+1) += vh_loc(1);
	      Vh(1)(3*j+2) += vh_loc(2);
	      
	      // then sigma_xz, sigma_yz, sigma_zz
	      tmp(0) = B(offset_elt+2);
	      tmp(1) = B(offset_elt+4);
	      tmp(2) = B(offset_elt+5);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(2)(3*j) += vh_loc(0);
	      Vh(2)(3*j+1) += vh_loc(1);
	      Vh(2)(3*j+2) += vh_loc(2);
	      
	      offset_elt = offset + 6*(2*nb_quad+j);
              
	      // retrieving sigma_xx, sigma_xy, sigma_xz
	      tmp(0) = B(offset_elt);
	      tmp(1) = B(offset_elt+1);
	      tmp(2) = B(offset_elt+2);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(0)(3*j) += vh_loc(0);
	      Vh(0)(3*j+1) += vh_loc(1);
	      Vh(0)(3*j+2) += vh_loc(2);
              
	      // then sigma_xy, sigma_yy, sigma_yz
	      tmp(0) = B(offset_elt+1);
	      tmp(1) = B(offset_elt+3);
	      tmp(2) = B(offset_elt+4);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(1)(3*j) += vh_loc(0);
	      Vh(1)(3*j+1) += vh_loc(1);
	      Vh(1)(3*j+2) += vh_loc(2);
	      
	      // then sigma_xz, sigma_yz, sigma_zz
	      tmp(0) = B(offset_elt+2);
	      tmp(1) = B(offset_elt+4);
	      tmp(2) = B(offset_elt+5);
              
	      Mlt(dfjm1, tmp, vh_loc);
	      Mlt(coef, vh_loc);
              
	      Vh(2)(3*j) += vh_loc(0);
	      Vh(2)(3*j+1) += vh_loc(1);
	      Vh(2)(3*j+2) += vh_loc(2);
	      
	    }
	}
    else
      {
	const TinyMatrix<Real_wp, General, 3, 3> & dfjm1 = var.Glob_DFjm1(i)(0);
	
	for (int j = 0; j < nb_quad; j++)
	  { 	  
	    coef = alpha*Fb.WeightsND(j);
            
	    int offset_elt = offset + 6*j;
	    // retrieving sigma_xx, sigma_xy, sigma_xz
	    tmp(0) = B(offset_elt);
	    tmp(1) = B(offset_elt+1);
	    tmp(2) = B(offset_elt+2);
            
	    Mlt(dfjm1, tmp, vh_loc);
	    Mlt(coef, vh_loc);
            
	    Vh(0)(3*j) = vh_loc(0);
	    Vh(0)(3*j+1) = vh_loc(1);
	    Vh(0)(3*j+2) = vh_loc(2);
	    
	    // then sigma_xy, sigma_yy, sigma_yz
	    tmp(0) = B(offset_elt+1);
	    tmp(1) = B(offset_elt+3);
	    tmp(2) = B(offset_elt+4);
            
	    Mlt(dfjm1, tmp, vh_loc);
	    Mlt(coef, vh_loc);
            
	    Vh(1)(3*j) = vh_loc(0);
	    Vh(1)(3*j+1) = vh_loc(1);
	    Vh(1)(3*j+2) = vh_loc(2);
            
	    // then sigma_xz, sigma_yz, sigma_zz
	    tmp(0) = B(offset_elt+2);
	    tmp(1) = B(offset_elt+4);
	    tmp(2) = B(offset_elt+5);
            
	    Mlt(dfjm1, tmp, vh_loc);
	    Mlt(coef, vh_loc);
            
	    Vh(2)(3*j) = vh_loc(0);
	    Vh(2)(3*j+1) = vh_loc(1);
	    Vh(2)(3*j+2) = vh_loc(2);
            
	    if (pml)
	      {
		offset_elt = offset + 6*(nb_quad+j);
                
		// retrieving sigma_xx, sigma_xy, sigma_xz
		tmp(0) = B(offset_elt);
		tmp(1) = B(offset_elt+1);
		tmp(2) = B(offset_elt+2);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(0)(3*j) += vh_loc(0);
		Vh(0)(3*j+1) += vh_loc(1);
		Vh(0)(3*j+2) += vh_loc(2);
                
		// then sigma_xy, sigma_yy, sigma_yz
		tmp(0) = B(offset_elt+1);
		tmp(1) = B(offset_elt+3);
		tmp(2) = B(offset_elt+4);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(1)(3*j) += vh_loc(0);
		Vh(1)(3*j+1) += vh_loc(1);
		Vh(1)(3*j+2) += vh_loc(2);
                
		// then sigma_xz, sigma_yz, sigma_zz
		tmp(0) = B(offset_elt+2);
		tmp(1) = B(offset_elt+4);
		tmp(2) = B(offset_elt+5);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(2)(3*j) += vh_loc(0);
		Vh(2)(3*j+1) += vh_loc(1);
		Vh(2)(3*j+2) += vh_loc(2);
                
		offset_elt = offset + 6*(2*nb_quad+j);
                
		// retrieving sigma_xx, sigma_xy, sigma_xz
		tmp(0) = B(offset_elt);
		tmp(1) = B(offset_elt+1);
		tmp(2) = B(offset_elt+2);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(0)(3*j) += vh_loc(0);
		Vh(0)(3*j+1) += vh_loc(1);
		Vh(0)(3*j+2) += vh_loc(2);
                
		// then sigma_xy, sigma_yy, sigma_yz
		tmp(0) = B(offset_elt+1);
		tmp(1) = B(offset_elt+3);
		tmp(2) = B(offset_elt+4);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(1)(3*j) += vh_loc(0);
		Vh(1)(3*j+1) += vh_loc(1);
		Vh(1)(3*j+2) += vh_loc(2);
                
		// then sigma_xz, sigma_yz, sigma_zz
		tmp(0) = B(offset_elt+2);
		tmp(1) = B(offset_elt+4);
		tmp(2) = B(offset_elt+5);
                
		Mlt(dfjm1, tmp, vh_loc);
		Mlt(coef, vh_loc);
                
		Vh(2)(3*j) += vh_loc(0);
		Vh(2)(3*j+1) += vh_loc(1);
		Vh(2)(3*j+2) += vh_loc(2);
                
	      }                
	  }
      }
    
    // integration against \nabla \varphi
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    
    if (pml)
      {
	Fb.ApplyRhSplit(Vh(0), Uh(0), Uh(1), Uh(2));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		int num_pml = mesh_num.GetDofPML(num_dof);
		C(num_dof) -= Uh(0)(j) + Uh(1)(j) + Uh(2)(j);
		C(num_pml+Nvol) -= Uh(0)(j) - Uh(1)(j);
		C(num_pml+Nvol+Nvol_pml) -= Uh(0)(j) - Uh(2)(j);
	      }
	  }
	
	Fb.ApplyRhSplit(Vh(1), Uh(0), Uh(1), Uh(2));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		int num_pml = mesh_num.GetDofPML(num_dof);
		C(num_dof+Nscal) -= Uh(0)(j) + Uh(1)(j) + Uh(2)(j);
		C(num_pml+Nvol+Nscal) -= Uh(0)(j) - Uh(1)(j);
		C(num_pml+Nvol+Nvol_pml+Nscal) -= Uh(0)(j) - Uh(2)(j);
	      }
	  }
	
	Fb.ApplyRhSplit(Vh(2), Uh(0), Uh(1), Uh(2));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		int num_pml = mesh_num.GetDofPML(num_dof);
		C(num_dof+2*Nscal) -= Uh(0)(j) + Uh(1)(j) + Uh(2)(j);
		C(num_pml+Nvol+2*Nscal) -= Uh(0)(j) - Uh(1)(j);
		C(num_pml+Nvol+Nvol_pml+2*Nscal) -= Uh(0)(j) - Uh(2)(j);
	      }
	  }
      }
    else
      {
	Fb.ApplyRh(Vh(0), Uh(0));
	Fb.ApplyRh(Vh(1), Uh(1));
	Fb.ApplyRh(Vh(2), Uh(2));
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		C(num_dof) -= Uh(0)(j);
		C(num_dof+Nscal) -= Uh(1)(j);
		C(num_dof+2*Nscal) -= Uh(2)(j);
	      }
	  }            
      }
  }
#endif


  //! Performs matrix vector product for stiffness matrix related to vectorial unknowns
  //! (an element only)
  /*!
    \param[in] alpha coefficient (C = C + alpha*A*B)
    \param[in] A finite element matrix
    \param[in] i element number
    \param[in] B vector to multiply
    \param[in] C result
    \param[in] dU_loc intermediary vector
    \param[in] Uh intermediary vector
    \param[in] Vh intermediary vector
    \param[in] Fb finite element
   */  
  template<class Complexe, class TypeEquation>
  void MltAdd_ElementElasH1_Vectorial(const Complexe& alpha,
                                      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                      int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                                      TinyVector<TinyVector<Complexe, 2>, 2>& dU_loc,
				      TinyVector<Vector<Complexe>, 2>& Uh,
                                      TinyVector<Vector<Complexe>, 2>& Vh,
				      const ElementReference<Dimension2, 1>& Fb)
  {
    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension2>& mesh = var.mesh;
    const MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
    TinyVector<TinyVector<Complexe, 2>, 2> epsilon, sigma, epsilonX, epsilonY,
      sigmaX, sigmaY, epsilonZ, sigmaZ;
    
    TinyVector<Complexe, 2> tmp, vh_loc;
    int dim_N = tmp.GetM();
    Complexe coef;
    
    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    bool affine = mesh.IsElementAffine(i);
    bool pml = var.InsidePML(i);
    int ref_domain = mesh.Element(i).GetReference();
    
    for (int m = 0; m < 2; m++)
      {
	Uh(m).Reallocate(nb_dof_elt);
	Uh(m).Fill(0);
	if (Vh(m).GetM() != dim_N*nb_quad)
	  {
	    Vh(m).Reallocate(dim_N*nb_quad);
	    Vh(m).Fill(0);
	  }
      }
    
    //int Nvol = mesh_num.GetNbDof();
    int Nscal = var.offset_dof_unknown(1);
    // we get values of u for each unknown and dof
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    for (int j = 0; j < nb_dof_elt; j++ )
      {
	int num_dof = Nodle(j);
	if (num_dof >= 0)
	  {
	    Uh(0)(j) = B(num_dof);
	    Uh(1)(j) = B(num_dof + Nscal);
	  }
      }
    
    // computation of gradient of each unknown of u
    Fb.ApplyRhTranspose(Uh(0), Vh(0));
    Fb.ApplyRhTranspose(Uh(1), Vh(1));
    
    int offset = var.GetOffsetDofV(i); 
    // loop on each quadrature point
    // geometry and physical properties
    if (!affine)
      for (int j = 0; j < nb_quad; j++)
	{ 
	  coef = alpha*Fb.WeightsND(j);
	  const TinyMatrix<Real_wp, General, 2, 2>& dfjm1 = var.Glob_DFjm1(i)(j);
	  
	  // application of DF_j^*-1 to obtain the real gradient
	  int offset_elt = offset + 3*j;
          
	  tmp(0) = Vh(0)(2*j);
	  tmp(1) = Vh(0)(2*j+1);
          
	  MltTrans(dfjm1, tmp, vh_loc);
          
	  epsilon(0)(0) = coef*vh_loc(0);
	  epsilon(0)(1) = coef*vh_loc(1);
          
	  tmp(0) = Vh(1)(2*j);
	  tmp(1) = Vh(1)(2*j+1);
          
	  MltTrans(dfjm1, tmp, vh_loc);
          
	  epsilon(1)(0) = coef*vh_loc(0);
	  epsilon(1)(1) = coef*vh_loc(1);
          
	  if (pml)
	    {
	      epsilonX(0)(0) = epsilon(0)(0);
	      epsilonX(1)(0) = epsilon(1)(0);
	      epsilonY(0)(1) = epsilon(0)(1);
	      epsilonY(1)(1) = epsilon(1)(1);
              
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonX, sigmaX);
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonY, sigmaY);
              
	      C(offset_elt) += sigmaX(0)(0);
	      C(offset_elt+1) += sigmaX(0)(1);
	      C(offset_elt+2) += sigmaX(1)(1);
              
	      offset_elt += 3*nb_quad;
	      C(offset_elt) += sigmaY(0)(0);
	      C(offset_elt+1) += sigmaY(0)(1);
	      C(offset_elt+2) += sigmaY(1)(1);
	    }
	  else
	    {
	      // application of tensor C to get sigma
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilon, sigma);
              
	      // adding result to C
	      C(offset_elt) += sigma(0)(0);
	      C(offset_elt+1) += sigma(0)(1);
	      C(offset_elt+2) += sigma(1)(1);
	    }
	}
    else
      {
	const TinyMatrix<Real_wp, General, 2, 2>& dfjm1 = var.Glob_DFjm1(i)(0);
	
	for (int j = 0; j < nb_quad; j++)
	  { 
	    coef = alpha*Fb.WeightsND(j);
            
	    // application of DF_j^*-1 to obtain the real gradient
	    int offset_elt = offset + 3*j;
            
	    tmp(0) = Vh(0)(2*j);
	    tmp(1) = Vh(0)(2*j+1);
            
	    MltTrans(dfjm1, tmp, vh_loc);
            
	    epsilon(0)(0) = coef*vh_loc(0);
	    epsilon(0)(1) = coef*vh_loc(1);
            
	    tmp(0) = Vh(1)(2*j);
	    tmp(1) = Vh(1)(2*j+1);
            
	    MltTrans(dfjm1, tmp, vh_loc);
            
	    epsilon(1)(0) = coef*vh_loc(0);
	    epsilon(1)(1) = coef*vh_loc(1);
            
	    // application of tensor C to get sigma
	    var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilon, sigma);
            
	    // adding result to C
	    if (pml)
	      {
		epsilonX(0)(0) = epsilon(0)(0);
		epsilonX(1)(0) = epsilon(1)(0);
		epsilonY(0)(1) = epsilon(0)(1);
		epsilonY(1)(1) = epsilon(1)(1);
                
		var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonX, sigmaX);
		var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonY, sigmaY);
                
		C(offset_elt) += sigmaX(0)(0);
		C(offset_elt+1) += sigmaX(0)(1);
		C(offset_elt+2) += sigmaX(1)(1);
                
		offset_elt += 3*nb_quad;
                
		C(offset_elt) += sigmaY(0)(0);
		C(offset_elt+1) += sigmaY(0)(1);
		C(offset_elt+2) += sigmaY(1)(1);
	      }
	    else
	      {
		C(offset_elt) += sigma(0)(0);
		C(offset_elt+1) += sigma(0)(1);
		C(offset_elt+2) += sigma(1)(1);
	      }
	  }
      }          
  }


#ifdef MONTJOIE_WITH_THREE_DIM
  //! Performs matrix vector product for stiffness matrix related to vectorial unknowns
  //! (an element only)
  /*!
    \param[in] alpha coefficient (C = C + alpha*A*B)
    \param[in] A finite element matrix
    \param[in] i element number
    \param[in] B vector to multiply
    \param[in] C result
    \param[in] dU_loc intermediary vector
    \param[in] Uh intermediary vector
    \param[in] Vh intermediary vector
    \param[in] Fb finite element
   */  
  template<class Complexe, class TypeEquation>
  void MltAdd_ElementElasH1_Vectorial(const Complexe& alpha,
                                      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                      int i, const Vector<Complexe>& B, Vector<Complexe>& C,
                                      TinyVector<TinyVector<Complexe, 3>, 3>& dU_loc,
				      TinyVector<Vector<Complexe>, 3>& Uh,
                                      TinyVector<Vector<Complexe>, 3>& Vh,
				      const ElementReference<Dimension3, 1>& Fb)
  {
    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    TinyVector<TinyVector<Complexe, 3>, 3> epsilon, sigma,
      epsilonX, epsilonY, sigmaX, sigmaY, epsilonZ, sigmaZ;
    TinyVector<Complexe, 3> tmp, vh_loc;
    int dim_N = tmp.GetM();
    Complexe coef;
    
    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    bool affine = mesh.IsElementAffine(i);
    bool pml = var.InsidePML(i);
    int ref_domain = mesh.Element(i).GetReference();
    
    for (int m = 0; m < 3; m++)
      {
	Uh(m).Reallocate(nb_dof_elt);
	Uh(m).Fill(0);
	if (Vh(m).GetM() != dim_N*nb_quad)
	  {
	    Vh(m).Reallocate(dim_N*nb_quad);
	    Vh(m).Fill(0);
	  }
      }
    
    //int Nvol = mesh_num.GetNbDof();
    int Nscal = var.offset_dof_unknown(1);
    // we get values of u for each unknown and dof
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    for (int j = 0; j < nb_dof_elt; j++ )
      {
	int num_dof = Nodle(j);
	if(num_dof >= 0)
	  {
	    Uh(0)(j) = B(num_dof);
	    Uh(1)(j) = B(num_dof + Nscal);
	    Uh(2)(j) = B(num_dof + 2*Nscal);
	  }
      }
    
    // computation of gradient of each unknown of u
    Fb.ApplyRhTranspose(Uh(0), Vh(0));
    Fb.ApplyRhTranspose(Uh(1), Vh(1));
    Fb.ApplyRhTranspose(Uh(2), Vh(2));
    
    int offset = var.GetOffsetDofV(i); 
    // loop on each quadrature point
    // geometry and physical properties
    if (!affine)
      for (int j = 0; j < nb_quad; j++)
	{ 
	  coef = alpha*Fb.WeightsND(j);
	  const TinyMatrix<Real_wp, General, 3, 3>& dfjm1 = var.Glob_DFjm1(i)(j);
	  
	  // application of DF_j^*-1 to obtain the real gradient
	  int offset_elt = offset + 6*j;
          
	  tmp(0) = Vh(0)(3*j);
	  tmp(1) = Vh(0)(3*j+1);
	  tmp(2) = Vh(0)(3*j+2);
          
	  MltTrans(dfjm1, tmp, vh_loc);
          
	  epsilon(0)(0) = coef*vh_loc(0);
	  epsilon(0)(1) = coef*vh_loc(1);
	  epsilon(0)(2) = coef*vh_loc(2);
          
	  tmp(0) = Vh(1)(3*j);
	  tmp(1) = Vh(1)(3*j+1);
	  tmp(2) = Vh(1)(3*j+2);
          
	  MltTrans(dfjm1, tmp, vh_loc);
          
	  epsilon(1)(0) = coef*vh_loc(0);
	  epsilon(1)(1) = coef*vh_loc(1);
	  epsilon(1)(2) = coef*vh_loc(2);
	  
	  tmp(0) = Vh(2)(3*j);
	  tmp(1) = Vh(2)(3*j+1);
	  tmp(2) = Vh(2)(3*j+2);
          
	  MltTrans(dfjm1, tmp, vh_loc);
          
	  epsilon(2)(0) = coef*vh_loc(0);
	  epsilon(2)(1) = coef*vh_loc(1);
	  epsilon(2)(2) = coef*vh_loc(2);
          
	  if (pml)
	    {
	      epsilonX(0)(0) = epsilon(0)(0);
	      epsilonX(1)(0) = epsilon(1)(0);
	      epsilonX(2)(0) = epsilon(2)(0);
	      epsilonY(0)(1) = epsilon(0)(1);
	      epsilonY(1)(1) = epsilon(1)(1);
	      epsilonY(2)(1) = epsilon(2)(1);
	      epsilonZ(0)(2) = epsilon(0)(2);
	      epsilonZ(1)(2) = epsilon(1)(2);
	      epsilonZ(2)(2) = epsilon(2)(2);
              
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonX, sigmaX);
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonY, sigmaY);
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonZ, sigmaZ);
              
	      C(offset_elt) += sigmaX(0)(0);
	      C(offset_elt+1) += sigmaX(0)(1);
	      C(offset_elt+2) += sigmaX(0)(2);
	      C(offset_elt+3) += sigmaX(1)(1);
	      C(offset_elt+4) += sigmaX(1)(2);
	      C(offset_elt+5) += sigmaX(2)(2);
              
	      offset_elt += 6*nb_quad;
	      C(offset_elt) += sigmaY(0)(0);
	      C(offset_elt+1) += sigmaY(0)(1);
	      C(offset_elt+2) += sigmaY(0)(2);
	      C(offset_elt+3) += sigmaY(1)(1);
	      C(offset_elt+4) += sigmaY(1)(2);
	      C(offset_elt+5) += sigmaY(2)(2);
	      
	      offset_elt += 6*nb_quad;
	      C(offset_elt) += sigmaZ(0)(0);
	      C(offset_elt+1) += sigmaZ(0)(1);
	      C(offset_elt+2) += sigmaZ(0)(2);
	      C(offset_elt+3) += sigmaZ(1)(1);
	      C(offset_elt+4) += sigmaZ(1)(2);
	      C(offset_elt+5) += sigmaZ(2)(2);
	    }
	  else
	    {
	      // application of tensor C to get sigma
	      var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilon, sigma);
              
	      // adding result to C
	      C(offset_elt) += sigma(0)(0);
	      C(offset_elt+1) += sigma(0)(1);
	      C(offset_elt+2) += sigma(0)(2);
	      C(offset_elt+3) += sigma(1)(1);
	      C(offset_elt+4) += sigma(1)(2);
	      C(offset_elt+5) += sigma(2)(2);
	    }
	}
    else
      {
	const TinyMatrix<Real_wp, General, 3, 3> & dfjm1 = var.Glob_DFjm1(i)(0);
            
	for (int j = 0; j < nb_quad; j++)
	  { 
	    coef = alpha*Fb.WeightsND(j);
	    int offset_elt = offset + 6*j;
            
	    tmp(0) = Vh(0)(3*j);
	    tmp(1) = Vh(0)(3*j+1);
	    tmp(2) = Vh(0)(3*j+2);
            
	    MltTrans(dfjm1, tmp, vh_loc);
            
	    epsilon(0)(0) = coef*vh_loc(0);
	    epsilon(0)(1) = coef*vh_loc(1);
	    epsilon(0)(2) = coef*vh_loc(2);
            
	    tmp(0) = Vh(1)(3*j);
	    tmp(1) = Vh(1)(3*j+1);
	    tmp(2) = Vh(1)(3*j+2);
            
	    MltTrans(dfjm1, tmp, vh_loc);
            
	    epsilon(1)(0) = coef*vh_loc(0);
	    epsilon(1)(1) = coef*vh_loc(1);
	    epsilon(1)(2) = coef*vh_loc(2);
            
	    tmp(0) = Vh(2)(3*j);
	    tmp(1) = Vh(2)(3*j+1);
	    tmp(2) = Vh(2)(3*j+2);
            
	    MltTrans(dfjm1, tmp, vh_loc);
            
	    epsilon(2)(0) = coef*vh_loc(0);
	    epsilon(2)(1) = coef*vh_loc(1);
	    epsilon(2)(2) = coef*vh_loc(2);
            
	    if (pml)
	      {
		epsilonX(0)(0) = epsilon(0)(0);
		epsilonX(1)(0) = epsilon(1)(0);
		epsilonX(2)(0) = epsilon(2)(0);
		epsilonY(0)(1) = epsilon(0)(1);
		epsilonY(1)(1) = epsilon(1)(1);
		epsilonY(2)(1) = epsilon(2)(1);
		epsilonZ(0)(2) = epsilon(0)(2);
		epsilonZ(1)(2) = epsilon(1)(2);
		epsilonZ(2)(2) = epsilon(2)(2);
                
		var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonX, sigmaX);
		var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonY, sigmaY);
		var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilonZ, sigmaZ);
                
		C(offset_elt) += sigmaX(0)(0);
		C(offset_elt+1) += sigmaX(0)(1);
		C(offset_elt+2) += sigmaX(0)(2);
		C(offset_elt+3) += sigmaX(1)(1);
		C(offset_elt+4) += sigmaX(1)(2);
		C(offset_elt+5) += sigmaX(2)(2);
                
		offset_elt += 6*nb_quad;
		C(offset_elt) += sigmaY(0)(0);
		C(offset_elt+1) += sigmaY(0)(1);
		C(offset_elt+2) += sigmaY(0)(2);
		C(offset_elt+3) += sigmaY(1)(1);
		C(offset_elt+4) += sigmaY(1)(2);
		C(offset_elt+5) += sigmaY(2)(2);
                
		offset_elt += 6*nb_quad;
		C(offset_elt) += sigmaZ(0)(0);
		C(offset_elt+1) += sigmaZ(0)(1);
		C(offset_elt+2) += sigmaZ(0)(2);
		C(offset_elt+3) += sigmaZ(1)(1);
		C(offset_elt+4) += sigmaZ(1)(2);
		C(offset_elt+5) += sigmaZ(2)(2);
	      }
	    else
	      {
		// application of tensor C to get sigma
		var.ref_tensorC(ref_domain).MltMatrix(var, i, j, epsilon, sigma);
                
		// adding result to C
		C(offset_elt) += sigma(0)(0);
		C(offset_elt+1) += sigma(0)(1);
		C(offset_elt+2) += sigma(0)(2);
		C(offset_elt+3) += sigma(1)(1);
		C(offset_elt+4) += sigma(1)(2);
		C(offset_elt+5) += sigma(2)(2);
	      }
	  }
      }          
  }
#endif
  

  //! matrix-vector product for nodal finite element and mixed formulation (with R_h)
  /*!
    \param[in] alpha coefficient
    \param[in] level the operation can be done only on a subset of elements
    \param[in] A black-box matrix
    \param[in] B vector to multiply
    \param[in] beta coefficient
    \param[in,out] C result vector
    alpha*A*B + beta*C -> C
   */
  template <class Complexe, class TypeEquation>
  void MltAdd_SquareElasHex_VectorialH1(const Complexe& alpha, int level,
					const FemMatrixFreeClass<Complexe, TypeEquation>& A,
					const Vector<Complexe>& B, const Complexe& beta, Vector<Complexe>& C)
  {    
    Complexe zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    if (beta == zero)
      C.Fill(zero);
    else if (beta != one)
      Mlt(beta, C);
    
    typedef typename TypeEquation::Dimension Dimension;
    
    typedef TinyVector<Complexe, Dimension::dim_N> R_N_Complexe;
    
    TinyVector<Vector<Complexe>, Dimension::dim_N> Uh, Vh;
    TinyVector<R_N_Complexe, Dimension::dim_N> dU_loc;
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      MltAdd_ElementElasH1_Vectorial(alpha, A, i, B, C, dU_loc,
				     Uh, Vh, A.var.GetReferenceElementH1(i));
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
    // B.Write("bh.dat"); C.Write("ch.dat");
    // int test_input; cout<<"vectorial "<<endl; cin>>test_input;    
  }
  
  
  //! matrix-vector product for nodal finite element and mixed formulation (with R_h^*)
  /*!
    \param[in] alpha coefficient
    \param[in] level the operation can be done only on a subset of elements
    \param[in] A black-box matrix
    \param[in] B vector to multiply
    \param[in] beta coefficient
    \param[in,out] C result vector
    alpha*A*B + beta*C -> C
   */
  template<class Complexe, class TypeEquation>
  void MltAdd_SquareElasHex_ScalarH1(const Complexe& alpha, int level,
				     const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				     const Vector<Complexe>& B, const Complexe& beta, Vector<Complexe>& C)
  {
    Complexe zero, one;
    SetComplexZero(zero); SetComplexOne(one);

    glob_chrono.Start(VirtualTimer::STIFFNESS);
    if (beta == zero)
      C.Fill(0);
    else if (beta != one)
      Mlt(beta, C);
    
    typedef typename TypeEquation::Dimension Dimension;
    TinyVector<Vector<Complexe>, Dimension::dim_N> Uh, Vh;
    
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      MltAdd_ElementElasH1_Scalar(alpha, A, i, B, C,
				  Uh, Vh, A.var.GetReferenceElementH1(i));
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
    
    A.var.ImposeNullDirichletCondition(C);    
  }
    
}

#define MONTJOIE_FILE_PROD_MAT_VECT_ELASTIC_CXX
#endif
