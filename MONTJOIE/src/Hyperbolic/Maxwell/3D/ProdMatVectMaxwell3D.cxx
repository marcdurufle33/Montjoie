#ifndef MONTJOIE_FILE_PRODMATVECT_MAXWELL3D_CXX

namespace Montjoie
{
  
#ifdef MONTJOIE_WITH_NODAL_HCURL
  //! matrix vector product by stiffness matrix for Nedelec's second family on hexahedra
  /*!
    \param[in] alpha coefficient
    \param[in] level matrix-vector product can be done for a subset of elements
    \param[in] var given problem
    \param[in] B vector to be multiplied
    \param[in] beta coefficient
    \param[in,out] C result vector
    alpha*A*B + beta*C -> C
  */
  template<class TypeEquation>
  void MltAdd_SquareHex_VectorialHcurl3D(const Real_wp& alpha, int level,
                                         HyperbolicProblem<TypeEquation>& var,
                                         const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C)
  {
    Real_wp zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (beta == zero)
      C.Zero();
    else if (beta != one)
      Mlt(beta, C);

    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var.var_harmonic;
    const Mesh<Dimension3>& mesh = vars.mesh;
    const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);

    VectReal_wp Bstar;
    int Nvol = mesh_num.GetNbDof();
    if (vars.GetNbGlobalEltPML() > 0)
      {
	Bstar.Reallocate(Nvol);
	for (int i = 0; i < Nvol; i++)
	  {
	    Bstar(i) = B(i);
	    int npml = mesh_num.GetDofPML(i);
	    if (npml >= 0)
	      Bstar(i) = B(Nvol + npml);
	  }
      }
    
    TinyVector<VectReal_wp, 1> Uh;
    VectReal_wp Vh;
    TinyVector<Real_wp, 3> h_tmp, vh_loc;    
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);

	int nb_dof_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();

	bool affine = mesh.IsElementAffine(i);
	bool pml_elt = vars.InsidePML(i);

	// we check that the element does not touch a PML layer
	bool close_pml_elt = false;
	if ( (!pml_elt) && (vars.GetNbGlobalEltPML() > 0))
	  for (int j = 0; j < nb_dof_elt; j++)
	    {
	      int num_dof = mesh_num.Element(i).GetNumberDof(j);
	      if (num_dof >= 0)
		{
		  int npml = mesh_num.GetDofPML(num_dof);
		  if (npml >= 0)
		    close_pml_elt = true;
		}
	    }
	
	// we retrieve values of E* on the element
	Uh(0).Reallocate(nb_dof_elt);

	if (close_pml_elt || pml_elt)
	  vars.GetLocalUnknownVector(Bstar, i, Uh);
	else
	  vars.GetLocalUnknownVector(B, i, Uh);

	// curl of E* on the reference element
	Vh.Reallocate(3*nb_points_quad);
	Fb.ApplyRhTranspose(Uh(0), Vh);
	
	// then application of jacobian matrix
	int offset = vars.GetOffsetDofV(i);
	for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	  {
	    CopyVector(Vh, j, h_tmp);

	    if (affine)
	      Mlt(vars.Glob_DFj(i)(0), h_tmp, vh_loc);
	    else
	      Mlt(vars.Glob_DFj(i)(j), h_tmp, vh_loc);
		    
	    Mlt(alpha*Fb.WeightsND(j), vh_loc); 
	    C(offset + 3*j) -= vh_loc(0);
	    C(offset + 3*j+1) -= vh_loc(1);
	    C(offset + 3*j+2) -= vh_loc(2);
	  }
      }
  }
  
  
  //! matrix vector product by stiffness matrix for Nedelec's second family on hexahedra
  /*!
    \param[in] alpha coefficient
    \param[in] level matrix-vector product can be done for a subset of elements
    \param[in] var given problem
    \param[in] B vector to be multiplied
    \param[in] beta coefficient
    \param[in,out] C result vector
    alpha*A*B + beta*C -> C
  */
  template <class TypeEquation>
  void MltAdd_SquareHex_ScalarHcurl3D(const Real_wp& alpha, int level,
                                      HyperbolicProblem<TypeEquation>& var,
                                      const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C)
  {
    typedef Real_wp Complexe;

    Complexe zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    
    if (beta == zero)
      C.Zero();
    else if (beta != one)
      Mlt(beta, C);
    
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var.var_harmonic;
    const Mesh<Dimension3>& mesh = vars.mesh;
    //const MeshNumbering<Dimension3>& mesh_num = vars.mesh_num;
    
    TinyVector<VectReal_wp, 1> Uh;
    VectReal_wp Vh;
    R3 h_tmp, vh_loc;
    
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);

	int nb_dof_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();

	bool affine = mesh.IsElementAffine(i);
	bool pml_elt = vars.InsidePML(i);

	int offset = vars.GetOffsetDofV(i);
	// for PML elements, the product is performed with H* instead of H
	if (pml_elt)
	  offset += 3*nb_points_quad;

	Vh.Reallocate(3*nb_points_quad);
	Uh(0).Reallocate(nb_dof_elt);
	for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	  {
	    vh_loc(0) = B(offset + 3*j);
	    vh_loc(1) = B(offset + 3*j+1);
	    vh_loc(2) = B(offset + 3*j+2);
	    
	    if (affine)
	      MltTrans(vars.Glob_DFj(i)(0), vh_loc, h_tmp);
	    else
	      MltTrans(vars.Glob_DFj(i)(j), vh_loc, h_tmp);
	    
	    Mlt(Fb.WeightsND(j), h_tmp); 
	    Vh(3*j) = h_tmp(0); Vh(3*j+1) = h_tmp(1); Vh(3*j+2) = h_tmp(2);
	  }
	
	Fb.ApplyRh(Vh, Uh(0));
	
	vars.AddLocalUnknownVector(alpha, Uh, i, C);
      }
  }  
#endif
  
  
  /**************
   * DG Methods *
   **************/
  
  
  template<class T>
  class MaxwellVariablesDG
  {
  public :
    static Vector<T> Ex, Ey, Ez, Ex_quad, Ey_quad, Ez_quad, dEx, dEy, dEz;
    static Vector<T> Hx, Hy, Hz, HxLoc, HyLoc, HzLoc, Ephix, Ephiy, Ephiz, Vx, Vy, Vz;
    static TinyVector<Vector<T>, 3> curlE;
  };

  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ex;
  
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ey;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ez;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ex_quad;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ey_quad;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ez_quad;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::dEx;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::dEy;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::dEz;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Hx;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Hy;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Hz;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::HxLoc;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::HyLoc;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::HzLoc;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ephix;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ephiy;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Ephiz;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Vx;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Vy;
  template<class T>
  Vector<T> MaxwellVariablesDG<T>::Vz;
  
  template<class T>
  TinyVector<Vector<T>, 3> MaxwellVariablesDG<T>::curlE;
  
  template<class Complexe, class TypeEquation>
  void 
  MltAdd_ElementRh1Maxwell(const EllipticProblem<TypeEquation>& var,
                           int num_elem, const Vector<Complexe>& B,
			   TinyVector<Vector<Complexe>, 3>& En_quad,
                           const Complexe& alpha, Vector<Complexe>& C,
			   const ElementReference<Dimension3, 1>& Fb)
  {
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    bool affine = mesh.IsElementAffine(num_elem);
    
    int nb_points_elt = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    if ((affine)&&(!Fb.UseQuadraturePointsForSh()))
      nb_points_quad = nb_points_elt;

    // getting values of B inside the element
    Vector<Complexe>& Ex = MaxwellVariablesDG<Real_wp>::Ex;
    Vector<Complexe>& Ey = MaxwellVariablesDG<Real_wp>::Ey;
    Vector<Complexe>& Ez = MaxwellVariablesDG<Real_wp>::Ez;
    Ex.Reallocate(nb_points_elt);
    Ey.Reallocate(nb_points_elt);
    Ez.Reallocate(nb_points_elt);
    Vector<Complexe>& Ex_quad = MaxwellVariablesDG<Real_wp>::Ex_quad;
    Vector<Complexe>& Ey_quad = MaxwellVariablesDG<Real_wp>::Ey_quad;
    Vector<Complexe>& Ez_quad = MaxwellVariablesDG<Real_wp>::Ez_quad;
    int size_du = nb_points_quad*3;
    Vector<Complexe>& dEx = MaxwellVariablesDG<Real_wp>::dEx;
    Vector<Complexe>& dEy = MaxwellVariablesDG<Real_wp>::dEy;
    Vector<Complexe>& dEz = MaxwellVariablesDG<Real_wp>::dEz;
    dEx.Reallocate(size_du);
    dEy.Reallocate(size_du);
    dEz.Reallocate(size_du);
    
    TinyVector<Vector<Complexe>, 3>& curlE = MaxwellVariablesDG<Real_wp>::curlE;
    curlE(0).Reallocate(nb_points_quad);
    curlE(1).Reallocate(nb_points_quad);
    curlE(2).Reallocate(nb_points_quad);
    
    int offset_x = mesh_num.Element(num_elem).GetNumberDof(0);
    int offset_y = offset_x + mesh_num.GetNbDof();
    int offset_z = offset_y + mesh_num.GetNbDof();
    if (Fb.UseQuadraturePointsForSh())
      {
	Fb.ApplyRhQuadratureTranspose(En_quad(0), dEx);
	Fb.ApplyRhQuadratureTranspose(En_quad(1), dEy);
	Fb.ApplyRhQuadratureTranspose(En_quad(2), dEz);
	Mlt(alpha, dEx); Mlt(alpha, dEy); Mlt(alpha, dEz);
      }
    else
      {
	for (int j = 0; j < nb_points_elt; j++)
	  {
	    Ex(j) = alpha*B(offset_x + j);
	    Ey(j) = alpha*B(offset_y + j);
	    Ez(j) = alpha*B(offset_z + j);
	  }
	
	// computation of grad(E)
	if (affine)
	  {
	    Fb.ApplyConstantRhTranspose(Ex, dEx);
	    Fb.ApplyConstantRhTranspose(Ey, dEy);
	    Fb.ApplyConstantRhTranspose(Ez, dEz);	
	  }
	else
	  {
	    Fb.ApplyRhTranspose(Ex, dEx);
	    Fb.ApplyRhTranspose(Ey, dEy);
	    Fb.ApplyRhTranspose(Ez, dEz);
            
            if (ElementReference_Base::use_warburton_trick)
              {
                Ex_quad.Reallocate(nb_points_quad);
                Ey_quad.Reallocate(nb_points_quad);
                Ez_quad.Reallocate(nb_points_quad);
                Fb.ApplyChTranspose(Ex, Ex_quad);
                Fb.ApplyChTranspose(Ey, Ey_quad);
                Fb.ApplyChTranspose(Ez, Ez_quad);
              }
	  }
	
      }
    
    TinyVector<Complexe, 3> vec_u, vecE, gradEx, gradEy, gradEz;
    Real_wp coef, invJacob;
    TinyMatrix<Real_wp, General, 3, 3> dfjm1;
    // application of transformation DF_i^{*-1}
    if (affine)
      {
        dfjm1 = var.Glob_DFjm1(num_elem)(0);
	
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    vec_u(0) = dEx(3*k); vec_u(1) = dEx(3*k+1); vec_u(2) = dEx(3*k+2);	    
	    MltTrans(dfjm1, vec_u, gradEx);
	    vec_u(0) = dEy(3*k); vec_u(1) = dEy(3*k+1); vec_u(2) = dEy(3*k+2);	    
	    MltTrans(dfjm1, vec_u, gradEy);
	    vec_u(0) = dEz(3*k); vec_u(1) = dEz(3*k+1); vec_u(2) = dEz(3*k+2);	    
	    MltTrans(dfjm1, vec_u, gradEz);
	    
	    curlE(0)(k) = -gradEz(1) + gradEy(2);
	    curlE(1)(k) = -gradEx(2) + gradEz(0);
	    curlE(2)(k) = -gradEy(0) + gradEx(1);
	  }
        
        if (Fb.UseQuadraturePointsForSh())
          for (int k = 0; k < nb_points_quad; k++)
            {
              coef = Fb.WeightsND(k);
              curlE(0)(k) *= coef;
              curlE(1)(k) *= coef;
              curlE(2)(k) *= coef;
            }
      }
    else
      for (int k = 0; k < nb_points_quad; k++)
	{
          dfjm1 = var.Glob_DFjm1(num_elem)(k);
          
	  // grad Ex
	  vec_u(0) = dEx(3*k); vec_u(1) = dEx(3*k+1); vec_u(2) = dEx(3*k+2);
	  
	  MltTrans(dfjm1, vec_u, gradEx);
	  
	  coef = Fb.WeightsND(k); Mlt(coef, gradEx);
	  
	  // grad Ey
	  vec_u(0) = dEy(3*k); vec_u(1) = dEy(3*k+1); vec_u(2) = dEy(3*k+2);
	  
	  MltTrans(dfjm1, vec_u, gradEy);
	  Mlt(coef, gradEy);
	  
	  // grad Ez
	  vec_u(0) = dEz(3*k); vec_u(1) = dEz(3*k+1); vec_u(2) = dEz(3*k+2);
	  
	  MltTrans(dfjm1, vec_u, gradEz);
	  Mlt(coef, gradEz);
	 
	  curlE(0)(k) = -gradEz(1) + gradEy(2);
	  curlE(1)(k) = -gradEx(2) + gradEz(0);
	  curlE(2)(k) = -gradEy(0) + gradEx(1);
          
          if (ElementReference_Base::use_warburton_trick)
            {
              invJacob = coef/var.Glob_jacobian(num_elem)(k);
              curlE(0)(k) *= invJacob;
              curlE(1)(k) *= invJacob;
              curlE(2)(k) *= invJacob;
              
              // we add 0.5/J_i grad(J_i) \times E \cdot phi
              if (Fb.UseQuadraturePointsForSh())
                {
                  vecE(0) = En_quad(0)(k);
                  vecE(1) = En_quad(1)(k);
                  vecE(2) = En_quad(2)(k);
                }
              else
                {
                  vecE(0) = Ex_quad(k);
                  vecE(1) = Ey_quad(k);
                  vecE(2) = Ez_quad(k);
                }
              
              TimesProd(mesh.Glob_GradJacobian(num_elem)(k), vecE, vec_u);
              curlE(0)(k) += coef*vec_u(0);
              curlE(1)(k) += coef*vec_u(1);
              curlE(2)(k) += coef*vec_u(2);
            }
	}
    
    // application of Ch    
    if (Fb.UseQuadraturePointsForSh())
      {
	En_quad(0) = curlE(0);
	En_quad(1) = curlE(1);
	En_quad(2) = curlE(2);
      }
    else 
      {
	if (affine)
	  {
	    Ex = curlE(0); Ey = curlE(1); Ez = curlE(2);
	  }
	else
	  {
	    Fb.ApplyCh(curlE(0), Ex);
	    Fb.ApplyCh(curlE(1), Ey);
	    Fb.ApplyCh(curlE(2), Ez);
	  }
    
	for (int k = 0; k < nb_points_elt; k++)
	  {
	    C(offset_x + k) += Ex(k);
	    C(offset_y + k) += Ey(k);
	    C(offset_z + k) += Ez(k);
	  }
      }
  }

  
  template<class Complexe, class TypeEquation>
  void
  MltAdd_ElementRh2Maxwell(const EllipticProblem<TypeEquation>& var,
                           int num_elem, const Vector<Complexe>& B,
			   TinyVector<Vector<Complexe>, 3>& Hn_quad,
                           const Complexe& alpha, Vector<Complexe>& C,
			   const ElementReference<Dimension3, 1>& Fb)
  {
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    bool affine = mesh.IsElementAffine(num_elem);

    // number of dofs inside the element
    int nb_points_elt = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    if ((affine)&&(!Fb.UseQuadraturePointsForSh()))
      nb_points_quad = nb_points_elt;
    
    // getting values of H
    Vector<Complexe>& HxLoc = MaxwellVariablesDG<Complexe>::HxLoc;
    Vector<Complexe>& HyLoc = MaxwellVariablesDG<Complexe>::HyLoc;
    Vector<Complexe>& HzLoc = MaxwellVariablesDG<Complexe>::HzLoc;
    Vector<Complexe>& Hx = MaxwellVariablesDG<Complexe>::Hx;
    Vector<Complexe>& Hy = MaxwellVariablesDG<Complexe>::Hy;
    Vector<Complexe>& Hz = MaxwellVariablesDG<Complexe>::Hz;
    Vector<Complexe>& Ephix = MaxwellVariablesDG<Complexe>::Ephix;
    Vector<Complexe>& Ephiy = MaxwellVariablesDG<Complexe>::Ephiy;
    Vector<Complexe>& Ephiz = MaxwellVariablesDG<Complexe>::Ephiz;
    HxLoc.Reallocate(nb_points_elt);
    HyLoc.Reallocate(nb_points_elt);
    HzLoc.Reallocate(nb_points_elt);
    Hx.Reallocate(nb_points_quad);
    Hy.Reallocate(nb_points_quad);
    Hz.Reallocate(nb_points_quad);
    
    Vector<Complexe>& Vx = MaxwellVariablesDG<Complexe>::Vx;
    Vector<Complexe>& Vy = MaxwellVariablesDG<Complexe>::Vy;
    Vector<Complexe>& Vz = MaxwellVariablesDG<Complexe>::Vz;
    Vx.Reallocate(3*nb_points_quad);
    Vy.Reallocate(3*nb_points_quad);
    Vz.Reallocate(3*nb_points_quad);
        
    if ((!affine) && (ElementReference_Base::use_warburton_trick))
      {
        Ephix.Reallocate(nb_points_quad);
        Ephiy.Reallocate(nb_points_quad);
        Ephiz.Reallocate(nb_points_quad);
      }
    
    int Nvol = mesh_num.GetNbDof();
    int offset_x = mesh_num.Element(num_elem).GetNumberDof(0);
    int offset_y = offset_x + Nvol;
    int offset_z = offset_y + Nvol;
    Real_wp coef, invJacob;
    if (Fb.UseQuadraturePointsForSh())
      {
	Hx = Hn_quad(0);
	Hy = Hn_quad(1);
	Hz = Hn_quad(2);
	Mlt(alpha, Hx); Mlt(alpha, Hy); Mlt(alpha, Hz);
        
        if (affine)
          for (int k = 0; k < nb_points_quad; k++)
            {
              coef = Fb.WeightsND(k);
              Hx(k) *= coef;
              Hy(k) *= coef;
              Hz(k) *= coef;
            }
      }
    else
      {
	for (int j = 0; j < nb_points_elt; j++)
	  {
	    HxLoc(j) = alpha*B(offset_x + j);
	    HyLoc(j) = alpha*B(offset_y + j);
	    HzLoc(j) = alpha*B(offset_z + j);
	  }
	
	if (affine)
	  {
	    Hx = HxLoc; Hy = HyLoc; Hz = HzLoc;
	  }
	else
	  {
	    Fb.ApplyChTranspose(HxLoc, Hx);
	    Fb.ApplyChTranspose(HyLoc, Hy);
	    Fb.ApplyChTranspose(HzLoc, Hz);
	  }
      }
    
    TinyVector<Complexe, 3> vecH, vecE;
    // application of transformation DF_i^{-1}
    TinyMatrix<Real_wp, General, 3, 3> dfjm1;
    if (affine)
      {
        dfjm1 = var.Glob_DFjm1(num_elem)(0);
	
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    vecH(0) = 0; vecH(1) = -Hz(k); vecH(2) = Hy(k);
	    Mlt(dfjm1, vecH, vecE);
	    
	    Vx(3*k) = vecE(0); Vx(3*k+1) = vecE(1); Vx(3*k+2) = vecE(2);
	    
	    vecH(0) = Hz(k); vecH(1) = 0; vecH(2) = -Hx(k);
	    Mlt(dfjm1, vecH, vecE);
	    
	    Vy(3*k) = vecE(0); Vy(3*k+1) = vecE(1); Vy(3*k+2) = vecE(2);
	    
	    vecH(0) = -Hy(k); vecH(1) = Hx(k); vecH(2) = 0;
	    Mlt(dfjm1, vecH, vecE);
	    
	    Vz(3*k) = vecE(0); Vz(3*k+1) = vecE(1); Vz(3*k+2) = vecE(2);
	  }
      }
    else
      for (int k = 0; k < nb_points_quad; k++)
	{
          dfjm1 = var.Glob_DFjm1(num_elem)(k);
	  vecH(0) = 0; vecH(1) = -Hz(k); vecH(2) = Hy(k);
	  
	  Mlt(dfjm1, vecH, vecE);
	  
	  coef = Fb.WeightsND(k); Mlt(coef, vecE);
	  Vx(3*k) = vecE(0); Vx(3*k+1) = vecE(1); Vx(3*k+2) = vecE(2);

	  vecH(0) = Hz(k); vecH(1) = 0; vecH(2) = -Hx(k);
	  
	  Mlt(dfjm1, vecH, vecE);
	  Mlt(coef, vecE);
	  Vy(3*k) = vecE(0); Vy(3*k+1) = vecE(1); Vy(3*k+2) = vecE(2);

	  vecH(0) = -Hy(k); vecH(1) = Hx(k); vecH(2) = 0;
	  
	  Mlt(dfjm1, vecH, vecE);
	  Mlt(coef, vecE);
	  Vz(3*k) = vecE(0); Vz(3*k+1) = vecE(1); Vz(3*k+2) = vecE(2);

          if (ElementReference_Base::use_warburton_trick)
            {
              invJacob = coef/var.Glob_jacobian(num_elem)(k);
              Vx(3*k) *= invJacob; Vx(3*k+1) *= invJacob; Vx(3*k+2) *= invJacob;
              Vy(3*k) *= invJacob; Vy(3*k+1) *= invJacob; Vy(3*k+2) *= invJacob;
              Vz(3*k) *= invJacob; Vz(3*k+1) *= invJacob; Vz(3*k+2) *= invJacob;
              
              // we add 0.5/J_i grad(J_i) \times H \cdot phi
              vecH(0) = Hx(k);
              vecH(1) = Hy(k);
              vecH(2) = Hz(k);
              
              TimesProd(mesh.Glob_GradJacobian(num_elem)(k), vecH, vecE);
              Ephix(k) = coef*vecE(0);
              Ephiy(k) = coef*vecE(1);
              Ephiz(k) = coef*vecE(2);
            }
	}
    
    // integration against grad phi
    if (Fb.UseQuadraturePointsForSh())
      {
	Fb.ApplyRhQuadrature(Vx, Hx);
	Fb.ApplyRhQuadrature(Vy, Hy);
	Fb.ApplyRhQuadrature(Vz, Hz);
	Hn_quad(0) = Hx; Hn_quad(1) = Hy; Hn_quad(2) = Hz;
        if ((!affine) && (ElementReference_Base::use_warburton_trick))
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                Hn_quad(0)(k) += Ephix(k);
                Hn_quad(1)(k) += Ephiy(k);
                Hn_quad(2)(k) += Ephiz(k);
              }
          }
      }
    else
      {
	if (affine)
	  {
	    Fb.ApplyConstantRh(Vx, HxLoc);
	    Fb.ApplyConstantRh(Vy, HyLoc);
	    Fb.ApplyConstantRh(Vz, HzLoc);
	  }
	else
	  {
 	    Fb.ApplyRh(Vx, HxLoc);
	    Fb.ApplyRh(Vy, HyLoc);
	    Fb.ApplyRh(Vz, HzLoc);

            if (ElementReference_Base::use_warburton_trick)
              {
                Hx.Reallocate(nb_points_elt);
                Hy.Reallocate(nb_points_elt);
                Hz.Reallocate(nb_points_elt);
                Fb.ApplyCh(Ephix, Hx);
                Fb.ApplyCh(Ephiy, Hy);
                Fb.ApplyCh(Ephiz, Hz);
                for (int k = 0; k < nb_points_elt; k++)
                  {
                    HxLoc(k) += Hx(k);
                    HyLoc(k) += Hy(k);
                    HzLoc(k) += Hz(k);
                  }
              }
	  }
	
	for (int k = 0; k < nb_points_elt; k++)
	  {
	    C(offset_x + k) += HxLoc(k);
	    C(offset_y + k) += HyLoc(k);
	    C(offset_z + k) += HzLoc(k);
	  }
      }
  }
  
  
  template<class Complexe, class TypeEquation>
  void
  MltAdd_ElementShMaxwell(const EllipticProblem<TypeEquation>& var,
                          int num_elem1, const Vector<Complexe>& B,
			  TinyVector<Vector<Complexe>, 3>& En_quad,
                          const TinyVector<Vector<Complexe>, 3>& extrapolE, 
			  const TinyVector<Vector<Complexe>, 3>& Eneighbor,
                          const TinyVector<Vector<Complexe>, 3>& extrapolH,
			  const TinyVector<Vector<Complexe>, 3>& Hneighbor, 
                          const Complexe& alpha, Vector<Complexe>& C,
			  const ElementReference<Dimension3, 1>& Fb)
  {
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    
    int nb_points_elt = Fb.GetNbDof();
    int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();
    
    TinyVector<Vector<Complexe>, 3> Vh;
        
    int offset_face1 = mesh_num.OffsetQuadElementNumber(num_elem1);
    bool affine = mesh.IsElementAffine(num_elem1);
    int Nvol = mesh_num.GetNbDof();
    int ref_domain = mesh.Element(num_elem1).GetReference();
    TinyVector<Complexe, 3> valE1, valE2, valH1, valH2, jumpE, jumpH, vec_v, fH;
    R3 normale; Real_wp dsj, coef;
    VectReal_wp Cloc;
    bool presence_penalization = (var.delta_penalization != Real_wp(0));
    
    for (int num_pos1_face = 0; num_pos1_face < Fb.GetNbBoundaries(); num_pos1_face++)
      {
	int num_face = mesh.Element(num_elem1).numBoundary(num_pos1_face);
	int rot1 = mesh.Element(num_elem1).GetOrientationBoundary(num_pos1_face);
        int ref = mesh.Boundary(num_face).GetReference();
        int rf = mesh_num.GetOrderQuadrature(num_face);
	int num_elem2 = mesh.Boundary(num_face).numElement(0);
	bool new_face = var.IsNewFace(num_elem1)(num_pos1_face);
	int cond = mesh.GetBoundaryCondition(ref);
	    
	int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        const VectReal_wp& PoidsFlux
          = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
        const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	  GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
        
        Vh(0).Reallocate(nb_points_face);
        Vh(1).Reallocate(nb_points_face);
        Vh(2).Reallocate(nb_points_face);
        
        int num_pos2_face = -1, krot = 0;
	int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0, rot2;
	if (mesh.Boundary(num_face).GetNbElements() == 2)
	  if (num_elem2 == num_elem1)
	    num_elem2 = mesh.Boundary(num_face).numElement(1);
	        
	int offset_face2 = mesh_num.OffsetQuadElementNumber(num_elem2);
	if (num_elem2 != num_elem1)
	  {
	    int nb = 0;
	    for (int k = 0; k < mesh.Element(num_elem2).GetNbBoundary(); k++)
	      {
		int nf = mesh.Element(num_elem2).numBoundary(k);
		if (nf == num_face)
		  {
		    num_pos2_face = k;
		    break;
		  }
		else
		  nb += mesh_num.GetNbPointsQuadratureBoundary(nf);
	      }
	    
	    if (num_pos2_face < 0)
	      {
		int nf2 = mesh_num.GetPeriodicBoundary(num_face);
		nb = 0;
		for (int k = 0; k < mesh.Element(num_elem2).GetNbBoundary(); k++)
		  {
		    int nf = mesh.Element(num_elem2).numBoundary(k);
		    if (nf == nf2)
		      {
			num_pos2_face = k;
			break;
		      }
		    else
		      nb += mesh_num.GetNbPointsQuadratureBoundary(nf);
		  }
	      }
	    
	    offset_face2 += nb;
	    
	    rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
            rot = mesh_num.GetRotationFace(rot1, rot2, nv);
          }
	else
	  offset_face2 = offset_face1;
	
        Real_wp imped = var.delta_penalization*var.coefficient_impedance_absorbing(ref_domain);
                    		  
	if (cond != BoundaryConditionEnum::LINE_INSIDE)
	  {
	    // face with a boundary condition
	    for (int k = 0; k < nb_points_face; k++)
	      {
		// internal term
		krot = FacesQuadRotation(rot, k);
                
		// value of E
		int num_point = offset_face1 + k;
		valE1(0) = -extrapolE(0)(num_point);
		valE1(1) = -extrapolE(1)(num_point);
		valE1(2) = -extrapolE(2)(num_point);
					    
		// normale
                if (new_face)
                  {
                    normale = var.Glob_normale(num_face)(k);
                    dsj = var.Glob_dsj(num_face)(k);
                  }
                else
                  {
                    normale = var.Glob_normale(num_face)(krot);
                    dsj = var.Glob_dsj(num_face)(krot);
                    Mlt(Real_wp(-1), normale);
                  }
		
		TimesProd(normale, valE1, fH);
		    
		// term from the boundary condition
		if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
		  {
		    // neumann
		    Mlt(2.0, fH);
		  }
		else if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
		  {
                    int offset_neighbor = var.GetOffsetNeighboringFace(num_face);
		    int num_dof2 = offset_neighbor + k;
		    valE2(0) = Eneighbor(0)(num_dof2);
		    valE2(1) = Eneighbor(1)(num_dof2);
		    valE2(2) = Eneighbor(2)(num_dof2);
		    
		    TimesProd(normale, valE2, vec_v);                    
		    fH += vec_v;
		    
		    if (presence_penalization)
		      {
                        // jump of H
			jumpH(0) = -extrapolH(0)(num_point) + Hneighbor(0)(num_dof2);
			jumpH(1) = -extrapolH(1)(num_point) + Hneighbor(1)(num_dof2);
			jumpH(2) = -extrapolH(2)(num_point) + Hneighbor(2)(num_dof2);
			Real_wp H_dot_n = DotProd(jumpH, normale);
                        
			fH(0) += imped*(jumpH(0) - normale(0)*H_dot_n);
			fH(1) += imped*(jumpH(1) - normale(1)*H_dot_n);
			fH(2) += imped*(jumpH(2) - normale(2)*H_dot_n);
		      }
		    
		  }
                else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
                  {
                    Real_wp H_dot_n = extrapolH(0)(num_point)*normale(0)
                      + extrapolH(1)(num_point)*normale(1) + extrapolH(2)(num_point)*normale(2);
                    Real_wp imped_abc = var.coefficient_impedance_absorbing(ref_domain);
                    
                    fH(0) += imped_abc*(extrapolH(0)(num_point) - normale(0)*H_dot_n);
                    fH(1) += imped_abc*(extrapolH(1)(num_point) - normale(1)*H_dot_n);
                    fH(2) += imped_abc*(extrapolH(2)(num_point) - normale(2)*H_dot_n);
                  }
		else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
                  fH.Zero();
                
                if (Fb.UseQuadratureFreeSh())
                  {
                    Vh(0)(k) = 0.5*dsj*fH(0);
                    Vh(1)(k) = 0.5*dsj*fH(1);
                    Vh(2)(k) = 0.5*dsj*fH(2);
                  }
                else
                  {
                    coef = dsj*PoidsFlux(k);
                    Vh(0)(k) = coef*fH(0);
                    Vh(1)(k) = coef*fH(1);
                    Vh(2)(k) = coef*fH(2);
                  }  
              }	    
	  }
	else if (mesh.Boundary(num_face).GetNbElements()==2)
	  {
	    // internal edge	    
	    for (int k = 0; k < nb_points_face; k++)
	      {
                krot = FacesQuadRotation(rot, k);
		int num_dof1 = offset_face1 + k;
		int num_dof2 = offset_face2 + krot;
				
		// jump -> (E2 - E1)
		jumpE(0) = -extrapolE(0)(num_dof1) + extrapolE(0)(num_dof2);
		jumpE(1) = -extrapolE(1)(num_dof1) + extrapolE(1)(num_dof2);
		jumpE(2) = -extrapolE(2)(num_dof1) + extrapolE(2)(num_dof2);
                
                if (new_face)
                  {
                    normale = var.Glob_normale(num_face)(k);
                    dsj = var.Glob_dsj(num_face)(k);
                  }
                else
                  {
                    normale = var.Glob_normale(num_face)(krot);
                    dsj = var.Glob_dsj(num_face)(krot);
                    Mlt(Real_wp(-1), normale);
                  }
                
		TimesProd(normale, jumpE, fH);
		                
                if (presence_penalization)
                  {
                    // jump of H
                    jumpH(0) = -extrapolH(0)(num_dof1) + extrapolH(0)(num_dof2);
                    jumpH(1) = -extrapolH(1)(num_dof1) + extrapolH(1)(num_dof2);
                    jumpH(2) = -extrapolH(2)(num_dof1) + extrapolH(2)(num_dof2);
                    Real_wp H_dot_n = DotProd(jumpH, normale);
                    
                    fH(0) += imped*(jumpH(0) - normale(0)*H_dot_n);
                    fH(1) += imped*(jumpH(1) - normale(1)*H_dot_n);
                    fH(2) += imped*(jumpH(2) - normale(2)*H_dot_n);
                  }
                
                if (Fb.UseQuadratureFreeSh())
                  {
                    Vh(0)(k) = 0.5*dsj*fH(0);
                    Vh(1)(k) = 0.5*dsj*fH(1);
                    Vh(2)(k) = 0.5*dsj*fH(2);
                  }
                else
                  {
                    coef = dsj*PoidsFlux(k);
                    Vh(0)(k) = coef*fH(0);
                    Vh(1)(k) = coef*fH(1);
                    Vh(2)(k) = coef*fH(2);
                  }
	      }
	  }
        
        if ((!affine)&&(ElementReference_Base::use_warburton_trick))
          {
            for (int k = 0; k < nb_points_face; k++)
              {
                Real_wp jacob = mesh.Glob_invSqrtJacobian(num_elem1)(offset_war + k);
                Vh(0)(k) *= jacob;
                Vh(1)(k) *= jacob;
                Vh(2)(k) *= jacob;
              }
            
            offset_war += nb_points_face;
          }
        
        int offset_x = mesh_num.Element(num_elem1).GetNumberDof(0);
        int offset_y = offset_x + Nvol;
        int offset_z = offset_y + Nvol;
        if (Fb.UseQuadraturePointsForSh())
          {
            // we add contribution to Un_quad
            Fb.ApplyShQuadrature(-alpha, num_pos1_face, Vh(0), En_quad(0), rf);
            Fb.ApplyShQuadrature(-alpha, num_pos1_face, Vh(1), En_quad(1), rf);
            Fb.ApplyShQuadrature(-alpha, num_pos1_face, Vh(2), En_quad(2), rf);	
          }
        else
          {
            Cloc.SetData(nb_points_elt, &C(offset_x)); 
            Fb.ApplySh(-alpha, num_pos1_face, Vh(0), Cloc, rf);
            Cloc.Nullify();
            
            Cloc.SetData(nb_points_elt, &C(offset_y)); 
            Fb.ApplySh(-alpha, num_pos1_face, Vh(1), Cloc, rf);
            Cloc.Nullify();
            
            Cloc.SetData(nb_points_elt, &C(offset_z)); 
            Fb.ApplySh(-alpha, num_pos1_face, Vh(2), Cloc, rf);
            Cloc.Nullify();
          }
        
        offset_face1 += nb_points_face;
      }
    
  }
				  
  
  template<class Complexe, class TypeEquation>
  void MltAdd_ElementShtMaxwell(const EllipticProblem<TypeEquation>& var,
                                int num_elem1, const Vector<Complexe>& B,
				TinyVector<Vector<Complexe>, 3>& Hn_quad,
				const TinyVector<Vector<Complexe>, 3>& extrapolE,
				const TinyVector<Vector<Complexe>, 3>& Eneighbor,
				const TinyVector<Vector<Complexe>, 3>& extrapolH,
				const TinyVector<Vector<Complexe>, 3>& Hneighbor, 
				const Complexe& alpha, Vector<Complexe>& C,
				const ElementReference<Dimension3, 1>& Fb)
  {
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
 
    int ref_domain = mesh.Element(num_elem1).GetReference();
    int nb_points_elt = Fb.GetNbDof();
    int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();
    
    TinyVector<Vector<Complexe>, 3> Vh;
    
    bool affine = mesh.IsElementAffine(num_elem1);
    int offset_face1 = mesh_num.OffsetQuadElementNumber(num_elem1);
    TinyVector<Complexe, 3> valH1, valH2, valE1, valE2, fE, vec_v, sumH, jumpE;
    R3 normale; Real_wp dsj, coef;
    bool presence_penalization = (var.alpha_penalization != Real_wp(0));
    
    for (int num_pos1_face = 0; num_pos1_face < Fb.GetNbBoundaries(); num_pos1_face++)
      {
	int num_face = mesh.Element(num_elem1).numBoundary(num_pos1_face);
	int rot1 = mesh.Element(num_elem1).GetOrientationBoundary(num_pos1_face);
	int rf = mesh_num.GetOrderQuadrature(num_face);
        int ref = mesh.Boundary(num_face).GetReference();
	int num_elem2 = mesh.Boundary(num_face).numElement(0);
	bool new_face = var.IsNewFace(num_elem1)(num_pos1_face);
	int cond = mesh.GetBoundaryCondition(ref);
	
        int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        const VectReal_wp& PoidsFlux = mesh_num.number_map
	  .GetFluxWeight(rf, mesh.Boundary(num_face));
        
	const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	  GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
        
        Vh(0).Reallocate(nb_points_face);
        Vh(1).Reallocate(nb_points_face);
        Vh(2).Reallocate(nb_points_face);
        
        int num_pos2_face = -1, krot = 0;
	int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0, rot2;
	if (mesh.Boundary(num_face).GetNbElements() == 2)
	  if (num_elem2 == num_elem1)
	    num_elem2 = mesh.Boundary(num_face).numElement(1);
	        
	int offset_face2 = mesh_num.OffsetQuadElementNumber(num_elem2);
	if (num_elem2 != num_elem1)
	  {
	    int nb = 0;
	    for (int k = 0; k < mesh.Element(num_elem2).GetNbBoundary(); k++)
	      {
		int nf = mesh.Element(num_elem2).numBoundary(k);
		if (nf == num_face)
		  {
		    num_pos2_face = k;
		    break;
		  }
		else
		  nb += mesh_num.GetNbPointsQuadratureBoundary(nf);
	      }
	    
	    if (num_pos2_face < 0)
	      {
		int nf2 = mesh_num.GetPeriodicBoundary(num_face);
		nb = 0;
		for (int k = 0; k < mesh.Element(num_elem2).GetNbBoundary(); k++)
		  {
		    int nf = mesh.Element(num_elem2).numBoundary(k);
		    if (nf == nf2)
		      {
			num_pos2_face = k;
			break;
		      }
		    else
		      nb += mesh_num.GetNbPointsQuadratureBoundary(nf);
		  }
	      }
	    
	    offset_face2 += nb;
	    
	    rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
            rot = mesh_num.GetRotationFace(rot1, rot2, nv);
          }
	else
	  offset_face2 = offset_face1;
	
        //DISP(offset_face1); DISP(num_elem1); DISP(num_elem2);
        //DISP(nb_points_face); DISP(rot1); DISP(rot2); DISP(rot);
        Real_wp imped = var.alpha_penalization/var.coefficient_impedance_absorbing(ref_domain);
                            
	if (cond != BoundaryConditionEnum::LINE_INSIDE)
	  {
	    // face with a boundary condition
	    for (int k = 0; k < nb_points_face; k++)
	      {
		// internal term
		krot = FacesQuadRotation(rot, k);
                
		// value of H
		int num_point = offset_face1 + k;
		valH1(0) = extrapolH(0)(num_point);
		valH1(1) = extrapolH(1)(num_point);
		valH1(2) = extrapolH(2)(num_point);
			  
		// normale
                if (new_face)
                  {
                    normale = var.Glob_normale(num_face)(k);
                    dsj = var.Glob_dsj(num_face)(k);
                  }
                else
                  {
                    normale = var.Glob_normale(num_face)(krot);
                    dsj = var.Glob_dsj(num_face)(krot);
                    Mlt(Real_wp(-1), normale);
                  }
		
		TimesProd(normale, valH1, fE);
		
		// term from the boundary condition
		if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
		  {
		    Mlt(2.0, fE);
		  }
		else if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
		  {
                    int offset_neighbor = var.GetOffsetNeighboringFace(num_face);
		    int num_dof2 = offset_neighbor + k;
		    valH2(0) = Hneighbor(0)(num_dof2);
		    valH2(1) = Hneighbor(1)(num_dof2);
		    valH2(2) = Hneighbor(2)(num_dof2);
		    
		    TimesProd(normale, valH2, vec_v);
		    fE += vec_v;
                    
		    if (presence_penalization)
		      {
			jumpE(0) = -extrapolE(0)(num_point) + Eneighbor(0)(num_dof2);
			jumpE(1) = -extrapolE(1)(num_point) + Eneighbor(1)(num_dof2);
			jumpE(2) = -extrapolE(2)(num_point) + Eneighbor(2)(num_dof2);
			
                        Real_wp E_dot_n = DotProd(jumpE, normale);
						
                        fE(0) -= imped*(jumpE(0) - normale(0)*E_dot_n);
			fE(1) -= imped*(jumpE(1) - normale(1)*E_dot_n);
			fE(2) -= imped*(jumpE(2) - normale(2)*E_dot_n);
		      }		    
		  }
                else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
                  {
                    Real_wp E_dot_n = extrapolE(0)(num_point)*normale(0)
                      + extrapolE(1)(num_point)*normale(1) + extrapolE(2)(num_point)*normale(2);
                    Real_wp imped_abc = 1.0/var.coefficient_impedance_absorbing(ref_domain);
                    
                    fE(0) -= imped_abc*(extrapolE(0)(num_point) - normale(0)*E_dot_n);
                    fE(1) -= imped_abc*(extrapolE(1)(num_point) - normale(1)*E_dot_n);
                    fE(2) -= imped_abc*(extrapolE(2)(num_point) - normale(2)*E_dot_n);
                  }
		else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
		  fE.Zero();
			  
		if (Fb.UseQuadratureFreeSh())
                  {
                    Vh(0)(k) = 0.5*dsj*fE(0);
                    Vh(1)(k) = 0.5*dsj*fE(1);
                    Vh(2)(k) = 0.5*dsj*fE(2);
                  }
                else
                  {
                    coef = dsj*PoidsFlux(k);
                    Vh(0)(k) = coef*fE(0);
                    Vh(1)(k) = coef*fE(1);
                    Vh(2)(k) = coef*fE(2);
                  }
	      }
	  }
	else if (mesh.Boundary(num_face).GetNbElements()==2)
	  {
	    // internal edge
	    for (int k = 0; k < nb_points_face; k++)
	      {
                krot = FacesQuadRotation(rot, k);
		int num_dof1 = offset_face1 + k;
		int num_dof2 = offset_face2 + krot;
		
                if (new_face)
                  {
                    normale = var.Glob_normale(num_face)(k);
                    dsj = var.Glob_dsj(num_face)(k);
                  }
                else
                  {
                    normale = var.Glob_normale(num_face)(krot);
                    dsj = var.Glob_dsj(num_face)(krot);
                    Mlt(Real_wp(-1), normale);
                  }
                
		// mean value (u1+u2)
		sumH(0) = extrapolH(0)(num_dof1) + extrapolH(0)(num_dof2);
		sumH(1) = extrapolH(1)(num_dof1) + extrapolH(1)(num_dof2);
		sumH(2) = extrapolH(2)(num_dof1) + extrapolH(2)(num_dof2);
			  
		TimesProd(normale, sumH, fE);
                
                if (presence_penalization)
                  {
                    jumpE(0) = -extrapolE(0)(num_dof1) + extrapolE(0)(num_dof2);
                    jumpE(1) = -extrapolE(1)(num_dof1) + extrapolE(1)(num_dof2);
                    jumpE(2) = -extrapolE(2)(num_dof1) + extrapolE(2)(num_dof2);
                    
                    Real_wp E_dot_n = DotProd(jumpE, normale);
                    
                    fE(0) -= imped*(jumpE(0) - normale(0)*E_dot_n);
                    fE(1) -= imped*(jumpE(1) - normale(1)*E_dot_n);
                    fE(2) -= imped*(jumpE(2) - normale(2)*E_dot_n);
                  }
                
                if (Fb.UseQuadratureFreeSh())
                  {
                    Vh(0)(k) = 0.5*dsj*fE(0);
                    Vh(1)(k) = 0.5*dsj*fE(1);
                    Vh(2)(k) = 0.5*dsj*fE(2);
                  }
                else
                  {
                    coef = dsj*PoidsFlux(k);
                    Vh(0)(k) = coef*fE(0);
                    Vh(1)(k) = coef*fE(1);
                    Vh(2)(k) = coef*fE(2);
                  }
	      }
	  }
        
        if ((!affine)&&(ElementReference_Base::use_warburton_trick))
          {
            for (int k = 0; k < nb_points_face; k++)
              {
                Real_wp jacob = mesh.Glob_invSqrtJacobian(num_elem1)(offset_war + k);
                Vh(0)(k) *= jacob;
                Vh(1)(k) *= jacob;
                Vh(2)(k) *= jacob;
              }
            offset_war += nb_points_face;
          }
        
        if (Fb.UseQuadraturePointsForSh())
          {
            // we add contribution to Un_quad
            Fb.ApplyShQuadrature(alpha, num_pos1_face, Vh(0), Hn_quad(0), rf);
            Fb.ApplyShQuadrature(alpha, num_pos1_face, Vh(1), Hn_quad(1), rf);
            Fb.ApplyShQuadrature(alpha, num_pos1_face, Vh(2), Hn_quad(2), rf);	
          }
        else
          {
            int pos = mesh_num.Element(num_elem1).GetNumberDof(0);
            VectReal_wp Cloc;
            Cloc.SetData(nb_points_elt, &C(pos)); 
            Fb.ApplySh(alpha, num_pos1_face, Vh(0), Cloc, rf);
            Cloc.Nullify();
            
            pos += mesh_num.GetNbDof();
            Cloc.SetData(nb_points_elt, &C(pos)); 
            Fb.ApplySh(alpha, num_pos1_face, Vh(1), Cloc, rf);
            Cloc.Nullify();
            
            pos += mesh_num.GetNbDof();
            Cloc.SetData(nb_points_elt, &C(pos)); 
            Fb.ApplySh(alpha, num_pos1_face, Vh(2), Cloc, rf);
            Cloc.Nullify();
          }
        
        offset_face1 += nb_points_face;
      }
  }
  
  
  template<class Complexe, class TypeEquation>
  void AssembleVectorElementMaxwellDG(const EllipticProblem<TypeEquation>& vars,
                                      int num_elem1,
                                      const TinyVector<Vector<Complexe>, 3>& Un_quad, Vector<Complexe>& C,
				      const ElementReference<Dimension3, 1>& Fb)
  {
    const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
    if (Fb.UseQuadraturePointsForSh())
      {
        //const Mesh<Dimension>& mesh = vars.mesh;
        int offset_x = mesh_num.Element(num_elem1).GetNumberDof(0);
        int offset_y = offset_x + mesh_num.GetNbDof();
        int offset_z = offset_y + mesh_num.GetNbDof();
        int nb_points_elt = Fb.GetNbDof();
        
        // then we add that to C
        VectReal_wp Cloc;
        Cloc.Reallocate(nb_points_elt); Cloc.Fill(0);
        Fb.ApplyCh(Un_quad(0), Cloc);
	for (int i = 0; i < nb_points_elt; i++)
          C(offset_x + i) += Cloc(i);
        
        Fb.ApplyCh(Un_quad(1), Cloc);
        for (int i = 0; i < nb_points_elt; i++)
          C(offset_y + i) += Cloc(i);
        
        Fb.ApplyCh(Un_quad(2), Cloc);
        for (int i = 0; i < nb_points_elt; i++)
          C(offset_z + i) += Cloc(i);
      }
  }
  
  
  template<class Complexe, class TypeEquation>
  void AssembleVectorMaxwellDG(const EllipticProblem<TypeEquation>& vars,
                               const Vector<TinyVector<Vector<Complexe>, 3> >& Un_quad, Vector<Complexe>& C)
  {
    const Mesh<Dimension3>& mesh = vars.mesh;

    TinyVector<Vector<Real_wp>, 3> Un;
    for (int num_elem = 0; num_elem < mesh.GetNbElt(); num_elem++)
      {
        Un = Un_quad(num_elem);
	AssembleVectorElementMaxwellDG(vars, num_elem, Un, C, vars.GetReferenceElementH1(num_elem));
      }
  }
  
  
  template<class Complexe, class TypeEquationTime>
  void GetExtrapolationMaxwell3D(HyperbolicProblem<TypeEquationTime>& var_time,
                                 const Vector<Complexe>& E, const Vector<Complexe>& H, int level, 
                                 TinyVector<Vector<Complexe>, 3>& extrapolE,
				 TinyVector<Vector<Complexe>, 3>& extrapolH,
                                 TinyVector<Vector<Complexe>, 3>& Eneighbor,
				 TinyVector<Vector<Complexe>, 3>& Hneighbor,
                                 Vector<TinyVector<Vector<Complexe>, 3> >& En_quad,
				 Vector<TinyVector<Vector<Complexe>, 3> >& Hn_quad, bool treatE, bool treatH)
  {
    typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
    const EllipticProblem<TypeEquation>& vars = var_time.var_harmonic;
    
    if (!vars.FormulationDG())
      return;
    
    const Mesh<Dimension3>& mesh = vars.mesh;
    const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
    int size_extrapol = mesh_num.OffsetQuadElementNumber(mesh.GetNbElt());
    int size_neighbor = vars.GetNbPointsQuadratureNeighbor();

    glob_chrono.Start(VirtualTimer::EXTRAPOL);
    
    if (treatE)
      if (En_quad.GetM() != mesh.GetNbElt())
        En_quad.Reallocate(mesh.GetNbElt());

    if (treatH)
      if (Hn_quad.GetM() != mesh.GetNbElt())
        Hn_quad.Reallocate(mesh.GetNbElt());
    
    for (int k = 0; k < 3; k++)
      {
	if (treatE)
          {
            extrapolE(k).Reallocate(size_extrapol);
            if (var_time.GetTimeSchemeType() != TimeSchemeEnum::LOCAL_PIPERNO_SCHEME)
              extrapolE(k).Fill(0);

            Eneighbor(k).Reallocate(size_neighbor);
            if (var_time.GetTimeSchemeType() != TimeSchemeEnum::LOCAL_PIPERNO_SCHEME)
              Eneighbor(k).Fill(0);
          }
        
        if (treatH)
          {
            extrapolH(k).Reallocate(size_extrapol);
            if (var_time.GetTimeSchemeType() != TimeSchemeEnum::LOCAL_PIPERNO_SCHEME)
              extrapolH(k).Fill(0);

            Hneighbor(k).Reallocate(size_neighbor);
            if (var_time.GetTimeSchemeType() != TimeSchemeEnum::LOCAL_PIPERNO_SCHEME)
              Hneighbor(k).Fill(0);
          }
      }

    TinyVector<VectReal_wp, 3> Eloc, Hloc;
    int Nvol = mesh_num.GetNbDof();
    
    // first step : extrapolation step
    TinyVector<Vector<Real_wp>, 3> En, Hn;
    MatrixVectorProductLevel& level_time_scheme = var_time.GetTimeLevelDistribution();
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
          // order of element
          int nb_dof_elt = mesh_num.GetNbLocalDof(num_elem);
          int offset = mesh_num.Element(num_elem).GetNumberDof(0);
          if (treatE)
            {
              for (int m = 0; m < 3; m++)
                {
                  Eloc(m).Reallocate(nb_dof_elt);
                  for (int j = 0; j < nb_dof_elt; j++)
                    Eloc(m)(j) = E(offset+j);
                  
                  offset += Nvol;
                }
              
              offset -= 3*Nvol;
            }
          
          if (treatH)
            {
              for (int m = 0; m < 3; m++)
                {
                  Hloc(m).Reallocate(nb_dof_elt);
                  for (int j = 0; j < nb_dof_elt; j++)
                    Hloc(m)(j) = H(offset+j);
                  
                  offset += Nvol;
                }
            }
          
          // type of element (triangular, quadrangular ...)
	  if (treatE)
	    GetExtrapolationU(vars, num_elem, Eloc, En,
			      extrapolE, vars.GetReferenceElementH1(num_elem),
			      Eneighbor, true);
	  
	  if (treatH)
	    GetExtrapolationU(vars, num_elem, Hloc, Hn,
			      extrapolH, vars.GetReferenceElementH1(num_elem),
			      Hneighbor, true);
	  
          if (treatE)
            En_quad(num_elem) = En;
          
          if (treatH)
            Hn_quad(num_elem) = Hn;
          
        }
    
    glob_chrono.Stop(VirtualTimer::EXTRAPOL);
    
#ifdef SELDON_WITH_MPI
    Vector<Vector<MPI_Request> > requestE(3), requestH(3);
    TinyVector<Vector<Vector<Real_wp> >, 3> EXneighbor, HXneighbor, EXsend, HXsend;
    TinyVector<Vector<Vector<int64_t> >, 3> EXneighbor_tmp, HXneighbor_tmp, EXsend_tmp, HXsend_tmp;
    
    // sending and receiving 
    for (int k = 0; k < 3; k++)
      {
        if (treatE)
          vars.ExchangeUfaceDomains(Eneighbor(k), EXsend(k), EXsend_tmp(k),
                                    EXneighbor(k), EXneighbor_tmp(k), requestE(k), 110+k);
        
        if (treatH)
          vars.ExchangeUfaceDomains(Hneighbor(k), HXsend(k), HXsend_tmp(k),
                                    HXneighbor(k), HXneighbor_tmp(k), requestH(k), 113+k);
      }
    
    // finalizing transfers
    for (int k = 0; k < 3; k++)
      {
        if (treatE)
          vars.GetUfaceDomains(Eneighbor(k), EXsend(k), EXsend_tmp(k),
                               EXneighbor(k), EXneighbor_tmp(k), requestE(k), 110+k);
        
        if (treatH)
          vars.GetUfaceDomains(Hneighbor(k), HXsend(k), HXsend_tmp(k),
                               HXneighbor(k), HXneighbor_tmp(k), requestH(k), 113+k);
      }
    
#endif
    
  }
  
  
  //! matrix vector product with the stiffness matrix for maxwell equation and DG elements
  //! part - \int H rot phi - \int_gamma n \times {H} phi
  /*!
    \param[in] alpha multiplication coefficient
    \param[in] level not used
    \param[in] var instationary problem associated
    \param[in] B vector to be multiplied
    \param[in] beta multiplication coefficient
    \param[out] C resulting vector 
  */
  template<class TypeEquation>
  void MltAddStiffnessScalar_MaxwellDG(const Real_wp& alpha, int level,
                                       HyperbolicProblem<TypeEquation>& var,
                                       const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C)
  {
    if (beta == Real_wp(0))
      C.Fill(0);
    else
      Mlt(beta, C);
    
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars
      = var.var_harmonic;
        
    ////////////////////////////////
    // term -\int_K v \cdot \nabla \varphi
    ///////////////////////////////
    
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    MatrixVectorProductLevel& level_time_scheme = var.GetTimeLevelDistribution();
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	MltAdd_ElementRh2Maxwell(vars, num_elem, B, var.evalHn_quad(num_elem),
				 alpha, C, vars.GetReferenceElementH1(num_elem));
      }
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);    
    
    // fluxes terms
    glob_chrono.Start(VirtualTimer::FLUX);    
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	MltAdd_ElementShtMaxwell(vars, num_elem, B, var.evalHn_quad(num_elem),
				 var.extrapolE, var.Eneighbor,
				 var.extrapolH, var.Hneighbor,
				 alpha, C, vars.GetReferenceElementH1(num_elem));
      }
    
    //C.Write("C4.dat");
    glob_chrono.Stop(VirtualTimer::FLUX);
    
    glob_chrono.Start(VirtualTimer::EXTRAPOL);
    
    AssembleVectorMaxwellDG(vars, var.evalHn_quad, C);
    
    glob_chrono.Stop(VirtualTimer::EXTRAPOL);
  }

  
  //! matrix vector product with the stiffness matrix for Maxwell equation and DG elements
  /*!
  // part \int_K \nabla \times E \varphi + \int_{\partial K} n \times [E] \cdot \varphi \f$
  \param[in] alpha multiplication coefficient
  \param[in] level not used
  \param[in] var instationary problem associated
  \param[in] B vector to be multiplied
  \param[in] beta multiplication coefficient
  \param[out] C resulting vector 
  */
  template<class TypeEquation>
  void MltAddStiffnessVectorial_MaxwellDG(const Real_wp& alpha, int level,
                                          HyperbolicProblem<TypeEquation>& var,
                                          const VectReal_wp& B, const Real_wp& beta,
                                          VectReal_wp& C)
  {
    if (beta == Real_wp(0))
      C.Fill(0);
    else
      Mlt(beta, C);
    
    // basic alias
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars
      = var.var_harmonic;
    
    ////////////////////////////////
    // term -\int_K \nabla u varphi
    ///////////////////////////////
    
    // loop on elements
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    MatrixVectorProductLevel& level_time_scheme = var.GetTimeLevelDistribution();
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	MltAdd_ElementRh1Maxwell(vars, num_elem, B, var.evalEn_quad(num_elem),
				 alpha, C, vars.GetReferenceElementH1(num_elem));
      }
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
    
    
    //////////////////////////////////////
    // term \int_{\partial K} [u] \varphi
    //////////////////////////////////////
    
        
    // loop on all elements
    glob_chrono.Start(VirtualTimer::FLUX);
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	MltAdd_ElementShMaxwell(vars, num_elem, B, var.evalEn_quad(num_elem),
				var.extrapolE, var.Eneighbor,
				var.extrapolH, var.Hneighbor,
				alpha, C, vars.GetReferenceElementH1(num_elem));
      }
    
    glob_chrono.Stop(VirtualTimer::FLUX);

    glob_chrono.Start(VirtualTimer::EXTRAPOL);

    AssembleVectorMaxwellDG(vars, var.evalEn_quad, C);
    
    glob_chrono.Stop(VirtualTimer::EXTRAPOL);
  }


  template<class TypeEquationTime, class Vector1, class Vector2, class Vector3, class Vector4>
  void GetExtrapolationMaxwellHDG(HyperbolicProblem<TypeEquationTime>& var, int level,
				  const VectReal_wp& E, const VectReal_wp& H,
				  Vector1& En_quad, Vector2& Hn_quad,
				  Vector3& extrapolE, Vector4& extrapolH, VectReal_wp& Lambda)
  {
    typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
    
    const EllipticProblem<TypeEquation>& vars = var.var_harmonic;
    const Mesh<Dimension3>& mesh = vars.mesh;
    const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
     
    MatrixVectorProductLevel& list_level = var.GetTimeLevelDistribution();
    list_level.SetLevel(level);
    
    /* if (En_quad.GetM() != mesh.GetNbElt())
      En_quad.Reallocate(mesh.GetNbElt());
    
    if (Hn_quad.GetM() != mesh.GetNbElt())
    Hn_quad.Reallocate(mesh.GetNbElt()); */
    
    int nb_dof_L = mesh_num.GetNbDof();
    Lambda.Reallocate(nb_dof_L);

    if (level < 0)
      Lambda.Zero();
    
    /* int size_extrapol = mesh_num.OffsetQuadElementNumber(mesh.GetNbElt());

    for (int k = 0; k < 3; k++)
      {
	extrapolE(k).Reallocate(size_extrapol);
	extrapolH(k).Reallocate(size_extrapol);
      }
    */
    
    VectReal_wp Eloc, Hloc, ProdLambda;
    TinyVector<VectReal_wp, 1> ProdLloc;
    for (int i0 = 0; i0 < list_level.GetNbElt(); i0++)
      {
	int i = list_level.GetElementNumber(i0);
	const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
	const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	
	// retrieving values of E and H
	int offset = vars.GetOffsetDofV(i);
	Eloc.Reallocate(Fb.GetNbDof());	    
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  Eloc(j) = E(offset + j);

	Hloc.Reallocate(Fb.GetNbDof());	    
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  Hloc(j) = H(offset + j);
	
	/* int ref_domain = mesh.Element(i).GetReference();
	   Real_wp tau = -vars.alpha_penalization;
	if (vars.upwind_fluxes)
	tau = -vars.alpha_penalization*vars.coefficient_impedance_absorbing(ref_domain); */
	
	// then loop over boundaries of the element
	ProdLloc(0).Reallocate(mesh_num.Element(i).GetNbDof());
	offset = 0;
	int offsetS = 0; TinyVector<Real_wp, 2> vec_us, vec_vs;
	for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
	    int num_face, ref_boundary, rf, rot; bool new_face; Real_wp phase, phase_conj;
	    vars.GetGeometryPhaseData(i, num_loc,
				      num_face, ref_boundary, rf, new_face, rot,
				      phase, phase_conj);
	    
	    const ElementReference<Dimension2, 2>& Fb_s = vars.GetSurfaceFiniteElementHcurl(num_face);
	    const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);
	    
	    if (rf != Fb.GetOrder())
	      {
		cout << "Variable order not implemented" << endl;
		abort();
	      }

	    bool dirichlet_boundary = false;
	    if ((ref_boundary != 0) &&
		(mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_DIRICHLET))
	      dirichlet_boundary = true;
	    
	    bool signF = true;
	    if ((num_loc == 0) || (num_loc == 2) || (num_loc == 4))
	      signF = false;

	    // forming - H \cdot n \times q - tau n \times E \cdot n \times q
	    ProdLambda.Reallocate(Fb_s.GetNbDof());
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
		
		kx = Fb_quad.GetXdofNumber(k);
		ky = Fb_quad.GetYdofNumber(k);
		
		vec_us.Init(Eloc(jx), Eloc(jy));
		
		Mlt(vars.Glob_matMass_elem(i).GetBhSurf(offsetS+k), vec_us, vec_vs);
			
		ProdLambda(kx) = -vec_vs(0);
		ProdLambda(ky) = -vec_vs(1);
		
		if (!dirichlet_boundary)
		  {
		    Real_wp poids = Fb_s.WeightsND(k);
		    if (signF)
		      poids = -poids;
		    
		    ProdLambda(ky) -= poids*Hloc(jx);
		    ProdLambda(kx) += poids*Hloc(jy);
		  }
	      }

	    int order = rf;
	    for (int j = 0; j < Fb_s.GetNbDof(); j++)
	      {
		if ((j >= 2*(order+1)) && (j < 4*(order+1)))
		  ProdLloc(0)(offset+j) = -ProdLambda(j);
		else
		  ProdLloc(0)(offset+j) = ProdLambda(j);
	      }
	    
	    offset += Fb_s.GetNbDof();
	    offsetS += Fb_s.GetNbPointsQuadratureInside();
	  }
	
	// values are added to the result Lambda
	mesh_num.number_map.AddLocalUnknownVector(mesh_num, Real_wp(1), ProdLloc,
						  i, Lambda);
      }  
  }
  

  template<class TypeEquation, class Vector1, class Vector2, class Vector3, class Vector4>
  void MltAddStiffnessMaxwellHDG(const Real_wp& alpha, int level, int level_vol,
				 HyperbolicProblem<TypeEquation>& var,
				 const VectReal_wp& En, const VectReal_wp& Hn, const VectReal_wp& Lambda_n,
				 const Vector1& En_quad, const Vector2& Hn_quad,
				 const Vector3& extrapolE, const Vector4& extrapolH,
				 const Real_wp& beta, VectReal_wp& ProdE, VectReal_wp& ProdH)
  {
    if (beta == Real_wp(0))
      {
	ProdE.Zero();
	ProdH.Zero();
      }
    else if (beta != Real_wp(1))
      {
	Mlt(beta, ProdE);
	Mlt(beta, ProdH);
      }
    
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var.var_harmonic;
    const Mesh<Dimension3>& mesh = vars.mesh;
    const MeshNumbering<Dimension3>& mesh_num = vars.GetMeshNumbering(0);
    
    MatrixVectorProductLevel& list_level = var.GetTimeLevelDistribution();
    list_level.SetLevel(level);

    MatrixVectorProductLevel level_inside(list_level);
    level_inside.SetLevel(level_vol);
    
    TinyVector<VectReal_wp, 1> Lloc;
    VectReal_wp Eloc, Hloc, curlH_quad, ProdEloc, Lambda;
    for (int i0 = 0; i0 < list_level.GetNbElt(); i0++)
      {
	int i = list_level.GetElementNumber(i0);
	
	int offset_E = vars.GetOffsetDofV(i);
	const ElementReference<Dimension3, 2>& Fb = vars.GetReferenceElementHcurl(i);
	
	int nb_dof_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	
	int nb_dof_L = mesh_num.Element(i).GetNbDof();
	Lloc(0).Reallocate(nb_dof_L); Lloc(0).Zero();
	Eloc.Reallocate(nb_dof_elt);
	Hloc.Reallocate(nb_dof_elt);
	ProdEloc.Reallocate(nb_dof_elt);
	curlH_quad.Reallocate(3*nb_points_quad);
	
	// we retrieve values of lambda on the element
	mesh_num.number_map.GetLocalUnknownVector(mesh_num, Lambda_n, i, Lloc);
	
	// we retrieve values of E and H on the element
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    Eloc(j) = En(offset_E + j);
	    Hloc(j) = Hn(offset_E + j);
	  }
	
	if (Fb.LumpedMassMatrix())
	  {
	    const HexahedronHcurlLobatto& Fb_hex = dynamic_cast<const HexahedronHcurlLobatto&>(Fb);
	    const Matrix<Real_wp, General, ArrayRowSparse>& R = Fb_hex.GetPermutedRh();
	    
	    // curl of H
	    //Fb.ApplyRhTranspose(Hloc, curlH_quad);
	    Mlt(SeldonTrans, R, Hloc, curlH_quad);
	    Mlt(R, Eloc, ProdEloc);
	    
	    for (int j = 0; j < nb_dof_elt; j++)
	      {
		ProdE(offset_E+j) += curlH_quad(j);
		ProdH(offset_E+j) -= ProdEloc(j);
	      }
	    
	    // loop over boundaries of the element
	    int offsetS = 0, offset_L = 0;
	    TinyVector<Real_wp, 2> vec_us, vec_vs;
	    for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
	      {
		int num_face, ref_boundary, rf, rot; bool new_face; Real_wp phase, phase_conj;
		vars.GetGeometryPhaseData(i, num_loc,
					  num_face, ref_boundary, rf, new_face, rot,
					  phase, phase_conj);
		
		const ElementReference<Dimension2, 2>& Fb_s = vars.GetSurfaceFiniteElementHcurl(num_face);
		const QuadrangleHcurlLobatto& Fb_quad = dynamic_cast<const QuadrangleHcurlLobatto&>(Fb_s);
		    
		int order = Fb_quad.GetOrder();
		Lambda.Reallocate(Fb_s.GetNbDof());
		for (int j = 0; j < Fb_s.GetNbDof(); j++)
		  {
		    if ((j >= 2*(order+1)) && (j < 4*(order+1)))
		      Lambda(j) = -Lloc(0)(offset_L+j);
		    else
		      Lambda(j) = Lloc(0)(offset_L+j);
		  }

		bool dirichlet_boundary = false;
		if ((ref_boundary != 0) &&
		    (mesh.GetBoundaryCondition(ref_boundary) == BoundaryConditionEnum::LINE_DIRICHLET))
		  dirichlet_boundary = true;
		
		bool signF = true;
		if ((num_loc == 0) || (num_loc == 2) || (num_loc == 4))
		  signF = false;
		
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
		    
		    kx = Fb_quad.GetXdofNumber(k);
		    ky = Fb_quad.GetYdofNumber(k);
		    
		    vec_us.Init(Eloc(jx)-Lambda(kx), Eloc(jy)-Lambda(ky));
		    Mlt(vars.Glob_matMass_elem(i).GetBhSurf(offsetS+k), vec_us, vec_vs);
			
		    ProdE(offset_E+jx) -= vec_vs(0);
		    ProdE(offset_E+jy) -= vec_vs(1);
		    
		    if (!dirichlet_boundary)
		      {
			Real_wp coef = Fb_s.WeightsND(k);
			if (signF)
			  coef = -coef;
			
			ProdH(offset_E+jx) -= coef*Lambda(ky);
			ProdH(offset_E+jy) += coef*Lambda(kx);
		      }
		  }
		
		offset_L += Fb_s.GetNbDof();
		offsetS += Fb_s.GetNbPointsQuadratureInside();
	      }
	  }
	else
	  {
	    cout << "Case not implemented" << endl;
	    abort();
	  }
      }
  }

} // end namespace

#define MONTJOIE_FILE_PRODMATVECT_MAXWELL3D_CXX
#endif
