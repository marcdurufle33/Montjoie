#ifndef MONTJOIE_FILE_PROD_MAT_VECT_ACOUSTIC_CXX

namespace Montjoie
{
  
  template<class TypeEquation, class Vector1,
           class Vector2, class Vector3, class T0, class Dimension>
  void MltAdd_ElementRh1(const EllipticProblem<TypeEquation>& var, int num_elem,
			 const Vector1& B, Vector2& Un_quad, Vector3& Vn_quad,
			 const T0& alpha, Vector1& C, const ElementReference<Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Complexe Complexe;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    
    bool affine = mesh.IsElementAffine(num_elem);
    
    int nb_points_elt = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
      nb_points_quad = nb_points_elt;
    
    int Nvol = mesh_num.GetNbDof();
    int Nvol_pml = mesh_num.GetNbDofPML();
    // getting values of B inside the element
    Vector1 Uloc(nb_points_elt), Uquad;
    int size_du = nb_points_quad*Dimension::dim_N;     
    Vector1 dU(size_du);
    dU.Fill(0);
    
    TinyVector<Vector1, Dimension::dim_N> V;
    for (int i = 0; i < Dimension::dim_N; i++)
      V(i).Reallocate(nb_points_quad);
	
    int offset_elt = mesh_num.Element(num_elem).GetNumberDof(0);
    
    // computation of grad(U)
    if (Fb.UseQuadraturePointsForSh())
      Fb.ApplyRhQuadratureTranspose(Un_quad(0), dU);
    else
      {
	for (int j = 0; j < nb_points_elt; j++)
	  Uloc(j) = B(offset_elt + j);
        
        if (var.InsidePML(num_elem))
          {
            // split PML, we compute ux+uy  (ux+uy+uz in 3-D)
            if (Dimension::dim_N == 2)
              for (int j = 0; j < nb_points_elt; j++)
                Uloc(j) += B(Nvol + mesh_num.GetDofPML(offset_elt + j));
            else
              for (int j = 0; j < nb_points_elt; j++)
                {
                  int num_dof = mesh_num.GetDofPML(offset_elt + j);
                  Uloc(j) += B(Nvol+num_dof) + B(Nvol+Nvol_pml+num_dof);
                }
          }
	
        if (Fb.UseQuadraturePointsForRh()
            || ( !affine && ElementReference_Base::use_warburton_trick) )
          {
            Uquad.Reallocate(nb_points_quad);
            Uquad.Fill(0);
            Fb.ApplyChTranspose(Uloc, Uquad);
          }
        
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
    
    TinyVector<Complexe, Dimension::dim_N> vec_u, grad_U; Real_wp coef;
    // application of transformation DF_i^{*-1}
    if (affine)
      {
	const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N>&
	  dfjm1 = var.Glob_DFjm1(num_elem)(0);
	
	for (int k = 0; k < nb_points_quad; k++)
	  {
	    CopyVector(dU, k, vec_u);
	    
	    MltTrans(dfjm1, vec_u, grad_U);
	    
	    CopyVector(grad_U, k, V);
	  }
	
	if (Fb.UseQuadraturePointsForRh())
	  for (int k = 0; k < nb_points_quad; k++)
	    {
	      MltVector(Fb.WeightsND(k), k, V);
	    }	
      }
    else
      for (int k = 0; k < nb_points_quad; k++)
	{
	  CopyVector(dU, k, vec_u);
	  
	  MltTrans(var.Glob_DFjm1(num_elem)(k), vec_u, grad_U);
	  
	  coef = Fb.WeightsND(k);
	  
	  // adding - \int 0.5 grad(J)/J \cdot psi u
	  if (ElementReference_Base::use_warburton_trick)
            {
	      Real_wp invJacob = coef/var.Glob_jacobian(num_elem)(k);
	      Mlt(coef*invJacob, grad_U);
	      
	      Real_wp val_u(0);
	      if (Fb.UseQuadraturePointsForSh())
		val_u = Un_quad(0)(k);
	      else
		val_u = Uquad(k);
	      
	      Add(-coef*val_u, mesh.Glob_GradJacobian(num_elem)(k), grad_U);
	      CopyVector(grad_U, k, V);
	    }
	  else
	    {
	      Mlt(coef, grad_U);
	      CopyVector(grad_U, k, V);
	    }
	}
    
    // application of Ch    
    if (Fb.UseQuadraturePointsForSh())
      {
	for (int n = 0; n < Dimension::dim_N; n++)
          Copy(V(n), Vn_quad(n));
      }
    else
      {
	for (int n = 0; n < Dimension::dim_N; n++)
	  {
	    if ( affine && !Fb.UseQuadraturePointsForRh())
	      Copy(V(n), Uloc);
	    else
	      Fb.ApplyCh(V(n), Uloc);
	    
	    for (int k = 0; k < nb_points_elt; k++)
	      C(offset_elt + k + n*Nvol) += alpha*Uloc(k);
	  }
      }
  }
  
  
#ifdef MONTJOIE_WITH_TWO_DIM
  template<class TypeEquation, class Vector1,
           class Vector2, class Vector3, class T0>
  void MltAdd_ElementRh2(const EllipticProblem<TypeEquation>& var,
                         int num_elem, const Vector1& B, Vector2& Un_quad,
                         Vector3& Vn_quad, const T0& alpha, Vector1& C,
			 const ElementReference<Dimension2, 1>& Fb)
  {
    typedef typename TypeEquation::Complexe Complexe;
    const Mesh<Dimension2>& mesh = var.mesh;
    const MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
    bool affine = mesh.IsElementAffine(num_elem);
    
    bool ortho = (var.OrthogonalElement(num_elem) == 0);
    bool ortho_z = (var.OrthogonalElement(num_elem) == 1);
    if (!var.InsidePML(num_elem))
      ortho = true;
    
    // number of dofs inside the element
    int nb_points_elt = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
      nb_points_quad = nb_points_elt;
    
    // getting values of V
    Vector1 Uphi(nb_points_quad), Uquad(nb_points_quad), Uh;
    Uphi.Fill(0); Uquad.Fill(0);
    TinyVector<Vector1, 2> dU, Vh;
    for (int n = 0; n < 2; n++)
      {
	dU(n).Reallocate(nb_points_elt);
	Vh(n).Reallocate(nb_points_quad);
	Vh(n).Fill(0);
      }
    
    int Nvol = mesh_num.GetNbDof();
    //int Nvol_pml = mesh_num.GetNbDofPML();
    int offset_elt = mesh_num.Element(num_elem).GetNumberDof(0);
    if (Fb.UseQuadraturePointsForSh())
      {
	for (int n = 0; n < 2; n++)
	  {
	    Vh(n) = Vn_quad(n);
	    if (affine)
	      for (int k = 0; k < nb_points_quad; k++)
		{
		  Real_wp coef = Fb.WeightsND(k);
		  Vh(n)(k) *= coef;
		}
	  }
      }
    else
      {
	for (int j = 0; j < nb_points_elt; j++)
	  {
	    for (int n = 0; n < 2; n++)
	      dU(n)(j) = B(offset_elt + j + n*Nvol);
	  }
        
        if (Fb.UseQuadraturePointsForRh() || (!affine))
          for (int n = 0; n < 2; n++)
            Fb.ApplyChTranspose(dU(n), Vh(n));
        else
          for (int n = 0; n < 2; n++)
            Copy(dU(n), Vh(n));        
      }

    TinyVector<Complexe, 2> vec_u, vec_v; Real_wp coef;
    // application of transformation DF_i^{-1}
    Vector1 Vloc(2*nb_points_quad), Vtilde;
    if (ortho_z)
      Vtilde.Reallocate(2*nb_points_quad);
        
    if (affine)
      {
	const TinyMatrix<Real_wp, General, 2, 2>&
	  dfjm1 = var.Glob_DFjm1(num_elem)(0);
	
        if (ortho)
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                CopyVector(Vh, k, vec_u);
                
                Mlt(dfjm1, vec_u, vec_v);
                Mlt(-1.0, vec_v);
                if (Fb.UseQuadraturePointsForRh())
                  vec_v *= Fb.WeightsND(k);
                
                CopyVector(vec_v, k, Vloc);
              }
          }
      }
    else
      {
        if (ortho)
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                CopyVector(Vh, k, vec_u);
                
                Mlt(var.Glob_DFjm1(num_elem)(k), vec_u, vec_v);
                
                coef = Fb.WeightsND(k);
                Mlt(-coef, vec_v);
                
                if (ElementReference_Base::use_warburton_trick)
                  {
                    Real_wp invJacob = coef/var.Glob_jacobian(num_elem)(k);
                    Mlt(invJacob, vec_v);
                    CopyVector(vec_v, k, Vloc);
                    
                    Uphi(k) = coef*DotProd(vec_u, mesh.Glob_GradJacobian(num_elem)(k));
                  }
                else
                  CopyVector(vec_v, k, Vloc);
              }
          }
      }
    
    // integration against grad phi
    Vector1 Uloc(nb_points_elt);
    Uloc.Fill(0);
    
    if (var.InsidePML(num_elem))
      {
        if (Fb.UseQuadraturePointsForRh())
          {
	    if (!ortho)
	      {
		abort();
	      }
	    
	    Vector1 Uy(nb_points_quad);
	    Uy.Fill(0);
	    Fb.ApplyRhQuadratureSplit(Vloc, Uquad, Uy);
	    if (Fb.UseQuadraturePointsForSh())
	      {
		Un_quad(0) = Uquad;
		Un_quad(1) = Uy;
	      }
	    else
	      {                
		Fb.ApplyCh(Uquad, Uloc);
		for (int k = 0; k < nb_points_elt; k++)
		  C(offset_elt + k) += alpha*Uloc(k);
		
		Fb.ApplyCh(Uy, Uloc);
		for (int k = 0; k < nb_points_elt; k++)
		  C(Nvol + mesh_num.GetDofPML(offset_elt + k)) += alpha*Uloc(k);
	      }                
          }
        else
          {	
	    if (!ortho)
	      {
		abort();
	      }
	    
	    Vector1 Uy(nb_points_elt); Uy.Fill(0);
	    if (affine)
	      Fb.ApplyConstantRhSplit(Vloc, Uloc, Uy);
	    else
	      Fb.ApplyRhSplit(Vloc, Uloc, Uy);
	    
	    for (int k = 0; k < nb_points_elt; k++)
	      C(offset_elt + k) += alpha*Uloc(k);
	    
	    for (int k = 0; k < nb_points_elt; k++)
	      C(Nvol + mesh_num.GetDofPML(offset_elt + k)) += alpha*Uy(k);
          }
      }
    else
      {
        if (Fb.UseQuadraturePointsForRh())
          {
            VectReal_wp Vtmp = Vloc;
            Fb.ApplyRhQuadrature(Vloc, Uquad);
            if (Fb.UseQuadraturePointsForSh())
              {
                Un_quad(0) = Uquad;
                if ((!affine) && (ElementReference_Base::use_warburton_trick))
                  {
                    for (int k = 0; k < nb_points_quad; k++)
                      Un_quad(0)(k) += Uphi(k);
                  }
              }
            else
              {
                if ((!affine) && (ElementReference_Base::use_warburton_trick))
                  {
                    for (int k = 0; k < nb_points_quad; k++)
                      Uquad(k) += Uphi(k);
                  }
                
                Fb.ApplyCh(Uquad, Uloc);
                
                for (int k = 0; k < nb_points_elt; k++)
                  C(offset_elt + k) += alpha*Uloc(k);
              }
          }
        else
          {	
            if (affine)
              Fb.ApplyConstantRh(Vloc, Uloc);
            else
              {
                Fb.ApplyRh(Vloc, Uloc);
                
                if (ElementReference_Base::use_warburton_trick)
                  {
                    Uh.Reallocate(nb_points_elt);
                    Uh.Fill(0);
                    Fb.ApplyCh(Uphi, Uh);
                    for (int k = 0; k < nb_points_elt; k++)
                      Uloc(k) += Uh(k);
                  }
              }
            
            for (int k = 0; k < nb_points_elt; k++)
              C(offset_elt + k) += alpha*Uloc(k);
          }
      }
  }
#endif  

#ifdef MONTJOIE_WITH_THREE_DIM
  template<class TypeEquation, class Vector1,
           class Vector2, class Vector3, class T0>
  void MltAdd_ElementRh2(const EllipticProblem<TypeEquation>& var,
                         int num_elem, const Vector1& B, Vector2& Un_quad,
                         Vector3& Vn_quad, const T0& alpha, Vector1& C,
			 const ElementReference<Dimension3, 1>& Fb)
  {
    typedef typename TypeEquation::Complexe Complexe;
    const Mesh<Dimension3>& mesh = var.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    bool affine = mesh.IsElementAffine(num_elem);
    
    bool ortho = (var.OrthogonalElement(num_elem) == 0);
    bool ortho_z = (var.OrthogonalElement(num_elem) == 1);
    if (!var.InsidePML(num_elem))
      ortho = true;
    
    // number of dofs inside the element
    int nb_points_elt = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
      nb_points_quad = nb_points_elt;
    
    // getting values of V
    Vector1 Uphi(nb_points_quad), Uquad(nb_points_quad), Uh;
    Uphi.Fill(0); Uquad.Fill(0);
    TinyVector<Vector1, 3> dU, Vh;
    for (int n = 0; n < 3; n++)
      {
	dU(n).Reallocate(nb_points_elt);
	Vh(n).Reallocate(nb_points_quad);
	Vh(n).Fill(0);
      }
    
    int Nvol = mesh_num.GetNbDof();
    int Nvol_pml = mesh_num.GetNbDofPML();
    int offset_elt = mesh_num.Element(num_elem).GetNumberDof(0);
    if (Fb.UseQuadraturePointsForSh())
      {
	for (int n = 0; n < 3; n++)
	  {
	    Vh(n) = Vn_quad(n);
	    if (affine)
	      for (int k = 0; k < nb_points_quad; k++)
		{
		  Real_wp coef = Fb.WeightsND(k);
		  Vh(n)(k) *= coef;
		}
	  }
      }
    else
      {
	for (int j = 0; j < nb_points_elt; j++)
	  {
	    for (int n = 0; n < 3; n++)
	      dU(n)(j) = B(offset_elt + j + n*Nvol);
	  }
        
        if (Fb.UseQuadraturePointsForRh() || (!affine))
          for (int n = 0; n < 3; n++)
            Fb.ApplyChTranspose(dU(n), Vh(n));
        else
          for (int n = 0; n < 3; n++)
            Copy(dU(n), Vh(n));        
      }

    TinyVector<Complexe, 3> vec_u, vec_v; Real_wp coef;
    // application of transformation DF_i^{-1}
    Vector1 Vloc(3*nb_points_quad), Vtilde, Vterce;
    if (ortho_z)
      Vtilde.Reallocate(3*nb_points_quad);
    else if (!ortho)
      {
        Vtilde.Reallocate(3*nb_points_quad);
        Vterce.Reallocate(3*nb_points_quad);
      }
    
    if (affine)
      {
	const TinyMatrix<Real_wp, General, 3, 3>& dfjm1 = var.Glob_DFjm1(num_elem)(0);
	
        if (ortho)
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                CopyVector(Vh, k, vec_u);
                
                Mlt(dfjm1, vec_u, vec_v);
                Mlt(-1.0, vec_v);
                if (Fb.UseQuadraturePointsForRh())
                  vec_v *= Fb.WeightsND(k);
                
                CopyVector(vec_v, k, Vloc);
              }
          }
        else if (ortho_z)
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                CopyVector(Vh, k, vec_u);
                
                Mlt(-1.0, vec_u);
                if (Fb.UseQuadraturePointsForRh())
                  vec_u *= Fb.WeightsND(k);
                
                vec_v(0) = dfjm1(0, 0)*vec_u(0);
                vec_v(1) = dfjm1(1, 1)*vec_u(1);
                vec_v(2) = dfjm1(2, 2)*vec_u(2);

                CopyVector(vec_v, k, Vloc);
                
                vec_v(0) = dfjm1(0, 1)*vec_u(1);
                vec_v(1) = dfjm1(1, 0)*vec_u(0);
                vec_v(2) = 0;
                
                CopyVector(vec_v, k, Vtilde);
              }
          }
        else
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                CopyVector(Vh, k, vec_u);
                
                Mlt(-1.0, vec_u);
                if (Fb.UseQuadraturePointsForRh())
                  vec_u *= Fb.WeightsND(k);
                
                vec_v(0) = dfjm1(0, 0)*vec_u(0);
                vec_v(1) = dfjm1(1, 1)*vec_u(1);
                vec_v(2) = dfjm1(2, 2)*vec_u(2);

                CopyVector(vec_v, k, Vloc);
                
                vec_v(0) = dfjm1(0, 1)*vec_u(1);
                vec_v(1) = dfjm1(1, 0)*vec_u(0);
                vec_v(2) = dfjm1(2, 0)*vec_u(0);
                
                CopyVector(vec_v, k, Vtilde);

                vec_v(0) = dfjm1(0, 2)*vec_u(2);
                vec_v(1) = dfjm1(1, 2)*vec_u(2);
                vec_v(2) = dfjm1(2, 1)*vec_u(1);
                
                CopyVector(vec_v, k, Vterce);
              }
          }
      }
    else
      {
        if (ortho)
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                CopyVector(Vh, k, vec_u);
                
                Mlt(var.Glob_DFjm1(num_elem)(k), vec_u, vec_v);
                
                coef = Fb.WeightsND(k);
                Mlt(-coef, vec_v);
                
                if (ElementReference_Base::use_warburton_trick)
                  {
                    Real_wp invJacob = coef/var.Glob_jacobian(num_elem)(k);
                    Mlt(invJacob, vec_v);
                    CopyVector(vec_v, k, Vloc);
                    
                    Uphi(k) = coef*DotProd(vec_u, mesh.Glob_GradJacobian(num_elem)(k));
                  }
                else
                  CopyVector(vec_v, k, Vloc);
              }
          }
        else if (ortho_z)
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                const TinyMatrix<Real_wp, General, 3, 3>&
                  dfjm1 = var.Glob_DFjm1(num_elem)(k);

                CopyVector(Vh, k, vec_u);
                
                Mlt(-Fb.WeightsND(k), vec_u);
                
                vec_v(0) = dfjm1(0, 0)*vec_u(0);
                vec_v(1) = dfjm1(1, 1)*vec_u(1);
                vec_v(2) = dfjm1(2, 2)*vec_u(2);

                CopyVector(vec_v, k, Vloc);
                
                vec_v(0) = dfjm1(0, 1)*vec_u(1);
                vec_v(1) = dfjm1(1, 0)*vec_u(0);
                vec_v(2) = 0;
                
                CopyVector(vec_v, k, Vtilde);
              }
          }
        else
          {
            for (int k = 0; k < nb_points_quad; k++)
              {
                const TinyMatrix<Real_wp, General, 3, 3>&
                  dfjm1 = var.Glob_DFjm1(num_elem)(k);

                CopyVector(Vh, k, vec_u);
                
                Mlt(-Fb.WeightsND(k), vec_u);
                
                vec_v(0) = dfjm1(0, 0)*vec_u(0);
                vec_v(1) = dfjm1(1, 1)*vec_u(1);
                vec_v(2) = dfjm1(2, 2)*vec_u(2);

                CopyVector(vec_v, k, Vloc);
                
                vec_v(0) = dfjm1(0, 1)*vec_u(1);
                vec_v(1) = dfjm1(1, 0)*vec_u(0);
                vec_v(2) = dfjm1(2, 0)*vec_u(0);
                
                CopyVector(vec_v, k, Vtilde);

                vec_v(0) = dfjm1(0, 2)*vec_u(2);
                vec_v(1) = dfjm1(1, 2)*vec_u(2);
                vec_v(2) = dfjm1(2, 1)*vec_u(1);
                
                CopyVector(vec_v, k, Vterce);
              }
          }
      }
    
    // integration against grad phi
    Vector1 Uloc(nb_points_elt);
    Uloc.Fill(0);
    
    if (var.InsidePML(num_elem))
      {
        if (Fb.UseQuadraturePointsForRh())
          {
	    Vector1 Uy(nb_points_quad), Uz(nb_points_quad);
	    Uy.Fill(0); Uz.Fill(0);
	    Fb.ApplyRhQuadratureSplit(Vloc, Uquad, Uy, Uz);
	    if (ortho_z)
	      {
		Vector1 Uquad2(nb_points_quad), Uy2(nb_points_quad), Uz2(nb_points_quad);
		Fb.ApplyRhQuadratureSplit(Vtilde, Uquad2, Uy2, Uz2);
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    Uquad(k) += Uy2(k);
		    Uy(k) += Uquad2(k);
		  }
	      }
	    else if (!ortho)
	      {
		Vector1 Uquad2(nb_points_quad), Uy2(nb_points_quad), Uz2(nb_points_quad);
		Fb.ApplyRhQuadratureSplit(Vtilde, Uquad2, Uy2, Uz2);
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    Uquad(k) += Uy2(k) + Uz2(k);
		    Uy(k) += Uquad2(k);
		  }
		
		Fb.ApplyRhQuadratureSplit(Vterce, Uquad2, Uy2, Uz2);
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    Uy(k) += Uz2(k);
		    Uz(k) += Uquad2(k) + Uy2(k);
		  }
	      }
	    
	    if (Fb.UseQuadraturePointsForSh())
	      {
		Un_quad(0) = Uquad;
		Un_quad(1) = Uy;
		Un_quad(2) = Uz;
	      }
	    else
	      {                
		Fb.ApplyCh(Uquad, Uloc);
		for (int k = 0; k < nb_points_elt; k++)
		  C(offset_elt + k) += alpha*Uloc(k);
		
		Fb.ApplyCh(Uy, Uloc);
		for (int k = 0; k < nb_points_elt; k++)
		  C(Nvol + mesh_num.GetDofPML(offset_elt + k)) += alpha*Uy(k); 
		
		Fb.ApplyCh(Uz, Uloc);
		for (int k = 0; k < nb_points_elt; k++)
		  C(Nvol + Nvol_pml + mesh_num.GetDofPML(offset_elt + k)) += alpha*Uz(k); 
	      }                           
          }
        else
          {	
	    Vector1 Uy(nb_points_elt), Uz(nb_points_elt);
	    Uy.Fill(0); Uz.Fill(0);
	    if (affine)
	      Fb.ApplyConstantRhSplit(Vloc, Uloc, Uy, Uz);
	    else
	      Fb.ApplyRhSplit(Vloc, Uloc, Uy, Uz);
	    
	    if (ortho_z)
	      {
		Vector1 Uquad2(nb_points_elt), Uy2(nb_points_elt), Uz2(nb_points_elt);
		if (affine)
		  Fb.ApplyConstantRhSplit(Vtilde, Uquad2, Uy2, Uz2);
		else
		  Fb.ApplyRhSplit(Vtilde, Uquad2, Uy2, Uz2);
		
		for (int k = 0; k < nb_points_elt; k++)
		  {
		    Uquad(k) += Uy2(k);
		    Uy(k) += Uquad2(k);
		  }
	      }
	    else if (!ortho)
	      {
		Vector1 Uquad2(nb_points_elt), Uy2(nb_points_elt), Uz2(nb_points_elt);
		if (affine)
		  Fb.ApplyConstantRhSplit(Vtilde, Uquad2, Uy2, Uz2);
		else
		  Fb.ApplyRhSplit(Vtilde, Uquad2, Uy2, Uz2);
		
		for (int k = 0; k < nb_points_elt; k++)
		  {
		    Uquad(k) += Uy2(k) + Uz2(k);
		    Uy(k) += Uquad2(k);
		  }
		
		if (affine)
		  Fb.ApplyConstantRhSplit(Vterce, Uquad2, Uy2, Uz2);
		else
		  Fb.ApplyRhSplit(Vterce, Uquad2, Uy2, Uz2);
		
		for (int k = 0; k < nb_points_elt; k++)
		  {
		    Uy(k) += Uz2(k);
		    Uz(k) += Uquad2(k) + Uy2(k);
		  }
	      }
                
	    for (int k = 0; k < nb_points_elt; k++)
	      C(offset_elt + k) += alpha*Uloc(k);
	    
	    for (int k = 0; k < nb_points_elt; k++)
	      C(Nvol + mesh_num.GetDofPML(offset_elt + k)) += alpha*Uy(k);
	    
	    for (int k = 0; k < nb_points_elt; k++)
	      C(Nvol + Nvol_pml + mesh_num.GetDofPML(offset_elt + k)) += alpha*Uz(k);
                
	  }
      }
    else
      {
        if (Fb.UseQuadraturePointsForRh())
          {
            VectReal_wp Vtmp = Vloc;
            Fb.ApplyRhQuadrature(Vloc, Uquad);
            if (Fb.UseQuadraturePointsForSh())
              {
                Un_quad(0) = Uquad;
                if ((!affine) && (ElementReference_Base::use_warburton_trick))
                  {
                    for (int k = 0; k < nb_points_quad; k++)
                      Un_quad(0)(k) += Uphi(k);
                  }
              }
            else
              {
                if ((!affine) && (ElementReference_Base::use_warburton_trick))
                  {
                    for (int k = 0; k < nb_points_quad; k++)
                      Uquad(k) += Uphi(k);
                  }
                
                Fb.ApplyCh(Uquad, Uloc);
                
                for (int k = 0; k < nb_points_elt; k++)
                  C(offset_elt + k) += alpha*Uloc(k);
              }
          }
        else
          {	
            if (affine)
              Fb.ApplyConstantRh(Vloc, Uloc);
            else
              {
                Fb.ApplyRh(Vloc, Uloc);
                
                if (ElementReference_Base::use_warburton_trick)
                  {
                    Uh.Reallocate(nb_points_elt);
                    Uh.Fill(0);
                    Fb.ApplyCh(Uphi, Uh);
                    for (int k = 0; k < nb_points_elt; k++)
                      Uloc(k) += Uh(k);
                  }
              }
            
            for (int k = 0; k < nb_points_elt; k++)
              C(offset_elt + k) += alpha*Uloc(k);
          }
      }
  }
#endif

  
  template<class TypeEquation, class Vector1, class Vector2,
	   class Vector5, class Vector6, class T0, class Dimension>
  void MltAdd_ElementSh(const EllipticProblem<TypeEquation>& var, int num_elem1,
			const Vector1& B, Vector2& Un_quad,
			const Vector5& extrapolU, const Vector5& Uneighbor,
			const Vector6& extrapolV, const Vector6& Vneighbor,
			const T0& alpha, Vector1& C, const ElementReference<Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Complexe Complexe;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    
    int nb_points_elt = Fb.GetNbDof();
    TinyVector< Vector1, Dimension::dim_N> Vh;
    
    int offset_face1 = mesh_num.OffsetQuadElementNumber(num_elem1);
    int Nvol = mesh_num.GetNbDof();
    TinyVector<Real_wp, Dimension::dim_N> normale;    
    TinyVector<Complexe, Dimension::dim_N> jumpV, valV1, valV2, fV_n;
    Real_wp valU1, valU2, fV, jumpU, dsj, v_dot_n;
    
    bool presence_penalization = (var.delta_penalization != Real_wp(0));
    int ref_domain = mesh.Element(num_elem1).GetReference();
    bool affine = mesh.IsElementAffine(num_elem1);
    int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();
    int offset_elt = mesh_num.Element(num_elem1).GetNumberDof(0);
    
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
	const VectReal_wp& PoidsFlux
          = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
	const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	  GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));

	for (int n = 0; n < Dimension::dim_N; n++)
          {
            Vh(n).Reallocate(nb_points_face);
            Vh(n).Fill(0);
          }
	          
        int num_pos2_face = -1, krot = 0;
	int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0, rot2, ref2 = -1;
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

	    ref2 = mesh.Element(num_elem2).GetReference();
	    rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
	    rot = mesh_num.GetRotationFace(rot1, rot2, nv);
          }
	else
	  offset_face2 = offset_face1;
	
	if (cond != BoundaryConditionEnum::LINE_INSIDE)
	  {
            // face with a boundary condition
	    for (int k = 0; k < nb_points_face; k++)
	      {
		// internal term
		int num_dof1 = offset_face1 + k;
                krot = FacesQuadRotation(rot, k);
		
		// value of u
		valU1 = extrapolU(0)(num_dof1);
					    
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
		
		fV = valU1;
		
		// term of the boundary condition
		if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
		  {
		    fV *= 2.0;
		  }
		else if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
		  {
                    int offset_neighbor = var.GetOffsetNeighboringFace(num_face);
		    ref2 = var.GetRefDomainNeighboringFace(num_face);
		    int num_dof2 = offset_neighbor + k;
		    valU2 = Uneighbor(0)(num_dof2);
		    fV -= valU2;
		    if (presence_penalization)
		      {
			// jump of V
			CopyVector(extrapolV, num_dof1, valV1);
			CopyVector(Vneighbor, num_dof2, valV2);
			jumpV = valV2 - valV1;
			v_dot_n = DotProd(jumpV, normale);
			Real_wp imped = var.delta_penalization;
			if (var.upwind_fluxes)
			  {
			    Real_wp Ym = var.coefficient_impedance_absorbing(ref_domain);
			    Real_wp Yp = var.coefficient_impedance_absorbing(ref2);	
			    Real_wp inv_Ybar = 1.0/(Ym + Yp);		       
			    imped = -2.0*inv_Ybar;
			    
			    if (ref_domain != ref2)
			      {
				Real_wp coef_v = (Yp - Ym)*inv_Ybar;
				fV -= (valU2 - valU1)*coef_v;
			      }
			  }
			
			fV += v_dot_n*imped;
		      }
		  }
		else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
		  {
                    // jump of V
		    CopyVector(extrapolV, num_dof1, valV1);
                    v_dot_n = DotProd(valV1, normale);
		    Real_wp imped = 1.0/var.coefficient_impedance_absorbing(ref_domain);
		    
		    fV += v_dot_n*imped;
                  }
		else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
		  {
		    // neumann
		    fV = 0;
		  }
		
		fV_n = normale; Mlt(dsj*PoidsFlux(k)*fV, fV_n);
                                
		CopyVector(fV_n, k, Vh);
	      }
	    
	  }
	else if (mesh.Boundary(num_face).GetNbElements()==2)
	  {
	    for (int k = 0; k < nb_points_face; k++)
	      {
                krot = FacesQuadRotation(rot, k);
		int num_dof1 = offset_face1 + k;
		int num_dof2 = krot + offset_face2;
		
		// jump -> (u2 - u1)
		jumpU = -extrapolU(0)(num_dof1) + extrapolU(0)(num_dof2);
                
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
		
		fV = -jumpU;
		if (presence_penalization)
		  {
                    // jump of V
		    CopyVector(extrapolV, num_dof1, valV1);
		    CopyVector(extrapolV, num_dof2, valV2);
		    jumpV = valV2 - valV1;
		    v_dot_n = DotProd(jumpV, normale);
		    
		    Real_wp imped = var.delta_penalization;
		    if (var.upwind_fluxes)
		      {
			Real_wp Ym = var.coefficient_impedance_absorbing(ref_domain);
			Real_wp Yp = var.coefficient_impedance_absorbing(ref2);	
			Real_wp inv_Ybar = 1.0/(Ym + Yp);		       
			imped = -2.0*inv_Ybar;
			
			if (ref_domain != ref2)
			  {
			    Real_wp coef_v = (Yp - Ym)*inv_Ybar;
			    fV -= jumpU*coef_v;
			  }
		      }
		    
		    fV += v_dot_n*imped;
		  }
		
		fV_n = normale; Mlt(PoidsFlux(k)*dsj*fV, fV_n);
                CopyVector(fV_n, k, Vh);
	      }
	  }
        
        offset_face1 += nb_points_face;
        
        if ((!affine)&&(ElementReference_Base::use_warburton_trick))
          {
            for (int k = 0; k < nb_points_face; k++)
              {
                Real_wp jacob = mesh.Glob_invSqrtJacobian(num_elem1)(offset_war + k);
                MltVector(jacob, k, Vh);
              }
            
            offset_war += nb_points_face;
          }

        if (Fb.UseQuadraturePointsForSh())
          {	
            for (int n = 0; n < Dimension::dim_N; n++)
              Fb.ApplyShQuadrature(-1.0, num_pos1_face, Vh(n), Un_quad(n), rf);
          }
        else
          {
            VectReal_wp Cloc;
            for (int n = 0; n < Dimension::dim_N; n++)
              {
                int pos = offset_elt + n*Nvol;
                Cloc.SetData(nb_points_elt, &C(pos)); 
                Fb.ApplySh(-alpha, num_pos1_face, Vh(n), Cloc, rf);
                Cloc.Nullify();
              }
          }
      }
    
    if (Fb.UseQuadraturePointsForSh())
      {
        VectReal_wp Cloc(nb_points_elt);
        for (int n = 0; n < Dimension::dim_N; n++)
          {
            int offset = offset_elt + n*Nvol;
            Fb.ApplyCh(Un_quad(n), Cloc);
            for (int k = 0; k < nb_points_elt; k++)
              C(offset + k) += alpha*Cloc(k);
          }
      }
  }
				  
  
  template<class TypeEquation, class Vector1, class Vector2,
	   class Vector5, class Vector6, class T0, class Dimension>
  void MltAdd_ElementSht(const EllipticProblem<TypeEquation>& var, int num_elem1,
			 const Vector1& B, Vector2& Un_quad,
			 const Vector5& extrapolU, const Vector5& Uneighbor, 
			 const Vector6& extrapolV, const Vector6& Vneighbor, 
			 const T0& alpha, Vector1& C, const ElementReference<Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Complexe Complexe;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    
    int nb_points_elt = Fb.GetNbDof();
    
    Vector1 Uh, Uy, Uz;
    int offset_face1 = mesh_num.OffsetQuadElementNumber(num_elem1);
    TinyVector<Real_wp, Dimension::dim_N> normale;
    TinyVector<Complexe, Dimension::dim_N> valV1, valV2, fU_n;
    Complexe valU1, fU, jumpU; Real_wp dsj;
    bool presence_penalization = (var.alpha_penalization != Real_wp(0));
    bool pml = var.InsidePML(num_elem1);
    int ref_domain = mesh.Element(num_elem1).GetReference();
    bool affine = mesh.IsElementAffine(num_elem1);
    int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();
    int pos = mesh_num.Element(num_elem1).GetNumberDof(0);
    int Nvol = mesh_num.GetNbDof();
    int Nvol_pml = mesh_num.GetNbDofPML();
    
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
	const VectReal_wp& PoidsFlux
          = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
	const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	  GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
        
        Uh.Reallocate(nb_points_face);
        Uh.Fill(0);
        if (pml)
          {
            Uy.Reallocate(nb_points_face);
            Uz.Reallocate(nb_points_face);
            Uy.Fill(0);
            Uz.Fill(0);
          }
        
        int num_pos2_face = -1, krot = 0;
	int nv = mesh.Boundary(num_face).GetNbVertices(), rot = 0, rot2, ref2 = -1;
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

	    ref2 = mesh.Element(num_elem2).GetReference();
	    rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
            rot = mesh_num.GetRotationFace(rot1, rot2, nv);
          }
	else
	  offset_face2 = offset_face1;
	
	if (cond != BoundaryConditionEnum::LINE_INSIDE)
	  {
	    // face with a boundary condition
	    for (int k = 0; k < nb_points_face; k++)
	      {
		// internal term
		int num_dof1 = offset_face1 + k;
                krot = FacesQuadRotation(rot, k);
                
		// value of V
		CopyVector(extrapolV, num_dof1, valV1);
		fU_n = valV1;
		
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
                
		// term of the boundary condition
		if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
		  {
		    Mlt(2.0, fU_n);
		  }
		else if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
		  {
                    int offset_neighbor = var.GetOffsetNeighboringFace(num_face);
		    ref2 = var.GetRefDomainNeighboringFace(num_face);
		    int num_dof2 = offset_neighbor + k;
		    CopyVector(Vneighbor, num_dof2, valV2);
		    fU_n += valV2;
		    
		    if (presence_penalization)
		      {
			jumpU = -extrapolU(0)(num_dof1) + Uneighbor(0)(num_dof2);
			
                        Real_wp imped = var.alpha_penalization;
			if (var.upwind_fluxes)
			  {
			    Real_wp Ym = var.coefficient_impedance_absorbing(ref_domain);
			    Real_wp Yp = var.coefficient_impedance_absorbing(ref2);
			    Real_wp Zm = 1.0/Ym;
			    Real_wp Zp = 1.0/Yp;
			    Real_wp inv_Zbar = 1.0/(Zm + Zp); 
			    imped = -2.0*inv_Zbar;
			    
			    if (ref2 != ref_domain)
			      {
				Real_wp coef_u = (Zp - Zm)*inv_Zbar;
				fU_n += coef_u*(valV2 - valV1);
			      }
			  }
			
			Add(-jumpU*imped, normale, fU_n);
		      }
		  }
		else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
		  {
		    valU1 = -extrapolU(0)(num_dof1);
                    
		    Real_wp imped = var.coefficient_impedance_absorbing(ref_domain);
                    Add(valU1*imped, normale, fU_n);
		  }
		else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
		  {
		    fU_n.Fill(0);
		  }
		
		if (pml)
                  {
                    Uh(k) = fU_n(0)*normale(0)*PoidsFlux(k)*dsj;
                    Uy(k) = fU_n(1)*normale(1)*PoidsFlux(k)*dsj;
                    if (Dimension::dim_N == 3)
                      Uz(k) = fU_n(2)*normale(2)*PoidsFlux(k)*dsj;
                  }
                else
                  {
                    fU = DotProd(fU_n, normale)*PoidsFlux(k)*dsj;
                    Uh(k) = fU;
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
		int num_dof2 = krot + offset_face2;
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
		
		// mean value (v1+v2)
		CopyVector(extrapolV, num_dof1, valV1);
		CopyVector(extrapolV, num_dof2, valV2);
		Add(valV1, valV2, fU_n);
                		
		if (presence_penalization)
		  {
		    jumpU = -extrapolU(0)(num_dof1) + extrapolU(0)(num_dof2);
		    
                    Real_wp imped = var.alpha_penalization;
		    if (var.upwind_fluxes)
		      {
			Real_wp Ym = var.coefficient_impedance_absorbing(ref_domain);
			Real_wp Yp = var.coefficient_impedance_absorbing(ref2);
			Real_wp Zm = 1.0/Ym;
			Real_wp Zp = 1.0/Yp;
			Real_wp inv_Zbar = 1.0/(Zm + Zp); 
			imped = -2.0*inv_Zbar;

			if (ref2 != ref_domain)
			  {
			    Real_wp coef_u = (Zp - Zm)*inv_Zbar;
			    fU_n += coef_u*(valV2 - valV1);
			  }
		      }
		    
		    Add(-jumpU*imped, normale, fU_n);
		  }
                
                if (pml)
                  {
                    Uh(k) = fU_n(0)*normale(0)*PoidsFlux(k)*dsj;
                    Uy(k) = fU_n(1)*normale(1)*PoidsFlux(k)*dsj;
                    if (Dimension::dim_N == 3)
                      Uz(k) = fU_n(2)*normale(2)*PoidsFlux(k)*dsj;
                  }
                else
                  {
                    fU = DotProd(fU_n, normale)*PoidsFlux(k)*dsj;		
                    Uh(k) = fU;
                  }

	      }
            
	  }
        
        offset_face1 += nb_points_face;

        if ((!affine)&&(ElementReference_Base::use_warburton_trick))
          {
            for (int k = 0; k < nb_points_face; k++)
              Uh(k) *= mesh.Glob_invSqrtJacobian(num_elem1)(offset_war + k);
            
            offset_war += nb_points_face;
          }
        
        if (Fb.UseQuadraturePointsForSh())
          {
            Fb.ApplyShQuadrature(1.0, num_pos1_face, Uh, Un_quad(0), rf);
            if (pml)
              {
                Fb.ApplyShQuadrature(1.0, num_pos1_face, Uy, Un_quad(1), rf);
                if (Dimension::dim_N == 3)
                  Fb.ApplyShQuadrature(1.0, num_pos1_face, Uz, Un_quad(2), rf);
              }
          }
        else
          {
            VectReal_wp Cloc;
            Cloc.SetData(nb_points_elt, &C(pos)); 
            Fb.ApplySh(alpha, num_pos1_face, Uh, Cloc, rf);            
            Cloc.Nullify();

            if (pml)
              {
                Cloc.SetData(nb_points_elt, &C(Nvol + mesh_num.GetDofPML(pos)));
                Fb.ApplySh(alpha, num_pos1_face, Uy, Cloc, rf);
                Cloc.Nullify();                
                
                if (Dimension::dim_N == 3)
                  {
                    Cloc.SetData(nb_points_elt, &C(Nvol + Nvol_pml + mesh_num.GetDofPML(pos)));
                    Fb.ApplySh(alpha, num_pos1_face, Uz, Cloc, rf);
                    Cloc.Nullify();                 
                  }
              }
          }
      }
    
    if (Fb.UseQuadraturePointsForSh())
      {
        VectReal_wp Cloc(nb_points_elt);
        Fb.ApplyCh(Un_quad(0), Cloc);
        for (int i = 0; i < nb_points_elt; i++)
          C(pos + i) += alpha*Cloc(i);
        
        if (pml)
          {
            Fb.ApplyCh(Un_quad(1), Cloc);
            int offset = Nvol + mesh_num.GetDofPML(pos);
            for (int i = 0; i < nb_points_elt; i++)
              C(offset + i) += alpha*Cloc(i);

            if (Dimension::dim_N == 3)
              {
                Fb.ApplyCh(Un_quad(2), Cloc);
                offset += Nvol_pml;
                for (int i = 0; i < nb_points_elt; i++)
                  C(offset + i) += alpha*Cloc(i);
              }
          }
      }


  }
  

  //! computation of unknowns on the quadrature points of element
  template<class TypeEquation, class Vector1,
	   class Vector2, int p, class Vector4, class Vector5, class Dimension>
  void GetExtrapolationUacous(const EllipticProblem<TypeEquation>& var, int num_elem,
                              Vector1& Uloc, TinyVector<Vector2, p>& Uquad,
                              Vector4& extrapolU, const ElementReference<Dimension, 1>& Fb,
                              Vector5& Uneighbor, bool compute_uquad)
  {
    if (var.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
      return;
    
    // number of dofs inside the element
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int offset_face = mesh_num.OffsetQuadElementNumber(num_elem);
    
    typedef typename TypeEquation::Complexe Complexe;
    Vector<Complexe> face_U;
    int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();

    if (var.InsidePML(num_elem))
      {
        // split pml : we replace ux by ux+uy+uz
        if (Dimension::dim_N == 3)
          for (int i = 0; i < Uloc(0).GetM(); i++)
            Uloc(0)(i) += Uloc(1)(i) + Uloc(2)(i);
        else
          for (int i = 0; i < Uloc(0).GetM(); i++)
            Uloc(0)(i) += Uloc(1)(i);
      }
    
    if (Fb.UseQuadraturePointsForSh())
      if (compute_uquad)
        {
          Uquad(0).Reallocate(Fb.GetNbPointsQuadratureInside());
          Uquad(0).Fill(0);
          Fb.ApplyChTranspose(Uloc(0), Uquad(0));
          if (var.InsidePML(num_elem))
            {
              Uquad(1).Reallocate(Fb.GetNbPointsQuadratureInside());
              Uquad(1).Fill(0);
              if (Dimension::dim_N == 3)
                {
                  Uquad(2).Reallocate(Fb.GetNbPointsQuadratureInside());
                  Uquad(2).Fill(0);
                }
            }
        }
    
    // incrementing number of points of faces at the interface to other domains
    for (int num_loc = 0; num_loc < mesh.Element(num_elem).GetNbBoundary(); num_loc++)
      {
	int num_face = mesh.Element(num_elem).numBoundary(num_loc);
	int ref = mesh.Boundary(num_face).GetReference();
	int cond = mesh.GetBoundaryCondition(ref);
	int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        int rf = mesh_num.GetOrderQuadrature(num_face);
        
        bool affine = mesh.IsElementAffine(num_elem);
	face_U.Reallocate(nb_points_face); FillZero(face_U);        
        if (Fb.UseQuadraturePointsForSh())
          Fb.ApplyShQuadratureTranspose(num_loc, Uquad(0), face_U, rf);
        else
          Fb.ApplyShTranspose(num_loc, Uloc(0), face_U, rf);
        
        if ((!affine)&&(ElementReference_Base::use_warburton_trick))
          {
            for (int k = 0; k < nb_points_face; k++)
              {
                Real_wp invSqrtJacob = mesh.Glob_invSqrtJacobian(num_elem)(offset_war + k);
                face_U(k) *= invSqrtJacob;
              }                
          }
        
        if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            int offset_neighbor = var.GetOffsetNeighboringFace(num_face);
            for (int j = 0; j < nb_points_face; j++)
              Uneighbor(0)(offset_neighbor+j) = face_U(j);
          }
        
        // storing this extrapolation on vector extrapolU
        for (int j = 0; j < nb_points_face; j++)
          extrapolU(0)(offset_face + j) = face_U(j);
	    
        offset_face += nb_points_face;
	offset_war += nb_points_face;
      }
  }

  
  template<class TypeEquationTime,
           class Vector1, class Vector2, class Vector3, class Vector4, class Vector5>
  void GetExtrapolationAcoustic(HyperbolicProblem<TypeEquationTime>& var_time,
				const Vector1& U, const Vector1& V, int level,
				Vector2& extrapolU, Vector3& extrapolV,
				Vector2& Uneighbor, Vector3& Vneighbor,
				Vector4& Un_quad, Vector5& Vn_quad, bool treatU, bool treatV)
  {
    typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
    const EllipticProblem<TypeEquation>& vars
      = var_time.var_harmonic;
    
    if (vars.FormulationDG() == ElementReference_Base::CONTINUOUS)
      return;

    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
        
    int size_extrapol = mesh_num.OffsetQuadElementNumber(mesh.GetNbElt());
    int size_neighbor = vars.GetNbPointsQuadratureNeighbor();

    glob_chrono.Start(VirtualTimer::EXTRAPOL);
    
    if (treatU)
      if (Un_quad.GetM() != mesh.GetNbElt())
        Un_quad.Reallocate(mesh.GetNbElt());

    if (treatV)
      if (Vn_quad.GetM() != mesh.GetNbElt())
        Vn_quad.Reallocate(mesh.GetNbElt());
    
    if (treatU)
      {
	extrapolU(0).Reallocate(size_extrapol);
	extrapolU(0).Fill(0);
	
	Uneighbor(0).Reallocate(size_neighbor);
	Uneighbor(0).Fill(0);
      }
        
    if (treatV)
      for (int k = 0; k < Dimension::dim_N; k++)
	{
	  extrapolV(k).Reallocate(size_extrapol);
	  extrapolV(k).Fill(0);
	  
	  Vneighbor(k).Reallocate(size_neighbor);
	  Vneighbor(k).Fill(0);
	}

    TinyVector<VectReal_wp, Dimension::dim_N> Uloc;
    TinyVector<VectReal_wp, Dimension::dim_N> Vloc;
    int Nvol = mesh_num.GetNbDof();
    int Nvol_pml = mesh_num.GetNbDofPML();
    
    // first step : extrapolation step
    MatrixVectorProductLevel& level_time_scheme = var_time.GetTimeLevelDistribution();
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	// order of element
	int nb_dof_elt = mesh_num.GetNbLocalDof(num_elem);
	int offset = mesh_num.Element(num_elem).GetNumberDof(0);
	if (treatU)
	  {
	    Uloc(0).Reallocate(nb_dof_elt);
	    for (int j = 0; j < nb_dof_elt; j++)
	      Uloc(0)(j) = U(offset+j);
	    
	    if (vars.InsidePML(num_elem))
	      {
		Uloc(1).Reallocate(nb_dof_elt);
		for (int j = 0; j < nb_dof_elt; j++)
		  Uloc(1)(j) = U(Nvol + mesh_num.GetDofPML(offset+j));
		
		if (Dimension::dim_N == 3)
		  {
		    Uloc(2).Reallocate(nb_dof_elt);
		    for (int j = 0; j < nb_dof_elt; j++)
		      Uloc(2)(j) = U(Nvol + Nvol_pml + mesh_num.GetDofPML(offset+j));
		  } 
	      }
	  }
	
	if (treatV)
	  {
	    for (int m = 0; m < Dimension::dim_N; m++)
	      {
		Vloc(m).Reallocate(nb_dof_elt);
		for (int j = 0; j < nb_dof_elt; j++)
		  Vloc(m)(j) = V(offset+j);
		
		offset += Nvol;
	      }
	  }
	
	if (treatU)
	  GetExtrapolationUacous(vars, num_elem, Uloc, Un_quad(num_elem),
				 extrapolU, vars.GetReferenceElementH1(num_elem),
				 Uneighbor, true);
	
	if (treatV)
	  GetExtrapolationU(vars, num_elem, Vloc, Vn_quad(num_elem),
			    extrapolV, vars.GetReferenceElementH1(num_elem),
			    Vneighbor, true);
      }
    
    glob_chrono.Stop(VirtualTimer::EXTRAPOL);

#ifdef SELDON_WITH_MPI
    Vector<MPI_Request> requestU;
    Vector<Vector<MPI_Request> > requestV(Dimension::dim_N);
    Vector<Vector<Real_wp> > UXneighbor, UXsend;
    Vector<Vector<int64_t> > UXneighbor_tmp, UXsend_tmp;
    TinyVector<Vector<Vector<Real_wp> >, Dimension::dim_N> VXneighbor, VXsend;
    TinyVector<Vector<Vector<int64_t> >, Dimension::dim_N> VXneighbor_tmp, VXsend_tmp;
    
    // sending and receiving 
    if (treatU)
      vars.ExchangeUfaceDomains(Uneighbor(0), UXsend, UXsend_tmp,
                                UXneighbor, UXneighbor_tmp, requestU, 110);
    
    if (treatV)
      for (int k = 0; k < Dimension::dim_N; k++)
	vars.ExchangeUfaceDomains(Vneighbor(k), VXsend(k), VXsend_tmp(k),
                                  VXneighbor(k), VXneighbor_tmp(k), requestV(k), 111+k);
    
    // finalizing transfers
    if (treatU)
      vars.GetUfaceDomains(Uneighbor(0), UXsend, UXsend_tmp,
                           UXneighbor, UXneighbor_tmp, requestU, 110);
        
    if (treatV)
      for (int k = 0; k < Dimension::dim_N; k++)
	vars.GetUfaceDomains(Vneighbor(k), VXsend(k), VXsend_tmp(k),
                             VXneighbor(k), VXneighbor_tmp(k), requestV(k), 111+k);
    
#endif
    
  }
  
  
  //! matrix vector product with the stiffness matrix for acoustic equation and DG elements
  /*!
    \f$ C = C + alpha R_h B \f$
    \f$ V = R_h U  = -\int_K \nabla u \varphi - \int_{\partial K} [u] \varphi \f$
    where [u] = 1/2(u2 - u1)
    \param[in] alpha multiplication coefficient
    \param[in] level not used
    \param[in] var instationary problem associated
    \param[in] B vector to be multiplied
    \param[in] beta multiplication coefficient
    \param[out] C resulting vector 
   */
  template<class TypeEquation>
  void MltAddStiffnessVectorial_AcousticDG(const Real_wp& alpha, int level,
                                           HyperbolicProblem<TypeEquation>& var,
                                           const VectReal_wp& B,
                                           const Real_wp& beta, VectReal_wp& C)
  {
    if (beta == Real_wp(0))
      C.Fill(0);
    else
      Mlt(beta, C);
    
    // basic alias
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>&
      vars = var.var_harmonic;
    
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
	MltAdd_ElementRh1(vars, num_elem, B, var.evalUn_quad(num_elem),
			  var.evalVn_quad(num_elem),
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
	MltAdd_ElementSh(vars, num_elem, B, var.evalVn_quad(num_elem),
			 var.extrapolU, var.Uneighbor,
			 var.extrapolV, var.Vneighbor,
			 alpha, C, vars.GetReferenceElementH1(num_elem));
      }
    
    glob_chrono.Stop(VirtualTimer::FLUX);
  }
  
  
  //! matrix vector product with the stiffness matrix for acoustic equation and DG elements
  /*!
    \f$ C = C + alpha*R_h^* B \f$
    \f$ U = R_h^* V  = -\int_K v \nabla \varphi + \int_{\partial K} {u} \varphi \f$
    where {u} = 1/2(u1 + u2)
    \param[in] alpha multiplication coefficient
    \param[in] level not used
    \param[in] var instationary problem associated
    \param[in] B vector to be multiplied
    \param[in] beta multiplication coefficient
    \param[out] C resulting vector 
   */
  template<class TypeEquation>
  void MltAddStiffnessScalar_AcousticDG(const Real_wp& alpha, int level,
                                        HyperbolicProblem<TypeEquation>& var,
					const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C)
  {
    if (beta == Real_wp(0))
      C.Zero();
    else if (beta != Real_wp(1))
      Mlt(beta, C);
    
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var.var_harmonic;
    
    ////////////////////////////////
    // term -\int_K v \cdot \nabla \varphi
    ///////////////////////////////
    
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    MatrixVectorProductLevel& level_time_scheme = var.GetTimeLevelDistribution();
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	MltAdd_ElementRh2(vars, num_elem, B, var.evalUn_quad(num_elem),
			  var.evalVn_quad(num_elem),
			  alpha, C, vars.GetReferenceElementH1(num_elem));
      }
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);    
    
    // fluxes terms
    glob_chrono.Start(VirtualTimer::FLUX);    
    level_time_scheme.SetLevel(level);
    for (int num_elem_loc = 0; num_elem_loc < level_time_scheme.GetNbElt(); num_elem_loc++)
      {
	int num_elem = level_time_scheme.GetElementNumber(num_elem_loc);
	MltAdd_ElementSht(vars, num_elem, B, var.evalUn_quad(num_elem),
			  var.extrapolU, var.Uneighbor,
			  var.extrapolV, var.Vneighbor,
			  alpha, C, vars.GetReferenceElementH1(num_elem));
      }

    glob_chrono.Stop(VirtualTimer::FLUX);
  }


  template<class TypeEquationTime, class Vector1, class Vector2, class Vector3, class Vector4>
  void GetExtrapolationAcousticHDG(HyperbolicProblem<TypeEquationTime>& var, int level,
				   const VectReal_wp& U, const VectReal_wp& V,
				   Vector1& Un_quad, Vector2& Vn_quad,
				   Vector3& extrapolU, Vector4& extrapolV, VectReal_wp& Lambda)
  {
    typedef typename TypeEquationTime::TypeEquationStationary TypeEquation;
    typedef typename TypeEquation::Dimension Dimension;
    
    const EllipticProblem<TypeEquation>& vars = var.var_harmonic;
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    typedef typename Dimension::DimensionBoundary DimensionB;
    
    MatrixVectorProductLevel& list_level = var.GetTimeLevelDistribution();
    list_level.SetLevel(level);
    
    if (Un_quad.GetM() != mesh.GetNbElt())
      Un_quad.Reallocate(mesh.GetNbElt());
    
    if (Vn_quad.GetM() != mesh.GetNbElt())
      Vn_quad.Reallocate(mesh.GetNbElt());
    
    int nb_dof_L = mesh_num.GetNbDof();
    Lambda.Reallocate(nb_dof_L);

    if (level < 0)
      Lambda.Zero();
    
    int size_extrapol = mesh_num.OffsetQuadElementNumber(mesh.GetNbElt());

    extrapolU(0).Reallocate(size_extrapol);
    //extrapolU(0).Zero();
    
    for (int k = 0; k < Dimension::dim_N; k++)
      {
	extrapolV(k).Reallocate(size_extrapol);
	//extrapolV(k).Zero();
      }
    
    VectReal_wp Uloc, face_U;
    TinyVector<VectReal_wp, 1> contrib_all;
    TinyVector<VectReal_wp, Dimension::dim_N> Vloc, face_V;
    for (int i0 = 0; i0 < list_level.GetNbElt(); i0++)
      {
	int i = list_level.GetElementNumber(i0);
	const ElementReference<Dimension, 1>& Fb = vars.GetReferenceElementH1(i);
	// retrieving values of u and v
	int Nvol = vars.GetNbMainUnknownDof();
	int offset = vars.GetOffsetDofV(i);
	Uloc.Reallocate(Fb.GetNbDof());	    
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  Uloc(j) = U(offset + j);
	
	for (int k = 0; k < Dimension::dim_N; k++)
	  {
	    Vloc(k).Reallocate(Fb.GetNbDof());	    
	    for (int j = 0; j < Fb.GetNbDof(); j++)
	      Vloc(k)(j) = V(offset + j);
	    
	    offset += Nvol;
	  }
	
	// we compute values of u and v on quadrature points
	if (Fb.UseQuadraturePointsForSh())
	  {
	    Fb.ApplyChTranspose(Uloc, Un_quad(i)(0));
	    for (int k = 0; k < Dimension::dim_N; k++)
	      Fb.ApplyChTranspose(Vloc(k), Vn_quad(i)(k));
	  }
	
	//int ref_domain = mesh.Element(i).GetReference();
	Real_wp tau = -vars.alpha_penalization;
	
	// then loop over boundaries of the element
	contrib_all(0).Reallocate(mesh_num.Element(i).GetNbDof());
	offset = 0;
	int offset_face = mesh_num.OffsetQuadElementNumber(i);
	for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
	    int num_face = mesh.Element(i).numBoundary(num_loc);
            if (vars.upwind_fluxes)
              tau = -vars.alpha_penalization*vars.Glob_CoefPenalDG(num_face);

	    int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    
	    // values of u and v on quadrature points are computed
	    face_U.Reallocate(nb_points_face);
	    for (int k = 0; k < Dimension::dim_N; k++)
	      face_V(k).Reallocate(nb_points_face);
	    
	    if (Fb.UseQuadraturePointsForSh())
	      {
		Fb.ApplyShQuadratureTranspose(num_loc, Un_quad(i)(0), face_U, rf);
		for (int k = 0; k < Dimension::dim_N; k++)
		  Fb.ApplyShQuadratureTranspose(num_loc, Vn_quad(i)(k), face_V(k), rf);
	      }
	    else
	      {
		Fb.ApplyShTranspose(num_loc, Uloc, face_U, rf);
		for (int k = 0; k < Dimension::dim_N; k++)
		  Fb.ApplyShTranspose(num_loc, Vloc(k), face_V(k), rf);
	      }

	    // these values are stored in extrapolU and extrapolV
	    Real_wp dsj; typename Dimension::R_N normale, vec_v;
	    for (int j = 0; j < nb_points_face; j++)
	      {
		extrapolU(0)(offset_face+j) = face_U(j);
		CopyVector(face_V, j, vec_v);
		CopyVector(vec_v, offset_face+j, extrapolV);
	      }
	    
	    bool new_face = vars.IsNewFace(i)(num_loc);
	    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
	    
	    const ElementReference<DimensionB, 1>& Fb_s = vars.GetSurfaceFiniteElementH1(num_face);
	    VectReal_wp feval(nb_points_face), contrib(Fb_s.GetNbDof());
	    
	    int num_elem2 = vars.mesh.Boundary(num_face).numElement(0);
	    if ((num_elem2 == i)&&(vars.mesh.Boundary(num_face).GetNbElements()==2))
	      num_elem2 = vars.mesh.Boundary(num_face).numElement(1);
	    
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
	    
	    // forming - v \cdot n - tau u
	    for (int k = 0; k < nb_points_face; k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  {
		    normale = vars.Glob_normale(num_face)(k);
		    dsj = vars.Glob_dsj(num_face)(k);
		  }
		else
		  {
		    normale = vars.Glob_normale(num_face)(krot);
		    dsj = vars.Glob_dsj(num_face)(krot);
		    normale = -normale;
		  }
		
		CopyVector(face_V, k, vec_v);
		Real_wp poids = dsj*Fb_s.WeightsND(k);
		if (dirichlet_boundary)
		  feval(k) = - poids * tau * face_U(k);
		else
		  feval(k) = -poids*(-DotProd(vec_v, normale) + tau * face_U(k));
	      }
	    
	    // integral againt basis functions of lambda
	    Fb_s.ApplyCh(feval, contrib);
	    for (int j = 0; j < contrib.GetM(); j++)
	      contrib_all(0)(offset + j) = contrib(j);
	    
	    offset += Fb_s.GetNbDof();
	    offset_face += nb_points_face;
	  }
	
	// values are added to the result extrapolU
	mesh_num.number_map.AddLocalUnknownVector(mesh_num, Real_wp(1), contrib_all,
						  i, Lambda);
      }  
  }
  

  template<class TypeEquation, class Vector1, class Vector2, class Vector3, class Vector4>
  void MltAddStiffnessAcousticHDG(const Real_wp& alpha, int level, int level_vol,
				  HyperbolicProblem<TypeEquation>& var,
				  const VectReal_wp& U, const VectReal_wp& V, const VectReal_wp& Lambda,
				  const Vector1& Un_quad, const Vector2& Vn_quad,
				  const Vector3& extrapolU, const Vector4& extrapolV,
				  const Real_wp& beta, VectReal_wp& ProdU, VectReal_wp& ProdV)
  {
    typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
    if (beta == Real_wp(0))
      {
	ProdU.Zero();
	ProdV.Zero();
      }
    else if (beta != Real_wp(1))
      {
	Mlt(beta, ProdU);
	Mlt(beta, ProdV);
      }
    
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var.var_harmonic;
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    typedef typename Dimension::DimensionBoundary DimensionB;
    
    MatrixVectorProductLevel& list_level = var.GetTimeLevelDistribution();
    list_level.SetLevel(level);

    MatrixVectorProductLevel level_inside(list_level);
    level_inside.SetLevel(level_vol);
    
    TinyVector<VectReal_wp, 1> Lloc;
    VectReal_wp contrib, Lambda_quad, face_U, Uloc, dU;
    VectReal_wp ProdUloc, ProdUquad, feval, Uquad;
    TinyVector<VectReal_wp, Dimension::dim_N> fevalV, ProdVloc, ProdVquad, face_V;
    int Nvol = vars.GetNbMainUnknownDof();
    for (int i0 = 0; i0 < list_level.GetNbElt(); i0++)
      {
	int i = list_level.GetElementNumber(i0);
	const ElementReference<Dimension, 1>& Fb = vars.GetReferenceElementH1(i);
	ProdUloc.Reallocate(Fb.GetNbDof()); ProdUloc.Zero();
	ProdUquad.Reallocate(Fb.GetNbPointsQuadratureInside());
	ProdUquad.Zero();
	for (int d = 0; d < Dimension::dim_N; d++)
	  {
	    ProdVloc(d).Reallocate(Fb.GetNbDof());
	    ProdVquad(d).Reallocate(Fb.GetNbPointsQuadratureInside());
	    ProdVloc(d).Zero(); ProdVquad(d).Zero();
	  }
	
	int offset_u = vars.GetOffsetDofV(i);

	Real_wp tau = -vars.alpha_penalization;

	// we retrieve values of lambda on the element
	mesh_num.number_map.GetLocalUnknownVector(mesh_num, Lambda, i, Lloc);

	bool treat_inside = false;
	if (level_inside.TreatElement(i))
	  treat_inside = true;

	// loop over boundaries of the element
	int offset = 0;
	int offset_face = mesh_num.OffsetQuadElementNumber(i);
	for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
	    int num_face = mesh.Element(i).numBoundary(num_loc);
            if (vars.upwind_fluxes)
              tau = -vars.alpha_penalization*vars.Glob_CoefPenalDG(num_face);

	    int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    const ElementReference<DimensionB, 1>& Fb_s = vars.GetSurfaceFiniteElementH1(num_face);
	    
	    contrib.Reallocate(Fb_s.GetNbDof());
	    Lambda_quad.Reallocate(nb_points_face);
	    
	    // values of Lambda on quadrature points are computed
	    for (int j = 0; j < contrib.GetM(); j++)
	      contrib(j) = Lloc(0)(offset + j);
	    
	    Fb_s.ApplyChTranspose(contrib, Lambda_quad);

	    Real_wp dsj; typename Dimension::R_N normale, vec_v;
	    
	    // values of u and v on quadrature points are retrieved
	    if (treat_inside)
	      {
		face_U.Reallocate(nb_points_face);
		for (int k = 0; k < Dimension::dim_N; k++)
		  face_V(k).Reallocate(nb_points_face);
		
		for (int j = 0; j < nb_points_face; j++)
		  {
		    face_U(j) = extrapolU(0)(offset_face + j);
		    CopyVector(extrapolV, offset_face+j, vec_v);
		    CopyVector(vec_v, j, face_V);
		  }
	      }
	    
	    bool new_face = vars.IsNewFace(i)(num_loc);
	    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
	    
	    int num_elem2 = vars.mesh.Boundary(num_face).numElement(0);
	    if ((num_elem2 == i)&&(vars.mesh.Boundary(num_face).GetNbElements()==2))
	      num_elem2 = vars.mesh.Boundary(num_face).numElement(1);
	    
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
	    
	    for (int k = 0; k < nb_points_face; k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  {
		    normale = vars.Glob_normale(num_face)(k);
		    dsj = vars.Glob_dsj(num_face)(k);
		  }
		else
		  {
		    normale = vars.Glob_normale(num_face)(krot);
		    dsj = vars.Glob_dsj(num_face)(krot);
		    normale = -normale;
		  }

		Real_wp poids = dsj*Fb_s.WeightsND(k), diff_u(0);
		if (treat_inside)
		  {
		    CopyVector(face_V, k, vec_v);
		    feval(k) = poids* (tau * (Lambda_quad(k)-face_U(k)) + DotProd(vec_v, normale));
		    diff_u = -face_U(k);
		  }
		else
		  feval(k) = poids*tau*Lambda_quad(k);
		
		if (!dirichlet_boundary)
		  diff_u += Lambda_quad(k);
		
		vec_v = poids*normale*diff_u;
		CopyVector(vec_v, k, fevalV);
	      }
	    
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
	    offset_face += nb_points_face;
	  }

	// now we treat volume integrals
	// part for ProdVloc
	int nb_points_elt = Fb.GetNbDof();
	int nb_points_quad = Fb.GetNbPointsQuadratureInside();
	bool affine = mesh.IsElementAffine(i);
	typename Dimension::R_N vec_v, grad_U, vec_u;
	if ((affine)&&(!Fb.UseQuadraturePointsForRh()))
	  nb_points_quad = nb_points_elt;

	if (treat_inside)
	  {
	    dU.Reallocate(nb_points_quad*Dimension::dim_N);
	    contrib.Reallocate(nb_points_elt);
	    feval.Reallocate(nb_points_quad);
	    for (int d = 0; d < Dimension::dim_N; d++)
	      fevalV(d).Reallocate(nb_points_quad);
	    
	    // gradient of u on reference element
	    if (Fb.UseQuadraturePointsForSh())
	      Fb.ApplyRhQuadratureTranspose(Un_quad(0)(i), dU);
	    else
	      {
		Uloc.Reallocate(nb_points_elt);
		for (int j = 0; j < nb_points_elt; j++)
		  Uloc(j) = U(offset_u + j);
		
		if (Fb.UseQuadraturePointsForRh())
		  {
		    Uquad.Reallocate(nb_points_quad);
		    Fb.ApplyChTranspose(Uloc, Uquad);
		    Fb.ApplyRhQuadratureTranspose(Uquad, dU);
		  }
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
		  dfjm1 = vars.Glob_DFjm1(i)(0);
		
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    CopyVector(dU, k, vec_v);		
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
		  MltTrans(vars.Glob_DFjm1(i)(k), vec_v, grad_U);	  
		  grad_U *= Fb.WeightsND(k);
		  CopyVector(grad_U, k, fevalV);
		}
	    
	    // integration agains basis functions
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
	    
	    // computation of v on quadrature points
	    if (Fb.UseQuadraturePointsForSh())
	      {
		for (int n = 0; n < Dimension::dim_N; n++)
		  {
		    fevalV(n) = Vn_quad(n)(i);
		    if (affine)
		      for (int k = 0; k < nb_points_quad; k++)
			{
			  Real_wp coef = Fb.WeightsND(k);
			  fevalV(n)(k) *= coef;
			}
		  }
	      }
	    else
	      {
		for (int n = 0; n < Dimension::dim_N; n++)
		  {
		    for (int j = 0; j < nb_points_elt; j++)
		      Uloc(j) = V(offset_u + j + n*Nvol);
		    
		    if (Fb.UseQuadraturePointsForRh() || (!affine))
		      Fb.ApplyChTranspose(Uloc, fevalV(n));
		    else
		      fevalV(n) = Uloc;
		  }
	      }
	    
	    // application of transformation DF_i^{-1}	
	    if (affine)
	      {
		const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N>&
		  dfjm1 = vars.Glob_DFjm1(i)(0);
		
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    CopyVector(fevalV, k, vec_u);
		    
		    Mlt(dfjm1, vec_u, vec_v);
		    vec_v = -vec_v;
		    if (Fb.UseQuadraturePointsForRh())
		      vec_v *= Fb.WeightsND(k);
		    
		    CopyVector(vec_v, k, dU);
		  }
	      }
	    else
	      {
		for (int k = 0; k < nb_points_quad; k++)
		  {
		    CopyVector(fevalV, k, vec_u);
		    
		    Mlt(vars.Glob_DFjm1(i)(k), vec_u, vec_v);
		    vec_v *= -Fb.WeightsND(k);
		    CopyVector(vec_v, k, dU);
		  }
	      }
	    
	    // integration against grad phi
	    if (Fb.UseQuadraturePointsForRh())
	      {
		Fb.ApplyRhQuadrature(dU, Uquad);	    
		if (Fb.UseQuadraturePointsForSh())
		  {
		    ProdUquad += Uquad;
		    Fb.ApplyCh(ProdUquad, ProdUloc);
		  }
		else
		  {
		    Fb.ApplyCh(Uquad, Uloc);
		    ProdUloc += Uloc;
		  }
	      }
	    else
	      {	
		if (affine)
		  Fb.ApplyConstantRh(dU, Uloc);
		else
		  Fb.ApplyRh(dU, Uloc);
		
		ProdUloc += Uloc;
	      }
	  }
	else
	  {
	    if (Fb.UseQuadraturePointsForSh())
	      {
		for (int d = 0; d < Dimension::dim_N; d++)
		  Fb.ApplyCh(ProdVquad(d), ProdVloc(d));
		
		Fb.ApplyCh(ProdUquad, ProdUloc);
	      }
	  }
	
	// contributions are added to ProdU and ProdV
	for (int j = 0; j < ProdUloc.GetM(); j++)
	  ProdU(offset_u + j) += alpha*ProdUloc(j);

	for (int d = 0; d < Dimension::dim_N; d++)
	  {
	    for (int j = 0; j < ProdVloc(d).GetM(); j++)
	      ProdV(offset_u + j) += alpha*ProdVloc(d)(j);

	    offset_u += Nvol;
	  }
      }
  }

  
  template<class TypeEquation, class T>
  void MltAddLambdaAcousticHDG(const T& alpha, int level,
			       HyperbolicProblem<TypeEquation>& var,
			       const Vector<T>& Lambda, Vector<T>& ProdUV)
  {
    typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
    const EllipticProblem<typename TypeEquation::TypeEquationStationary>& vars = var.var_harmonic;
    const Mesh<Dimension>& mesh = vars.mesh;
    const MeshNumbering<Dimension>& mesh_num = vars.GetMeshNumbering(0);
    typedef typename Dimension::DimensionBoundary DimensionB;

    MatrixVectorProductLevel& list_level = var.GetTimeLevelDistribution();
    list_level.SetLevel(level);
    
    TinyVector<Vector<T>, 1> Lloc;
    Vector<T> contrib, feval, ProdUloc, ProdUquad, Lambda_quad;
    TinyVector<Vector<T>, Dimension::dim_N> fevalV, ProdVloc, ProdVquad;
    int Nvol = vars.GetNbMainUnknownDof();
    typename Dimension::R_N normale;
    TinyVector<T, Dimension::dim_N> vec_v; Real_wp dsj;
    for (int i0 = 0; i0 < list_level.GetNbElt(); i0++)
      {
	int i = list_level.GetElementNumber(i0);
	const ElementReference<Dimension, 1>& Fb = vars.GetReferenceElementH1(i);
	ProdUloc.Reallocate(Fb.GetNbDof()); ProdUloc.Zero();
	if (Fb.UseQuadraturePointsForSh())
	  ProdUquad.Reallocate(Fb.GetNbPointsQuadratureInside());

	ProdUquad.Zero();
	for (int d = 0; d < Dimension::dim_N; d++)
	  {
	    ProdVloc(d).Reallocate(Fb.GetNbDof());
	    if (Fb.UseQuadraturePointsForSh())
	      ProdVquad(d).Reallocate(Fb.GetNbPointsQuadratureInside());
	    
	    ProdVloc(d).Zero(); ProdVquad(d).Zero();
	  }
	
	int offset_u = vars.GetOffsetDofV(i);
	Real_wp tau = -vars.alpha_penalization;
        
	// we retrieve values of lambda on the element
	mesh_num.number_map.GetLocalUnknownVector(mesh_num, Lambda, i, Lloc);
	
	// loop over boundaries of the element
	int offset = 0;
	int offset_face = mesh_num.OffsetQuadElementNumber(i);
	for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
	    int num_face = mesh.Element(i).numBoundary(num_loc);
            if (vars.upwind_fluxes)
              tau = -vars.alpha_penalization*vars.Glob_CoefPenalDG(num_face);

	    int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    const ElementReference<DimensionB, 1>& Fb_s = vars.GetSurfaceFiniteElementH1(num_face);
	    
	    contrib.Reallocate(Fb_s.GetNbDof());
	    Lambda_quad.Reallocate(nb_points_face);
	    
	    // values of Lambda on quadrature points are computed
	    for (int j = 0; j < contrib.GetM(); j++)
	      contrib(j) = Lloc(0)(offset + j);

	    Fb_s.ApplyChTranspose(contrib, Lambda_quad);

	    bool new_face = vars.IsNewFace(i)(num_loc);
	    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
	    
	    int num_elem2 = vars.mesh.Boundary(num_face).numElement(0);
	    if ((num_elem2 == i)&&(vars.mesh.Boundary(num_face).GetNbElements()==2))
	      num_elem2 = vars.mesh.Boundary(num_face).numElement(1);
	    
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
	    
	    for (int k = 0; k < nb_points_face; k++)
	      {
		int krot = FacesQuadRotation(rot, k);
		
		if (new_face)
		  {
		    normale = vars.Glob_normale(num_face)(k);
		    dsj = vars.Glob_dsj(num_face)(k);
		  }
		else
		  {
		    normale = vars.Glob_normale(num_face)(krot);
		    dsj = vars.Glob_dsj(num_face)(krot);
		    normale = -normale;
		  }

		Real_wp poids = dsj*Fb_s.WeightsND(k);
		feval(k) = poids* tau * Lambda_quad(k);
		T diff_u(0);
		if (!dirichlet_boundary)
		  diff_u += Lambda_quad(k);
		
		vec_v = poids*normale*diff_u;
		CopyVector(vec_v, k, fevalV);
	      }
	    
	    if (Fb.UseQuadraturePointsForSh())
	      {
		Fb.ApplyShQuadrature(T(1), num_loc, feval, ProdUquad, rf);
		for (int d = 0; d < Dimension::dim_N; d++)
		  Fb.ApplyShQuadrature(T(1), num_loc, fevalV(d), ProdVquad(d), rf);
	      }
	    else
	      {
		Fb.ApplySh(T(1), num_loc, feval, ProdUloc, rf);
		for (int d = 0; d < Dimension::dim_N; d++)
		  Fb.ApplySh(T(1), num_loc, fevalV(d), ProdVloc(d), rf);
	      }
	    
	    offset += Fb_s.GetNbDof();
	    offset_face += nb_points_face;
	  }

	if (Fb.UseQuadraturePointsForSh())
	  {
	    Fb.ApplyCh(ProdUquad, ProdUloc);
	    for (int d = 0; d < Dimension::dim_N; d++)
	      Fb.ApplyCh(ProdVquad(d), ProdVloc(d));
	  }

	// contributions are added to ProdU and ProdV
	for (int j = 0; j < ProdUloc.GetM(); j++)
	  ProdUV(offset_u + j) += alpha*ProdUloc(j);

	offset_u += Nvol;
	for (int d = 0; d < Dimension::dim_N; d++)
	  {
	    for (int j = 0; j < ProdVloc(d).GetM(); j++)
	      ProdUV(offset_u + j) += alpha*ProdVloc(d)(j);

	    offset_u += Nvol;
	  }
      }

  }
  
}

#define MONTJOIE_FILE_PROD_MAT_VECT_ACOUSTIC_CXX
#endif



