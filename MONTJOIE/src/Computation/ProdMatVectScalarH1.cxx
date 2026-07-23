#ifndef MONTJOIE_FILE_PROD_MAT_VECT_SCALAR_H1_CXX

#ifdef MONTJOIE_WITH_META_PROGRAMMING
//#include "OptProdMatVectNodalH1.cxx"
#endif

namespace Montjoie
{
  
  //! computation of unknowns on the quadrature points of element
  template<class Complexe, class TypeEquation, class T1>
  void GetExtrapolationSecondOrder_Scalar(const EllipticProblem<TypeEquation>& var,
                                          int num_elem, const Vector<T1>& B,
					  ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra,
                                          const GlobalGenericMatrix<Complexe>& nat_mat,
					  const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int offset_face = mesh_num.OffsetQuadElementNumber(num_elem);
    bool affine = mesh.IsElementAffine(num_elem);    
    
    // values of u and grad u on quadrature points of face
    Vector<T1> gradU, valU;
    
    T1 du_dn;
    typename Dimension::R_N normale;
    typename Dimension::MatrixN_N dfjm1;
    TinyVector<T1, Dimension::dim_N> grad_v, vec_v;
    TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;

    TinyVector<Vector<T1>, 1>& Uloc = var_extra.Uloc;
    TinyVector<Vector<T1>, 1>& Uquad = var_extra.GetUnQuad(num_elem, Fb);
    
    // we get values of u on the element
    var.GetLocalUnknownVector(B, num_elem, var_extra.Uloc);
    
    int num_point = Fb.GetNbPointsQuadratureInside();  
    if (Fb.UseQuadraturePointsForSh())
      {
        // computation of u on quadrature points if needed
        Uquad(0).Reallocate(num_point);
        Fb.ApplyChTranspose(Uloc(0), Uquad(0));
      }
    
    int ref_d = mesh.Element(num_elem).GetReference();
    int num_loc2, num_elem2, rot = 0;
    int ref, cond, offset_face2; Complexe phase, phase_conj; bool new_face;
    // loop over boundaries of the element
    for (int num_loc = 0; num_loc < mesh.Element(num_elem).GetNbBoundary(); num_loc++)
      {
	int num_face = mesh.Element(num_elem).numBoundary(num_loc);
        int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        if (var.FaceHasToBeConsideredForBoundaryIntegral(num_face))
          {            
            // the face is involved in the fluxes terms
            int rf = mesh_num.GetOrderQuadrature(num_face);        
            var_extra.GetFaceOrientationAndNumber(var, num_face,
                                                  num_elem, num_loc, offset_face,
                                                  ref, cond, new_face, num_elem2, num_loc2,
                                                  offset_face2, rot, phase, phase_conj);
            
            const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
              GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
            
            // computation of u on quadrature points of face
            valU.Reallocate(nb_points_face);            
            if (Fb.UseQuadraturePointsForSh())
              Fb.ApplyShQuadratureTranspose(num_loc, Uquad(0), valU, rf);
            else
              Fb.ApplyShTranspose(num_loc, Uloc(0), valU, rf);
            
            // computation of gradient of u on quadrature points of face
            gradU.Reallocate(Dimension::dim_N*nb_points_face);
            if (Fb.UseQuadraturePointsForSh())
              Fb.ApplyNablaShQuadratureTranspose(num_loc, Uquad(0), gradU, rf);
            else
              Fb.ApplyNablaShTranspose(num_loc, Uloc(0), gradU, rf);

	    
            int offset_neighbor = 0;
            bool neighbor_face = false;
	    if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
	      {
                neighbor_face = true;
                offset_neighbor = var.GetOffsetNeighboringFace(num_face);
              }
            
            // loop over quadrature points
            for (int j = 0; j < nb_points_face; j++)
              {
                int jrot = FacesQuadRotation(rot, j);
                
                if (new_face)
                  normale = var.Glob_normale(num_face)(j);
                else
                  {
                    normale = var.Glob_normale(num_face)(jrot);
                    Mlt(Real_wp(-1), normale);
                  }
                
                // storing u
                var_extra.extrapolU(0)(offset_face+j) = valU(j);
                
                // storing mu du/dn
                if (affine)
                  {
                    dfjm1 = var.Glob_DFjm1(num_elem)(0);
                    Mlt(1.0/var.Glob_jacobian(num_elem)(0), dfjm1);
                  }
                else
                  dfjm1 = var.Glob_DFjm1(num_elem)(num_point+j);
                
                CopyVector(gradU, j, grad_v);
                MltTrans(dfjm1, grad_v, vec_v);
                
                mu = var.GetMu_QuadraturePoint(ref_d, num_elem, num_point+j);
                Mlt(mu, vec_v, grad_v);
                
                du_dn = DotProd(normale, grad_v);
                
                var_extra.extrapolCgradU(0)(offset_face+j) = du_dn; 
                if (neighbor_face)
                  {
                    var_extra.Uneighbor(0)(offset_neighbor+j) = valU(j);
                    var_extra.CgradUNeighbor(0)(offset_neighbor+j) = du_dn;
                  }
              }                      
          } 
        
        offset_face += nb_points_face;
        num_point += nb_points_face;
      }
  }
  
  
  //! computation of du/dn for SIPG formulation of Helmholtz equation
  //! values of du/dn on quadrature points are exchanged between processors
  template<class Complexe, class TypeEquation, class T1>
  void GetExtrapolationAndExchange_Scalar(const GlobalGenericMatrix<Complexe>& nat_mat,
                                          const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                          const Vector<T1>& B,
					  ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra)
  {
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = A.var.mesh;
    const MeshNumbering<Dimension>& mesh_num = A.var.GetMeshNumbering(0);
    
    int size_extrapol = 0;
    int size_neighbor = 0;

    if (A.var.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      //if (true)
      {        
        // for DG formulation, we will store U on quadrature points of each boundary
        // and exchange U on interface (Uneighbor)
        size_extrapol = mesh_num.OffsetQuadElementNumber(mesh.GetNbElt());
        size_neighbor = A.var.GetNbPointsQuadratureNeighbor();
        
        var_extra.extrapolU(0).Reallocate(size_extrapol);
        var_extra.Uneighbor(0).Reallocate(size_neighbor);
        
        var_extra.extrapolCgradU(0).Reallocate(size_extrapol);
        var_extra.CgradUNeighbor(0).Reallocate(size_neighbor);
      }
    
    // loop over all elements of the mesh to compute u on quadrature points of faces
    for (int i = 0; i < mesh.GetNbElt(); i++)
      GetExtrapolationSecondOrder_Scalar(A.var, i, B, var_extra, nat_mat,
					 A.var.GetReferenceElementH1(i));
    
    
#ifdef SELDON_WITH_MPI    
    // exchanging datas between processors
    if (A.var.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      {
        var_extra.ReallocateExchangeVector(2);
        
        if (size_neighbor > 0)
          {
            // sending and receiving u
            A.var.ExchangeUfaceDomains(var_extra.Uneighbor(0), var_extra.Xsend(0), var_extra.Xsend_tmp(0),
                                       var_extra.Xneighbor(0), var_extra.Xneighbor_tmp(0), var_extra.request(0), 90);
            
            // sending and receiving mu du/dn
            A.var.ExchangeUfaceDomains(var_extra.CgradUNeighbor(0), var_extra.Xsend(1), var_extra.Xsend_tmp(1),
                                       var_extra.Xneighbor(1), var_extra.Xneighbor_tmp(1), var_extra.request(1), 91);
          }
      }
#endif
    
  }
  
  
  //! exchange of values of u and du/dn are completed and numerical fluxes are finalized
  template<class Complexe, class TypeEquation, class T0, class T1>
  void CompleteExchangeAndFlux_Scalar(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                                      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                                      Vector<T1>& C, ExtrapolVariablesProductFEM<T1,
                                      TypeEquation>& var_extra)
  {
    
#ifdef SELDON_WITH_MPI
    if (A.var.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      //if (true)
      {
        int size_neighbor = A.var.GetNbPointsQuadratureNeighbor();
        if (size_neighbor > 0)
          {
            // finalizing transfers if needed
            A.var.GetUfaceDomains(var_extra.Uneighbor(0), var_extra.Xsend(0), var_extra.Xsend_tmp(0),
                                  var_extra.Xneighbor(0), var_extra.Xneighbor_tmp(0), var_extra.request(0), 90);
            
            A.var.GetUfaceDomains(var_extra.CgradUNeighbor(0), var_extra.Xsend(1), var_extra.Xsend_tmp(1),
                                  var_extra.Xneighbor(1), var_extra.Xneighbor_tmp(1), var_extra.request(1), 91);
          }
      }
#endif
    
    // 
    
  }
  
#ifdef MONTJOIE_WITH_META_PROGRAMMING
  //! to do
  template<class T0, class TypeEquation,
	   class Complexe, class Vector1, int r>
  void MltAdd_ElementH1_Scalar(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			       const SeldonTranspose& trans,
			       const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			       int i, const Vector1& B, Vector1& C, const ElementReference<Dimension, 1>& Fb,
			       const TinyVector<int, r>& rvec,
			       const ElementReference<Dimension, 1>& Fbr, GhostIf<false>& optim_implemented)
  {
    // order for which optimization not implemented, we use normal function
  }

  //! to do
  template<class T0, class TypeEquation,
	   class Complexe, class Vector1, int r>
  void MltAdd_ElementH1_Scalar(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			       const SeldonTranspose& trans,
			       const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			       int i, const Vector1& B, Vector1& C, const ElementReference<Dimension, 1>& Fb,
			       const TinyVector<int, r>& rvec,
			       const ElementReference<Dimension, 1>& Fbr, GhostIf<true>& optim_implemented)
  {
    // order for which optimization is implemented
    MltAdd_Optimized(alpha, nat_mat, trans, A, i, B, C, Fb, Fbr);
  }
#endif


  //! matrix-vector product for an element and for scalar equation (Helmholtz equation)
  template<class T0, class TypeEquation, class Complexe, class T1>
  void MltAdd_ElementH1_Scalar(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			       const SeldonTranspose& trans,
			       const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			       int i, const Vector<T1>& B, Vector<T1>& C,
                               ExtrapolVariablesProductFEM<T1, TypeEquation>& var_extra,
			       const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
#ifdef MONTJOIE_WITH_META_PROGRAMMING
    // if the required order is implemented in an optimized version ...
    /*
      int r = Fb.GetOrder();
    switch (r)
      {
      case 1 :
	{
	  TinyVector<int, 1> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r1, optim);
	  return;
	}
      case 2 :
	{
	  TinyVector<int, 2> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r2, optim);
	  return;
	}
      case 3 :
	{
	  TinyVector<int, 3> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r3, optim);
	  return;
	}
      case 4 :
	{
	  TinyVector<int, 4> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r4, optim);
	  return;
	}
      case 5 :
	{
	  TinyVector<int, 5> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r5, optim);
	  return;
	}
      case 6 :
	{
	  TinyVector<int, 6> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r6, optim);
	  return;
	}
      case 7 :
	{
	  TinyVector<int, 7> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r7, optim);
	  return;
	}
      case 8 :
	{
	  TinyVector<int, 8> rvec;
	  MltAdd_ElementH1_Scalar(alpha, nat_mat, A, i, B, C, Fb, rvec, Fb.elt_r8, optim);
	  return;
	}
      }
    */
#endif
    
    // otherwise, we use the non-optimized version
    typedef typename TypeEquation::Dimension Dimension;
    const EllipticProblem<TypeEquation>& var = A.var;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    
    TinyVector<T1, Dimension::dim_N> tmp, vh_loc;
    int dim_N = tmp.GetM();
    bool variable = var.UseNumericalIntegration(i);
    bool affine = mesh.IsElementAffine(i);
    
    // for Helmholtz and Laplace equation, mass coefficient is related to \int phi_i phi_j
    // while stiffness coefficient is related to \int grad phi_i grad phi_j
    bool presence_mass = (nat_mat.GetCoefMass()!= Complexe(0));
    bool presence_stiff = (nat_mat.GetCoefStiffness() != Complexe(0));
    bool presence_damping = (nat_mat.GetCoefDamping() != Complexe(0));
    
    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad_elt = Fb.GetNbPointsQuadratureInside();
    int ref = mesh.Element(i).GetReference();
        
    TinyVector<Vector<T1>, 1>& Uh_dof = var_extra.Uloc;
    TinyVector<Vector<T1>, 1>& Uh_quad = var_extra.GetUnQuad(i, Fb);

    // we get values of u on the element
    var.GetLocalUnknownVector(B, i, Uh_dof);
    
    // interpolation on quadrature points (if necessary)
    TinyVector<Vector<T1>, 1> Prod_Uh;
    Vector<T1> Vh_quad, ProdUh_quad;
    ProdUh_quad.Reallocate(nb_quad_elt);
    Vh_quad.Reallocate(dim_N*nb_quad_elt);
    Prod_Uh(0).Reallocate(nb_dof_elt);
    Prod_Uh(0).Zero(); ProdUh_quad.Zero(); Vh_quad.Zero();
    T1 cone; SetComplexOne(cone);

    // u is evaluated on quadrature points
    if (!Fb.UseQuadraturePointsForSh())
      {
        //if ((presence_mass) || (Fb.UseQuadraturePointsForRh()))
          {
            Uh_quad(0).Reallocate(nb_quad_elt);
            Fb.ApplyChTranspose(Uh_dof(0), Uh_quad(0));
          }
      }

    if (presence_stiff || presence_damping)
      {
        // gradient of u
        if (Fb.UseQuadraturePointsForRh())
          Fb.ApplyRhQuadratureTranspose(Uh_quad(0), Vh_quad);
        else
          Fb.ApplyRhTranspose(Uh_dof(0), Vh_quad);
      }

    if (presence_mass || presence_damping)
      {
        if (variable)
          {
	    if (trans.NoTrans())
	      for (int j = 0; j < nb_quad_elt; j++)
		{
		  ProdUh_quad(j) = Uh_quad(0)(j)
		    * TypeEquation::GetMassCoefficient(A.var, i, j, nat_mat, ref);
		  
		  if (presence_damping)
                    {
                      CopyVector(Vh_quad, j, tmp);
                      ProdUh_quad(j) += DotProd(TypeEquation::GetNablaU_Coefficient(A.var, i, j, nat_mat, ref), tmp);					  
                    }
		}
	    else if (trans.Trans())
	      for (int j = 0; j < nb_quad_elt; j++)
		{
		  ProdUh_quad(j) = Uh_quad(0)(j)
		    * TypeEquation::GetMassCoefficient(A.var, i, j, nat_mat, ref);
		  
		  if (presence_damping)
                    {
                      CopyVector(Vh_quad, j, tmp);
                      ProdUh_quad(j) += DotProd(TypeEquation::GetNablaPhi_Coefficient(A.var, i, j, nat_mat, ref), tmp);
                    }					  
		}	     
	    else
	      {
		cout << "Not implemented" << endl;
		abort();
	      }
          }
        else
          {
	    Complexe mass = TypeEquation::GetMassCoefficient(A.var, i, 0, nat_mat, ref);
	    TinyVector<Complexe, Dimension::dim_N>
	      grad = TypeEquation::GetNablaU_Coefficient(A.var, i, 0, nat_mat, ref);
	    
            if (Fb.UseQuadraturePointsForRh())
              {
                if (presence_damping)
                  for (int j = 0; j < nb_quad_elt; j++)
                    {
                      CopyVector(Vh_quad, j, tmp);
                      ProdUh_quad(j) = Fb.WeightsND(j)*
                        (mass*Uh_quad(0)(j) + DotProd(grad, tmp));
                    }
                else
                  for (int j = 0; j < nb_quad_elt; j++)
                    ProdUh_quad(j) = Fb.WeightsND(j)*mass*Uh_quad(0)(j);
              }
            else
              {
		if (!grad.IsZero())
		  {
		    cout << "Not implemented" << endl;
		    abort();
		  }
		
                Copy(Uh_dof(0), Prod_Uh(0));
                Fb.MltMassMatrix(Prod_Uh(0));
                Mlt(mass, Prod_Uh(0));                
              }
          }
      }
    
    if (presence_stiff)
      {
        if (!variable)
          {
            TinyMatrix<Complexe, Symmetric,
		       Dimension::dim_N, Dimension::dim_N> matC_dfj
              = TypeEquation::GetStiffCoefficient(A.var, i, 0, nat_mat, ref);

	    TinyVector<Complexe, Dimension::dim_N>
	      grad = TypeEquation::GetNablaPhi_Coefficient(A.var, i, 0, nat_mat, ref);
	    
	    if (trans.Trans())
	      grad = TypeEquation::GetNablaU_Coefficient(A.var, i, 0, nat_mat, ref);	      

	    // we apply geometrical transformation J_i DF_i^{-1} \mu DF_i^{*-1}
            for (int j = 0; j < nb_quad_elt; j++)
              { 
                CopyVector(Vh_quad, j, tmp);
                
                Mlt(matC_dfj, tmp, vh_loc);
		Add(Uh_quad(0)(j), grad, vh_loc);
		
                Mlt(Fb.WeightsND(j), vh_loc);
                
                CopyVector(vh_loc, j, Vh_quad);
              }
          }
        else
          {            
            // we apply geometrical transformation J_i DF_i^{-1} \mu DF_i^{*-1}
            if (trans.NoTrans())
	      for (int j = 0; j < nb_quad_elt; j++)
		{ 
		  CopyVector(Vh_quad, j, tmp);
		  
		  Mlt(TypeEquation::GetStiffCoefficient(A.var, i, j, nat_mat, ref), tmp, vh_loc);
		  Add(Uh_quad(0)(j), TypeEquation::GetNablaPhi_Coefficient(A.var, i, j, nat_mat, ref), vh_loc);
		  
		  CopyVector(vh_loc, j, Vh_quad);
		}
	    else
	      for (int j = 0; j < nb_quad_elt; j++)
		{ 
		  CopyVector(Vh_quad, j, tmp);
		  
		  Mlt(TypeEquation::GetStiffCoefficient(A.var, i, j, nat_mat, ref), tmp, vh_loc);
		  Add(Uh_quad(0)(j), TypeEquation::GetNablaU_Coefficient(A.var, i, j, nat_mat, ref), vh_loc);
		  
		  CopyVector(vh_loc, j, Vh_quad);
		}
          }

        // integration againt grad(phi)
        if (Fb.UseQuadraturePointsForRh())
          {
            Fb.ApplyRhQuadrature(Vh_quad, Uh_quad(0));
	    ProdUh_quad += Uh_quad(0);
	  }
        else
          {
            Fb.ApplyRh(Vh_quad, Uh_dof(0));
	    Prod_Uh(0) += Uh_dof(0);
	  }

        // then flux terms
        Vector<T1> V, dV;
        int offset_face1 = mesh_num.OffsetQuadElementNumber(i);
        int num_elem2, ref_line, cond, num_pos2_face, offset_face2, rot;
        bool new_face; Complexe phase, phase_conj;
        typename Dimension::R_N normale; Real_wp dsj;
        typename Dimension::MatrixN_N dfjm1;
        TinyMatrix<Complexe, Symmetric, Dimension::dim_N, Dimension::dim_N> mu;
        TinyVector<T1, Dimension::dim_N> vec_u, vec_v;
        T1 u1, u2, du1_dn, du2_dn, jump_u, phi_loc;
        int offset_quad = Fb.GetNbPointsQuadratureInside();
        for (int num_pos1_face = 0; num_pos1_face < Fb.GetNbBoundaries();
             num_pos1_face++)
          {
            int num_face = mesh.Element(i).numBoundary(num_pos1_face);
            int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
            if (var.FaceHasToBeConsideredForBoundaryIntegral(num_face))
              {                
                int rf = mesh_num.GetOrderQuadrature(num_face);
                const VectReal_wp& PoidsFlux = mesh_num.number_map
		  .GetFluxWeight(rf, mesh.Boundary(num_face));
                
		const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
                  GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
                
                V.Reallocate(nb_points_face);
                dV.Reallocate(nb_points_face*Dimension::dim_N);
                
                FillZero(V);
                
                var_extra.GetFaceOrientationAndNumber(var, num_face, i,
                                                      num_pos1_face, offset_face1,
                                                      ref_line, cond, new_face,
                                                      num_elem2, num_pos2_face,
                                                      offset_face2, rot, phase, phase_conj);
                
                bool neighbor_face = false;
                int offset_neighbor = 0;
                if (cond != 0)
                  {
                    if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
                      {
                        neighbor_face = true;                        
                        offset_neighbor = var.GetOffsetNeighboringFace(num_face);                        
                      }
                  }
                
                if ( ((cond == 0) && (mesh.Boundary(num_face).GetNbElements() == 2))
                     || neighbor_face)
                  {
                    // case of an internal boundary, with two adjacent elements	
                    for (int k = 0; k < nb_points_face; k++)
                      { 
                        int krot = FacesQuadRotation(rot, k);
                        int num_dof1 = k + offset_face1;
                        
                        int num_dof2 = krot + offset_face2;
                        if (neighbor_face)
                          num_dof2 = k + offset_neighbor;

                        int num_point = offset_quad + k;
                        
                        // outward normale
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
                        
                        // values of u1 and u2
                        u1 = var_extra.extrapolU(0)(num_dof1);
                        if (neighbor_face)
                          u2 = var_extra.Uneighbor(0)(num_dof2);
                        else
                          u2 = var_extra.extrapolU(0)(num_dof2);
                        
                        u2 *= phase;
                        jump_u = u2 - u1;
                        
                        // part due to C : - {mu grad u n} phi + [u] mu grad phi n
                        if (affine)
                          {
                            dfjm1 = var.Glob_DFjm1(i)(0);
                            Mlt(1.0/var.Glob_jacobian(i)(0), dfjm1);
                          }
                        else
                          dfjm1 = var.Glob_DFjm1(i)(num_point);
                        
                        mu = var.GetMu_QuadraturePoint(ref, i, num_point);
                        du1_dn = var_extra.extrapolCgradU(0)(num_dof1);
                        if (neighbor_face)
                          du2_dn = var_extra.CgradUNeighbor(0)(num_dof2);
                        else
                          du2_dn = var_extra.extrapolCgradU(0)(num_dof2);
                        
                        du2_dn *= phase;
                        // penalty term P [u] phi
                        phi_loc = var.alpha_penalization*jump_u;
                        if (var.automatic_choice_penalization)
                          phi_loc *= var.Glob_CoefPenalDG(num_face);
                        
                        phi_loc += - du1_dn + du2_dn;
                                                
                        Mlt(mu, normale, vec_u);
                        Mlt(dfjm1, vec_u, vec_v);
                        vec_v *= jump_u*dsj*PoidsFlux(k);
                        CopyVector(vec_v, k, dV);
                        
                        // mulltiplication by weight of integration
                        phi_loc *= dsj*PoidsFlux(k);
                        V(k) = phi_loc;
                      }
                    
                    if (Fb.UseQuadraturePointsForSh())
                      Fb.ApplyNablaShQuadrature(nat_mat.GetCoefStiffness(),
                                                num_pos1_face, dV, ProdUh_quad, rf);
                    else
                      Fb.ApplyNablaSh(nat_mat.GetCoefStiffness(),
                                      num_pos1_face, dV, Prod_Uh(0), rf);
                    
                    // adding the result to ProdUh
                    if (Fb.UseQuadraturePointsForSh())
                      Fb.ApplyShQuadrature(nat_mat.GetCoefStiffness(),
                                           num_pos1_face, V, ProdUh_quad, rf);
                    else
                      Fb.ApplySh(nat_mat.GetCoefStiffness(), num_pos1_face, V, Prod_Uh(0), rf);
                  }
              }
            
            offset_face1 += nb_points_face;
            offset_quad += nb_points_face;
          }
	
        // projection on dofs
        if ((presence_mass && variable) || (Fb.UseQuadraturePointsForRh())
            || (Fb.UseQuadraturePointsForSh()))
          {
            Fb.ApplyCh(ProdUh_quad, Uh_dof(0));
            Add(cone, Uh_dof(0), Prod_Uh(0));
          }
      }
    else
      {
        if (presence_mass)
          {
            if ( (variable) || (Fb.UseQuadraturePointsForRh()))
              Fb.ApplyCh(ProdUh_quad, Prod_Uh(0));
          }
      }

    // we add contribution to C
    var.AddLocalUnknownVector(alpha, Prod_Uh, i, C);
    
  }
  
  
  //! generic matrix-vector product for nodal elements
  //! and scalar equation (Helmholtz equation)
  /*
    this algorithm separates the geometry \f$ DF_i^{*-1} \f$ from the 
    spatial derivatives. This allows low-storage and fast-algorithm
    in \f$ O(r^4) \f$ instead of \f$ O(r^6) \f$ if the full matrix is stored, where
    r is the order of approximation
  */
  template <class T0, class Complexe, class T1, class TypeEquation>
  void MltAddScalarH1(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
		      const SeldonTranspose& trans, int level,
		      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
		      const Vector<T1>& B2, const T0& beta, Vector<T1>& C2, bool assemble)
  {
    typedef typename TypeEquation::Dimension Dimension;

    Vector<T1> B, C;
    
    A.ApplyRightScaling(B2, C2, B, C);    
    
    if (beta == T0(0))
      C2.Fill(0);
    else if (beta != T0(1))
      Mlt(beta, C2);
    
    const Mesh<Dimension>& mesh = A.var.mesh;
    //const MeshNumbering<Dimension>& mesh_num = A.var.mesh_num;
    
    ExtrapolVariablesProductFEM<T1, TypeEquation> var_extra;
    var_extra.Un_quad.Reallocate(mesh.GetNbElt());

    glob_chrono.Start(VirtualTimer::EXTRAPOL);
    
    // computation of du/dn on quadrature points of faces for SIPG method
    GetExtrapolationAndExchange_Scalar(nat_mat, A, B, var_extra); 
    
    // we finalize exchanges and we add terms due to neighbor elements
    CompleteExchangeAndFlux_Scalar(alpha, nat_mat, A, C, var_extra);
    
    glob_chrono.Stop(VirtualTimer::EXTRAPOL);
    
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    
    // loop over all elements of the mesh
    for (int i = 0; i < mesh.GetNbElt(); i++)
      MltAdd_ElementH1_Scalar(alpha, nat_mat, trans, A, i, B, C, var_extra,
			      A.var.GetReferenceElementH1(i));

    A.AddExtraBoundaryTerms(alpha, B, C);
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
    
    A.ApplyLeftScaling(B2, C2, B, C);    
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C2);
  }

}

#define MONTJOIE_FILE_PROD_MAT_VECT_SCALAR_H1_CXX
#endif
