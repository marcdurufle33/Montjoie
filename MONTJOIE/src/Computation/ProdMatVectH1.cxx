#ifndef MONTJOIE_FILE_PROD_MAT_VECT_H1_CXX

#include "ProdMatVectScalarH1.cxx"

namespace Montjoie
{
    
  //! computation of unknowns on the quadrature points of element
  /*!
    \param[in] var considered problem
    \param[in] num_elem element number
    \param[in] Uloc solution u on the element
    \param[inout] Uquad solution u on quadrature points of the element
    \param[inout] extrapolU array storing values of u on quadrature points of the faces
    \param[in] Fb finite element associated with the considered element
    \param[inout] Uneighbor array storing values u on quadrature points of "neighboring" faces
                            (i.e. faces close to another processor)
    \param[in] compute_uquad if true the values of u
                             on quadrature points of the element are computed
      (if Fb.UseQuadraturePointsForSh() is true), otherwise these values
      are assumed to be already available in Uquad
   */
  template<class TypeEquation, class Vector1, class Vector2, int p, class Vector4, class Vector5>
  void GetExtrapolationU(const EllipticProblem<TypeEquation>& var, int num_elem,
			 const Vector1& Uloc, TinyVector<Vector2, p>& Uquad, Vector4& extrapolU,
                         const ElementReference<typename TypeEquation::Dimension, 1>& Fb,
			 Vector5& Uneighbor, bool compute_uquad)
  {
    typedef typename TypeEquation::Dimension Dimension;
    if (var.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
      return;
    
    // number of dofs inside the element
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int offset_face = mesh_num.OffsetQuadElementNumber(num_elem);
    
    typedef typename TypeEquation::Complexe Complexe;
    Vector<Complexe> face_U;
    int offset_war = Fb.GetNbPointsNodalElt() + Fb.GetNbPointsQuadratureInside();

    if (Fb.UseQuadraturePointsForSh())
      if (compute_uquad)
        {
          for (int n = 0; n < Uquad.GetM(); n++)
            {
              Uquad(n).Reallocate(Fb.GetNbPointsQuadratureInside());
              Fb.ApplyChTranspose(Uloc(n), Uquad(n));
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
	face_U.Reallocate(nb_points_face);
	for (int n = 0; n < Uquad.GetM(); n++)
	  {
            if (Fb.UseQuadraturePointsForSh())
              Fb.ApplyShQuadratureTranspose(num_loc, Uquad(n), face_U, rf);
            else
	      Fb.ApplyShTranspose(num_loc, Uloc(n), face_U, rf);
	    
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
                  Uneighbor(n)(offset_neighbor+j) = face_U(j);
              }
            
	    // storing this extrapolation on vector extrapolU
	    for (int j = 0; j < nb_points_face; j++)
	      extrapolU(n)(offset_face + j) = face_U(j);
	    
	  }
	
	offset_face += nb_points_face;
	offset_war += nb_points_face;
      }
  }
  
  
  //! computation of unknowns on the quadrature points of element
  /*!
    \param[in] var considered problem
    \param[in] num_elem element number
    \param[in] B components of the solution on all degrees of freedom of the mesh
    \param[inout] var_extra object storing values of the solution
                           on quadrature points of elements and faces
    \param[in] nat_mat mass, stiffness and damping coefficients
    \param[in] Fb finite element of the element
   */
  template<class TypeEquation, class Complexe, class Vector1>
  void GetExtrapolationSecondOrder(const EllipticProblem<TypeEquation>& var,
                                   int num_elem, const Vector1& B,
				   ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
                                   const GlobalGenericMatrix<Complexe>& nat_mat,
				   const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    int offset_face = mesh_num.OffsetQuadElementNumber(num_elem);
    bool affine = mesh.IsElementAffine(num_elem);    
    
    // values of u and grad u on quadrature points of face
    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns> valU, gradU;
    
    // small vectors
    TinyVector<Complexe, TypeEquation::nb_unknowns> vec_u, vec_w;
    typename Dimension::R_N normale;
    typename Dimension::MatrixN_N dfjm1;
    TinyVector<TinyVector<Complexe, Dimension::dim_N>, TypeEquation::nb_unknowns> vec_v, vec_Cv;
    TinyVector<Complexe, Dimension::dim_N> tau, tau_normale, grad_v;

    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>& Uloc = var_extra.Uloc;
    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>&
      Uquad = var_extra.GetUnQuad(num_elem, Fb);

    // we get values of u on the element
    var.GetLocalUnknownVector(B, num_elem, var_extra.Uloc);
    
    int num_point = Fb.GetNbPointsQuadratureInside();  
    if (Fb.UseQuadraturePointsForSh())
      for (int n = 0; n < TypeEquation::nb_unknowns; n++)
        {
          Uquad(n).Reallocate(num_point);
          Fb.ApplyChTranspose(Uloc(n), Uquad(n));
        }
    
    int i1 = num_elem - mesh.GetNbElt() + var.GetNbEltPML();
    bool pml_element = (var.IsComplexProblem() && var.InsidePML(num_elem));
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
            int rf = mesh_num.GetOrderQuadrature(num_face);        
            var_extra.GetFaceOrientationAndNumber(var, num_face,
                                                  num_elem, num_loc, offset_face,
                                                  ref, cond, new_face,
                                                  num_elem2, num_loc2, offset_face2,
                                                  rot, phase, phase_conj);
            
            const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
              GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
            
            int offset_neighbor = 0;
            if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
	      offset_neighbor = var.GetOffsetNeighboringFace(num_face);
	    
            for (int n = 0; n < Uquad.GetM(); n++)
              {
                // computation of u on quadrature points of face
                valU(n).Reallocate(nb_points_face);            
                if (Fb.UseQuadraturePointsForSh())
                  Fb.ApplyShQuadratureTranspose(num_loc, Uquad(n), valU(n), rf);
                else
                  Fb.ApplyShTranspose(num_loc, Uloc(n), valU(n), rf);
                
                if (!TypeEquation::FirstOrderFormulation)
                  {
                    // computation of gradient of u on quadrature points of face
                    gradU(n).Reallocate(Dimension::dim_N*nb_points_face);
                    if (Fb.UseQuadraturePointsForSh())
                      Fb.ApplyNablaShQuadratureTranspose(num_loc, Uquad(n), gradU(n), rf);
                    else
                      Fb.ApplyNablaShTranspose(num_loc, Uloc(n), gradU(n), rf);
                  }
              }
            
            // loop over quadrature points
            for (int j = 0; j < nb_points_face; j++)
              {
                int jrot = FacesQuadRotation(rot, j);
                
                CopyVector(valU, j, vec_u);
                if (new_face)
                  normale = var.Glob_normale(num_face)(j);
                else
                  {
                    normale = var.Glob_normale(num_face)(jrot);
                    Mlt(Real_wp(-1), normale);
                  }
                
                if (pml_element)
                  for (int m = 0; m < Dimension::dim_N; m++)
                    {
                      tau(m) = var.GetTauPML(i1, num_point+j, m);
                      tau_normale(m) = normale(m)*tau(m);
                    }
                else
                  tau_normale = normale;
                
                // storing u
                if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
                  CopyVector(vec_u, offset_neighbor+j, var_extra.Uneighbor);
                
                CopyVector(vec_u, offset_face+j, var_extra.extrapolU);
                
                // if matrices Di are discontinuous, we store D u n
                if (TypeEquation::DiscontinuousDiMatrix)
                  {            
                    TypeEquation::ApplyGradientFctTest(var, num_elem, num_point+j, nat_mat,
                                                       ref_d, vec_u, vec_v);
                    
                    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                      vec_w(m) = DotProd(tau_normale, vec_v(m));
                    
                    if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
                      CopyVector(vec_w, offset_neighbor+j, var_extra.DunNeighbor);
                    
                    CopyVector(vec_w, offset_face+j, var_extra.extrapolDun);
                  }
                
                // for second-order formulation, we store C grad u n
                if (!TypeEquation::FirstOrderFormulation)
                  {
                    if (affine)
                      {
                        dfjm1 = var.Glob_DFjm1(num_elem)(0);
                        Mlt(1.0/var.Glob_jacobian(num_elem)(0), dfjm1);
                      }
                    else
                      dfjm1 = var.Glob_DFjm1(num_elem)(num_point+j);
                    
                    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                      {
                        CopyVector(gradU(m), j, grad_v);
                        MltTrans(dfjm1, grad_v, vec_v(m));
                      }
                    
                    TypeEquation::ApplyTensorStiffness(var, num_elem, num_point+j,
                                                       nat_mat, ref_d, vec_v, vec_Cv);
                    
                    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                      vec_w(m) = DotProd(tau_normale, vec_Cv(m));
                    
                    if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
                      CopyVector(vec_w, offset_neighbor+j, var_extra.CgradUNeighbor);
                    
                    CopyVector(vec_w, offset_face+j, var_extra.extrapolCgradU);
                  }	    
              }
            
          }
        
        offset_face += nb_points_face;
        num_point += nb_points_face;   
      }
  }

  
  //! computation of u on quadrature points of the faces and exchange with neighboring processors
  /*!
    \param[in] nat_mat mass, stiffness and damping coefficients
    \param[in] A finite element matrix
    \param[in] B solution u on all degrees of freedom of the mesh
    \param[inout] var_extra object containing values u on quadrature points of elements and faces
   */
  template<class Complexe, class TypeEquation, class T1>
  void GetExtrapolationAndExchange(const GlobalGenericMatrix<Complexe>& nat_mat,
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
        
        for (int p = 0; p < TypeEquation::nb_unknowns; p++)
          {
            var_extra.extrapolU(p).Reallocate(size_extrapol);
            var_extra.Uneighbor(p).Reallocate(size_neighbor);
            if (TypeEquation::DiscontinuousDiMatrix)
              {
                var_extra.extrapolDun(p).Reallocate(size_extrapol);
                var_extra.DunNeighbor(p).Reallocate(size_neighbor);            
              }

            if (!TypeEquation::FirstOrderFormulation)
              {
                var_extra.extrapolCgradU(p).Reallocate(size_extrapol);
                var_extra.CgradUNeighbor(p).Reallocate(size_neighbor);            
              }
          }
      }
    
    // loop over all elements of the mesh to compute u on quadrature points of faces
    for (int i = 0; i < mesh.GetNbElt(); i++)
      GetExtrapolationSecondOrder(A.var, i, B, var_extra, nat_mat,
				  A.var.GetReferenceElementH1(i));
    
#ifdef SELDON_WITH_MPI
    int nb_exchange = 1;
    if (TypeEquation::DiscontinuousDiMatrix)
      nb_exchange++;
    
    if (!TypeEquation::FirstOrderFormulation)
      nb_exchange++;
    
    var_extra.ReallocateExchangeVector(nb_exchange*TypeEquation::nb_unknowns);
    
    // exchanging datas between processors
    if (A.var.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
      //if (true)
      {
        // sending and receiving u
        for (int n = 0; n < TypeEquation::nb_unknowns; n++)
          A.var.ExchangeUfaceDomains(var_extra.Uneighbor(n), var_extra.Xsend(n), var_extra.Xsend_tmp(n),
                                     var_extra.Xneighbor(n), var_extra.Xneighbor_tmp(n), var_extra.request(n), 90+n);
        
        int p = TypeEquation::nb_unknowns;
        
        if (TypeEquation::DiscontinuousDiMatrix)
          {
            // sending and receiving D u n
            for (int n = 0; n < TypeEquation::nb_unknowns; n++)
              A.var.ExchangeUfaceDomains(var_extra.DunNeighbor(n), var_extra.Xsend(n+p), var_extra.Xsend_tmp(n+p),
                                         var_extra.Xneighbor(n+p), var_extra.Xneighbor_tmp(n+p), var_extra.request(n+p), 90+n+p);
            
            p += TypeEquation::nb_unknowns;
          }


        if (!TypeEquation::FirstOrderFormulation)
          {
            // sending and receiving C grad u n
            for (int n = 0; n < TypeEquation::nb_unknowns; n++)
              A.var.ExchangeUfaceDomains(var_extra.CgradUNeighbor(n), var_extra.Xsend(n+p), var_extra.Xsend_tmp(n+p),
                                         var_extra.Xneighbor(n+p), var_extra.Xneighbor_tmp(n+p), var_extra.request(n+p), 90+n+p);
          }
        
      }
#endif
    
  }
  
  
  //! adds local flux associated with an element
  /*!
    \param[in] alpha we add alpha \int_{\partial K} f ...
    \param[in] nat_mat mass, damping and stiffness coefficients
    \param[in] A finite element matrix
    \param[in] i element number
    \param[in] B solution u on all degrees of freedom
    \param[inout] C vector to which the flux is added
    \param[inout] var_extra object containing solution of u on quadrature points
    \param[in] Fb finite element
   */
  template<class Complexe, class TypeEquation, class T0, class Vector1>
  void AddLocalFluxNeighbor(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                            const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                            int i, const Vector1& B, Vector1& C,
                            ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
			    const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = A.var.mesh;
    const MeshNumbering<Dimension>& mesh_num = A.var.GetMeshNumbering(0);
    const EllipticProblem<TypeEquation>& var = A.var;
    int i1 = i - mesh.GetNbElt() + var.GetNbEltPML();
    int nb_dof_elt = Fb.GetNbDof();
    bool pml_element = (var.IsComplexProblem() && var.InsidePML(i));
    int offset_face1 = mesh_num.OffsetQuadElementNumber(i);    
    int offset_quad = Fb.GetNbPointsQuadratureInside();
    bool affine = mesh.IsElementAffine(i);
    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns> ProdUh, V, dV;
    Complexe cone(1.0);

    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>& Uh = var_extra.Uloc;
    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>& Uh_quad
      = var_extra.GetUnQuad(i, Fb);
    
    bool neighbor_elt = false;
    for (int num_pos1_face = 0; num_pos1_face < Fb.GetNbBoundaries(); num_pos1_face++)
      {
	int num_face = mesh.Element(i).numBoundary(num_pos1_face);
        int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        if (var.FaceHasToBeConsideredForBoundaryIntegral(num_face))
          {           
            int ref_line = mesh.Boundary(num_face).GetReference();
            int cond = mesh.GetBoundaryCondition(ref_line);
            if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
              {
                if (!neighbor_elt)
                  {
                    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                      ProdUh(m).Reallocate(nb_dof_elt);
                      
                    FillZero(ProdUh);
                    
                    // we get values of u on the element
                    var.GetLocalUnknownVector(B, i, Uh);

                    neighbor_elt = true;
                  }
                
                for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                  {
                    V(p).Reallocate(nb_points_face);
                    if (!TypeEquation::FirstOrderFormulation)
                      dV(p).Reallocate(nb_points_face*Dimension::dim_N);
                  }
                
                FillZero(V);
                
                int ref = mesh.Element(i).GetReference();
                int ref2 = var.GetRefDomainNeighboringFace(num_face);
                int offset_neighbor = var.GetOffsetNeighboringFace(num_face);
                bool new_face = var.IsNewFace(i)(num_pos1_face);
                typename Dimension::R_N normale; typename Dimension::MatrixN_N dfjm1;
                Real_wp dsj; TinyVector<Complexe, Dimension::dim_N> tau, tau_normale, tmp;
                TinyVector<Complexe, TypeEquation::nb_unknowns> u1, u2, jump_u, du1_n, du2_n,
                  phi_loc, U_loc, psi_loc;
                
                TinyVector<TinyVector<Complexe, Dimension::dim_N>, TypeEquation::nb_unknowns>
                  dU_loc, sigma_loc, tau_loc;
                
                int rf = mesh_num.GetOrderQuadrature(num_face);
                const VectReal_wp& PoidsFlux
                  = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
                const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
                  GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
                
                Complexe phase(1);
                int rot = var.GetRotationNeighboringFace(num_face);
                
                // loop over quadrature points
                for (int k = 0; k < nb_points_face; k++)
                  { 
                    int krot = FacesQuadRotation(rot, k);
                    int num_dof1 = k + offset_face1;
                    int num_dof2 = k + offset_neighbor;
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
                    
                    if (pml_element)
                      for (int m = 0; m < Dimension::dim_N; m++)
                        {
                          tau(m) = var.GetTauPML(i1, num_point, m);
                          tau_normale(m) = normale(m)*tau(m);
                        }
                    else
                      tau_normale = normale;
                    
                    // values of u1 and u2
                    CopyVector(var_extra.extrapolU, num_dof1, u1);
                    CopyVector(var_extra.Uneighbor, num_dof2, u2);
                    
                    u2 *= phase;
                    jump_u = u2 - u1;
                    
                    // part due to D : - {D u n} phi
                    if (TypeEquation::DiscontinuousDiMatrix)
                      {
                        CopyVector(var_extra.extrapolDun, num_dof1, du1_n);
                        CopyVector(var_extra.DunNeighbor, num_dof2, du2_n);
                        
                        du2_n *= phase;
                        phi_loc = - du1_n + du2_n;
                      }
                    else
                      {
                        U_loc = u1 + u2;
                        TypeEquation::ApplyGradientFctTest(var, i, num_point, nat_mat,
                                                           ref, U_loc, tau_loc);
                        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                          phi_loc(m) = - DotProd(tau_loc(m), tau_normale);
                      }
                    
                    // part due to E : E [u] n phi
                    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                      {
                        dU_loc(m) = tau_normale;
                        dU_loc(m) *= jump_u(m);
                      }
                    
                    TypeEquation::ApplyGradientUnknown(var, i, num_point, nat_mat,
                                                       ref, dU_loc, psi_loc);
                    phi_loc += psi_loc;
                    
                    // part due to penalization : P [u] 
                    TypeEquation::MltPenalDG(normale, jump_u, psi_loc, i, num_point, num_face,
                                             nat_mat, ref, ref2, var, Fb);
                    phi_loc += psi_loc;
                    
                    // part due to C : - {C grad u n} phi + [u] C grad phi n
                    if (!TypeEquation::FirstOrderFormulation)
                      {
                        if (affine)
                          {
                            dfjm1 = var.Glob_DFjm1(i)(0);
                            Mlt(1.0/var.Glob_jacobian(i)(0), dfjm1);
                          }
                        else
                          dfjm1 = var.Glob_DFjm1(i)(num_point);
                        
                        CopyVector(var_extra.extrapolCgradU, num_dof1, du1_n);
                        CopyVector(var_extra.CgradUNeighbor, num_dof2, du2_n);
                        
                        du2_n *= phase;
                        phi_loc += - du1_n + du2_n;
                        
                        if (TypeEquation::TensorStiffnessSymmetric)
                          {
                            TypeEquation::ApplyTensorStiffness(var, i, num_point, nat_mat, ref,
                                                               dU_loc, sigma_loc);
                            
                            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                              {
                                Mlt(dfjm1, sigma_loc(m), tmp);
                                tmp *= dsj*PoidsFlux(k);
                                CopyVector(tmp, k, dV(m));
                              }
                          }
                        else
                          {
                            cout << "Not implemented " << endl;
                            abort();
                          }
                      }
                    
                    // mulltiplication by weight of integration
                    phi_loc *= dsj*PoidsFlux(k);
                    CopyVector(phi_loc, k, V);
                  }
                
                if (!TypeEquation::FirstOrderFormulation)
                  {
                    if (Fb.UseQuadraturePointsForSh())
                      for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                        Fb.ApplyNablaShQuadrature(cone, num_pos1_face, dV(m), Uh_quad(m), rf);
                    else
                      for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                        Fb.ApplyNablaSh(cone, num_pos1_face, dV(m), ProdUh(m), rf);
                  }              
                
                // adding the result to ProdUh
                if (Fb.UseQuadraturePointsForSh())
                  for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                    Fb.ApplyShQuadrature(cone, num_pos1_face, V(m), Uh_quad(m), rf);
                else
                  for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                    Fb.ApplySh(cone, num_pos1_face, V(m), ProdUh(m), rf);
              }
          }
        
        offset_face1 += nb_points_face;
        offset_quad += nb_points_face;
      }
    
    // adding the vector to C
    if (neighbor_elt)
      {
        if (Fb.UseQuadraturePointsForSh())
          for (int m = 0; m < TypeEquation::nb_unknowns; m++)
            Fb.ApplyCh(Uh_quad(m), ProdUh(m));
        
        var.AddLocalUnknownVector(alpha, ProdUh, i, C);
      }
    
  }
  
  //! completes transfers between processors and computation of fluxes
  /*!
    \param[in] alpha we add alpha \int_{\partial K} f ...
    \param[in] nat_mat mass, damping and stiffness coefficients
    \param[in] A finite element matrix
    \param[in] B solution u on all degrees of freedom
    \param[inout] C vector to which the fluxes are added
    \param[inout] var_extra object containing solution of u on quadrature points
    */
  template<class Complexe, class TypeEquation, class T0, class Vector1>
  void CompleteExchangeAndFlux(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
                               const FemMatrixFreeClass<Complexe, TypeEquation>& A,
                               const Vector1& B, Vector1& C,
                               ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra)
  {
    
#ifdef SELDON_WITH_MPI
    if (A.var.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
    //if (true)
      {
        // finalizing transfers if needed
        for (int n = 0; n < TypeEquation::nb_unknowns; n++)
          A.var.GetUfaceDomains(var_extra.Uneighbor(n), var_extra.Xsend(n), var_extra.Xsend_tmp(n),
                                var_extra.Xneighbor(n), var_extra.Xneighbor_tmp(n),
                                var_extra.request(n), 90+n);

        int p = TypeEquation::nb_unknowns;
        
        if (TypeEquation::DiscontinuousDiMatrix)
          {
            for (int n = 0; n < TypeEquation::nb_unknowns; n++)
              A.var.GetUfaceDomains(var_extra.DunNeighbor(n), var_extra.Xsend(n+p), var_extra.Xsend_tmp(n+p),
                                    var_extra.Xneighbor(n+p), var_extra.Xneighbor_tmp(n+p),
                                    var_extra.request(n+p), 90+n+p);
            
            p += TypeEquation::nb_unknowns;
          }

        if (!TypeEquation::FirstOrderFormulation)
          {
            for (int n = 0; n < TypeEquation::nb_unknowns; n++)
              A.var.GetUfaceDomains(var_extra.CgradUNeighbor(n), var_extra.Xsend(n+p), var_extra.Xsend_tmp(n+p),
                                    var_extra.Xneighbor(n+p), var_extra.Xneighbor_tmp(n+p),
                                    var_extra.request(n+p), 90+n+p);
          }
      }
#endif
    
    // adding fluxes terms on boundaries at the interface between processors
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      AddLocalFluxNeighbor(alpha, nat_mat, A, i, B, C,
			   var_extra, A.var.GetReferenceElementH1(i));
  }
  
  
  //! Adds volume integrals associated with an element to a vector
  /*!
    \param[in] alpha we add alpha \int_{K} f ...
    \param[in] nat_mat mass, damping and stiffness coefficients
    \param[in] A finite element matrix
    \param[in] i element number
    \param[in] B solution u on all degrees of freedom
    \param[inout] C vector to which the integrals are added
    \param[inout] var_extra object containing solution of u on quadrature points
    \param[in] Fb finite element
    */  
  template<class T0, class TypeEquation, class Complexe, class Vector1>
  void MltAdd_ElementH1(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
			const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			int i, const Vector1& B, Vector1& C,
                        ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra,
			const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = A.var.mesh;
    const MeshNumbering<Dimension>& mesh_num = A.var.GetMeshNumbering(0);
    const EllipticProblem<TypeEquation>& var = A.var;
    int i1 = i - mesh.GetNbElt() + var.GetNbEltPML();
    int nb_dof_elt = Fb.GetNbDof();
    int nb_points_quad = Fb.GetNbPointsQuadratureInside();
    //int nb_points_nodal = Fb.GetNbPointsNodalElt();
    bool variable = var.UseNumericalIntegration(i);
    bool affine = mesh.IsElementAffine(i);
    int ref = mesh.Element(i).GetReference();
    
    TinyVector<Complexe, Dimension::dim_N> tmp, tau, vh_loc;
    TinyVector<Complexe, TypeEquation::nb_unknowns> U_loc, phi_loc, psi_loc;
    TinyVector<TinyVector<Complexe, Dimension::dim_N>,
      TypeEquation::nb_unknowns> dU_loc, sigma_loc, tau_loc;
    
    bool pml_element = (var.IsComplexProblem() && var.InsidePML(i));
    bool first_order = TypeEquation::FirstOrderFormulation;
    Real_wp coef, coef_stiff;
    Complexe cone, czero;
    SetComplexOne(cone); SetComplexZero(czero);
    
    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns> ProdUh, V, dV;

    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>& Uh = var_extra.Uloc;
    TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>& Uh_quad
      = var_extra.GetUnQuad(i, Fb);
    
    // we get values of u on the element
    var.GetLocalUnknownVector(B, i, Uh);

    if ((variable)||(Fb.UseQuadraturePointsForRh())
        ||(Fb.UseQuadraturePointsForSh())||(!first_order))
      {
        TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>
          ProdUh_quad, gradUh_quad;
        
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          {
            ProdUh_quad(m).Reallocate(nb_points_quad);
            gradUh_quad(m).Reallocate(nb_points_quad*Dimension::dim_N);
            ProdUh(m).Reallocate(nb_dof_elt);
            ProdUh(m).Fill(0);
            ProdUh_quad(m).Fill(0);
          }   

        // computation of u on quadrature points
        if (!Fb.UseQuadraturePointsForSh())
          for (int m = 0; m < TypeEquation::nb_unknowns; m++)
            {
              Uh_quad(m).Reallocate(nb_points_quad);
              Fb.ApplyChTranspose(Uh(m), Uh_quad(m));
            }
        
        // computation of needed derivatives
        if (Fb.UseQuadraturePointsForRh())
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              if (var_extra.unknown_to_derive(m))
                Fb.ApplyRhQuadratureTranspose(Uh_quad(m), gradUh_quad(m));
          }
        else
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              if (var_extra.unknown_to_derive(m))
                Fb.ApplyRhTranspose(Uh(m), gradUh_quad(m));
          }
        
        // loop on each quadrature point
        // applying geometry transformations and physical coefficients
        for (int j = 0; j < nb_points_quad; j++)
          { 
            if (pml_element)
              for (int k = 0; k < Dimension::dim_N; k++)
                tau(k) = var.GetTauPML(i1, j, k);
            
            int jmat = j;
            if (affine)
              jmat = 0;
            
            const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N>&
              dfjm1 = var.Glob_DFjm1(i)(jmat);
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              {
                U_loc(m) = Uh_quad(m)(j);
                if (var_extra.unknown_to_derive(m))
                  {
                    CopyVector(gradUh_quad(m), j, tmp);
                    		
                    // we apply DF_i^{-t}
                    MltTrans(dfjm1, tmp, dU_loc(m));
                    
                    if (pml_element)
                      dU_loc(m) = dU_loc(m)*tau;
                  }
                else
                  dU_loc(m).Fill(czero);
              }
            
            // application of physical properties
            if (!first_order)
              TypeEquation::ApplyTensorStiffness(var, i, j, nat_mat, ref, dU_loc, sigma_loc);
            
            TypeEquation::ApplyGradientUnknown(var, i, j, nat_mat, ref, dU_loc, phi_loc);
            TypeEquation::ApplyGradientFctTest(var, i, j, nat_mat, ref, U_loc, tau_loc);
            TypeEquation::ApplyTensorMass(var, i, j, nat_mat, ref, U_loc, psi_loc);
            
            if (pml_element)
              for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                {
		  tau_loc(m) = tau_loc(m)*tau;
		  sigma_loc(m) = sigma_loc(m)*tau;
                }
            
            Real_wp jacob;
            if (affine)
              jacob = var.Glob_jacobian(i)(0);
            else
              jacob = var.Glob_jacobian(i)(j)/Fb.WeightsND(j);
            
            coef = 1.0/jacob;
            coef_stiff = Fb.WeightsND(j);
            if (!first_order)
              {
                // second order terms du/dx_i dphi/dx_j
                for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                  {
                    ProdUh_quad(m)(j) += coef_stiff*(phi_loc(m) + psi_loc(m)*jacob);
                    
                    Mlt(coef, sigma_loc(m));
                    sigma_loc(m) += tau_loc(m);
                    
                    // we apply DF_i^{-1}
                    Mlt(dfjm1, sigma_loc(m), vh_loc);
                    
                    Mlt(coef_stiff, vh_loc);
                    CopyVector(vh_loc, j, gradUh_quad(m));
                  }
              }
            else
              {
                // only first order terms
                for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                  {
                    ProdUh_quad(m)(j) += coef_stiff*(phi_loc(m) + psi_loc(m)*jacob);
                    
                    // we apply DF_i^{-1}
                    Mlt(dfjm1, tau_loc(m), vh_loc);
                    
                    Mlt(coef_stiff, vh_loc);
                    CopyVector(vh_loc, j, gradUh_quad(m));
                  }
              }
          }
        
        // integration against grad(phi) for each unknown
        if (Fb.UseQuadraturePointsForRh())
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              if (var_extra.fct_test_to_derive(m))
                {
                  Fb.ApplyRhQuadrature(gradUh_quad(m), Uh_quad(m));
                  Add(Real_wp(1), Uh_quad(m), ProdUh_quad(m));
                }
          }
        else
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              if (var_extra.fct_test_to_derive(m))
                {
                  Fb.ApplyRh(gradUh_quad(m), Uh(m));
                  Add(Real_wp(1), Uh(m), ProdUh(m));                  
                }
          }
        
        // integration against basis functions
        if (Fb.UseQuadraturePointsForSh())
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              Copy(ProdUh_quad(m), Uh_quad(m));
          }
        else
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              {
                Fb.ApplyCh(ProdUh_quad(m), Uh(m));
                Add(Real_wp(1), Uh(m), ProdUh(m));
                Uh_quad(m).Clear();
              }
          }
      }
    else
      {
        TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>
          ProdDh, gradUh;
        
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          {
            ProdDh(m).Reallocate(nb_dof_elt);
            gradUh(m).Reallocate(nb_dof_elt*Dimension::dim_N);
            ProdUh(m).Reallocate(nb_dof_elt);
            ProdUh(m).Fill(0);
            ProdDh(m).Fill(0);
          }
        
        // affine element and constant physical coefficients
        const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N>&
          dfjm1 = var.Glob_DFjm1(i)(0);
	
        Real_wp jacob = var.Glob_jacobian(i)(0);
        
        for (int k = 0; k < nb_dof_elt; k++)
          {
            CopyVector(Uh, k, U_loc);
            TypeEquation::ApplyTensorMass(var, i, 0, nat_mat, ref, U_loc, psi_loc);
            TypeEquation::ApplyGradientFctTest(var, i, 0, nat_mat, ref, U_loc, tau_loc);
            psi_loc *= jacob;
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              {
                Mlt(dfjm1, tau_loc(m), vh_loc);
                CopyVector(vh_loc, k, gradUh(m));
              }
            
            CopyVector(psi_loc, k, ProdUh);
          }
        
        // multiplication by mass matrix
	for (int p = 0; p < TypeEquation::nb_unknowns; p++)
          {
            Fb.MltMassMatrix(ProdUh(p));
            if (var_extra.fct_test_to_derive(p))
              {
                Fb.ApplyConstantRh(gradUh(p), ProdDh(p));
                Add(Real_wp(1), ProdDh(p), ProdUh(p));
              }
            
            if (var_extra.unknown_to_derive(p))
              Fb.ApplyConstantRhTranspose(Uh(p), gradUh(p));            
          }
        
        // constant stiffness matrices
        for (int j = 0; j < nb_dof_elt; j++)
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              {                
                if (var_extra.unknown_to_derive(m))
                  {
                    CopyVector(gradUh(m), j, tmp);
                    
                    // we apply DF_i^{-t}
                    MltTrans(dfjm1, tmp, dU_loc(m));
                  }
                else
                  dU_loc(m).Fill(czero);
              }
            
            TypeEquation::ApplyGradientUnknown(var, i, j, nat_mat, ref, dU_loc, phi_loc);
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              ProdUh(m)(j) += phi_loc(m);
          }
      }
    
    
    // offset for quadrature points on face 1
    int offset_face1 = mesh_num.OffsetQuadElementNumber(i);
    int num_elem2, ref_line, cond, num_pos2_face, offset_face2, rot;
    bool new_face; Complexe phase, phase_conj;
    typename Dimension::R_N normale; Real_wp dsj;
    typename Dimension::MatrixN_N dfjm1;
    TinyVector<Complexe, Dimension::dim_N> tau_normale;
    TinyVector<Complexe, TypeEquation::nb_unknowns> u1, u2, jump_u, du1_n, du2_n;
    int offset_quad = Fb.GetNbPointsQuadratureInside();
    bool neighbor_elt = false;
    for (int num_pos1_face = 0; num_pos1_face < Fb.GetNbBoundaries(); num_pos1_face++)
      {
	int num_face = mesh.Element(i).numBoundary(num_pos1_face);
        int nb_points_face = mesh_num.GetNbPointsQuadratureBoundary(num_face);
        if (var.FaceHasToBeConsideredForBoundaryIntegral(num_face))
          {            
            int rf = mesh_num.GetOrderQuadrature(num_face);
            const VectReal_wp& PoidsFlux
              = mesh_num.number_map.GetFluxWeight(rf, mesh.Boundary(num_face));
            const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
              GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
            
            for (int p = 0; p < TypeEquation::nb_unknowns; p++)
              {
                V(p).Reallocate(nb_points_face);
                if (!TypeEquation::FirstOrderFormulation)
                  dV(p).Reallocate(nb_points_face*Dimension::dim_N);
              }
            
            FillZero(V);
            
            var_extra.GetFaceOrientationAndNumber(var, num_face,
                                                  i, num_pos1_face, offset_face1,
                                                  ref_line, cond, new_face,
                                                  num_elem2, num_pos2_face,
                                                  offset_face2, rot, phase, phase_conj);
            
            if (cond != 0)
              {
                if (cond == BoundaryConditionEnum::LINE_NEIGHBOR)
                  {
                    neighbor_elt = true;
                  }
                else if (TypeEquation::FirstOrderFormulation)
                  {
                    // case of a boundary condition (Dirichlet, Neumann ...)
                    for (int k = 0; k < nb_points_face; k++)
                      {
                        int num_dof1 = k + offset_face1;
                        int num_point = offset_quad + k;
                        
                        // outward normale
                        int krot = FacesQuadRotation(rot, k);
                        
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
                        
                        if (pml_element)
                          for (int m = 0; m < Dimension::dim_N; m++)
                            {
                              tau(m) = var.GetTauPML(i1, num_point, m);
                              tau_normale(m) = normale(m)*tau(m);
                            }
                        else
                          tau_normale = normale;
                        
                        // values on the face are retrieved
                        CopyVector(var_extra.extrapolU, num_dof1, U_loc);
                        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                          {
                            dU_loc(m) = tau_normale;
                            dU_loc(m) *= -U_loc(m);
                          }
                        
                        // we compute phi_loc = - D u1 n - E u1 n
                        TypeEquation::ApplyGradientFctTest(var, i, num_point,
                                                           nat_mat, ref, U_loc, tau_loc);
                        TypeEquation::ApplyGradientUnknown(var, i, num_point,
                                                           nat_mat, ref, dU_loc, phi_loc);
                        
                        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                          phi_loc(m) -= DotProd(tau_loc(m), tau_normale);
                        
                        // and we add boundary condition term
                        TypeEquation::MltNabc(normale, ref_line, U_loc, psi_loc, i, num_point,
                                              nat_mat, ref, var, Fb);
                        
                        phi_loc += psi_loc;
                        
                        // mulltiplication by weight of integration
                        phi_loc *= dsj*PoidsFlux(k);
                        CopyVector(phi_loc, k, V);
                      }
                  }
              }
            else if (mesh.Boundary(num_face).GetNbElements() == 2)
              {                
                int ref2 = mesh.Element(num_elem2).GetReference();
                // case of an internal boundary, with two adjacent elements	
                for (int k = 0; k < nb_points_face; k++)
                  { 
                    int krot = FacesQuadRotation(rot, k);
                    int num_dof1 = k + offset_face1;
                    int num_dof2 = krot + offset_face2;
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
                    
                    if (pml_element)
                      for (int m = 0; m < Dimension::dim_N; m++)
                        {
                          tau(m) = var.GetTauPML(i1, num_point, m);
                          tau_normale(m) = normale(m)*tau(m);
                        }
                    else
                      tau_normale = normale;
                    
                    // values of u1 and u2
                    CopyVector(var_extra.extrapolU, num_dof1, u1);
                    CopyVector(var_extra.extrapolU, num_dof2, u2);
                    
                    u2 *= phase;
                    jump_u = u2 - u1;
                    
                    // part due to D : - {D u n} phi
                    if (TypeEquation::DiscontinuousDiMatrix)
                      {
                        CopyVector(var_extra.extrapolDun, num_dof1, du1_n);
                        CopyVector(var_extra.extrapolDun, num_dof2, du2_n);
                        
                        du2_n *= phase;
                        phi_loc = - du1_n + du2_n;
                      }
                    else
                      {
                        U_loc = u1 + u2;
                        TypeEquation::ApplyGradientFctTest(var, i, num_point, nat_mat,
                                                           ref, U_loc, tau_loc);
                        
                        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                          phi_loc(m) = - DotProd(tau_loc(m), tau_normale);
                      }
                    
                    // part due to E : E [u] n phi
                    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                      {
                        dU_loc(m) = tau_normale;
                        dU_loc(m) *= jump_u(m);
                      }
                    
                    TypeEquation::ApplyGradientUnknown(var, i, num_point, nat_mat,
                                                       ref, dU_loc, psi_loc);
                    phi_loc += psi_loc;
                    
                    // part due to penalization : P [u] 
                    TypeEquation::MltPenalDG(normale, jump_u, psi_loc, i, num_point, num_face,
					     nat_mat, ref, ref2, var, Fb);
                    
                    phi_loc += psi_loc;
                    
                    // part due to C : - {C grad u n} phi + [u] C grad phi n
                    if (!TypeEquation::FirstOrderFormulation)
                      {
                        if (affine)
                          {
                            dfjm1 = var.Glob_DFjm1(i)(0);
                            Mlt(Real_wp(1)/var.Glob_jacobian(i)(0), dfjm1);
                          }
                        else
                          dfjm1 = var.Glob_DFjm1(i)(num_point);
                        
                        CopyVector(var_extra.extrapolCgradU, num_dof1, du1_n);
                        CopyVector(var_extra.extrapolCgradU, num_dof2, du2_n);
                        
                        du2_n *= phase;
                        phi_loc += - du1_n + du2_n;
                        
                        if (TypeEquation::TensorStiffnessSymmetric)
                          {
                            TypeEquation::ApplyTensorStiffness(var, i, num_point, nat_mat, ref,
                                                               dU_loc, sigma_loc);
                            
                            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                              {
                                Mlt(dfjm1, sigma_loc(m), tmp);
                                tmp *= dsj*PoidsFlux(k);
                                CopyVector(tmp, k, dV(m));
                              }
                          }
                        else
                          {
                            cout << "Not implemented " << endl;
                            abort();
                          }
                      }
                    
                    // mulltiplication by weight of integration
                    phi_loc *= dsj*PoidsFlux(k);
                    CopyVector(phi_loc, k, V);
                  }
                
                if (!TypeEquation::FirstOrderFormulation)
                  {
                    if (Fb.UseQuadraturePointsForSh())
                      for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                        Fb.ApplyNablaShQuadrature(cone, num_pos1_face, dV(m), Uh_quad(m), rf);
                    else
                      for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                        Fb.ApplyNablaSh(cone, num_pos1_face, dV(m), ProdUh(m), rf);
                  }
              }
            
            // adding the result to ProdUh
            if (Fb.UseQuadraturePointsForSh())
              for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                Fb.ApplyShQuadrature(cone, num_pos1_face, V(m), Uh_quad(m), rf);
            else
              for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                Fb.ApplySh(cone, num_pos1_face, V(m), ProdUh(m), rf);
          }

        offset_face1 += nb_points_face;
        offset_quad += nb_points_face;
      }
    
    // adding the vector to C
    if (Fb.UseQuadraturePointsForSh())
      {
        if (!neighbor_elt)
          {
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              Fb.ApplyCh(Uh_quad(m), ProdUh(m));
            
            var.AddLocalUnknownVector(alpha, ProdUh, i, C);
          }
      }
    else
      var.AddLocalUnknownVector(alpha, ProdUh, i, C);
    
  }
    
  
  //! generic matrix-vector product for tensorized nodal elements (quad or hex) for any equation
  /*
    \param[in] alpha coefficient
    \param[in] beta coefficient
    \param[in] nat_mat mass, damping and stiffness coefficients
    \param[in] phi ghost argument to know the number of unknowns
    \param[in] trans SeldonTrans or SeldonNoTrans
    \param[in] level matrix-vector product can be performed for a subset of elements
    \param[in] A finite element matrix
    \param[in] B2 vector that will be multiplied
    \param[inout] C2 result
    \param[in] assemble if true the result is assembled (components on overlapped dofs are summed)
    Computation of C2 = beta C2 + alpha *A * B2
    this algorithm separates the geometry \f$ DF_i^{*-1} \f$ from the 
    spatial derivatives. This allows low-storage and fast-algorithm
    in \f$ O(r^4) \f$ instead of \f$ O(r^6) \f$ if the full matrix is stored, where
    r is the order of approximation
    NB : we suppose that quadrature and dof points are the same (mass lumping)
  */
  template <class T0, class Complexe, class Vector1, class TypeEquation>
  void MltAddVectorH1(const T0& alpha, const GlobalGenericMatrix<Complexe>& nat_mat,
		      const SeldonTranspose& trans, int level,
		      const FemMatrixFreeClass<Complexe, TypeEquation>& A,
		      const Vector1& B2, const T0& beta, Vector1& C2, bool assemble)
  {
    typedef typename TypeEquation::Dimension Dimension;

    Vector1 B, C;
    A.ApplyRightScaling(B2, C2, B, C);    

    if (beta == T0(0))
      C2.Fill(0);
    else
      Mlt(beta, C2);
    
    const Mesh<Dimension>& mesh = A.var.mesh;
    //const MeshNumbering<Dimension>& mesh_num = A.var.mesh_num;
    ExtrapolVariablesProductFEM<Complexe, TypeEquation>& var_extra
      = const_cast<FemMatrixFreeClass<Complexe, TypeEquation>& >(A).GetExtrapolVariables();
    
    var_extra.Un_quad.Reallocate(mesh.GetNbElt());
    TypeEquation::GetNeededDerivative(A.var, nat_mat, var_extra.unknown_to_derive,
                                      var_extra.fct_test_to_derive);
    
    glob_chrono.Start(VirtualTimer::EXTRAPOL);
    
    // computation of u, C grad u n on quadrature points of each face
    GetExtrapolationAndExchange(nat_mat, A, B, var_extra); 
    
    glob_chrono.Stop(VirtualTimer::EXTRAPOL);
    
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    // loop over all elements of the mesh
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      MltAdd_ElementH1(alpha, nat_mat, A, i, B, C, var_extra, A.var.GetReferenceElementH1(i));
    
    A.AddExtraBoundaryTerms(Complexe(alpha), B, C);
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
    
    // we finalize exchanges and we add terms due to neighbor elements
    CompleteExchangeAndFlux(alpha, nat_mat, A, B, C, var_extra);
    
    A.ApplyLeftScaling(B2, C2, B, C);    
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C2);
  }
  

  template<class Complexe, class Vector2>
  void ExtractApplyDF_PML(const Vector<Complexe>&, int nb_unknowns, int j,
                          bool ortho_z, Dimension2& dim, int offset_elt,
                          const Matrix2_2& dfjm1, const Complexe& coef,
                          Vector2& Vh, Vector2& Vtilde, Vector2& Vterce)
  {
    abort();
  }


  template<class Complexe, class Vector2>
  void ExtractApplyDF_PML(const Vector<Complexe>& B, int nb_unknowns, int j,
                          bool ortho_z, Dimension3& dim, int offset_elt,
                          const Matrix3_3& dfjm1, const Complexe& coef,
                          Vector2& Vh, Vector2& Vtilde, Vector2& Vterce)
  {
    TinyVector<Complexe, 3> tmp, vh_loc;

    if (ortho_z)
      {
        for (int m = 0; m < nb_unknowns; m++)
          {
            ExtractVector(B, offset_elt + m*3, tmp);
                        
            vh_loc(0) = dfjm1(0, 0)*tmp(0);
            vh_loc(1) = dfjm1(1, 1)*tmp(1);
            vh_loc(2) = dfjm1(2, 2)*tmp(2);
            Mlt(coef, vh_loc);
            
            ExtractVector(vh_loc, 3*j, Vh(m));
            
            vh_loc(1) = dfjm1(1, 0)*tmp(0);
            vh_loc(0) = dfjm1(0, 1)*tmp(1);
            vh_loc(2) = 0;
            Mlt(coef, vh_loc);
            
            ExtractVector(vh_loc, 3*j, Vtilde(m));
          }
      }
    else
      {
        for (int m = 0; m < nb_unknowns; m++)
          {
            ExtractVector(B, offset_elt + m*3, tmp);
            
            vh_loc(0) = dfjm1(0, 0)*tmp(0);
            vh_loc(1) = dfjm1(1, 1)*tmp(1);
            vh_loc(2) = dfjm1(2, 2)*tmp(2);
            Mlt(coef, vh_loc);
            
            ExtractVector(vh_loc, 3*j, Vh(m));
            
            vh_loc(0) = dfjm1(0, 1)*tmp(1);
            vh_loc(1) = dfjm1(1, 0)*tmp(0);
            vh_loc(2) = dfjm1(2, 0)*tmp(0);
            Mlt(coef, vh_loc);
            
            ExtractVector(vh_loc, 3*j, Vtilde(m));
            
            vh_loc(0) = dfjm1(0, 2)*tmp(2);
            vh_loc(1) = dfjm1(1, 2)*tmp(2);
            vh_loc(2) = dfjm1(2, 1)*tmp(1);
            Mlt(coef, vh_loc);
            
            ExtractVector(vh_loc, 3*j, Vterce(m));
          }
      }
  }


#ifdef MONTJOIE_WITH_TWO_DIM
  template<class Vector2, class Complexe, class TypeEquation>
  void ApplyRhSplitElement(Vector2& Vh, Vector2& Uh, const IVect& Nodle,
			   int nb_dof_elt, Vector2& Vtilde, Vector2& Vterce, int nb_unknowns, int Nvol, int Npml,
			   bool ortho, bool ortho_z, Vector<Complexe>& C,
			   const EllipticProblem<TypeEquation>& var, const ElementReference<Dimension2, 1>& Fb)
  {
    if (! ortho)
      {
	abort();
      }

    const MeshNumbering<Dimension2>& mesh_num = var.GetMeshNumbering(0);
    VectReal_wp Uy(nb_dof_elt); Uy.Fill(0);
    for (int m = 0; m < nb_unknowns; m++)
      {
	Fb.ApplyRhSplit(Vh(m), Uh(m), Uy);
	for (int j = 0; j < nb_dof_elt; j++)
	  {
	    int num_dof = Nodle(j);
	    if (num_dof >= 0)
	      {
		int num_pml = Nvol*nb_unknowns + mesh_num.GetDofPML(num_dof);
		C(num_dof + m*Nvol) -= Uh(m)(j) + Uy(j);
		C(num_pml + m*Npml) -= Uh(m)(j) - Uy(j);
	      }
	  }
      }
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  template<class Vector2, class Complexe, class TypeEquation>
  void ApplyRhSplitElement(Vector2& Vh, Vector2& Uh, const IVect& Nodle,
			   int nb_dof_elt, Vector2& Vtilde, Vector2& Vterce, int nb_unknowns, int Nvol, int Npml,
			   bool ortho, bool ortho_z, Vector<Complexe>& C,
			   const EllipticProblem<TypeEquation>& var, const ElementReference<Dimension3, 1>& Fb)
  {
    Complexe vx, vy, vz;
    VectReal_wp Uy(nb_dof_elt), Uz(nb_dof_elt);
    Uy.Fill(0); Uz.Fill(0);
    const MeshNumbering<Dimension3>& mesh_num = var.GetMeshNumbering(0);
    if (ortho)
      {
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    Fb.ApplyRhSplit(Vh(m), Uh(m), Uy, Uz);
	    for (int j = 0; j < nb_dof_elt; j++)
	      {
		int num_dof = Nodle(j);
		if (num_dof >= 0)
		  {
		    int num_pml = Nvol*nb_unknowns + mesh_num.GetDofPML(num_dof);
		    int num_pml2 = num_pml + Npml*nb_unknowns;
		    C(num_dof + m*Nvol) -= Uh(m)(j) + Uy(j) + Uz(j);
                    C(num_pml + m*Npml) -= Uh(m)(j) - Uy(j);
                    C(num_pml2 + m*Npml) -= Uh(m)(j) - Uz(j);
		  }
	      }                                        
	  }
      }
    else if (ortho_z)
      {
	VectReal_wp Ux2(nb_dof_elt), Uy2(nb_dof_elt), Uz2(nb_dof_elt);
	Ux2.Fill(0); Uy2.Fill(0); Uz2.Fill(0);
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    Fb.ApplyRhSplit(Vh(m), Uh(m), Uy, Uz);
	    Fb.ApplyRhSplit(Vtilde(m), Ux2, Uy2, Uz2);
	    for (int j = 0; j < nb_dof_elt; j++)
	      {
		int num_dof = Nodle(j);
		if (num_dof >= 0)
		  {
		    int num_pml = Nvol*nb_unknowns + mesh_num.GetDofPML(num_dof);
		    int num_pml2 = num_pml + Npml*nb_unknowns;
		    vx = Uh(m)(j) + Uy2(j);
		    vy = Ux2(j) + Uy(j);
		    vz = Uz(j);
		    C(num_dof + m*Nvol) -= vx + vy + vz;
		    C(num_pml + m*Npml) -= vx - vy;
		    C(num_pml2 + m*Npml) -= vx - vz;
		  }
	      }                                        
	  }
      }
    else
      {
	VectReal_wp Ux2(nb_dof_elt), Uy2(nb_dof_elt), Uz2(nb_dof_elt);
	VectReal_wp Ux3(nb_dof_elt), Uy3(nb_dof_elt), Uz3(nb_dof_elt);
	Ux2.Fill(0); Uy2.Fill(0); Uz2.Fill(0);
	Ux3.Fill(0); Uy3.Fill(0); Uz3.Fill(0);
	for (int m = 0; m < nb_unknowns; m++)
	  {
	    Fb.ApplyRhSplit(Vh(m), Uh(m), Uy, Uz);
	    Fb.ApplyRhSplit(Vtilde(m), Ux2, Uy2, Uz2);
	    Fb.ApplyRhSplit(Vterce(m), Ux3, Uy3, Uz3);
	    for (int j = 0; j < nb_dof_elt; j++)
	      {
		int num_dof = Nodle(j);
		if (num_dof >= 0)
		  {
		    int num_pml = Nvol*nb_unknowns + mesh_num.GetDofPML(num_dof);
		    int num_pml2 = num_pml + Npml*nb_unknowns;
		    vx = Uh(m)(j) + Uy2(j) + Uz2(j);
		    vy = Ux2(j) + Uy(j) + Uz3(j);
		    vz = Ux3(j) + Uy3(j) + Uz(j);
		    C(num_dof + m*Nvol) -= vx + vy + vz;
		    C(num_pml + m*Npml) -= vx - vy;
		    C(num_pml2 + m*Npml) -= vx - vz;
		  }
	      }                                        
	  }
      }
  }
#endif
  
  
  //! Performs matrix vector product for stiffness matrix related to scalar unknowns
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
  template<class T0, class Complexe, class TypeEquation,
	   class Vector, class Vector2>
  void MltAdd_ElementH1_Scalar(const T0& alpha,
                               const FemMatrixFreeClass<Complexe, TypeEquation>& A,
			       int i, const Vector& B, Vector& C, Vector2& Uh,
			       Vector2& Vh, const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    const EllipticProblem<TypeEquation>& var = A.var;
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    TinyVector<Complexe, Dimension::dim_N> tmp, vh_loc;
    Vector2 Vtilde, Vterce;
    int dim_N = tmp.GetM();
    Complexe coef;
    bool affine = mesh.IsElementAffine(i);
    bool ortho = (var.OrthogonalElement(i) == 0);
    bool ortho_z = (var.OrthogonalElement(i) == 1);
    if (!var.InsidePML(i))
      ortho = true;

    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    
    int nb_unknowns = Uh.GetM();
    int Nvol = mesh_num.GetNbDof();
    for (int m = 0; m < nb_unknowns; m++)
      {
	if (Uh(m).GetM() != nb_dof_elt)
	  {
	    Uh(m).Reallocate(nb_dof_elt);
	    Uh(m).Fill(0);
	  }
	
	Vh(m).Reallocate(dim_N*nb_quad);
        if (ortho_z)
          Vtilde(m).Reallocate(dim_N*nb_quad);
        else if (!ortho)
          {
            Vtilde(m).Reallocate(dim_N*nb_quad);
            Vterce(m).Reallocate(dim_N*nb_quad);
          }
      }
    
    int offset = var.GetOffsetDofV(i);
    Dimension dim;
    if (!affine)
      for (int j = 0; j < nb_quad; j++)
	{ 	  
	  coef = alpha*Fb.WeightsND(j);
	  const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> &
	    dfjm1 = var.Glob_DFjm1(i)(j);
	  
          int offset_elt = offset + nb_unknowns*dim_N*j;
	  if (ortho)
            {
              for (int m = 0; m < nb_unknowns; m++)
                {
                  ExtractVector(B, offset_elt + m*Dimension::dim_N, tmp);
                  
                  Mlt(dfjm1, tmp, vh_loc);
                  Mlt(coef, vh_loc);
                  
                  ExtractVector(vh_loc, Dimension::dim_N*j, Vh(m));
                }
            }
          else 
            ExtractApplyDF_PML(B, nb_unknowns, j, ortho_z, dim, offset_elt,
                               dfjm1, coef, Vh, Vtilde, Vterce);

	}
    else
      {
	const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> &
          dfjm1 = var.Glob_DFjm1(i)(0);
        
	for (int j = 0; j < nb_quad; j++)
	  { 
	    coef = alpha*Fb.WeightsND(j);
            int offset_elt = offset + nb_unknowns*dim_N*j;
	    if (ortho)
              {
                for (int m = 0; m < nb_unknowns; m++)
                  {
                    ExtractVector(B, offset_elt + m*Dimension::dim_N, tmp);
                    
                    Mlt(dfjm1, tmp, vh_loc);
                    Mlt(coef, vh_loc);
                    
                    ExtractVector(vh_loc, Dimension::dim_N*j, Vh(m));
                  }
              }
            else
              ExtractApplyDF_PML(B, nb_unknowns, j, ortho_z, dim, offset_elt,
                                 dfjm1, coef, Vh, Vtilde, Vterce);
          }
      }
    
    // integration against \nabla \varphi
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    if (var.InsidePML(i))
      {
        int Npml = mesh_num.GetNbDofPML();
	ApplyRhSplitElement(Vh, Uh, Nodle, nb_dof_elt, Vtilde, Vterce,
			    nb_unknowns, Nvol, Npml, ortho, ortho_z, C, var, Fb);
	
      }
    else
      {
        for (int m = 0; m < nb_unknowns; m++)
          Fb.ApplyRh(Vh(m), Uh(m));
        
        // adding interactions to final vector
        for (int j = 0; j < nb_dof_elt; j++)
          {
            int num_dof = Nodle(j);
            if (num_dof >= 0)
              for (int m = 0; m < nb_unknowns; m++)
                C(num_dof + m*Nvol) -= Uh(m)(j);
          }
      }
  }
  

  //! Performs matrix vector product for stiffness matrix related to vectorial unknowns
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
  template<class T0, class Complexe, class TypeEquation,
	   class Vector, class Vector2, class Vector3>
  void MltAdd_ElementH1_Vectorial(const T0& alpha,
                                  const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				  int i, const Vector& B, Vector& C, Vector2& dU_loc, Vector3& Uh,
				  Vector3& Vh, const ElementReference<typename TypeEquation::Dimension, 1>& Fb)
  {
    const EllipticProblem<TypeEquation>& var = A.var;
    typedef typename TypeEquation::Dimension Dimension;
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    TinyVector<Complexe, Dimension::dim_N> tmp, vh_loc;
    int dim_N = tmp.GetM();
    Complexe coef;
    
    int nb_dof_elt = Fb.GetNbDof();
    int nb_quad = Fb.GetNbPointsQuadratureInside();
    bool affine = mesh.IsElementAffine(i);
    
    int nb_unknowns = Uh.GetM();
    for (int m = 0; m < nb_unknowns; m++)
      {
	Uh(m).Reallocate(nb_dof_elt);
	Uh(m).Fill(0);
	if (Vh(m).GetM() != dim_N*nb_quad)
	  {
	    Vh(m).Reallocate(dim_N*nb_quad);
	    Vh(m).Fill(0);
	  }
      }
    
    int Nvol = mesh_num.GetNbDof();
    // we get values of u for each unknown and dof
    const IVect& Nodle = mesh_num.Element(i).GetNodle();
    for (int j = 0; j < nb_dof_elt; j++ )
      {
	int num_dof = Nodle(j);
	if(num_dof >= 0)
	  for (int m = 0; m < nb_unknowns; m++)
	    Uh(m)(j) = B(num_dof + m*Nvol);
      }
    
    // computation of gradient of each unknown of u
    for (int m = 0; m < nb_unknowns; m++)
      Fb.ApplyRhTranspose(Uh(m), Vh(m));
    
    int offset = var.GetOffsetDofV(i); 
    // loop on each quadrature point
    // geometry and physical properties
    if (!affine)
      for (int j = 0; j < nb_quad; j++)
	{ 
	  coef = alpha*Fb.WeightsND(j);
	  const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> &
	    dfjm1 = var.Glob_DFjm1(i)(j);
	  
          int offset_elt = offset + nb_unknowns*dim_N*j;
	  for (int m = 0; m < nb_unknowns; m++)
	    {
              ExtractVector(Vh(m), dim_N*j, tmp);
              
	      // we apply DF_i^{*-1}
	      MltTrans(dfjm1, tmp, dU_loc(m));
	      Mlt(coef, dU_loc(m));
              
              AddVector(dU_loc(m), offset_elt+m*Dimension::dim_N, C);
	    }
	}
    else
      {
	const TinyMatrix<Real_wp, General, Dimension::dim_N, Dimension::dim_N> &
	  dfjm1 = var.Glob_DFjm1(i)(0);
	
	for (int j = 0; j < nb_quad; j++)
	  { 
            int offset_elt = offset + nb_unknowns*dim_N*j;
            coef = alpha*Fb.WeightsND(j);            
	    for (int m = 0; m < nb_unknowns; m++)
	      {
                ExtractVector(Vh(m), dim_N*j, tmp);
		
		// we apply DF_i^{*-1}
		MltTrans(dfjm1, tmp, dU_loc(m));
		Mlt(coef, dU_loc(m));
                
                AddVector(dU_loc(m), offset_elt+m*Dimension::dim_N, C);
	      }
	  }
      }
  }
  

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
  template <class T0, class Vector1, class Complexe, class TypeEquation>
  void MltAdd_SquareHex_VectorialH1(const T0& alpha, int level,
				    const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				    const Vector1& B, const T0& beta, Vector1& C)
  {    
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    if (beta == T0(0))
      C.Fill(0);
    else if (beta != T0(1))
      Mlt(beta, C);
    
    typedef typename TypeEquation::Dimension Dimension;    
    typedef TinyVector<Complexe, Dimension::dim_N> R_N_Complexe;
    
    TinyVector<Vector1, TypeEquation::nb_unknowns> Uh, Vh;
    TinyVector<R_N_Complexe, TypeEquation::nb_unknowns> dU_loc;
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      MltAdd_ElementH1_Vectorial(alpha, A, i, B, C, dU_loc,
				 Uh, Vh, A.var.GetReferenceElementH1(i));
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
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
  template <class T0, class Vector1, class Complexe, class TypeEquation>
  void MltAdd_SquareHex_ScalarH1(const T0& alpha, int level,
				 const FemMatrixFreeClass<Complexe, TypeEquation>& A,
				 const Vector1& B, const T0& beta, Vector1& C)
  {
    glob_chrono.Start(VirtualTimer::STIFFNESS);
    if (beta == T0(0))
      C.Fill(0);
    else if (beta != T0(1))
      Mlt(beta, C);
    
    TinyVector<Vector1, TypeEquation::nb_unknowns> Uh, Vh;
    
    for (int i = 0; i < A.var.mesh.GetNbElt(); i++)
      MltAdd_ElementH1_Scalar(alpha, A, i, B, C,
			      Uh, Vh, A.var.GetReferenceElementH1(i));
    
    glob_chrono.Stop(VirtualTimer::STIFFNESS);
    
    if (!A.DirichletDofIgnored())
      A.var.ImposeNullDirichletCondition(C);    
  }  

} // namespace Montjoie

#define MONTJOIE_FILE_PROD_MAT_VECT_H1_CXX
#endif
