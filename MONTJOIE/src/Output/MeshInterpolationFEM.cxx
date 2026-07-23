#ifndef MONTJOIE_FILE_MESH_INTERPOLATION_FEM_CXX

namespace Montjoie
{

  template<class T>
  void MeshInterpolationFEM_Base::SetContribU(const Vector<Vector<TinyVector<T, 1> > >& u,
                                              int nb_points_quad, const Vector<int>& offset, Vector<Vector<T> >& u_boundary)
  {
    for (int j = 0; j < u.GetM(); j++)
      for (int k = 0; k < nb_points_quad; k++)
        u_boundary(offset(j))(k) = u(j)(k)(0);
  }


  template<class T>
  void MeshInterpolationFEM_Base::SetContribU(const Vector<Vector<TinyVector<T, 2> > >& u,
                                              int nb_points_quad, const Vector<int>& offset, Vector<Vector<T> >& u_boundary)
  {
    for (int j = 0; j < u.GetM(); j++)
      for (int k = 0; k < nb_points_quad; k++)
        {
          u_boundary(offset(j))(k) = u(j)(k)(0);
          u_boundary(offset(j)+1)(k) = u(j)(k)(1);
        }
  }


  template<class T>
  void MeshInterpolationFEM_Base::SetContribU(const Vector<Vector<TinyVector<T, 3> > >& u,
                                              int nb_points_quad, const Vector<int>& offset, Vector<Vector<T> >& u_boundary)
  {
    for (int j = 0; j < u.GetM(); j++)
      for (int k = 0; k < nb_points_quad; k++)
        {
          u_boundary(offset(j))(k) = u(j)(k)(0);
          u_boundary(offset(j)+1)(k) = u(j)(k)(1);
          u_boundary(offset(j)+2)(k) = u(j)(k)(2);
        }
  }

    
  //! internal function
  template<>
  void MeshInterpolationFEM<Dimension2>::
  StoreGeometricDataSurface(const DistributedProblem<Dimension2>& var, const Mesh<Dimension2>& mesh,
                            int i1, int num_elem, int num_loc)
  {
    // points and jacobian matrices
    SetPoints<Dimension2> pts;
    SetMatrices<Dimension2> mat;
    VectR2 s;    
    
    // computation of transformations Fi and DFi
    const ElementGeomReference<Dimension2>& Fb = var.GetShapeElement(num_elem);
    mesh.GetVerticesElement(num_elem, s);
    Fb.FjElemNodal(s, pts, mesh, num_elem);
    Fb.DFjElemNodal(s, pts, mat, mesh, num_elem);
    
    // storing DF_i and DF_i^-1 on nodal points of the boundary
    int nb_node = Fb.GetNbNodalBoundary(num_loc);
    this->PointsNodalBoundary(i1).Reallocate(nb_node);
    this->DFiNodalBoundary(i1).Reallocate(nb_node);
    this->invDFiNodalBoundary(i1).Reallocate(nb_node);
    Matrix2_2 mat_dfj, dfjm1;
    for (int j = 0; j < nb_node; j++)
      {
        int node = Fb.GetNodalNumber(num_loc, j);
        mat_dfj = mat.GetPointNodal(node);
        GetInverse(mat_dfj, dfjm1);
	this->PointsNodalBoundary(i1)(j) = pts.GetPointNodal(node);
        this->DFiNodalBoundary(i1)(j) = mat_dfj;
        this->invDFiNodalBoundary(i1)(j) = dfjm1;
      }
    
    // loop on quadrature points
    int nb_points = step_subdiv.GetM();
    this->PointsBoundary(i1).Reallocate(nb_points);
    this->NormaleBoundary(i1).Reallocate(nb_points);
    this->WeightsBoundary(i1).Reallocate(nb_points);
    this->DsBoundary(i1).Reallocate(nb_points);
    this->Nodle(i1).Reallocate(nb_points);
    this->DFiBoundary(i1).Reallocate(nb_points);
    this->invDFiBoundary(i1).Reallocate(nb_points);
    if (i1 == 0)
      this->Nodle(i1)(0) = 0;
    else
      this->Nodle(i1)(0) = this->Nodle(i1-1)(0) + this->Nodle(i1-1).GetM();
    
    R2 point_loc, point_glob, normale_fj; Real_wp dsj;
    for (int j = 0; j < nb_points; j++)
      {
        Fb.GetLocalCoordOnBoundary(num_loc, step_subdiv(j), point_loc);
        
        Fb.Fj(s, pts, point_loc, point_glob, mesh, num_elem);
        Fb.DFj(s, pts, point_loc, mat_dfj, mesh, num_elem);
        GetInverse(mat_dfj, dfjm1);
        Fb.GetNormale(dfjm1, normale_fj, dsj, num_loc);
            	    
        this->Nodle(i1)(j) = this->Nodle(i1)(0) + j;
        this->PointsBoundary(i1)(j) = point_glob;
        this->WeightsBoundary(i1)(j) = dsj*this->weights_surf(0)(j);
        this->DFiBoundary(i1)(j) = mat_dfj;
        this->invDFiBoundary(i1)(j) = dfjm1;
        this->DsBoundary(i1)(j) = dsj;
        this->NormaleBoundary(i1)(j) = normale_fj;
      }
  }

  
  //! extracts a surface mesh and computes normales, jacobian matrices, etc
  template<>
  void MeshInterpolationFEM<Dimension2>::
  ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension2>& mesh,
                     Mesh<Dimension2>& mesh_subdiv,
		     const DistributedProblem<Dimension2>& var, int ref_domain)
  {
    SurfacicMesh<Dimension2> mesh_surf;
    
    // extraction of the surfacic mesh
    mesh.GetBoundaryMesh(1, mesh_surf, ref_cond, ref_domain);
    
    // some local variables
    int num_elem, num_loc;
    
    // subdivision of the initial mesh
    mesh_subdiv = mesh_surf;

    // computation of the normale, projection of points
    int N = mesh_surf.ListeBoundaries.GetM();
    this->NumElement.Reallocate(N);
    this->LocalPositionBoundary.Reallocate(N);        
    this->Nodle.Reallocate(N);
    this->TypeBoundary.Reallocate(N); this->TypeBoundary.Fill(0);
    
    this->NormaleBoundary.Reallocate(N);
    this->PointsBoundary.Reallocate(N);
    this->WeightsBoundary.Reallocate(N);
    this->DsBoundary.Reallocate(N);
    this->PointsNodalBoundary.Reallocate(N);
    this->DFiNodalBoundary.Reallocate(N);
    this->invDFiNodalBoundary.Reallocate(N);
    this->DFiBoundary.Reallocate(N);
    this->invDFiBoundary.Reallocate(N);
    
    // loop on each face
    for (int i1 = 0; i1 < N; i1++)
      {
        //int i = mesh_surf.ListeBoundaries(i1);
        //num_boundary = i;
        num_elem = mesh_surf.NumElement(i1);
        
        this->NumElement(i1) = num_elem;
        num_loc = mesh_surf.NumLocalBoundary(i1);
        this->LocalPositionBoundary(i1) = num_loc;
        
        
        this->StoreGeometricDataSurface(var, mesh, i1, num_elem, num_loc);
      }
  }
  

#ifdef MONTJOIE_WITH_THREE_DIM
  //! internal function
  template<>
  void MeshInterpolationFEM<Dimension3>::
  StoreGeometricDataSurface(const DistributedProblem<Dimension3>& var, const Mesh<Dimension3>& mesh,
                            int i1, int num_elem, int num_loc)
  {
    // points and jacobian matrices
    SetPoints<Dimension3> pts;
    SetMatrices<Dimension3> mat;
    VectR3 s;    
    
    const ElementGeomReference<Dimension3>& Fb = var.GetShapeElement(num_elem);
    mesh.GetVerticesElement(num_elem, s);    
    Fb.FjElemNodal(s, pts, mesh, num_elem);
    Fb.DFjElemNodal(s, pts, mat, mesh, num_elem);
    
    // storing DF_i and DF_i^-1 on nodal points of the boundary
    int nb_node = Fb.GetNbNodalBoundary(num_loc);
    this->PointsNodalBoundary(i1).Reallocate(nb_node);
    this->DFiNodalBoundary(i1).Reallocate(nb_node);
    this->invDFiNodalBoundary(i1).Reallocate(nb_node);
    Matrix3_3 mat_dfj, dfjm1;
    for (int j = 0; j < nb_node; j++)
      {
        int node = mesh.GetNodalNumber(num_elem, num_loc, j);
        mat_dfj = mat.GetPointNodal(node);
        GetInverse(mat_dfj, dfjm1);
        this->PointsNodalBoundary(i1)(j) = pts.GetPointNodal(node);
	this->DFiNodalBoundary(i1)(j) = mat_dfj;
        this->invDFiNodalBoundary(i1)(j) = dfjm1;
      }
    
    // loop on quadrature points
    int nb_points = this->points_div_surf(this->TypeBoundary(i1)).GetM();
    this->PointsBoundary(i1).Reallocate(nb_points);
    this->NormaleBoundary(i1).Reallocate(nb_points);
    this->WeightsBoundary(i1).Reallocate(nb_points);
    this->DsBoundary(i1).Reallocate(nb_points);
    this->Nodle(i1).Reallocate(nb_points);
    this->DFiBoundary(i1).Reallocate(nb_points);
    this->invDFiBoundary(i1).Reallocate(nb_points);
    if (i1 == 0)
      this->Nodle(i1)(0) = 0;
    else
      this->Nodle(i1)(0) = this->Nodle(i1-1)(0) + this->Nodle(i1-1).GetM();
    
    R3 point_loc, point_glob, normale_fj; Real_wp dsj;
    for (int j = 0; j < nb_points; j++)
      {
        Fb.GetLocalCoordOnBoundary(num_loc, this->points_div_surf(this->TypeBoundary(i1))(j), point_loc);
        
        Fb.Fj(s, pts, point_loc, point_glob, mesh, num_elem);
        Fb.DFj(s, pts, point_loc, mat_dfj, mesh, num_elem);
        GetInverse(mat_dfj, dfjm1);
        Fb.GetNormale(dfjm1, normale_fj, dsj, num_loc);
	
        this->Nodle(i1)(j) = this->Nodle(i1)(0) + j;
        this->PointsBoundary(i1)(j) = point_glob;
        this->WeightsBoundary(i1)(j) = dsj*this->weights_surf(this->TypeBoundary(i1))(j);
        this->DsBoundary(i1)(j) = dsj;
        this->DFiBoundary(i1)(j) = mat_dfj;
        this->invDFiBoundary(i1)(j) = dfjm1;
        this->NormaleBoundary(i1)(j) = normale_fj;
      }
  }
  
  
  //! extracts a surface mesh and computes normales, jacobian matrices, etc
  template<>
  void MeshInterpolationFEM<Dimension3>::
  ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension3>& mesh,
		     Mesh<Dimension3>& mesh_subdiv,
		     const DistributedProblem<Dimension3>& var, int ref_domain)
  {
    SurfacicMesh<Dimension3> mesh_surf;
    
    // extraction of the surfacic mesh
    mesh.GetBoundaryMesh(1, mesh_surf, ref_cond, ref_domain);
    
    // some local variables
    int num_elem, num_loc;
    
    // subdivision of the initial mesh
    mesh_subdiv = mesh_surf;

    // computation of the normale, projection of points
    int N = mesh_surf.ListeBoundaries.GetM();
    this->NumElement.Reallocate(N);
    this->LocalPositionBoundary.Reallocate(N);    
    this->NormaleBoundary.Reallocate(N);
    this->Nodle.Reallocate(N);
    this->TypeBoundary.Reallocate(N); this->TypeBoundary.Fill(0);
    this->PointsBoundary.Reallocate(N);
    this->WeightsBoundary.Reallocate(N);
    this->DsBoundary.Reallocate(N);
    this->PointsNodalBoundary.Reallocate(N);
    this->DFiNodalBoundary.Reallocate(N);
    this->invDFiNodalBoundary.Reallocate(N);
    this->DFiBoundary.Reallocate(N);
    this->invDFiBoundary.Reallocate(N);
    
    // loop on each face
    for (int i1 = 0; i1 < N; i1++)
      {
	int i = mesh_surf.ListeBoundaries(i1);
        num_elem = mesh_surf.NumElement(i1);
        this->TypeBoundary(i1) = mesh.BoundaryRef(i).GetHybridType();
        NumElement(i1) = num_elem;
        num_loc = mesh_surf.NumLocalBoundary(i1);
        LocalPositionBoundary(i1) = num_loc;
        
	// computation of transformations Fi and DFi        
        this->StoreGeometricDataSurface(var, mesh, i1, num_elem, num_loc);
      }
  }
#endif


  //! computation of u and du/dn (for Helmholtz) E x n and H x n (for Maxwell's equations)
  /*!
    \param[in] var considered problem
    \param[in] U0 components of u (or E) on degrees of freedom of the mesh
    \param[out] trace_En values of E x n (or u) on quadrature points of the surface
    \param[out] trace_Hn values of H x n (or du/dn) on quadrature points of the surface
   */
  template<class Dimension> template<class Complexe>
  void MeshInterpolationFEM<Dimension>::
  ComputeEnHnOnBoundary(const VarProblem<Dimension>& var,
                        const Vector<Complexe>& U0, Vector<Complexe>& trace_En,
                        Vector<Complexe>& trace_Hn, bool assemble, bool compute_H) const
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    int nb_boundaries_ref = this->NumElement.GetM();
    
    trace_En.Reallocate(this->nb_points_all_quadrature*var.nb_components_en);
    if (compute_H)
      trace_Hn.Reallocate(this->nb_points_all_quadrature*var.nb_components_hn);

    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    
    int nb_mesh_num = var.GetNbMeshNumberings();
    IVect nb_unknowns_per_mesh(nb_mesh_num), type_u(nb_mesh_num);
    IVect nb_comp_mesh(nb_mesh_num), nb_grad_mesh(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      {
        int nb_u = var.GetNbComponentsUnknown(nm);
        int nb_du = var.GetNbComponentsGradient(nm);
        nb_comp_mesh(nm) = nb_u; nb_grad_mesh(nm) = nb_du;
        if (nm == 0)
          type_u(nm) = var.type_element;
        else
          type_u(nm) = var.other_type_element(nm-1);
      }

    nb_unknowns_per_mesh.Zero();
    int nb_comp_all = 0, nb_grad_all = 0, offset = 0;
    Vector<bool> MeshUsed(nb_mesh_num); MeshUsed.Fill(false);
    Vector<int> nb_comp_per_mesh(nb_mesh_num), offset_dof(nb_mesh_num);
    for (int i = 0; i < var.nb_unknowns; i++)
      {
        int nm = var.mesh_num_unknown(i);
        if (!MeshUsed(nm))
          {
            offset_dof(nm) = offset;
            MeshUsed(nm) = true;
          }

        nb_unknowns_per_mesh(nm)++;
        nb_comp_all += nb_comp_mesh(nm);
        nb_grad_all += nb_grad_mesh(nm);
        offset += var.GetNbDofUnknown(i);
      }
    
    Vector<IVect> num_comp_unknown(nb_mesh_num), num_grad_unknown(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      {
        num_comp_unknown(i).Reallocate(nb_unknowns_per_mesh(i)*nb_comp_mesh(i));
        num_grad_unknown(i).Reallocate(nb_unknowns_per_mesh(i)*nb_grad_mesh(i));
      }
    
    nb_unknowns_per_mesh.Zero();
    int offset_comp = 0, offset_grad = 0;
    for (int i = 0; i < var.nb_unknowns; i++)
      {
        int nm = var.mesh_num_unknown(i);
        int p = nb_unknowns_per_mesh(nm);
        nb_unknowns_per_mesh(nm)++;
        num_comp_unknown(nm)(p) = offset_comp;
        num_grad_unknown(nm)(p) = offset_grad;
        offset_comp += nb_comp_mesh(nm);
        offset_grad += nb_grad_mesh(nm);
      }
    
    Vector<Vector<Vector<Complexe> > > uloc(nb_mesh_num);
    for (int nm = 0; nm < nb_mesh_num; nm++)
      uloc(nm).Reallocate(nb_unknowns_per_mesh(nm));

    Vector<Complexe> Un; Vector<Vector<Complexe> > u_boundary, grad_boundary;
    Vector<Vector<TinyVector<Complexe, 1> > > u_h1;
    Vector<Vector<TinyVector<Complexe, Dimension::dim_N> > > grad_h1;

    Vector<Vector<TinyVector<Complexe, Dimension::dim_N> > > u_hdiv;
    Vector<Vector<TinyVector<Complexe, 1> > > div_hdiv;

    enum {nb_u_hcurl = FiniteElementHcurl<Dimension>::nb_components_u,
          nb_du_hcurl = FiniteElementHcurl<Dimension>::nb_components_grad};
    
    Vector<Vector<TinyVector<Complexe, nb_u_hcurl> > > u_hcurl;
    Vector<Vector<TinyVector<Complexe, nb_du_hcurl> > > curl_hcurl;
    
    Vector<Vector<Complexe> > u_quadrature(nb_comp_all), grad_quadrature(nb_grad_all);
    
    bool compute_grad = compute_H;
    if ((var.FormulationDG() != ElementReference_Base::CONTINUOUS) && (var.FirstOrderFormulationDG()))
      compute_grad = false;
    
    Vector<Complexe> En_quad, Hn_quad;
    u_boundary.Reallocate(offset_comp);
    grad_boundary.Reallocate(offset_grad);

    // loop 
    int offset_E = 0, offset_H = 0;
    for (int i = 0; i < nb_boundaries_ref; i++)
      {
	int num_elem = this->NumElement(i);
	int num_loc = this->LocalPositionBoundary(i);
	int num_face = var.mesh.Element(num_elem).numBoundary(num_loc);
        int type_elt = var.mesh.Boundary(num_face).GetHybridType();
        int nb_points_quad = this->GetNbPointsQuadrature(i);
        
        for (int k = 0; k < u_boundary.GetM(); k++)
          u_boundary(k).Reallocate(nb_points_quad);

        for (int k = 0; k < grad_boundary.GetM(); k++)
          grad_boundary(k).Reallocate(nb_points_quad);
        
	for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(nm);
            int size_Un = nb_unknowns_per_mesh(nm)*mesh_num.GetNbDof();
            Un.SetData(size_Un, const_cast<Complexe*>(&U0(offset_dof(nm))));
            mesh_num.number_map.GetLocalUnknownVector(mesh_num, Un, num_elem, uloc(nm));
            Un.Nullify();
          }
        
	// we store DFi and DFi^-1 in MatricesElem
	MatricesElem.SetPointNodalBoundary(this->DFiNodalBoundary(i));
	MatricesElem.SetInversePointNodalBoundary(this->invDFiNodalBoundary(i));
    
	// computation of u and grad(u)
        for (int nm = 0; nm < nb_mesh_num; nm++)
          {
            switch (type_u(nm))
              {
              case 1:
                var.GetReferenceElementH1(num_elem).
                  ComputeValueBoundary(MatricesElem, uloc(nm), u_h1, var.mesh,
                                       num_elem, num_loc);

                this->SetContribU(u_h1, nb_points_quad, num_comp_unknown(nm), u_boundary);
                break;
              case 2:
                var.GetReferenceElementHcurl(num_elem).
                  ComputeValueBoundary(MatricesElem, uloc(nm), u_hcurl, var.mesh,
                                       num_elem, num_loc);
                
                this->SetContribU(u_hcurl, nb_points_quad, num_comp_unknown(nm), u_boundary);
                break;
              case 3:
                var.GetReferenceElementHdiv(num_elem).
                  ComputeValueBoundary(MatricesElem, uloc(nm), u_hdiv, var.mesh,
                                       num_elem, num_loc);

                this->SetContribU(u_hdiv, nb_points_quad, num_comp_unknown(nm), u_boundary);
                break;
              }            
          }
        
	for (int m = 0; m < nb_comp_all; m++)
	  {
	    u_quadrature(m).Reallocate(nb_points_quad);
	    this->ProjectSurface(u_boundary(m), u_quadrature(m), type_elt);
	  }
	
	if (compute_grad)
	  {
            for (int nm = 0; nm < nb_mesh_num; nm++)
              {
                switch (type_u(nm))
                  {
                  case 1:
                    var.GetReferenceElementH1(num_elem).
                      ComputeGradientBoundary(MatricesElem, uloc(nm), grad_h1,
                                              var.mesh, num_elem, num_loc);
                    
                    this->SetContribU(grad_h1, nb_points_quad, num_grad_unknown(nm), grad_boundary);
                    break;
                  case 2:
                    var.GetReferenceElementHcurl(num_elem).
                      ComputeGradientBoundary(MatricesElem, uloc(nm), curl_hcurl,
                                              var.mesh, num_elem, num_loc);
                    
                    this->SetContribU(curl_hcurl, nb_points_quad, num_grad_unknown(nm), grad_boundary);
                    break;
                  case 3:
                    var.GetReferenceElementHdiv(num_elem).
                      ComputeGradientBoundary(MatricesElem, uloc(nm), div_hdiv,
                                              var.mesh, num_elem, num_loc);
                    
                    this->SetContribU(div_hdiv, nb_points_quad, num_grad_unknown(nm), grad_boundary);
                    break;
                  }
              }
            
	    for (int m = 0; m < nb_grad_all; m++)
	      {
		grad_quadrature(m).Reallocate(nb_points_quad);
		this->ProjectSurface(grad_boundary(m), grad_quadrature(m), type_elt);
	      }
	  }
        
	var.ComputeEnHnQuadrature(u_quadrature, grad_quadrature, num_elem,
				  this->PointsQuadrature(i), this->NormaleQuadrature(i),
				  compute_H, En_quad, Hn_quad);
    
	for (int j = 0; j < nb_points_quad*var.nb_components_en; j++)
	  trace_En(offset_E + j) = En_quad(j);
	
	if (compute_H)
	  for (int j = 0; j < nb_points_quad*var.nb_components_hn; j++)
	    trace_Hn(offset_H + j) = Hn_quad(j);
	
        offset_E += nb_points_quad*var.nb_components_en;
        offset_H += nb_points_quad*var.nb_components_hn;
      }
    
    // for parallel execution, we gather values of En/Hn on other processors
    if (nb_proc > 1 && assemble)
      {
        
#ifdef SELDON_WITH_MPI
        Vector<Complexe> EnProc = trace_En;
        Vector<Complexe> HnProc = trace_Hn;
        Vector<int64_t> xtmp;

        int Nall = this->AllPoints.GetM();
        trace_En.Reallocate(Nall*var.nb_components_en);
        if (compute_H)
          trace_Hn.Reallocate(Nall*var.nb_components_hn);
        
        offset_E = this->offset_quadrature_per_proc(this->rank_processor)*var.nb_components_en;
        offset_H = this->offset_quadrature_per_proc(this->rank_processor)*var.nb_components_hn;
        for (int j = 0; j < this->nb_points_quadrature_per_proc(this->rank_processor)*var.nb_components_en; j++)
	  trace_En(offset_E+j) = EnProc(j);
	
	if (compute_H)
	  for (int j = 0; j < this->nb_points_quadrature_per_proc(this->rank_processor)*var.nb_components_hn; j++)
	    trace_Hn(offset_H+j) = HnProc(j);
	
        for (int i = 0; i < nb_proc; i++)
          if (this->nb_points_quadrature_per_proc(i) > 0)
            {
	      MpiBcast(var.comm_group_mode, &trace_En(this->offset_quadrature_per_proc(i)*var.nb_components_en),
		       xtmp, var.nb_components_en*this->nb_points_quadrature_per_proc(i), i);
              
              if (compute_H)
		MpiBcast(var.comm_group_mode, &trace_Hn(this->offset_quadrature_per_proc(i)*var.nb_components_hn),
			 xtmp, var.nb_components_hn*this->nb_points_quadrature_per_proc(i), i);
            }        
#endif
      }
    
  }

  
  //! on calcule rhs = rhs + \int_\Gamma g_source phi ds
  template<class Dimension> template<class T>
  void MeshInterpolationFEM<Dimension>
  ::AddSourceBoundary(const VarProblem<Dimension>& var,
                      const Vector<T>& g_source, Vector<T>& rhs) const
  {
    if (var.GetNbMeshNumberings() > 1)
      {
        cout << "Not implemented for several numberings" << endl;
        abort();
      }
    
    Vector<Vector<T> > contrib(1), scal_g(1);
    int nb = 0;
    for (int i = 0; i < this->GetNbBoundary(); i++)
      {        
	int num_elem = this->GetElementNumberOfSurface(i);
        // we retrieve number of dofs
	IVect Nodle = var.GetDofNumberOnElement(num_elem);
        
        int num_loc = this->GetLocalPositionOfSurface(i);
        const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0); 
	int nb_dof_loc = mesh_num.GetNbLocalDof(num_elem);
	int nb_points_quadrature = this->GetNbPointsQuadrature(i);
        scal_g(0).Reallocate(nb_points_quadrature);
	for (int j = 0; j < nb_points_quadrature; j++)
	  {
            scal_g(0)(j) = g_source(nb);
            nb++;
	  }

        int type_elt = var.type_element;
	switch (type_elt)
          {
          case 1:
            var.GetReferenceElementH1(num_elem).
              ComputeGaussIntegralSurface(this->Dfjm1Quadrature(i),
                                          this->WeightsQuadrature(i),
                                          this->DsQuadrature(i),
                                          scal_g, contrib, mesh_num, num_elem, num_loc);
            break;
          case 2:
            var.GetReferenceElementHcurl(num_elem).
              ComputeGaussIntegralSurface(this->Dfjm1Quadrature(i),
                                          this->WeightsQuadrature(i),
                                          this->DsQuadrature(i),
                                          scal_g, contrib, mesh_num, num_elem, num_loc);
            break;
          case 3:
            var.GetReferenceElementHdiv(num_elem).
              ComputeGaussIntegralSurface(this->Dfjm1Quadrature(i),
                                          this->WeightsQuadrature(i),
                                          this->DsQuadrature(i),
                                          scal_g, contrib, mesh_num, num_elem, num_loc);
            break;
          }
        
        for (int j = 0; j < nb_dof_loc; j++)
          {
            int num_dof_loc = Nodle(j);
            int num_dof = num_dof_loc;
            if (num_dof_loc >= 0)
              rhs(num_dof) += contrib(0)(j);
          }
      }
  }

} // namespace Montjoie

#define MONTJOIE_FILE_MESH_INTERPOLATION_FEM_CXX
#endif
