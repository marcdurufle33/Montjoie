#ifndef MONTJOIE_FILE_VAR_GEOMETRY_PROBLEM_CXX

namespace Montjoie
{
  
  /**********************
   * VarGeometryProblem *
   **********************/
  
  
  //! Sets attributes to default values
  template<class Dimension>
  void VarGeometryProblem<Dimension>::InitDefaultValues()
  {
    write_quadrature_points = false;
    write_quad_points_pml = false;
    
    kwave.Fill(0);
    kwave(0) = this->omega;
    reference_infinity = -1;
    
    this->all_mesh_num.Reallocate(1);
    this->all_mesh_num(0) = &this->mesh_num;
  }


  //! Destructor
  template<class Dimension>
  VarGeometryProblem<Dimension>::~VarGeometryProblem()
  {
    if (this->other_mesh_num.GetM() > 0)
      for (int i = 0; i < this->other_mesh_num.GetM(); i++)
        delete this->other_mesh_num(i);

    this->other_mesh_num.Clear();
  }

  
  //! returns the reference of the infinite homogeneous media
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetReferenceInfinity() const
  {
    return reference_infinity;
  }
  
  
  //! updates the wave vector is omega has changed
  template<class Dimension>
  void VarGeometryProblem<Dimension>::UpdateWaveVector()
  {
    Real_wp coef = this->omega / Norm2(this->kwave);
    kwave *= coef;
  }


  //! Inits polarization vector
  template<class Dimension>
  void VarGeometryProblem<Dimension>::InitPolarization()
  {
    this->polarization.Reallocate(this->GetNbComponentsAll());
    this->polarization.Fill(0);
    if (this->GetNbComponentsAll() <= 0)
      {
        cout << "Null number of components" << endl;
        abort();
      }

    this->polarization(0) = 1.0;

    this->polarization_grad.Reallocate(this->GetNbComponentsGradientAll());
    this->polarization_grad.Zero();    
  }

  
  //! Fills kx, ky from wave vector
  template<>
  void VarGeometryProblem<Dimension2>::FillWaveVectorComponents(Real_wp& kx, Real_wp& ky, Real_wp& kz)
  {
    kx = this->kwave(0);
    ky = this->kwave(1);
  }
  
  
  //! Fills kx, ky from wave vector
  template<>
  void VarGeometryProblem<Dimension3>::FillWaveVectorComponents(Real_wp& kx, Real_wp& ky, Real_wp& kz)
  {
    kx = this->kwave(0);
    ky = this->kwave(1);
    kz = this->kwave(2);
  }
  

  //! returns the number of components for the given type
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsType(int type) const
  {
    switch (type)
      {
      case 1: return 1;
      case 2: return FiniteElementHcurl<Dimension>::nb_components_u;
      case 3: return Dimension::dim_N;
      }

    return 0;
  }


  //! returns the number of components for gradient for the given type
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsGradType(int type) const
  {
    switch (type)
      {
      case 1: return Dimension::dim_N;
      case 2: return FiniteElementHcurl<Dimension>::nb_components_grad;
      case 3: return 1;
      }

    return 0;
  }

  //! returns the number of components for the mesh numbering n
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsUnknown(int n) const
  {
    if (n == 0)
      return GetNbComponentsType(this->type_element);

    return this->GetNbComponentsType(this->other_type_element(n-1));
  }


  //! returns the number of components (for the gradient) for the mesh numbering n
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsGradient(int n) const
  {
    if (n == 0)
      return GetNbComponentsGradType(this->type_element);

    return this->GetNbComponentsGradType(this->other_type_element(n-1));
  }


  //! returns the total number of components
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsAll(int nb_vec) const
  {
    if (nb_vec <= 0)
      nb_vec = this->nb_unknowns;
    
    int nb_comp = 0;
    for (int n = 0; n < nb_vec; n++)
      {
        int nm = this->mesh_num_unknown(n%this->nb_unknowns);
        nb_comp += this->GetNbComponentsUnknown(nm);
      }

    return nb_comp;
  }

  
  //! returns the total number of components (for derivatives of unknowns)
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsGradientAll(int nb_vec) const
  {
    if (nb_vec <= 0)
      nb_vec = this->nb_unknowns;
    
    int nb_comp = 0;
    for (int n = 0; n < nb_vec; n++)
      {
        int nm = this->mesh_num_unknown(n%this->nb_unknowns);
        nb_comp += this->GetNbComponentsGradient(nm);
      }
    
    return nb_comp;
  }


  //! returns the total number of components (for hessian of unknowns)
  template<class Dimension>
  int VarGeometryProblem<Dimension>::GetNbComponentsHessianAll(int nb_vec) const
  {
    // assuming that there are only scalar unknowns
    if (nb_vec <= 0)
      nb_vec = this->nb_unknowns;
    
    int nb_comp = nb_vec*((Dimension::dim_N+1)*Dimension::dim_N)/2;
    return nb_comp;
  }


  //! returns the n-th mesh numbering
  template<class Dimension>
  const MeshNumbering<Dimension>& VarGeometryProblem<Dimension>::GetMeshNumbering(int n) const
  {
    if (n == 0)
      return this->mesh_num;
    else
      return *this->other_mesh_num(n-1);
  }


  //! returns the n-th mesh numbering
  template<class Dimension>
  MeshNumbering<Dimension>& VarGeometryProblem<Dimension>::GetMeshNumbering(int n)
  {
    if (n == 0)
      return this->mesh_num;
    else
      return *this->other_mesh_num(n-1);
  }


  //! modifies the object with a line of the data file
  template<class Dimension>
  void VarGeometryProblem<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarProblem_Base::SetInputData(description_field, parameters);
    mesh.SetInputData(description_field, parameters);
    mesh_num.SetInputData(description_field, parameters);
    for (int i = 0; i < this->other_mesh_num.GetM(); i++)
      this->other_mesh_num(i)->SetInputData(description_field, parameters);

    if ((description_field.size() > 4) && (description_field[0] == 'U'))
      {
        // case of a keyword specific to a numbering
        size_t pos = description_field.find("_");
        bool is_number_provided = true;
        if ((pos <= 1) || (pos == string::npos))
          is_number_provided = false;

        if (is_number_provided)
          for (size_t k = 1; k < pos; k++)
            if (!isdigit(description_field[k]))
              is_number_provided = false;
        
        if (is_number_provided)
          {
            int num = to_num<int>(description_field.substr(1, pos-1));
            string keyword = description_field.substr(pos+1);
            if (num <= this->other_mesh_num.GetM())
              {
                if (num == 0)
                  mesh_num.SetInputData(keyword, parameters);
                else
                  {
                    if (keyword == "TypeElement")
                      this->name_other_elements(num-1) = parameters(0);
                    else if (keyword == "Discontinuous")
                      {
                        if (parameters(0) == "YES")
                          this->other_dg_formulation(num-1) = ElementReference_Base::DISCONTINUOUS;
                      }
                    else
                      this->other_mesh_num(num-1)->SetInputData(keyword, parameters);
                  }
              }
            else
              {
                cout << "Unknown number = " << num << " whereas the mesh contains " << all_mesh_num.GetM()+1 << " numberings " << endl;
                cout << "Provide NumberAdditionalNumberings first with a correct number" << endl;
                abort();
              }
          }
      }
    
    if (!description_field.compare("IncidentAngle"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarHarmonic" << endl;
            cout << "IncidentAngle needs more parameters, for instance :" << endl;
            cout << "IncidentAngle = teta" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        // in 2-D teta is the incident angle of the plane wave
        // in 3-D, the wave vector is oriented along
        // (sin(teta)*cos(phi) ; sin(teta)*sin(phi) ; cos(teta))
        Real_wp teta = 0.0, phi = 0.0;
        if (Dimension::dim_N >= 2)
          teta = to_num<Real_wp>(parameters(0))*pi_wp/180;
        
        if (Dimension::dim_N == 3)
          if (parameters.GetM() > 1)
            phi = to_num<Real_wp>(parameters(1))*pi_wp/180;
        
        Real_wp c0 = this->GetVelocityOfInfinity();
        SetIncidentAngle(this->omega/c0, this->kwave, teta, phi);
      }
    else if (description_field == "NumberAdditionalNumberings")
      {
        if (this->other_mesh_num.GetM() > 0)
          for (int i = 0; i < this->other_mesh_num.GetM(); i++)
            delete this->other_mesh_num(i);
        
        // how many other numberings must be constructed ?
        int p = to_num<int>(parameters(0));
        this->other_mesh_num.Reallocate(p);
        this->name_other_elements.Reallocate(p);
        this->other_dg_formulation.Reallocate(p);
        this->other_dg_formulation.Fill(this->dg_formulation);
        for (int i = 0; i < p; i++)
          this->other_mesh_num(i) = new MeshNumbering<Dimension>(mesh);

        this->all_mesh_num.Reallocate(1+p);
        this->all_mesh_num(0) = &this->mesh_num;
        for (int i = 0; i < p; i++)
          this->all_mesh_num(i+1) = this->other_mesh_num(i);
      }
    else if (!description_field.compare("WriteQuadraturePoints"))
      {
        if (parameters(0) == "YES")
          {
            this->write_quadrature_points = true;
            if (parameters.GetM() > 1 && parameters(1) == "NO_PML_POINTS")
              this->write_quad_points_pml = false;
            else
              this->write_quad_points_pml = true;
          }
        else
          this->write_quadrature_points = false;
      }
    else if (!description_field.compare("ReferenceInfinity"))
      {
        reference_infinity = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("Polarization"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarHarmonic" << endl;
            cout << "Polarization needs more parameters, for instance :" << endl;
            cout << "Polarization = Ex Ey Ez" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        this->polarization.Zero();
        for (int i = 0; i < min(polarization.GetM(), parameters.GetM()); i++)
          this->polarization(i) = to_num<Real_wp>(parameters(i));
      }
    else if (!description_field.compare("OriginePhase"))
      {
        if (parameters.GetM() < Dimension::dim_N)
          {
            cout << "In SetInputData of VarHarmonic" << endl;
            cout << "OriginePhase needs more parameters, for instance :" << endl;
            cout << "OriginePhase = xO yO zO" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        for (int k = 0; k < this->origine_phase.GetM(); k++)
          to_num(parameters(k), this->origine_phase(k));
      }
    else if (!description_field.compare("UseWarburtonTrick"))
      {
        if (parameters.GetM() <= 0)
          {
            cout << "In SetInputData of VarHarmonic" << endl;
            cout << "UseWarburtonTrick needs more parameters, for instance :" << endl;
            cout << "UseWarburtonTrick = YES" << endl;
            cout << "Current parameters are : " << endl << parameters << endl;
            abort();
          }
        
        ElementReference_Base::use_warburton_trick = false;
        if (!parameters(0).compare("YES"))
          ElementReference_Base::use_warburton_trick = true;
      }
  }
  
  
  //! returns true if the face num_face is involved in a boundary integral
  template<class Dimension>
  bool VarGeometryProblem<Dimension>::
  FaceHasToBeConsideredForBoundaryIntegral(int num_face) const
  {
    if (this->FormulationDG() == ElementReference_Base::CONTINUOUS)
      {
        if (num_face >= this->mesh.GetNbBoundaryRef())
          return false;
        
        if (this->mesh_num.GetFormulationForPeriodicCondition() != this->mesh_num.WEAK_PERIODIC)
          return false;
        
        if (this->mesh_num.GetPeriodicBoundary(num_face) < 0)
          return false;
      }
    
    return true;
  }
  
  
  //! adds memory used by the current object
  template<class Dimension>
  void VarGeometryProblem<Dimension>::GetMemoryUsed(map<string, size_t>& var) const
  {
    VarProblem_Base::GetMemoryUsed(var);
    
    size_t n = mesh.GetMemorySize();
    var["Mesh"] = n;
    
      n = mesh_num.GetMemorySize();
      var["Mesh Numbering"] = n;

      n = Seldon::GetMemorySize(Glob_PointsQuadrature);
      var["Quadrature Points"] = n;

      // geometric quantities
      n = Seldon::GetMemorySize(Glob_jacobian) + Seldon::GetMemorySize(Glob_decomp_jacobian)
        + Seldon::GetMemorySize(Glob_normale) + Seldon::GetMemorySize(Glob_dsj)
        + OrthogonalElement.GetMemorySize() + Seldon::GetMemorySize(Glob_DFjm1)
        + Seldon::GetMemorySize(IsNewFace) + polarization.GetMemorySize()
        + this->Glob_CoefPenalDG.GetMemorySize() + sizeof(*this);

      var["JacobianMatrix"] = n;
    }


  //! writes nodal points of the mesh in a file understood by write_index
  template<class Dimension>
  void VarGeometryProblem<Dimension>::WriteNodalPointsMesh()
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    if (nb_proc == 1)
      {
        Vector<typename Dimension::R_N> PosNodes; Vector<IVect> Nodle;
        this->mesh.GetNodesCurvedMesh(PosNodes, Nodle, this->mesh.GetGeometryOrder(), true);
        
        ofstream file_out("nodal_points.dat");
        file_out.precision(15);
        file_out << this->mesh.GetNbElt() << " " << PosNodes.GetM() << '\n';
        
        // on ecrit les points
        for (int i = 0; i < PosNodes.GetM(); i++)          
          {
            PrintNoBrackets(file_out, PosNodes(i));
            file_out << '\n';
          }
        
        // on ecrit la numerotation pour chaque element
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          {
            file_out << Nodle(i).GetM() << '\n';
            int n = -1;
            for (int j = 0; j < Nodle(i).GetM()-1; j++)
              {
                n = Nodle(i)(j);
                file_out << n << " ";
              }
            
            n = Nodle(i)(Nodle(i).GetM()-1);
            file_out << n << '\n';
          }
        
        file_out.close();
      }
  }
  
  
  //! localize points of the boundary on an element
  template<class Dimension>
  void VarGeometryProblem<Dimension>
  ::LocalizePointsBoundaryElement(const VectR_N& Points, int num_elem, int num_loc,
                                  FjInverseProblem<Dimension>& inverseFj,
                                  typename DimensionBoundary::VectR_N& pts, IVect& num) const
  {
    TinyVector<R_N, 2> enveloppe; R_N point_loc;
    typename DimensionBoundary::R_N xloc;
    
    VectR_N& s = inverseFj.GetVertices();
    SetPoints<Dimension>& PointsElem = inverseFj.GetSetPoints();
    
    // bounding rectangle
    mesh.GetBoundingBox(num_elem, s, PointsElem, enveloppe);
    const ElementGeomReference<Dimension>& Fb_geom = mesh.GetReferenceElement(num_elem);
    
    // loop over points
    pts.Clear(); num.Clear();
    for (int k = 0; k < Points.GetM(); k++)
      if (PointInsideBoundingBox(Points(k), enveloppe))
        {	  
          bool test_inside = inverseFj.Solve(Points(k), point_loc);
          
          if (test_inside)
            {
              Real_wp distance_boundary_elt
                = Fb_geom.GetDistanceToBoundary(point_loc, num_loc);
              
              if (abs(distance_boundary_elt) <= R_N::threshold)
                {
                  Fb_geom.ProjectPointOnBoundary(point_loc);
                  Fb_geom.GetLocalCoordOnBoundary(num_loc, point_loc, xloc);
                  for (int k2 = 0; k2 < pts.GetM(); k2++)
                    if (Norm2(xloc-pts(k2)) <= 1e-10)
                      {
                        cout << "Points should be unique" << endl;
                        abort();
                      }
                  
                  pts.PushBack(xloc);
                  num.PushBack(k);
                }
            }
        }
  }
  
  
  //! returns true if a numerical integration is necessary for element i
  template<class Dimension>
  bool VarGeometryProblem<Dimension>::UseNumericalIntegration(int i) const
  { 
    // numerical integration for non-affine elements
    // (curved or quadrilaterals for example)
    if (this->mesh.Element(i).IsCurved())
      return true;
    
    if (!this->mesh.IsElementAffine(i))
      return true;
    
    // for PML elements
    if (this->InsidePML(i))
      return true;
    
    // or varying media
    int ref = this->mesh.Element(i).GetReference();
    return this->IsVaryingMedia(ref);
  }
  
  
  //! returns the jacobian multiplied by the weight
  template<class Dimension>
  Real_wp VarGeometryProblem<Dimension>
  ::GetWeightedJacobian(int num_elem, int j, bool affine,
                        const ElementGeomReference<Dimension>& Fb) const
  {
    if (affine)
      return this->Glob_jacobian(num_elem)(0)*Fb.WeightsND(j);
    else
      return this->Glob_jacobian(num_elem)(j);
  }
  
  template<class Dimension>
  Real_wp VarGeometryProblem<Dimension>
  ::GetSurfaceWeightedJacobian(int num_elem, int num_loc, int num_edge, int j,
                               const ElementGeomReference<Dimension>& Fb) const
  {
    if (this->IsNewFace(num_elem)(num_loc))
      {
        return Fb.WeightsQuadratureBoundary(j, num_loc)*this->Glob_dsj(num_edge)(j);
      }
    else
      {
        int num_elem2 = this->mesh.Boundary(num_edge).numElement(0);
        if (num_elem2 == num_elem)
          num_elem2 = this->mesh.Boundary(num_edge).numElement(1);
        
        int rot1 = mesh.Element(num_elem).GetOrientationBoundary(num_loc);
        int num_pos2_face = mesh.Element(num_elem2).GetPositionBoundary(num_edge);
        int rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
        int nv = mesh.Boundary(num_edge).GetNbVertices();          
        int rot = this->mesh_num.GetRotationFace(rot1, rot2, nv);

        int rf = this->mesh_num.GetOrderQuadrature(num_edge);
        const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
          GetRotationQuadraturePoints(rf, mesh.Boundary(num_edge));
        
        int jrot = FacesQuadRotation(rot, j);
        return Fb.WeightsQuadratureBoundary(j, num_loc)*this->Glob_dsj(num_edge)(jrot);
      }
  }
  
  template<class Dimension>
  void VarGeometryProblem<Dimension>
  ::GetInverseJacobianMatrix(int num_elem, int j, bool affine,
                             const ElementGeomReference<Dimension>& Fb,
                             MatrixN_N& dfjm1, Real_wp& jacob, Real_wp& jacob_weighted) const
  {
    if (affine)
      {
        jacob = this->Glob_jacobian(num_elem)(0);
        dfjm1 = this->Glob_DFjm1(num_elem)(0);                
        Mlt(1.0/jacob, dfjm1);
        jacob_weighted = jacob*Fb.WeightsND(j);
      }
    else
      {
        jacob = this->Glob_jacobian(num_elem)(j)/Fb.WeightsND(j);
        dfjm1 = this->Glob_DFjm1(num_elem)(j);
        Mlt(1.0/jacob, dfjm1);
        jacob_weighted = this->Glob_jacobian(num_elem)(j);
      }
  }
  
  
  //! fills MatricesElem with Glob_DFjm1
  template<class Dimension>
  void VarGeometryProblem<Dimension>
  ::FillQuadratureJacobian(int i, int N, const ElementGeomReference<Dimension>& Fb,
                           const Vector<R_N>& s, SetPoints<Dimension>& PointsElem,
                           SetMatrices<Dimension>& MatricesElem) const
  {
    if (this->Glob_DFjm1.GetM() > 0)
      {
        MatricesElem.ReallocatePointsQuadrature(N);
        MatrixN_N dfjm1, mat_dfj;
        if (this->Glob_DFjm1(i).GetM() == 1)
          {
            // affine case
            dfjm1 = this->Glob_DFjm1(i)(0);
            Mlt(1.0/this->Glob_jacobian(i)(0), dfjm1);
            GetInverse(dfjm1, mat_dfj);
            MatricesElem.FillQuadrature(mat_dfj);
          }
        else
          {
            for (int j = 0; j < N; j++)
              {
                dfjm1 = this->Glob_DFjm1(i)(j);
                Mlt(Fb.WeightsND(j)/this->Glob_jacobian(i)(j), dfjm1);
                GetInverse(dfjm1, mat_dfj);
                MatricesElem.SetPointQuadrature(j, mat_dfj);
              }
          }
      }
    else
      {
        Fb.FjElemQuadrature(s, PointsElem, this->mesh, i);
        Fb.DFjElemQuadrature(s, PointsElem, MatricesElem, this->mesh, i);
      }
  }
  
  
  //! computation of geometrical quantities such as DF_i for an element of the mesh  
  /*!
i : element number
nb_points_quadrature : number of quadrature points inside the element
*/
  template<class Dimension>
  void VarGeometryProblem<Dimension>::
  ComputeLocalMassMatrix(int i, int nb_points_quadrature, bool linear_sparse_mass,
                         SetPoints<Dimension>& PointsElem,
                         SetMatrices<Dimension>& MatricesElem,
                         IVect& OrderFace, const ElementGeomReference<Dimension>& Fb)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    int nb_elt = this->mesh.GetNbElt();
    Real_wp jacob(0);    
    VectR_N s; MatrixN_N dfjm1, dfj;
    
    // number of nodal points
    int nb_points_nodal = Fb.GetNbPointsNodalElt();
    
    // vertices of the element i
    this->mesh.GetVerticesElement(i, s);
    
    // if the mass matrix for "straight elements" (=non-curved) is linear
    // the jacobian can be decomposed (e.g. in Q1 basis for hexahedron) 
    if (linear_sparse_mass && ! this->mesh.Element(i).IsCurved())
      Fb.ComputeCoefJacobian(s, Glob_decomp_jacobian(i));
    
    
    /***************************************************
     * Computation of geometric coefficient : Ji, DF_i *
     ***************************************************/
    
    
    bool affine = this->mesh.IsElementAffine(i);
    if (true)
      {
        // we compute F_i(\xi_j) where \xi_j are quadrature points on unit element
        Fb.FjElem(s, PointsElem, this->mesh, i);
        // we compute DF_i(\xi_j) where \xi_j are quadrature points on unit element
        Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, i);
        
        // allocation of arrays
        if (affine)
          Glob_jacobian(i).Reallocate(1);
        else
          Glob_jacobian(i).Reallocate(nb_points_quadrature);

        if (affine)
          Glob_DFjm1(i).Reallocate(1);
        else
          Glob_DFjm1(i).Reallocate(nb_points_quadrature);

        // loop on quadrature points
        for (int j = 0; j < Glob_jacobian(i).GetM(); j++)
        {
          // jacobian = J_i = det(DF_i)
          jacob = Det(MatricesElem.GetPointQuadrature(j));

          if (affine)
            Glob_jacobian(i)(j) = jacob; 
          else
            Glob_jacobian(i)(j) = jacob*Fb.WeightsND(j); 
          
          // detection of degenerated elements
          if (jacob < 0)
            {
              DISP(i); DISP(j); DISP(PointsElem.GetPointQuadrature(j));
              DISP(Fb.PointsND(j)); DISP(Fb.GetNbPointsQuadrature()); DISP(jacob);
              DISP(s); 
              for (int ne_loc = 0; ne_loc < Fb.GetNbEdges(); ne_loc++)
                {
                  int ne = this->mesh.Element(i).numEdge(ne_loc);                    
                  DISP(ne_loc); DISP(ne);DISP(this->mesh.GetEdge(ne).GetReference());
                  if (this->mesh.GetEdge(ne).GetReference() != 0)
                    for (int k = 0; k < Fb.GetGeometryOrder()-1; k++)
                      {
                        DISP(this->mesh.GetPointInsideEdge(ne, k));
                        DISP(Norm2(this->mesh.GetPointInsideEdge(ne, k)));
                      }
                }
              
              for (int num_loc = 0; num_loc < Fb.GetNbBoundaries(); num_loc++)
                {
                  int nf = this->mesh.Element(i).numBoundary(num_loc);                    
                  DISP(num_loc); DISP(nf);DISP(this->mesh.Boundary(nf).GetReference());
                  //for (int k = 0; k < this->mesh.GetNbPointsInsideFace(nf); k++)
                  //DISP(this->mesh.GetPointInsideFace(nf, k));
                }
              
              cout<<"Negative jacobian, probably a degenerated element "<<endl;
              
              WriteElementMesh(this->mesh, Fb, PointsElem, 
                               string("elem")+to_str(rank_proc)+".mesh", i);
              
              //this->mesh.Write(string("test")+to_str(rank_proc)+".mesh");
              cout << "waiting"  << endl; int test_input; cin >> test_input;
              break;
              //abort();
            }            
        }
        
        // we compute and store J_i DF_i^{-1} in Glob_DFjm1
        for (int j = 0; j < Glob_DFjm1(i).GetM(); j++)
          {
            jacob = Det(MatricesElem.GetPointQuadrature(j));
            GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);
            Glob_DFjm1(i)(j) = dfjm1; Mlt(jacob, Glob_DFjm1(i)(j));
          }
        
        // coefficient r(r+1) mu / (2 h_min) for the penalty term
        // the coefficient is multiplied by mu (stiffness coefficient)
        // the factor 1/2 is already contained in PoidsFlux
        // h_min is the minimal size of the element
        int r = this->mesh_num.GetOrderElement(i);
        Real_wp mu_max
          = this->GetCoefficientPenaltyStiffness(this->mesh.Element(i).GetReference());
        
        Real_wp h_min = Fb.GetMinimalSize(s);
        Glob_CoefPenalDG(i) = r*(r+1)*mu_max/h_min;
      }
    
    
    /*****************************
     * Storing quadrature points *
     *****************************/
    
    
    int nb_boundaries = Fb.GetNbBoundaries();
    R_N normale, pt_loc, pt_glob; Real_wp dsj;    
    int nb_pts_quad_faces = 0;
    for (int j = 0; j < nb_boundaries; j++)
      nb_pts_quad_faces += this->mesh_num
        .GetNbPointsQuadratureBoundary(this->mesh.Element(i).numBoundary(j));
    
    bool second_order_dg = !this->FirstOrderFormulationDG();
    int nb_points_all = nb_points_quadrature+nb_pts_quad_faces;
    this->Glob_PointsQuadrature(i).Reallocate(nb_points_all);
    if ((!affine) && (second_order_dg))
      {
        this->Glob_jacobian(i).Resize(nb_points_all);
        this->Glob_DFjm1(i).Resize(nb_points_all);
      }
    
    VectR_N AllPoints(nb_points_all);
    for (int j = 0; j < nb_points_quadrature; j++)
      {
        AllPoints(j) = PointsElem.GetPointQuadrature(j);
        this->Glob_PointsQuadrature(i)(j) = PointsElem.GetPointQuadrature(j);
      }
    
    int num_point = nb_points_quadrature;
    Vector<Vector<MatrixN_N> > MatDFjBoundary(nb_boundaries);
    // for PML, we need to know tau and primitive of tau on quadrature points of the faces
    for (int j = 0; j < nb_boundaries; j++)
      {
        int num_face = this->mesh.Element(i).numBoundary(j);
        int rf = this->mesh_num.GetOrderQuadrature(num_face);
        int nb_pts_face = this->mesh_num.GetNbPointsQuadratureBoundary(num_face);
        
        MatDFjBoundary(j).Reallocate(nb_pts_face);
        for (int k = 0; k < nb_pts_face; k++)
          {
            Fb.GetLocalCoordOnBoundary(j, this->mesh_num.number_map.
                                       GetQuadraturePoint(rf, k, this->mesh.Boundary(num_face)),
                                       pt_loc);
            
            Fb.Fj(s, PointsElem, pt_loc, pt_glob, this->mesh, i);
            Fb.DFj(s, PointsElem, pt_loc, dfjm1, this->mesh, i);
            MatDFjBoundary(j)(k) = dfjm1;
            
            AllPoints(num_point) = pt_glob;
            this->Glob_PointsQuadrature(i)(num_point) = pt_glob;
            
            if ((!affine) && (second_order_dg))
              {
                Fb.DFj(s, PointsElem, pt_loc, dfj, this->mesh, i);
                jacob = Det(dfj); GetInverse(dfj, dfjm1);
                this->Glob_jacobian(i)(num_point) = jacob;
                this->Glob_DFjm1(i)(num_point) = dfjm1;
              }
            
            num_point++;
          }
      }
    
    
    /***********************************
     * Computation of PML coefficients *
     ***********************************/
    
    int i1 = i - nb_elt + var_boundary.GetNbEltPML();
    if (this->InsidePML(i))
      this->ComputeStoreCoefficientsPML(i1, i, AllPoints);
    
    
    /*********************************************
     * Computation of Normales on the boundaries *
     *********************************************/
    
    
    IsNewFace(i).Reallocate(nb_boundaries);
    IsNewFace(i).Fill(false);
    
    //R_N sum; sum.Fill(0);
    for (int j = 0; j < nb_boundaries; j++)
      {
        int num_face = this->mesh.Element(i).numBoundary(j);
        int rf = this->mesh_num.GetOrderQuadrature(num_face);
        if (rf > OrderFace(num_face))
          {
            this->IsNewFace(i)(j) = true;
            OrderFace(num_face) = rf;
            int nb_pts_face = this->mesh_num.GetNbPointsQuadratureBoundary(num_face);
            Glob_normale(num_face).Reallocate(nb_pts_face);
            Glob_dsj(num_face).Reallocate(nb_pts_face);
            
            for (int k = 0; k < nb_pts_face; k++)
              {
                GetInverse(MatDFjBoundary(j)(k), dfjm1);                
                Fb.GetNormale(dfjm1, normale, dsj, j);
                
                Glob_normale(num_face)(k) = normale;
                Glob_dsj(num_face)(k) = dsj;
              }            
          }
      }
    
    
    /*******************
     * Warburton trick *
     *******************/
    
    
    if ((ElementReference_Base::use_warburton_trick) && (!affine))
      {
        // phi_i is replaced by phi_i/sqrt(J_i)
        // then, we need to store 1/sqrt(J_i) and grad(J_i)
        int nb = nb_points_nodal + nb_points_quadrature + nb_pts_quad_faces;
        this->mesh.Glob_invSqrtJacobian(i).Reallocate(nb);
        
        // in this case, we store 0.5 grad(J_i)/J_i and 1/sqrt(J_i)
        // 1/sqrt(J_i) is computed for quadrature points of the boundaries
        // and on quadrature and nodal points
        VectReal_wp JacobNodal(nb_points_nodal);
        for(int j = 0; j < nb_points_nodal; j++)
          {  
            jacob = Det(MatricesElem.GetPointNodal(j));
            JacobNodal(j) = jacob;
            this->mesh.Glob_invSqrtJacobian(i)(j) = 1.0/sqrt(jacob);
            //file_out << this->mesh.Glob_invSqrtJacobian(i)(j) << endl;
          }
        
        VectReal_wp JacobDof(Fb.GetNbPointsDof());
        for (int j = 0; j < Fb.GetNbPointsDof(); j++)
          {
            jacob = Det(MatricesElem.GetPointDof(j));
            JacobDof(j) = jacob;
          }
        
        VectReal_wp JacobQuadrature(nb_points_quadrature);
        for (int j = 0; j < nb_points_quadrature; j++)
          {
            jacob = Det(MatricesElem.GetPointQuadrature(j));
            JacobQuadrature(j) = jacob;
            this->mesh.Glob_invSqrtJacobian(i)(nb_points_nodal+j) = 1.0/sqrt(jacob);
          }
        
        nb = nb_points_nodal + nb_points_quadrature;
        for (int num_loc = 0; num_loc < nb_boundaries; num_loc++)
          for(int j = 0; j < MatDFjBoundary(num_loc).GetM(); j++)
            {  
              // jacobian = J_i = det(DF_i)
              jacob = Det(MatDFjBoundary(num_loc)(j));
              this->mesh.Glob_invSqrtJacobian(i)(nb) = 1.0/sqrt(jacob);
              nb++;
            }
        
        VectReal_wp GradJacobQuad(Dimension::dim_N*nb_points_quadrature);
        this->ComputeReferenceGradientElement(i, JacobDof, JacobQuadrature, GradJacobQuad);
        
        // we store 1/(2 Ji) grad(Ji) since that term appears when
        // derivating 1/sqrt(Ji)
        R_N grad_hat;
        this->mesh.Glob_GradJacobian(i).Reallocate(nb_points_quadrature);
        for (int j = 0; j < nb_points_quadrature; j++)
          {
            for (int k = 0; k < Dimension::dim_N; k++)
              grad_hat(k) = GradJacobQuad(j*Dimension::dim_N+k);
            
            GetInverse(MatricesElem.GetPointQuadrature(j), dfjm1);
            MltTrans(dfjm1, grad_hat, this->mesh.Glob_GradJacobian(i)(j));
            jacob = JacobQuadrature(j);
            this->mesh.Glob_GradJacobian(i)(j) *= 0.5/jacob;
          } 
      }

    
    /********************************
     * Checking boundary conditions *
     ********************************/
    
    
    // on regarde si sur chaque arete, on a bien une condition aux limites
    for (int j = 0; j < nb_boundaries; j++)
      {
        int num_face = this->mesh.Element(i).numBoundary(j);
        int ref = this->mesh.Boundary(num_face).GetReference();
        if (this->mesh.GetBoundaryCondition(ref) != BoundaryConditionEnum::LINE_INSIDE)
          {
            // arete ou face avec condition aux limites
          }
        else if (this->mesh.Boundary(num_face).GetNbElements() == 2)
          {
            // arete ou face completement interieure (deux elements de part et d'autre)
          }
        else
          {
            if (this->exit_if_no_boundary_condition)
              {
                cout << "Element " << i << " contains an isolated edge/face with "
                     << "no boundary condition " << endl;
                DISP(j); DISP(ref);
                DISP(this->mesh.GetBoundaryCondition(ref));
                DISP(this->mesh.Element(i));
                DISP(num_face); DISP(this->mesh.BoundaryRef(num_face));
                if (nb_proc == 1)
                  this->mesh.Write("test.mesh");
                
                abort();
              }
          }
      }
  }
  
  
  //! Computation of geometrical values : J_i, DF_i, normales
  /*!
    computation of Glob_jacobian(i,j) = \f$ \omega_k J_i(\xi_k) \f$
    computation of Glob_DFjm1(i,j) = \f$ J_i(\xi_k) DF_i(\xi_k) \f$
    if dfjm1_to_compute is set to true
    computation of ds n  (weighed normal)
    computation of tau and derivative (damping factor in PML layers)    
    */
  template<class Dimension>
  void VarGeometryProblem<Dimension>
  ::ComputeMassMatrix(bool compute_rho, bool delete_points)
  {    
#ifdef SELDON_WITH_MPI
      int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
      int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
      int rank_proc(0), nb_proc(1);
#endif

      // allocation of arrays
      this->AllocateMassMatrices();    
      
      // jacobian
      int nb_elt = this->mesh.GetNbElt();
      Glob_jacobian.Reallocate(nb_elt);
      if (ElementReference_Base::use_warburton_trick)
        {
          this->mesh.Glob_invSqrtJacobian.Reallocate(nb_elt);
          this->mesh.Glob_GradJacobian.Reallocate(nb_elt);
        }
      
      Glob_decomp_jacobian.Reallocate(nb_elt);
      
      // parfois on doit stocker DF_i^{-1}
      Glob_DFjm1.Reallocate(nb_elt);
      var_boundary.AllocateTauPML();

      // quadrature points
      this->Glob_PointsQuadrature.Reallocate(nb_elt);

      if (rank_proc == 0)
        if (this->print_level >= 1)
          cout << rank_proc << " We begin to compute the mass matrix - order "
               << this->mesh_num.GetOrder() << endl;
      
      IVect OrderFace(this->mesh.GetNbBoundary()); OrderFace.Fill(-1);
      Glob_normale.Reallocate(this->mesh.GetNbBoundary());
      Glob_dsj.Reallocate(this->mesh.GetNbBoundary());
      Glob_CoefPenalDG.Reallocate(this->mesh.GetNbElt());
      this->IsNewFace.Reallocate(nb_elt);
      
      this->ComputeGeometryQuantity(OrderFace);
      
      this->OrthogonalElement.Reallocate(nb_elt);
      this->OrthogonalElement.Fill(0);
      for (int i = 0; i < nb_elt; i++)
        if (this->InsidePML(i))
          {
            int ortho_elt = 0;
            for (int j = 0; j < this->Glob_DFjm1(i).GetM(); j++)
              {
                const typename Dimension::MatrixN_N& dfjm1 = this->Glob_DFjm1(i)(j);
                if (Dimension::dim_N == 2)
                  {
                    if ((abs(dfjm1(0, 1)) > R2::threshold) || (abs(dfjm1(1, 0)) > R2::threshold) )
                      ortho_elt = max(ortho_elt, 1);
                  }
                else
                  {
                    if ((abs(dfjm1(0, 2)) > R2::threshold) || (abs(dfjm1(1, 2)) > R2::threshold) 
                        || (abs(dfjm1(2, 0)) > R2::threshold) || (abs(dfjm1(2, 1)) > R2::threshold) )
                      ortho_elt = max(ortho_elt, 2);
                    
                    if ( (abs(dfjm1(0, 1)) > R2::threshold)  || (abs(dfjm1(1, 0)) > R2::threshold) )
                      ortho_elt = max(ortho_elt, 1);
                  }
              }
            
            this->OrthogonalElement(i) = ortho_elt;
          }
      
      // physical indices are evaluated if needed
      if (compute_rho)
        this->ComputePhysicalCoefficients();        
      
      this->ComputeArraySpecificEquation();
      
      if (write_quadrature_points)
        {
          string file_name = "quadrature_points";
          if (nb_proc > 1)
            file_name += "_P" + to_str(rank_proc) + ".dat";
          else
            file_name += ".dat";
          
          ofstream file_out(file_name.data());
          file_out.precision(15);
          for (int i = 0; i < this->Glob_PointsQuadrature.GetM(); i++)
            if (!this->InsidePML(i) || write_quad_points_pml)
              {
                for (int j = 0; j < this->GetNbPointsQuadratureInside(i); j++)
                  {
                    PrintNoBrackets(file_out, this->Glob_PointsQuadrature(i)(j));
                    file_out << '\n';
                  }
                
              }
          file_out.close();
        }
      
      if (delete_points)
        this->Glob_PointsQuadrature.Clear();
      
      if (!this->ComputeDFjm1())
        Glob_DFjm1.Clear();
      
      // for each face, we take the maximal penalty coefficient
      VectReal_wp coef_elt(Glob_CoefPenalDG);
      Glob_CoefPenalDG.Reallocate(this->mesh.GetNbBoundary());
      for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
        {
          int num_elem = this->mesh.Boundary(i).numElement(0);
          int num_elem2 = -1;
          if (this->mesh.Boundary(i).GetNbElements() == 2)
            num_elem2 = this->mesh.Boundary(i).numElement(1);
          
          if (num_elem >= 0)
            this->Glob_CoefPenalDG(i) = coef_elt(num_elem);
          
          if (num_elem2 >= 0)
            this->Glob_CoefPenalDG(i) = max(this->Glob_CoefPenalDG(i), coef_elt(num_elem2));
        }
      
      this->ComputeTauCoefficient();
  }
  
  
  //! we deallocate arrays, which contain mass matrices
  template<class Dimension>
  void VarGeometryProblem<Dimension>::ClearMassMatrix()
  {
    Glob_jacobian.Clear(); Glob_DFjm1.Clear();
    Glob_normale.Clear(); Glob_dsj.Clear();
    IsNewFace.Clear(); this->mesh.Clear();
  }
  
  
  //! computes the order for each element
  template<class Dimension>
  void VarGeometryProblem<Dimension>::ComputeVariableOrder()
  {
    // if a variable order is required
    // we compute the order of approximation
    // for each edge, face and element for the global mesh
    // this information will be distributed to other processors
    Vector<Real_wp> velocity;
    this->GetVelocityOnElements(velocity, this->mesh);
    
    // smoothing step
    /* Vector<Real_wp> v0_vertex(this->mesh.GetNbVertices());
       v0_vertex.Fill(1e300);
       for (int i = 0; i < this->mesh.GetNbElt(); i++)
       for (int j = 0; j < this->mesh.Element(i).GetNbVertices(); j++)
       {
       int nv = this->mesh.Element(i).numVertex(j);
       v0_vertex(nv) = min(v0_vertex(nv), velocity(i));
       }
       
       for (int i = 0; i < this->mesh.GetNbElt(); i++)
       for (int j = 0; j < this->mesh.Element(i).GetNbVertices(); j++)
       {
       int nv = this->mesh.Element(i).numVertex(j);
       velocity(i) = min(v0_vertex(nv), velocity(i));
       }
    */
    
    // adaptation en fonction de la vitesse et de la frequence
    for (int i = 0; i < velocity.GetM(); i++)
      velocity(i) = this->GetFrequency() / velocity(i);
    
    this->mesh_num.SetCoefficientVariableOrder(velocity);
    this->mesh_num.ComputeVariableOrder(this->FormulationDG());

    for (int i = 0; i < other_mesh_num.GetM(); i++)
      {
        this->other_mesh_num(i)->SetCoefficientVariableOrder(velocity);
        this->other_mesh_num(i)->ComputeVariableOrder(this->FormulationDG());
      }    

    if (this->print_level >= 2)
      this->mesh.WriteOrder("order.mesh", this->mesh_num);    
  }
  
  
  //! returns mesh associated with file_name
  template<class Dimension>
  Mesh<Dimension>& VarGeometryProblem<Dimension>::
  GetPhysicalCoefficientMesh(int i, const string& file_name, int r, bool same_mesh,
                             Vector<string>& all_names, Vector<Mesh<Dimension> >& all_mesh)
  {
    bool only_quad = (this->mesh.GetNbTriangles() == 0);
#ifdef SELDON_WITH_MPI
    bool only_quad_loc(only_quad);
    MPI_Allreduce(&only_quad_loc, &only_quad, 1,
                  MPI_CHAR, MPI_LAND, var_comm.comm_group_mode);
#endif
    
    // checking if the mesh is already constructed
    for (int j = 0; j < i; j++)
      {
        if (!file_name.compare(all_names(j)))
          {
            // we have already read this mesh
            if (r == all_mesh(j).GetGeometryOrder())
              {
                // same order, we return the old mesh
                return all_mesh(j);
              }
          }
      }
    
    // now we have to create a new mesh
    all_names(i) = file_name;
    
    // checking if a same mesh exists with different order
    for (int j = 0; j < i; j++)
      {
        if (!file_name.compare(all_names(j)))
          {
            // we copy the mesh
            all_mesh(i) = all_mesh(j);
            all_mesh(i).SetGeometryOrder(r, only_quad);
            return all_mesh(i);
          }
      }
    
    // last case, reading of the new mesh
    all_mesh(i).SetGeometryOrder(r, only_quad);
    if (same_mesh)
      all_mesh(i).CopyCurves(this->mesh);
    
    all_mesh(i).Read(file_name);
    all_mesh(i).RearrangeElements();
    return all_mesh(i);
  }
  
  
  //! returns projector associated with order r
  template<class Dimension>
  Vector<FiniteElementInterpolator>& VarGeometryProblem<Dimension>
  ::GetPhysicalCoefInterp(int r, const Mesh<Dimension>& meshb,
                          const TinyVector<IVect, 4>& order_mesh,
                          Vector<Vector<FiniteElementInterpolator> >& all_interp,
                          Vector<int> & order_interp)
  {
    int j = -1; bool index_found = false;
    for (int i = 0; i < order_interp.GetM(); i++)
      {
        int ri = order_interp(i);
        if (r == ri)
          {
            // we have already computed this projector
            return all_interp(i);
          }        
        
        if (!index_found)
          if (ri < 0)
            {
              // we found a new slot for the projector
              j = i;
              index_found = true;
            }
      }
    
    // computing the new projector
    Vector<FiniteElementInterpolator>& proj = all_interp(j);
    
    order_interp(j) = r;
    int order_max = 0;
    for (int p = 0; p < 4; p++)
      for (int i = 0; i < order_mesh(p).GetM(); i++)
        order_max = max(order_max, order_mesh(p)(i));
    
    proj.Reallocate(order_max+1);
    TinyVector<typename Dimension::VectR_N, 4> points_div;
    for (int rm = 1; rm <= order_max; rm++)
      {
        bool order_used = false;
        points_div(0).Clear(); points_div(1).Clear();
        points_div(2).Clear(); points_div(3).Clear();
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          if (this->mesh_num.GetOrderElement(i) == rm)	    
            {
              order_used = true;
              int type_elt = this->mesh.Element(i).GetHybridType();
              if (points_div(type_elt).GetM() <= 0)
                points_div(type_elt) = this->PointsQuadInsideND(i);
            }
        
        if (order_used)
          proj(rm).InitProjection(meshb.GetReferenceElement(), points_div);
      }
    
    return proj;
  }
  
  
  //! returns projector associated with order r
  template<class Dimension>
  Vector<FiniteElementInterpolator>& VarGeometryProblem<Dimension>::
  GetPhysicalCoefInterpSurf(int r, const Mesh<Dimension>& meshb,
                            const TinyVector<IVect, 4>& order_mesh,
                            Vector<Vector<FiniteElementInterpolator> >& all_interp,
                            Vector<int> & order_interp)
  {
    int j = -1; bool index_found = false;
    for (int i = 0; i < order_interp.GetM(); i++)
      {
        int ri = order_interp(i);
        if (r == ri)
          {
            // we have already computed this projector
            return all_interp(i);
          }        
        
        if (!index_found)
          if (ri < 0)
            {
              // we found a new slot for the projector
              j = i;
              index_found = true;
            }
      }
    
    // computing the new projector
    Vector<FiniteElementInterpolator>& proj = all_interp(j);
    
    order_interp(j) = r;
    int order_max = 0;
    for (int i = 0; i < order_mesh(0).GetM(); i++)
      order_max = max(order_max, order_mesh(0)(i));
    
    for (int i = 0; i < order_mesh(1).GetM(); i++)
      order_max = max(order_max, order_mesh(1)(i));
    
    proj.Reallocate(order_max+1);
    TinyVector<typename DimensionBoundary::VectR_N, 2> points_div;
    for (int rm = 1; rm <= order_max; rm++)
      {
        bool order_used = false;
        points_div(0).Clear();
        points_div(1).Clear();
        for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
          if (this->mesh_num.GetOrderQuadrature(i) == rm)
            {
              order_used = true;
              int type_elt = this->mesh.Boundary(i).GetHybridType();
              if (points_div(type_elt).GetM() <= 0)
                points_div(type_elt) = this->mesh_num.GetReferenceQuadrature(rm, type_elt);
            }

        Vector<const ElementGeomReference<DimensionBoundary>* > ref_elt;
        meshb.GetReferenceElementSurface(ref_elt);
        if (order_used)
          proj(rm).InitProjection(ref_elt, points_div);
      }
    
    return proj;
  }
  
  
  //! computation of a physical index : rho, mu, epsilon, etc
  /*!
    \param[in] ind number of the variable index to evaluate
    \param[in] ref reference associated with this index
    \param[in] nb_elt_on_ref number of elements of reference ref
    \param[in] num_index number associated with the physical index 
    \param[in] num_component component number of the physical index
    \param[in] compute_grad if true the gradient of the index is also computed
    \param[in] compute_hess if true the hessian of the index is also computed
    \param[in] PointsQuadrature quadrature points for which the index will be evaluated
    \param[inout] fct_rho definition of the variable index, in output it will store the values
    of the physical index at the quadrature points (and gradient if asked)
    \param[inout] all_mesh different meshes previously stored (during the computation of previous indexes)
    \param[inout] all_names names of the different meshes
    \param[inout] all_interp projectors that have been previously computed
    \param[inout] order_interp orders associated with projectors stored in all_interp
    \param[inout] all_interp_surf surface projectors that have been previously computed
    \param[inout] order_interp_surf orders associated with projectors stored in all_interp_surf
    \param[in] var_file_name name of the current mesh (with its path)
    */
  template<class Dimension> template<class T>
  void VarGeometryProblem<Dimension>::
  ComputePhysicalCoefficientsRef(int ind, int ref, int nb_elt_on_ref,
                                 int num_index, int num_component,
                                 bool compute_grad, bool compute_hess, Vector<VectR_N>& PointsQuadrature,
                                 PhysicalVaryingMedia<Dimension, T>& fct_rho,
                                 Vector<Mesh<Dimension> >& all_mesh,
                                 Vector<string>& all_names,
                                 Vector<Vector<FiniteElementInterpolator> >& all_interp,
                                 Vector<int>& order_interp,
                                 Vector<Vector<FiniteElementInterpolator> >& all_interp_surf,
                                 Vector<int>& order_interp_surf, const string& var_file_name)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif

    // orders present in the mesh (for elements)
    TinyVector<IVect, 4> order_mesh;
    this->mesh_num.GetOrder(order_mesh);
    
    // orders present in the mesh (for integration on surfaces)
    TinyVector<IVect, 4> order_quad;
    this->mesh_num.GetOrderQuadrature(order_quad);
    
    fct_rho.Reallocate(nb_elt_on_ref, compute_grad, compute_hess);
    Mesh<Dimension>* mesh_index = NULL;
    if (fct_rho.GetMediaType() == fct_rho.MESH)
      {
        // case where the index is defined on a mesh
        PhysicalMeshMedia<Dimension, T>& index = fct_rho.index_mesh;
        string file_name = index.GetMeshName();
        int r = index.GetOrder();
        bool same_mesh = fct_rho.IsSameMesh();
        if (!file_name.compare(var_file_name))
          same_mesh = true;

        // we retrieve the mesh
        Mesh<Dimension>& meshb = this->GetPhysicalCoefficientMesh(ind, file_name, r, same_mesh,
                                                                  all_names, all_mesh);
        
        mesh_index = &meshb;
        // if needed the gradient is computed
        // we use here meshb because the order of meshb may be different from this->mesh
        index.CheckMesh(meshb);
        if (compute_grad)
          index.ComputeGradient(meshb, compute_hess);
        
        if (same_mesh)
          {
            Vector<FiniteElementInterpolator>& proj 
              = this->GetPhysicalCoefInterp(r, meshb, order_mesh,
                                            all_interp, order_interp);
            
            Vector<FiniteElementInterpolator>& proj_surf
              = this->GetPhysicalCoefInterpSurf(r, meshb, order_quad,
                                                all_interp_surf, order_interp_surf);
            
            // restriction on element of considered processor
            if (var_comm.GetNbProcPerMode() > 1)
              if (index.GetNbElt() != this->mesh.GetNbElt())
                index.RemoveGlobalElements(this->mesh.GlobElementNumber_Subdomain);
            
            Vector<T> rho_loc, rho_quad;
            Vector<TinyVector<T, Dimension::dim_N> > grad_rho, grad_rho_quad;
            for (int i = 0; i < this->mesh.GetNbElt(); i++)
              if (this->mesh.Element(i).GetReference() == ref)
                {
                  int type_elt = this->mesh.Element(i).GetHybridType();           
                  int rm = this->mesh_num.GetOrderElement(i);
                  
                  // we retrieve the values of rho on the element
                  int nb_nodes = index.GetNbNodes(i);
                  rho_loc.Reallocate(nb_nodes);
                  for (int j = 0; j < nb_nodes; j++)
                    rho_loc(j) = index.GetValue(i, j);
                  
                  // projection on quadrature points of the volume
                  proj(rm).ProjectScalar(rho_loc, rho_quad, type_elt);
                  fct_rho.ReallocateOnElement(ElementRho(i), PointsQuadrature(i).GetM(),
                                              compute_grad, compute_hess);
                  
                  for (int j = 0; j < rho_quad.GetM(); j++)
                    fct_rho.SetCoefficient(ElementRho(i), j, rho_quad(j));
                  
                  // then quadrature points of the boundaries
                  int offset = rho_quad.GetM();
                  for (int num_loc = 0; num_loc < this->mesh.Element(i).GetNbBoundary(); num_loc++)
                    {
                      // we use here meshb because the order of meshb may be different from this->mesh
                      IVect node_surf = meshb.GetLocalNodalNumber(type_elt, num_loc);
                      Vector<T> rho_node(node_surf.GetM()), rho_surf;
                      for (int j = 0; j < rho_node.GetM(); j++)
                        rho_node(j) = rho_loc(node_surf(j));
                      
                      int num_face = this->mesh.Element(i).numBoundary(num_loc);
                      int rf = this->mesh_num.GetOrderQuadrature(num_face);
                      int type_elt = this->mesh.Boundary(num_face).GetHybridType();
                      proj_surf(rf).ProjectScalar(rho_node, rho_surf, type_elt);
                      for (int j = 0; j < rho_surf.GetM(); j++)
                        fct_rho.SetCoefficient(ElementRho(i), offset+j, rho_surf(j));
                      
                      offset += rho_surf.GetM();
                    }
                  
                  // if gradient is needed, projection of gradient on quadrature points
                  if (compute_grad)
                    {
                      int nb_calcul = 1;
                      if (compute_hess)
                        nb_calcul += Dimension::dim_N;
                      
                      for (int num_calcul = 0; num_calcul < nb_calcul; num_calcul++)
                        {
                          grad_rho.Reallocate(nb_nodes);
                          if (num_calcul == 0)
                            for (int j = 0; j < nb_nodes; j++)
                              grad_rho(j) = index.GetGradient(i, j);
                          else
                            for (int j = 0; j < nb_nodes; j++)
                              grad_rho(j) = index.GetHessian(i, j, num_calcul-1);
                          
                          // projection on quadrature points of the volume
                          proj(rm).Project(grad_rho, grad_rho_quad, type_elt);
                          if (num_calcul == 0)
                            for (int j = 0; j < rho_quad.GetM(); j++)
                              fct_rho.SetCoefGradient(ElementRho(i), j, grad_rho_quad(j));
                          else
                            for (int j = 0; j < rho_quad.GetM(); j++)
                              fct_rho.SetCoefHessian(ElementRho(i), j, num_calcul-1, grad_rho_quad(j));
                          
                          // then quadrature points of the boundaries
                          offset = rho_quad.GetM();
                          for (int num_loc = 0; num_loc < this->mesh.Element(i).GetNbBoundary(); num_loc++)
                            {
                              IVect node_surf = meshb.GetLocalNodalNumber(type_elt, num_loc);
                              Vector<TinyVector<T, Dimension::dim_N> >
                                grad_rho_node(node_surf.GetM()), grad_rho_surf;
                              
                              for (int j = 0; j < grad_rho_node.GetM(); j++)
                                grad_rho_node(j) = grad_rho(node_surf(j));
                              
                              int num_face = this->mesh.Element(i).numBoundary(num_loc);
                              int rf = this->mesh_num.GetOrderQuadrature(num_face);
                              int type_elt = this->mesh.Boundary(num_face).GetHybridType();
                              proj_surf(rf).Project(grad_rho_node, grad_rho_surf, type_elt);
                              if (num_calcul == 0)
                                for (int j = 0; j < grad_rho_surf.GetM(); j++)
                                  fct_rho.SetCoefGradient(ElementRho(i), offset+j, grad_rho_surf(j));
                              else
                                for (int j = 0; j < grad_rho_surf.GetM(); j++)
                                  fct_rho.SetCoefHessian(ElementRho(i), offset+j, num_calcul-1, grad_rho_surf(j));
                              
                              offset += grad_rho_surf.GetM();
                            }
                        }
                    }
                }
          }
        else
          {
            cout << "case currently not treated " << endl;
            abort();
          }            
      }
    else
      {
        // allocating arrays
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          if (this->mesh.Element(i).GetReference() == ref)
            {
              int nb_quad = PointsQuadrature(i).GetM();
              fct_rho.ReallocateOnElement(ElementRho(i), nb_quad, compute_grad, compute_hess);
            }
        
        // computation of values 
        switch (fct_rho.GetMediaType())
          {
          case PhysicalVaryingMedia<Dimension, T>::SINUSOIDE :
            fct_rho.index_sinus.
              ComputeIndex(this->mesh, ElementRho, ref,
                           PointsQuadrature, compute_grad, compute_hess,
                           fct_rho.GetValue(), fct_rho.GetGradient(), fct_rho.GetHessian());
            break;
          case PhysicalVaryingMedia<Dimension, T>::REGULAR_GRID :
            fct_rho.index_regular.
              ComputeIndex(this->mesh, ElementRho, ref,
                           PointsQuadrature, compute_grad, compute_hess,
                           fct_rho.GetValue(), fct_rho.GetGradient(), fct_rho.GetHessian());
            break;
          case PhysicalVaryingMedia<Dimension, T>::RADIAL :
            fct_rho.index_radial.
              ComputeIndex(this->mesh, ElementRho, ref,
                           PointsQuadrature, compute_grad, compute_hess,
                           fct_rho.GetValue(), fct_rho.GetGradient(), fct_rho.GetHessian());
            break;
          case PhysicalVaryingMedia<Dimension, T>::USER :
            ComputeVariableUserIndex(*this, ElementRho, ref,
                                     num_index, num_component,
                                     PointsQuadrature, compute_grad, compute_hess,
                                     fct_rho.GetValue(), fct_rho.GetGradient(),
                                     fct_rho.GetHessian(), fct_rho.GetOffset(),
                                     fct_rho.GetAmplitude());
            break;
          }
      }
    
    if ((this->print_level < 7) || (nb_proc > 1))
      return;
    
    // for print levels greater than 7, we try to interpolate the physical index
    // on the points of the interpolation grid, and write them as standard output files
    for (int n = 0; n < var_output.var_grid.GetM(); n++)
      {
        GridInterpolationFull<Dimension>& var_gr = var_output.var_grid(n);
        GridInterpolation<Dimension>& var_interp = var_output.all_points_display;
        const IVect& list_points = var_gr.GetPointNumber();
        const Vector<IVect>& list_points_proc = var_gr.GetPointNumberAllProc();
        
        int Nglob = var_gr.GetNbGlobalPoints();
        Vector<T> interp_index(Nglob);
        interp_index.Zero();
        T val; VectReal_wp phi;
        for (int i1 = 0; i1 < list_points.GetM(); i1++)
          {
            int i = list_points(i1);
            if (i < 0)
              continue;
            
            R_N pt_glob = var_interp.GetGlobalCoordinate(i);
            int num_elem = var_interp.GetElementNumber(i);
            R_N pt_loc = var_interp.GetLocalCoordinate(i);
            int ref_elem = this->mesh.Element(num_elem).GetReference();
            SetComplexZero(val);
            if (ref == ref_elem)              
              switch (fct_rho.GetMediaType())
                {
                case PhysicalVaryingMedia<Dimension, T>::SINUSOIDE :
                  fct_rho.index_sinus.ComputeIndexAtPoint(pt_glob, val);
                  break;
                case PhysicalVaryingMedia<Dimension, T>::REGULAR_GRID :
                  fct_rho.index_regular.ComputeIndexAtPoint(pt_glob, val);
                  break;
                case PhysicalVaryingMedia<Dimension, T>::RADIAL :
                  fct_rho.index_radial.ComputeIndexAtPoint(pt_glob, val);
                  break;
                case PhysicalVaryingMedia<Dimension, T>::MESH :
                  {
                    mesh_index->ComputeValuesPhiNodalRef(num_elem, pt_loc, phi);
                    fct_rho.index_mesh.ComputeIndexAtPoint(num_elem, phi, val);
                    break;
                  }
                }
            
            int iglob = list_points_proc(0)(i1);
            interp_index(iglob) = val;
          }
        
        // the index is written in the output file
        string file_name = this->GetPhysicalIndexName(num_index) + "_ref" + to_str(ref)
          + "_G" + to_str(n) + "_U" + to_str(num_component) + ".dat";
        
        cout << "Writing index " << num_index << " in " << file_name << endl;
        WriteMatlab(interp_index, var_gr, file_name, var_output.OutputWrittenInDoublePrecision());
      }
  }
  
  
  //! computation of physical indices : rho, mu, epsilon, etc
  template<class Dimension>
  void VarGeometryProblem<Dimension>::ComputePhysicalCoefficients()
  {
    // first step : we update ElementRho
    int nb_elt = this->mesh.GetNbElt();
    ElementRho.Reallocate(nb_elt); ElementRho.Zero();
    IVect NbElement_ByRef(this->GetNbPhysicalIndices()); NbElement_ByRef.Zero();
    for (int i = 0; i < nb_elt; i++)
      {
        int ref = this->mesh.Element(i).GetReference();
        ElementRho(i) = NbElement_ByRef(ref);
        NbElement_ByRef(ref)++;
      }
    
    // we retrieve all the varying indices of the problem
    Vector<PhysicalVaryingMedia<Dimension, Complex_wp>* > rho_complex;
    Vector<PhysicalVaryingMedia<Dimension, Real_wp>* > rho_real;
    IVect num_index, num_component, num_ref;
    Vector<bool> compute_grad, compute_hess;
    this->GetVaryingIndices(rho_complex, rho_real, num_ref,
                            num_index, num_component, compute_grad, compute_hess);
    
    Vector<Mesh<Dimension> > all_mesh(rho_complex.GetM() + rho_real.GetM());
    Vector<string> all_names(rho_complex.GetM() + rho_real.GetM());
    Vector<Vector<FiniteElementInterpolator> >
      all_interp(rho_complex.GetM() + rho_real.GetM());
    
    Vector<Vector<FiniteElementInterpolator> >
      all_interp_surf(rho_complex.GetM() + rho_real.GetM());
    
    Vector<int> order_interp(all_interp.GetM()); 
    Vector<int> order_interp_surf(all_interp_surf.GetM()); 
    order_interp.Fill(-1);
    order_interp_surf.Fill(-1);
    
    string var_file_name("NONE");
    if (this->mesh_data.GetM() > 0)
      var_file_name = this->mesh.GetPathName() + this->mesh_data(0)(0);
    
    // then we compute the values of indices on quadrature points
    for (int ind = 0; ind < rho_complex.GetM(); ind++)
      {
        PhysicalVaryingMedia<Dimension, Complex_wp>& fct_rho = *rho_complex(ind);
        this->ComputePhysicalCoefficientsRef(ind, num_ref(ind), NbElement_ByRef(num_ref(ind)),
                                             num_index(ind), num_component(ind),
                                             compute_grad(ind), compute_hess(ind), this->Glob_PointsQuadrature,
                                             fct_rho, all_mesh, all_names,
                                             all_interp, order_interp, all_interp_surf,
                                             order_interp_surf, var_file_name);
      }
    
    // same stuff for real indices
    for (int ind = rho_complex.GetM(); ind < rho_complex.GetM()+rho_real.GetM(); ind++)
      {
        PhysicalVaryingMedia<Dimension, Real_wp>& fct_rho = *rho_real(ind-rho_complex.GetM());
        this->ComputePhysicalCoefficientsRef(ind, num_ref(ind), NbElement_ByRef(num_ref(ind)),
                                             num_index(ind), num_component(ind),
                                             compute_grad(ind), compute_hess(ind), this->Glob_PointsQuadrature,
                                             fct_rho, all_mesh, all_names,
                                             all_interp, order_interp, all_interp_surf,
                                             order_interp_surf, var_file_name);
      }
    
    this->FinalizeComputationVaryingIndices();    
  }
  
  
  //! Checks and modifies if necessary the initial mesh
  /*!
    This method is called right after the construction of the mesh
  */
  template<class Dimension>
  void VarGeometryProblem<Dimension>::CheckInputMesh()
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    if ((nb_proc == 1) && this->mesh.GetNbElt() <= 0)
      {
        cout << "The mesh contains no element" << endl;
        cout << "Input = " << this->mesh_data(0) << endl;
        abort();
      }
  }
  
  
  //! construction of the mesh and finite element
  /*!
    \param[in] name_elt name of the finite element to use
    \param[in] split_mesh if true, we split the mesh into several parts
    and distribute them to other processors
    */
  template<class Dimension>
  void VarGeometryProblem<Dimension>
  ::ComputeMeshAndFiniteElement(const string& name_elt, bool split_mesh)
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_comm.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_comm.comm_group_mode, &rank_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    this->finite_element_name = name_elt;
    if (nb_proc <= 1)
      split_mesh = false;

    if ((this->dg_formulation != ElementReference_Base::CONTINUOUS)
        && (this->other_mesh_num.GetM() > 0))
      {
        cout << "Unknowns with different orders only implemented for continuous elements" << endl;
        abort();
      }
    
    ////////////////////
    // MESH TREATMENT //
    
    if ((rank_proc == 0)||(!split_mesh))
      {
        // reads the mesh in the master node
        // type of the mesh, triangular, quadrilateral or hybrid ?
        // method ConstructMesh is implemented in Mesh2D.cxx (2-D) or Mesh3D.cxx (3-D)
        if (this->print_level >= 2)
          cout << "Constructing mesh ... " << endl;
        
        if (this->mesh_data.GetM() > 0)
          this->mesh.ConstructMesh(0, this->mesh_data(0));
        else
          if (this->print_level >= 1)
            cout << "Warning : no mesh has been provided" << endl;

        this->CheckInputMesh();

        if (this->print_level >= 2)
          {
            cout << "Mesh successfully constructed " << endl;
            cout << "It contains " << this->mesh.GetNbElt() << " elements " << endl;
            cout << "It contains " << this->mesh.GetNbVertices() << " vertices " << endl;
            cout << "It contains " << this->mesh.GetNbEdges() << " edges " << endl;
            cout << "It contains " << this->mesh.GetNbBoundary() << " edges/faces " << endl;
          }
        
        // computing order for each edge, face, element if variable
        if ((nb_proc == 1) || (!split_mesh))
          {
            if (this->mesh_num.IsOrderVariable())
              this->ComputeVariableOrder();
            
            if (this->print_level >= 2)
              cout << "Constructing finite element ... " << endl;
            
            // constructing finite elements
            this->ConstructFiniteElement(name_elt);
            
            if (this->print_level >= 2)
              cout << "Finite elements successfully constructed " << endl;
          }
        
        
#ifdef MONTJOIE_WITH_TRANSMISSION        
        this->PartMeshTransmission();
#endif

        if ((nb_proc == 1) || (!split_mesh))
          { 
            // finding modes (for periodic and cyclic conditions)
            var_boundary.InitCyclicDomain();
            
            // finding pml
            var_boundary.FindElementsInsidePML();
            
            if (this->print_level >= 2)
              cout << "Numbering mesh ... " << endl;
            
            // numbering of the mesh, (see file Mesh2D.cxx or Mesh3D.cxx)
            this->mesh_num.NumberMesh();
            for (int i = 0; i < this->other_mesh_num.GetM(); i++)
              this->other_mesh_num(i)->NumberMesh();
            
#ifdef MONTJOIE_WITH_TRANSMISSION
            {
              IVect MinimalProc;
              this->TreatTransmission(MinimalProc);
            }
#endif
            
            {
              IVect MinimalProc;
              this->TreatGibc(MinimalProc);
            }
            
            Real_wp step_h = this->mesh.GetMeshSize();
            if (this->print_level >= 0)
              {
                cout << "The mesh contains " << this->mesh.GetNbElt() << " elements " << endl;
                cout << "The mesh contains " << this->mesh_num.GetNbDof()
                     << " degrees of freedom " << endl;
                
                cout << "The mesh contains " << this->mesh_num.GetNbDofPML()
                     << " dofs inside the PML " << endl;

                for (int i = 0; i < this->other_mesh_num.GetM(); i++)
                  cout << "The mesh contains " << this->other_mesh_num(i)->GetNbDof()
                       << " degrees of freedom for other numbering " << i << endl;
                
                cout << "Mesh size : " << step_h << endl;
              }
            
          }
        
      }
    
    if (split_mesh)
      {
        // finding modes (for periodic and cyclic conditions)
        var_boundary.InitCyclicDomain();
        
        this->SplitMeshForParallelComputation(name_elt);

        if (rank_proc != 0)
          this->mesh.print_level = 0;
      }
    else
      {
        // computation of dofs for the initial formulation
        this->ComputeNumberOfDofs();

        var_comm.nodl_mesh.Reallocate(this->GetNbMeshNumberings());
        var_comm.nodl_pml.Reallocate(this->GetNbMeshNumberings());
        var_comm.nodl_mesh(0) = this->mesh_num.GetNbDof();
        var_comm.nodl_pml(0) = this->mesh_num.GetNbDofPML();
        for (int i = 0; i < this->other_mesh_num.GetM(); i++)
          {
            var_comm.nodl_mesh(i+1) = this->other_mesh_num(i)->GetNbDof();
            var_comm.nodl_pml(i+1) = this->other_mesh_num(i)->GetNbDofPML();
          }
      }

    // Dirichlet condition
    if ((rank_proc == 0)||(this->print_level >= 10))
      if (this->print_level >= 2)
        cout << rank_proc << " Finding dofs with Dirichlet condition ... "<<endl;
    
    var_boundary.TreatDirichletCondition();
    
    if ((rank_proc == 0) && (this->print_level >= 2))
      cout << "Dirichlet done " << endl;
    

    // MESH TREATMENT //
    ////////////////////
    
    
    if (this->print_level >= 12)
      {
        if (rank_proc == 0)
          cout << "Checking continuity of basis functions ... " << endl;
        
        // checking continuity (for H1/Hcurl) in a "debug" mode
        if (this->FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          this->CheckContinuity();
        
        if (rank_proc == 0)
          cout << "Basis functions are continuous " << endl;
      }
    
    /* Vector<R_N> PointsDof(this->mesh_num.GetNbDof());
       SetPoints<Dimension> PointsElem;
       Vector<R_N> s;
         for (int i = 0; i < this->mesh.GetNbElt(); i++)
         {
         int r = this->mesh_num.GetOrderElement(i);
         this->mesh.GetVerticesElement(i, s);

         this->GetReferenceElement(i).FjElemDof(s, PointsElem, this->mesh, i);
         for (int j = 0; j < this->mesh_num.Element(i).GetNbDof(); j++)
         {
         int jglob = this->mesh_num.Element(i).GetNumberDof(j);
         PointsDof(jglob) = PointsElem.GetPointDof(j);
         }

      //cout << "Element " << i << endl;
      //cout << "Dof numbers : " << endl << this->mesh_num.Element(i).GetNodle() << endl;
      //cout << "Negative dofs : " << endl 
      //<< this->mesh_num.Element(i).GetNegativeDofNumber() << endl;
      }
      */
      /* Vector<int> permut(PointsDof.GetM());
         permut.Fill();      
         Sort(PointsDof, permut);

         int nb_u = ElementReference<Dimension, type>::nb_components_u;
         int r = this->mesh_num.GetOrder();
         string name_file = "points_dof_N"+to_str(nb_u)+"_R"+to_str(r)+".dat";
         ofstream file_out(name_file.data());
         file_out.precision(7);
         for (int i = 0; i < PointsDof.GetM(); i++)
         {
         PrintNoBrackets(file_out, PointsDof(i));
         file_out << '\n';
         }

         file_out.close();
         permut.WriteText("numbers_N"+to_str(nb_u)+"_R"+to_str(r)+".dat");

         file_out.open("points_mesh.dat");
         for (int i = 0; i < this->mesh.GetNbVertices(); i++)
         {
         PrintNoBrackets(file_out, this->mesh.Vertex(i));
         file_out << '\n';
         }

         file_out.close();

         file_out.open("num_quad.data");
         for (int i = 0; i < this->mesh.GetNbElt(); i++)
         {
         for (int j = 0; j < this->mesh.Element(i).GetNbVertices(); j++)
         file_out << this->mesh.Element(i).numVertex(j)+1 << " ";

         file_out << '\n';
         }

         file_out.close();
         */
      // computation of the output grids and meshes        
    if (var_output.grid_to_be_computed)
      {
        if (this->var_chrono.NameExists("MeshGeneration"))
          this->var_chrono.Stop("MeshGeneration");

        this->var_chrono.SetMessage("InterpolationGrid", "compute the interpolation grid");
        this->var_chrono.Start("InterpolationGrid");

        if ((rank_proc == 0) && (this->print_level >= 2))
          cout << "Computation of interpolation grid ... " << endl;

        var_output.InitVarGrid();
        var_output.ComputeVarGrid();

        if ((rank_proc == 0) && (this->print_level >= 2))
          cout << "Grid successfully computed " << endl;

        this->var_chrono.Stop("InterpolationGrid");

        if (this->var_chrono.NameExists("MeshGeneration"))
          this->var_chrono.Start("MeshGeneration");
      }
  }
  
  
  //! copies input parameters of another similar problem
  template<class Dimension>
  void VarGeometryProblem<Dimension>::CopyInputData(const VarProblem_Base& var_base)
  {
    VarProblem_Base::CopyInputData(var_base);
    
    const VarGeometryProblem<Dimension>& var = static_cast<const VarGeometryProblem<Dimension>& >(var_base);
    
    this->mesh.CopyInputData(var.mesh);
    this->mesh_num.CopyInputData(var.mesh_num);
    for (int k = 0; k < this->other_mesh_num.GetM(); k++)
      this->other_mesh_num(k)->CopyInputData(*var.other_mesh_num(k));
    
    kwave = var.kwave;
    polarization = var.polarization;
    polarization_grad = var.polarization_grad;
    origine_phase = var.origine_phase;
  }
  
  
  //! checks if continuity is ensured between elements
  template<class Dimension>
  void VarGeometryProblem<Dimension>::
  CheckLocalContinuity(int i, int num_elem, int num_elem2,
                       SetPoints<Dimension>& PointsElem, SetPoints<Dimension>& PointsElem2,
                       SetMatrices<Dimension>& MatricesElem,
                       SetMatrices<Dimension>& MatricesElem2,
                       const ElementReference<Dimension, 3>& Fb,
                       const ElementReference<Dimension, 3>& Fb2)
  {
#ifdef MONTJOIE_WITH_NODAL_HDIV
    
    // we check continuity of basis functions across boundaries
    Vector<R_N> s, s2; R_N tmp;
    
    this->mesh.GetVerticesElement(num_elem, s);
    this->mesh.GetVerticesElement(num_elem2, s2);
    
    int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
    int num_loc2 = -1, rot = 0;
    this->mesh_num.GetBoundaryRotation(i, num_elem, num_elem2, num_loc, num_loc2, rot);
    
    Fb.FjElem(s, PointsElem, this->mesh, num_elem);
    Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, num_elem);
    
    Fb.FjSurfaceElem(s, PointsElem, this->mesh, num_elem, num_loc);
    Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, this->mesh, num_elem, num_loc);
    
    Fb2.FjElem(s2, PointsElem2, this->mesh, num_elem2);
    Fb2.DFjElem(s2, PointsElem2, MatricesElem2, this->mesh, num_elem2);
    
    Fb2.FjSurfaceElem(s2, PointsElem2, this->mesh, num_elem2, num_loc2);
    Fb2.DFjSurfaceElem(s2, PointsElem2, MatricesElem2, this->mesh, num_elem2, num_loc2);
    
    int nb_quad = Fb.GetNbQuadBoundary(num_loc);
    for (int k = 0; k < nb_quad; k++)
      {
        int nk = Fb.GetQuadNumber(num_loc, k);
        int k2 = -1;
        for (int kp = 0; kp < nb_quad; kp++)
          {
            int nkp = Fb2.GetQuadNumber(num_loc2, kp);
            if (PointsElem.GetPointQuadrature(nk) == PointsElem2.GetPointQuadrature(nkp))
              k2 = kp;
          }
        
        int nk2 = Fb2.GetQuadNumber(num_loc2, k2);
        const Vector<R_N>& pts = Fb.PointsND();
        const Vector<R_N>& pts2 = Fb2.PointsND();

        // we check if normales are opposite
        R_N normale, normale2;
        normale = MatricesElem.GetNormaleQuadratureBoundary(k);
        normale2 = MatricesElem2.GetNormaleQuadratureBoundary(k2);
        Mlt(1.0/Norm2(normale), normale);
        Mlt(-1.0/Norm2(normale2), normale2);
        if (normale != normale2)
          {
            cout << "Normale are not the same "<<endl;
            abort();
          }
        
        // now we compute values of phi on the quadrature point
        typename Dimension::MatrixN_N dfj, dfjm1;
        Vector<TinyVector<Real_wp, Dimension::dim_N> > phi, phi2;
        dfj = MatricesElem.GetPointQuadrature(nk);
        GetInverse(dfj, dfjm1);
        Fb.ComputeValuesPhi(pts(nk), phi, dfjm1, this->mesh_num, num_elem);
        
        dfj = MatricesElem2.GetPointQuadrature(nk2);
        GetInverse(dfj, dfjm1);
        Fb2.ComputeValuesPhi(pts2(nk2), phi2, dfjm1, this->mesh_num, num_elem2);
        
        VectReal_wp phi1_n(phi.GetM()), phi2_n(phi2.GetM());
        int nb_dof = this->mesh_num.GetNbLocalDof(num_elem);
        int nb_dof2 = this->mesh_num.GetNbLocalDof(num_elem2);		
        for (int j = 0; j < nb_dof; j++)
          phi1_n(j) = DotProd(phi(j), normale);
        
        for (int j = 0; j < nb_dof2; j++)
          phi2_n(j) = DotProd(phi2(j), normale);
        
        // we check continuity of phi \cdot n
        for (int j = 0; j < nb_dof; j++)
          {
            int j2 = -1;
            for (int jp = 0; jp < nb_dof2; jp++)
              if (this->mesh_num.Element(num_elem).GetNumberDof(j)
                  == this->mesh_num.Element(num_elem2).GetNumberDof(jp))
                j2 = jp;
            
            if (j2 != -1)
              {
                if (abs(phi1_n(j) - phi2_n(j2)) > 1e-5)
                  {
                    cout << "Tangential continuity not ensured "<< endl;
                    DISP(phi(j)); DISP(phi2(j2));
                    DISP(phi1_n(j)); DISP(phi2_n(j2)); DISP(j); DISP(j2); DISP(rot);
                    DISP(num_loc); DISP(num_loc2);
                    abort();
                  }
              }
          }
      }
#endif
  }
  
  
  //! checks if continuity is ensured between elements  
  template<>
  void VarGeometryProblem<Dimension2>::
  CheckLocalContinuity(int i, int num_elem, int num_elem2,
                       SetPoints<Dimension2>& PointsElem, SetPoints<Dimension2>& PointsElem2,
                       SetMatrices<Dimension2>& MatricesElem,
                       SetMatrices<Dimension2>& MatricesElem2,
                       const ElementReference<Dimension2, 2>& Fb,
                       const ElementReference<Dimension2, 2>& Fb2)
  {
#ifdef MONTJOIE_WITH_NODAL_HCURL
    // we check continuity of basis functions across boundaries
    Vector<R2> s, s2; R2 tmp;
    
    this->mesh.GetVerticesElement(num_elem, s);
    this->mesh.GetVerticesElement(num_elem2, s2);
    int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
    int num_loc2 = -1, rot = 0;
    this->mesh_num.GetBoundaryRotation(i, num_elem, num_elem2, num_loc, num_loc2, rot);
    
    Fb.FjElem(s, PointsElem, this->mesh, num_elem);
    Fb.DFjElem(s, PointsElem, MatricesElem, this->mesh, num_elem);
    
    Fb.FjSurfaceElem(s, PointsElem, this->mesh, num_elem, num_loc);
    Fb.DFjSurfaceElem(s, PointsElem, MatricesElem, this->mesh, num_elem, num_loc);
    
    Fb2.FjElem(s2, PointsElem2, this->mesh, num_elem2);
    Fb2.DFjElem(s2, PointsElem2, MatricesElem2, this->mesh, num_elem2);
    
    Fb2.FjSurfaceElem(s2, PointsElem2, this->mesh, num_elem2, num_loc2);
    Fb2.DFjSurfaceElem(s2, PointsElem2, MatricesElem2, this->mesh, num_elem2, num_loc2);
    
    int nb_quad = Fb.GetNbQuadBoundary(num_loc);
    for (int k = 0; k < nb_quad; k++)
      {
        int nk = Fb.GetQuadNumber(num_loc, k);
        int k2 = -1;
        for (int kp = 0; kp < nb_quad; kp++)
          {
            int nkp = Fb2.GetQuadNumber(num_loc2, kp);
            if (PointsElem.GetPointQuadrature(nk) == PointsElem2.GetPointQuadrature(nkp))
              k2 = kp;
          }
        
        int nk2 = Fb2.GetQuadNumber(num_loc2, k2);
        const Vector<R2>& pts = Fb.PointsND();
        const Vector<R2>& pts2 = Fb2.PointsND();
        
        // we check if normales are opposite
        R2 normale, normale2;
        normale = MatricesElem.GetNormaleQuadratureBoundary(k);
        normale2 = MatricesElem2.GetNormaleQuadratureBoundary(k2);
        Mlt(1.0/Norm2(normale), normale);
        Mlt(-1.0/Norm2(normale2), normale2);
        if (normale != normale2)
          {
            cout << "Normale are not the same "<<endl;
            abort();
          }
        
        // now we compute values of phi on the quadrature point
        Matrix2_2 dfjm1;
        Vector<TinyVector<Real_wp, 2> > phi, phi2;
        GetInverse(MatricesElem.GetPointQuadrature(nk), dfjm1);
        Fb.ComputeValuesPhi(pts(nk), phi, dfjm1, this->mesh_num, num_elem);
        
        GetInverse(MatricesElem2.GetPointQuadrature(nk2), dfjm1);
        Fb2.ComputeValuesPhi(pts2(nk2), phi2, dfjm1, this->mesh_num, num_elem2);
        
        VectReal_wp phi1_n(phi.GetM()), phi2_n(phi2.GetM());
        int nb_dof = mesh_num.GetNbLocalDof(num_elem);
        int nb_dof2 = mesh_num.GetNbLocalDof(num_elem2);		
        for (int j = 0; j < nb_dof; j++)
          phi1_n(j) = phi(j)(0)*normale(1) - phi(j)(1)*normale(0);
        
        for (int j = 0; j < nb_dof2; j++)
          phi2_n(j) = phi2(j)(0)*normale(1) - phi2(j)(1)*normale(0);
        
        // we check continuity of phi \times n
        for (int j = 0; j < nb_dof; j++)
          {
            int j2 = -1;
            for (int jp = 0; jp < nb_dof2; jp++)
              if (this->mesh_num.Element(num_elem).GetNumberDof(j)
                  == this->mesh_num.Element(num_elem2).GetNumberDof(jp))
                j2 = jp;

            if (j2 != -1)
              {
                if (abs(phi1_n(j) - phi2_n(j2)) > 1e-5)
                  {
                    cout << "Tangential continuity not ensured "<< endl;
                    DISP(phi(j)); DISP(phi2(j2));
                    DISP(phi1_n(j)); DISP(phi2_n(j2)); DISP(j); DISP(j2); DISP(rot);
                    abort();
                  }
              }
          }
      }
    
#endif
  }
  

  //! checks if continuity is ensured between elements  
  template<>
  void VarGeometryProblem<Dimension3>::
  CheckLocalContinuity(int i, int num_elem, int num_elem2,
                       SetPoints<Dimension3>& PointsElem, SetPoints<Dimension3>& PointsElem2,
                       SetMatrices<Dimension3>& MatricesElem,
                       SetMatrices<Dimension3>& MatricesElem2,
                       const ElementReference<Dimension3, 2>& Fb,
                       const ElementReference<Dimension3, 2>& Fb2)
  {
#ifdef MONTJOIE_WITH_NODAL_HCURL
#ifdef MONTJOIE_WITH_THREE_DIM
    Vector<R3> s, s2; R3 tmp; Matrix3_3 mat_dfj1, mat_dfj2, invDFj1, invDFj2;
    Real_wp dsj1, dsj2;
    this->mesh.GetVerticesElement(num_elem, s);
    this->mesh.GetVerticesElement(num_elem2, s2);
    
    int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);	    
    int num_loc2 = this->mesh.Element(num_elem2).GetPositionBoundary(i);
    int rot1 = this->mesh.Element(num_elem).GetOrientationBoundary(num_loc);
    int rot2 = this->mesh.Element(num_elem2).GetOrientationBoundary(num_loc2);
    int rot = MeshNumbering<Dimension3>::
      GetRotationFace(rot1, rot2, this->mesh.Boundary(i).GetNbVertices());
    
    Fb.FjElemNodal(s, PointsElem, this->mesh, num_elem);
    Fb2.FjElemNodal(s2, PointsElem2, this->mesh, num_elem2);
    
    if (this->FormulationDG() == ElementReference_Base::HDG)
      {
        VarProblem<Dimension3>& var_p = dynamic_cast<VarProblem<Dimension3>& >(*this);
        int num_face = this->mesh.Element(num_elem).numBoundary(num_loc);
        const ElementReference<Dimension2, 2>& Fb_s = dynamic_cast<const ElementReference<Dimension2, 2>& >(var_p.GetSurfaceFiniteElement(num_face));
        
        int N = Fb_s.GetNbPointsNodalElt();
        Vector<R3> Pts1(N), Pts2(N);
        for (int k = 0; k < N; k++)
          {
            Fb.GetLocalCoordOnBoundary(num_loc, Fb_s.PointsNodalND(k), tmp);
            Fb.Fj(s, PointsElem, tmp, Pts1(k), this->mesh, num_elem);
            
            Fb2.GetLocalCoordOnBoundary(num_loc2, Fb_s.PointsNodalND(k), tmp);
            Fb2.Fj(s2, PointsElem2, tmp, Pts2(k), this->mesh, num_elem2);
          }
        
        R3 normale, normale2;
        VectR3 phi1, phi2; VectR2 phi_surf1, phi_surf2;
        VectR3 phi_vol1, phi_vol2;
        for (int k = 0; k < N; k++)
          {
            int k2 = -1;
            for (int kp = 0; kp < Pts2.GetM(); kp++)
              if (Pts1(k) == Pts2(kp))
                k2 = kp;
            
            if (k2 == -1)
              {
                cout << "Impossible" << endl;
                DISP(Pts1); DISP(Pts2); DISP(k); DISP(Pts1(k));
                abort();
              }
            
            Fb.GetLocalCoordOnBoundary(num_loc, Fb_s.PointsNodalND(k), tmp);
            Fb.DFj(s, PointsElem, tmp, mat_dfj1, this->mesh, num_elem);
            
            Fb2.GetLocalCoordOnBoundary(num_loc2, Fb_s.PointsNodalND(k2), tmp);
            Fb2.DFj(s2, PointsElem2, tmp, mat_dfj2, this->mesh, num_elem2);
            
            GetInverse(mat_dfj1, invDFj1);
            GetInverse(mat_dfj2, invDFj2);
            
            Fb.GetNormale(invDFj1, normale, dsj1, num_loc);
            Fb2.GetNormale(invDFj2, normale2, dsj2, num_loc2);
            Mlt(1.0/Norm2(normale), normale);
            Mlt(-1.0/Norm2(normale2), normale2);
            if (normale != normale2)
              {
                cout << "Normale are not the same "<<endl;
                abort();
              }
            
            Fb_s.ComputeValuesPhiRef(Fb_s.PointsNodalND(k), phi_surf1);
            Fb_s.ComputeValuesPhiRef(Fb_s.PointsNodalND(k2), phi_surf2);
            //DISP(k); DISP(k2); DISP(Fb_s.PointsNodalND(k)); DISP(Fb_s.PointsNodalND(k2));
            //DISP(phi_surf1); DISP(phi_surf2);
            
            IVect num_dof1(Fb_s.GetNbDof()), num_dof2(Fb_s.GetNbDof());
            int offset = 0;
            for (int j = 0; j < num_loc; j++)
              {
                int nf = this->mesh.Element(num_elem).numBoundary(j);
                offset += var_p.GetSurfaceFiniteElement(nf).GetNbDof();
              }
            
            int offset2 = 0;
            for (int j = 0; j < num_loc2; j++)
              {
                int nf = this->mesh.Element(num_elem2).numBoundary(j);
                offset2 += var_p.GetSurfaceFiniteElement(nf).GetNbDof();
              }
            
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              num_dof1(j) = this->mesh_num.Element(num_elem).GetNumberDof(offset+j);
            
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              num_dof2(j) = this->mesh_num.Element(num_elem2).GetNumberDof(offset2+j);
            
            //DISP(offset); DISP(offset2); DISP(num_dof1); DISP(num_dof2);
            phi_vol1.Reallocate(this->mesh_num.Element(num_elem).GetNbDof());
            phi_vol2.Reallocate(this->mesh_num.Element(num_elem2).GetNbDof());
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              {
                tmp = Fb.GetTangentialVector(num_loc, phi_surf1(j));
                MltTrans(invDFj1, tmp, phi_vol1(offset+j));
                
                tmp = Fb2.GetTangentialVector(num_loc2, phi_surf2(j));
                MltTrans(invDFj2, tmp, phi_vol2(offset2+j));
              }
            
            this->mesh_num.number_map.ModifyLocalUnknownVector(this->mesh_num, phi_vol1, num_elem);
            this->mesh_num.number_map.ModifyLocalUnknownVector(this->mesh_num, phi_vol2, num_elem2);
            
            phi1.Reallocate(Fb_s.GetNbDof());
            phi2.Reallocate(Fb_s.GetNbDof());
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              {
                phi1(j) = phi_vol1(offset+j);
                phi2(j) = phi_vol2(offset2+j);
              }
            
            //DISP(phi1); DISP(phi2);
            
            Sort(num_dof1, phi1);
            Sort(num_dof2, phi2);
            
            //DISP(phi1); DISP(phi2);
            TinyMatrix<Real_wp, Symmetric, 3, 3> ProjTang;
            GetTangentialProjector(normale, ProjTang);
            R3 vec_u, vec_v;
            //DISP(normale); DISP(ProjTang);
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              {
                Mlt(ProjTang, phi1(j), vec_u);
                Mlt(ProjTang, phi2(j), vec_v);
                if (vec_u.Distance(vec_v) > 1e-10)
                  {
                    cout << "tangential component of basis functions not continuous" << endl;
                    DISP(j); DISP(vec_u); DISP(vec_v);
                    abort();
                  }
              }
          }
        
        return;
      }
    
    int N = Fb.GetNbPointsNodalElt();
    Vector<R3> Pts1(N), Pts2(Fb2.GetNbPointsNodalElt());
    for (int k = 0; k < N; k++)
      Fb.Fj(s, PointsElem, Fb.PointsNodalND(k), Pts1(k), this->mesh, num_elem);
    
    for (int k = 0; k < Pts2.GetM(); k++)
      Fb2.Fj(s2, PointsElem2, Fb2.PointsNodalND(k), Pts2(k), this->mesh, num_elem2);
    
    R3 normale, normale2;
    for (int k = 0; k < N; k++)
      {
        int k2 = -1;
        for (int kp = 0; kp < Pts2.GetM(); kp++)
          if (Pts1(k) == Pts2(kp))
            k2 = kp;
        
        // apex of a pyramid is avoided
        if ((Fb.GetNbVertices() == 5) && (k == 4))
          k2 = -1;
        
        // apex of a pyramid is avoided
        if ((Fb2.GetNbVertices() == 5) && (k2 == 4))
          k2 = -1;
        
        if (k2 != -1)
          {
            Fb.DFj(s, PointsElem, Fb.PointsNodalND(k), mat_dfj1, this->mesh, num_elem);
            Fb2.DFj(s2, PointsElem2, Fb2.PointsNodalND(k2), mat_dfj2, this->mesh, num_elem2);
            
            GetInverse(mat_dfj1, invDFj1);
            GetInverse(mat_dfj2, invDFj2);
            
            Fb.GetNormale(invDFj1, normale, dsj1, num_loc);
            Fb2.GetNormale(invDFj2, normale2, dsj2, num_loc2);
            Mlt(1.0/Norm2(normale), normale);
            Mlt(-1.0/Norm2(normale2), normale2);
            if (normale != normale2)
              {
                cout << "Normale are not the same "<<endl;
                abort();
              }
            
            // now we compute values of phi on the quadrature point
            Vector<TinyVector<Real_wp, 3> > phi, phi2;
            Fb.ComputeValuesPhi(Fb.PointsNodalND(k), phi, invDFj1,
                                this->mesh_num, num_elem);
            
            Fb2.ComputeValuesPhi(Fb2.PointsNodalND(k2), phi2, invDFj2,
                                 this->mesh_num, num_elem2);
            
            VectR3 phib_1 = phi, phib_2 = phi2;
            R3 phi1_n, phi2_n;
            int nb_dof = Fb.GetNbDof();
            int nb_dof2 = Fb2.GetNbDof();
            for (int j = 0; j < nb_dof; j++)
              {
                TimesProd(phi(j), normale, phi1_n);
                phi(j) = phi1_n;
              }
            
            for (int j = 0; j < nb_dof2; j++)
              {
                TimesProd(phi2(j), normale, phi2_n);
                phi2(j) = phi2_n;
              }
            
            // we check continuity of phi \times n
            for (int j = 0; j < nb_dof; j++)
              if (this->mesh_num.Element(num_elem).GetNumberDof(j) != -1)
                {
                  int j2 = -1;
                  for (int jp = 0; jp < nb_dof2; jp++)
                    if (this->mesh_num.Element(num_elem).GetNumberDof(j)
                        == this->mesh_num.Element(num_elem2).GetNumberDof(jp))
                      j2 = jp;
                  
                  if (j2 != -1)
                    {
                      if (phi(j).Distance(phi2(j2)) > 1e-5)
                        {
                          cout << "Tangential continuity not ensured "<< endl;
                          DISP(num_loc); DISP(num_loc2); DISP(num_elem); DISP(num_elem2);
                          DISP(this->mesh.Element(num_elem).GetOrientationFace(num_loc));
                          DISP(this->mesh.Element(num_elem2).GetOrientationFace(num_loc2));
                          DISP(phi(j)); DISP(phi2(j2));
                          DISP(phi(j+1)); DISP(phi2(j2+1));
                          DISP(j); DISP(j2); DISP(rot);
                          abort();
                        }
                    }
                }
          }
      }    
#endif
#endif
  }

    //! checks if continuity is ensured between elements
  template<class Dimension>
  void VarGeometryProblem<Dimension>::
  CheckLocalContinuity(int i, int num_elem, int num_elem2,
                       SetPoints<Dimension>& PointsElem, SetPoints<Dimension>& PointsElem2,
                       SetMatrices<Dimension>& MatricesElem,
                       SetMatrices<Dimension>& MatricesElem2,
                       const ElementReference<Dimension, 1>& Fb,
                       const ElementReference<Dimension, 1>& Fb2)
  {
#ifdef MONTJOIE_WITH_NODAL_H1
    Vector<R_N> s, s2; R_N tmp; MatrixN_N dfjm1;
    this->mesh.GetVerticesElement(num_elem, s);
    this->mesh.GetVerticesElement(num_elem2, s2);
    
    int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);	    
    int num_loc2 = this->mesh.Element(num_elem2).GetPositionBoundary(i);
    int rot1 = this->mesh.Element(num_elem).GetOrientationBoundary(num_loc);
    int rot2 = this->mesh.Element(num_elem2).GetOrientationBoundary(num_loc2);
    int rot = this->mesh_num.GetRotationFace(rot1, rot2, mesh.Boundary(i).GetNbVertices());
    
    Fb.FjElemNodal(s, PointsElem, this->mesh, num_elem);
    Fb2.FjElemNodal(s2, PointsElem2, this->mesh, num_elem2);
    
    if (this->FormulationDG() == ElementReference_Base::HDG)
      {
        VarProblem<Dimension>& var_p = dynamic_cast<VarProblem<Dimension>& >(*this);
        typedef typename Dimension::DimensionBoundary DimensionB;
        int num_face = this->mesh.Element(num_elem).numBoundary(num_loc);
        const ElementReference<DimensionB, 1>& Fb_s = var_p.GetSurfaceFiniteElementH1(num_face);
        int N = Fb_s.GetNbPointsNodalElt();
        Vector<R_N> Pts1(N), Pts2(N);
        for (int k = 0; k < N; k++)
          {
            Fb.GetLocalCoordOnBoundary(num_loc, Fb_s.PointsNodalND(k), tmp);
            Fb.Fj(s, PointsElem, tmp, Pts1(k), this->mesh, num_elem);
            
            Fb2.GetLocalCoordOnBoundary(num_loc2, Fb_s.PointsNodalND(k), tmp);
            Fb2.Fj(s2, PointsElem2, tmp, Pts2(k), this->mesh, num_elem2);
          }
        
        VectReal_wp phi1, phi2;
        for (int k = 0; k < N; k++)
          {
            int k2 = -1;
            for (int kp = 0; kp < Pts2.GetM(); kp++)
              if (Pts1(k) == Pts2(kp))
                k2 = kp;
            
            if (k2 == -1)
              {
                cout << "Impossible" << endl;
                DISP(Pts1); DISP(Pts2); DISP(k); DISP(Pts1(k));
                abort();
              }
            
            Fb_s.ComputeValuesPhiRef(Fb_s.PointsNodalND(k), phi1);
            Fb_s.ComputeValuesPhiRef(Fb_s.PointsNodalND(k2), phi2);
            
            IVect num_dof1(Fb_s.GetNbDof()), num_dof2(Fb_s.GetNbDof());
            int offset = 0;
            for (int j = 0; j < num_loc; j++)
              {
                int nf = this->mesh.Element(num_elem).numBoundary(j);
                offset += var_p.GetSurfaceFiniteElement(nf).GetNbDof();
              }
            
            int offset2 = 0;
            for (int j = 0; j < num_loc2; j++)
              {
                int nf = this->mesh.Element(num_elem2).numBoundary(j);
                offset2 += var_p.GetSurfaceFiniteElement(nf).GetNbDof();
              }
            
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              num_dof1(j) = this->mesh_num.Element(num_elem).GetNumberDof(offset+j);
            
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              num_dof2(j) = this->mesh_num.Element(num_elem2).GetNumberDof(offset2+j);
            
            Sort(num_dof1, phi1);
            Sort(num_dof2, phi2);
            for (int j = 0; j < Fb_s.GetNbDof(); j++)
              if (abs(phi1(j) - phi2(j)) > 1e-10)
                {
                  cout << "basis functions not continuous" << endl;
                  DISP(j); DISP(phi1(j)); DISP(phi2(j));
                  abort();
                }
          }
        
        return;
      }
    
    int N = Fb.GetNbPointsNodalElt();
    Vector<R_N> Pts1(N), Pts2(Fb2.GetNbPointsNodalElt());
    for (int k = 0; k < N; k++)
      Fb.Fj(s, PointsElem, Fb.PointsNodalND(k), Pts1(k), this->mesh, num_elem);
    
    for (int k = 0; k < Pts2.GetM(); k++)
      Fb2.Fj(s2, PointsElem2, Fb2.PointsNodalND(k), Pts2(k), this->mesh, num_elem2);
    
    Vector<TinyVector<Real_wp, 1> > phi, phi2;
    for (int k = 0; k < N; k++)
      {
        int k2 = -1;
        for (int kp = 0; kp < Pts2.GetM(); kp++)
          if (Pts1(k) == Pts2(kp))
            k2 = kp;
        
        if (k2 != -1)
          {
            Fb.ComputeValuesPhi(Fb.PointsNodalND(k), phi, dfjm1,
                                this->mesh_num, num_elem);
            
            Fb2.ComputeValuesPhi(Fb2.PointsNodalND(k2), phi2, dfjm1,
                                 this->mesh_num, num_elem2);
            
            // we check continuity
            int nb_dof = Fb.GetNbDof();
            int nb_dof2 = Fb2.GetNbDof();
            for (int j = 0; j < nb_dof; j++)
              if (this->mesh_num.Element(num_elem).GetNumberDof(j) != -1)
                {
                  int j2 = -1;
                  for (int jp = 0; jp < nb_dof2; jp++)
                    if (this->mesh_num.Element(num_elem).GetNumberDof(j)
                        == this->mesh_num.Element(num_elem2).GetNumberDof(jp))
                      j2 = jp;
                  
                  if (j2 != -1)
                    {
                      if (abs(phi(j)(0)-phi2(j2)(0)) > 1e-12)
                        {
                          DISP(j); DISP(j2); DISP(rot);
                          DISP(phi(j)(0)); DISP(phi2(j2)(0));
                          abort();
                        }
                    }
                }
          }
      }
#endif 
  }

}

#define MONTJOIE_FILE_VAR_GEOMETRY_PROBLEM_CXX
#endif
