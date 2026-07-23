#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_CXX

namespace Montjoie
{
  
  //! Sets values of attributes to default values
  template<class Dimension>
  void VarTransmission_Base<Dimension>::InitDefaultValues()
  {
    OrdreTransmission = 1;
    R0 = 1.0;
    nb_sections = 1;
    RefOmegaMoins = 1;
    RefOmegaPlus = 2;
    RefGammaPlus = -1;
    RefGammaMoins = -1;
    
    Decalage = 0.5;
    delta = 0.0;
    
    ref_omega_layer = 3;
    exact_solution = false;
    duplicate_extremity = true;
    gamma_is_plate = false;
    gamma_is_sphere = false;
    nb_ddl_volume_on_surface = 0;
  }
  

  template<class Dimension>
  size_t VarTransmission_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = 10*sizeof(int) + 3*sizeof(Real_wp) + 4*sizeof(bool);
    taille += transmission_references.GetMemorySize() + mesh_plate.GetMemorySize()
      + mesh_surf.GetMemorySize() + mesh_plate_num.GetMemorySize()
      + NumFaceMinus.GetMemorySize() + NumFacePlus.GetMemorySize()
      + DdlVolMinus.GetMemorySize() + DdlVolPlus.GetMemorySize() + DdlSurfLambda.GetMemorySize()
      + Seldon::GetMemorySize(NodleSurf) + Seldon::GetMemorySize(NodleAcous)
      + Seldon::GetMemorySize(NodleMeca) + Seldon::GetMemorySize(NodleVolMinus)
      + Seldon::GetMemorySize(NodleVolPlus) + ProcVolMinus.GetMemorySize()
      + ProcVolPlus.GetMemorySize() + MatchingNumber_Subdomain.GetMemorySize()
      + Seldon::GetMemorySize(MatchingDofOrig_Subdomain) + Seldon::GetMemorySize(JacobianAcous)
      + DirichletAcous.GetMemorySize() + OldNumFaceSurf.GetMemorySize()
      + OldNumVertexSurf.GetMemorySize() + IndexVert.GetMemorySize();

    taille += Seldon::GetMemorySize(NormaleNodes) + Seldon::GetMemorySize(PosNodes)
      + Seldon::GetMemorySize(DsNodes) + Seldon::GetMemorySize(MatDfjNodes)
      + nb_pts_quad_local.GetMemorySize() + coef_dfj_quad_local.GetMemorySize();
    
    return taille;
  }
  
  
  //! modifying transmissions conditions with a line of the data file
  template<class Dimension>
  void VarTransmission_Base<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("OrdreTransmission") ) 
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "OrdreTransmission needs more parameters, for instance :" << endl;
	    cout << "OrdreTransmission = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	OrdreTransmission = to_num<int>(parameters(0) );
      }    
    else if ( !description_field.compare("R0") )
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "R0 needs more parameters, for instance :" << endl;
	    cout << "R0 = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	R0 = to_num<Real_wp>(parameters(0) );
      }
    else if ( !description_field.compare("NbSectionsTransmission") )
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "NbSectionsTransmission needs more parameters, for instance :" << endl;
	    cout << "NbSectionsTransmission = N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	nb_sections = to_num<int>(parameters(0) );
      }
    else if ( !description_field.compare("DuplicateExtremity") )
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "DuplicateExtremity needs more parameters, for instance :" << endl;
	    cout << "DuplicateExtremity = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          duplicate_extremity = true;
        else
          duplicate_extremity = false;
      }
    else if ( !description_field.compare("TransmissionSurface") )
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "TransmissionSurface needs more parameters, for instance :" << endl;
	    cout << "TransmissionSurface = Sphere" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("Sphere"))
          gamma_is_sphere = true;
        else
          gamma_is_sphere = false;
      }
    else if (!description_field.compare("TransmissionExactSolution") )
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "TransmissionExactSolution needs more parameters, for instance :" << endl;
	    cout << "TransmissionExactSolution = YES ref" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        exact_solution = false;
        if (!parameters(0).compare("YES"))
          {
            exact_solution = true;
            ref_omega_layer = to_num<int>(parameters(1));
          }
      }
    else if (!description_field.compare("OmegaMinusPlus"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "OmegaMinusPlus needs more parameters, for instance :" << endl;
	    cout << "OmegaMinusPlus = ref_minus ref_plus" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        RefOmegaMoins =  to_num<int>(parameters(0));
        RefOmegaPlus =  to_num<int>(parameters(1));
      }
    else if (!description_field.compare("GammaPlate"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Base" << endl;
	    cout << "GammaPlate needs more parameters, for instance :" << endl;
	    cout << "GammaPlate = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          gamma_is_plate = true;
        else
          gamma_is_plate = false;
      }
    
  }
  
  
  //! sets column numbers for impedance boundary condition
  template<class Dimension>
  void VarTransmission_Base<Dimension>::SetModifiedColNumbers(int num)
  {
    int Nvol = var_problem.offset_dof_unknown(1);
    int Nsurf = this->DdlVolPlus.GetM();    
    int Nvol_mesh = var_problem.GetMeshNumbering(0).GetNbDof();
    int offset;
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif
    
    if (nb_proc == 1)
      {	
        if (num == 0)
          {
	    for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  var_boundary.NewColumnNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                    = this->DdlVolMinus(i) + m*Nvol;
		  
                  var_boundary.NewColumnNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                    = this->DdlVolPlus(i) + m*Nvol;
		}
          }
        else
          {
            for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  var_boundary.NewColumnNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                    = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		  
                  var_boundary.NewColumnNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                    = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		}
          }
      }
    else
      {
        int Nvol_glob = var_problem.GetNbGlobalUnknownDof(1);
        if (num == 0)
          {
	    for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  if (ProcVolPlus(i) == rank_proc)
		    {
                      if (this->ProcVolMinus(i) == rank_proc)
                        offset = m*Nvol;
                      else
                        offset = m*Nvol_glob;
                      
		      var_boundary.NewColumnNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = this->DdlVolMinus(i) + offset;
		      
                      var_boundary.ProcColumnNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = this->ProcVolMinus(i);
		    }
		  
		  if (ProcVolMinus(i) == rank_proc)
		    {
                      if (this->ProcVolPlus(i) == rank_proc)
                        offset = m*Nvol;
                      else
                        offset = m*Nvol_glob;
                      
		      var_boundary.NewColumnNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = this->DdlVolPlus(i) + offset;
		      
                      var_boundary.ProcColumnNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = this->ProcVolPlus(i);
		    }
		}
          }
        else
          {
	    for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  if (ProcVolPlus(i) == rank_proc)
		    {
		      var_boundary.NewColumnNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		      
                      var_boundary.ProcColumnNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = rank_proc;
		    }
		  
		  if (ProcVolMinus(i) == rank_proc)
		    {
		      var_boundary.NewColumnNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		      
                      var_boundary.ProcColumnNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = rank_proc;
		    }
		}
          }        
      }
  }
  

  //! sets row numbers for impedance boundary condition  
  template<class Dimension>
  void VarTransmission_Base<Dimension>::SetModifiedRowNumbers(int num)
  {
    int Nvol = var_problem.offset_dof_unknown(1);
    int Nvol_mesh = var_problem.GetMeshNumbering(0).GetNbDof();
    int Nsurf = this->DdlVolPlus.GetM();
    int offset;

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int rank_proc(0), nb_proc(1);
#endif

    if (nb_proc == 1)
      {
        if (num == 0)
          {
            for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  var_boundary.NewRowNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                    = this->DdlVolMinus(i) + m*Nvol;
		  
                  var_boundary.NewRowNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                    = this->DdlVolPlus(i) + m*Nvol;
		}
          }
        else
          {
	    for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  var_boundary.NewRowNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                    = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
                  
		  var_boundary.NewRowNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                    = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		}
          }
      }
    else
      {
        int Nvol_glob = var_problem.GetNbGlobalUnknownDof(0);
        if (num == 0)
          {
	    for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  if (ProcVolPlus(i) == rank_proc)
		    {
                      if (this->ProcVolMinus(i) == rank_proc)
                        offset = m*Nvol;
                      else
                        offset = m*Nvol_glob;
                      
		      var_boundary.NewRowNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = this->DdlVolMinus(i) + offset;
		      
                      var_boundary.ProcRowNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = this->ProcVolMinus(i);
		    }
		  
		  if (ProcVolMinus(i) == rank_proc)
		    {
                      if (this->ProcVolPlus(i) == rank_proc)
                        offset = m*Nvol;
                      else
                        offset = m*Nvol_glob;
                      
		      var_boundary.NewRowNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = this->DdlVolPlus(i) + offset;
		      
                      var_boundary.ProcRowNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = this->ProcVolPlus(i);
		    }
		}
          }
        else
          {
	    for (int m = 0; m < var_problem.nb_unknowns; m++)
	      for (int i = 0; i < DdlVolPlus.GetM(); i++)
		{
		  if (ProcVolPlus(i) == rank_proc)
		    {
		      var_boundary.NewRowNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		      
                      var_boundary.ProcRowNumbers_Impedance(this->DdlVolPlus(i) + m*Nvol)
                        = rank_proc;
		    }
		  
		  if (ProcVolMinus(i) == rank_proc)
		    {
		      var_boundary.NewRowNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = m*Nvol + Nvol_mesh + (num-1)*Nsurf + i;
		      
                      var_boundary.ProcRowNumbers_Impedance(this->DdlVolMinus(i) + m*Nvol)
                        = rank_proc;
		    }
		}
          }        
      }
  }
  
  
  //! separation of mesh into two parts  
  template<>
  void VarTransmission_Base<Dimension2>::PartMeshTransmission()
  {
    if (transmission_references.GetM() <= 0)
      return;
    
    if (var_problem.print_level >= 2)
      cout << "Parting mesh because of transmission condition ... " << endl;
    
    Mesh<Dimension2>& mesh = var_problem.mesh;
    
    bool store_normale_gamma = false;
    
    int ref_plate = 1;
    IVect ref_cond(mesh.GetNbReferences()); ref_cond.Fill(0);
    for (int i = 0; i < transmission_references.GetM(); i++)
      ref_cond(transmission_references(i)) = ref_plate;
  
    // we are extracting the surface mesh
    IVect NumEdge, NumVertices, NumElement, NumLocalEdge;
    mesh.GetBoundaryMesh(ref_plate, mesh_surf, NumEdge, NumVertices, 
                         NumElement, NumLocalEdge, ref_cond);
  
    OldNumFaceSurf = NumEdge;
    OldNumVertexSurf = NumVertices;
    mesh_surf.FindConnectivity();
    
    // detection of vertices to duplicate
    Vector<bool> PointToDuplicate(mesh_surf.GetNbVertices());
    PointToDuplicate.Fill(true);

    IndexVert.Reallocate(mesh.GetNbVertices()); IndexVert.Fill();
    int nb_old_vertices = mesh.GetNbVertices();
    int nb_vertices = nb_old_vertices;
    for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
      if (PointToDuplicate(i))
        {
          IndexVert(NumVertices(i)) = nb_vertices;
          nb_vertices++;
        }
    
    // references of domains must be different
    if (RefOmegaMoins == RefOmegaPlus)
      {
        cout << "Use different references for Omega- and Omega+ " << endl;
        abort();
      }
    
    int r = mesh_surf.GetGeometryOrder();
    int nb_vert_surf = mesh_surf.GetNbVertices();
    int nb_edge_surf = mesh_surf.GetNbBoundaryRef();
    // allocating arrays that will contain normale associated with vertices,
    // nodal points of edges, and nodal points of faces
    VectR2 Normales(nb_vert_surf) ;
    Matrix<R2> NormaleEdge;
    if (r > 1)
      NormaleEdge.Reallocate(mesh.GetNbBoundaryRef(), r-1);
    
    // computation of normales for each vertex, edge and face
    // normales are outward to domain Omega-
    VectR2 S ;
    SetPoints<Dimension2> Points ;
    SetMatrices<Dimension2> Matrices ;
    Matrix2_2 mat_dfj, dfjm1; R2 normale_fj; Real_wp dsj;
    if (store_normale_gamma)
      {
        NormaleNodes.Reallocate(nb_edge_surf);
        DsNodes.Reallocate(nb_edge_surf);
      }
    
    MatDfjNodes.Reallocate(mesh.GetNbElt());
    for (int i = 0; i < nb_edge_surf; i++)
      {
        int Ne = NumElement(i) ;
        int NumLoc =  NumLocalEdge(i) ;
        mesh.GetVerticesElement(Ne, S) ;
        int nb_vert = mesh.Element(Ne).GetNbVertices();
        mesh.FjElemNodal(S, Points, mesh, Ne) ;
        mesh.DFjElemNodal(S, Points, Matrices, mesh, Ne) ;
        bool same_rot = mesh.Element(Ne).GetOrientationEdge(NumLoc);

	// storing normales, ds
        if (store_normale_gamma)
          {
            NormaleNodes(i).Reallocate(r+1);
            DsNodes(i).Reallocate(r+1);
          }
	
        for (int j = 0; j <= r; j++)
          {
            int node = -1;
            if (j <= 1)
              {
                if (same_rot)
                  node = (NumLoc+j)%nb_vert;
                else
                  node = (NumLoc+1-j)%nb_vert;
              }
            else
              {
                if (same_rot)
                  node = nb_vert + NumLoc*(r-1) + j - 2;
                else
                  node = nb_vert + NumLoc*(r-1) + r - j;
              }
	    
            mat_dfj = Matrices.GetPointNodal(node);
            GetInverse(mat_dfj, dfjm1);
            mesh.GetNormale(dfjm1, normale_fj, dsj, Ne, NumLoc);
            
            if (RefOmegaMoins > RefOmegaPlus)
              Mlt(-1.0, normale_fj);
            
            if (j < 2)
              Normales(mesh_surf.BoundaryRef(i).numVertex(j)) = normale_fj;
            else
              NormaleEdge(i, j-2) = normale_fj;
	    
	    // storing normale, ds
	    int jloc = j-1;
	    if (j == 0)
	      jloc = 0;
	    else if (j == 1)
	      jloc = r;
	    
            if (store_normale_gamma)
              {
                NormaleNodes(i)(jloc) = normale_fj;
                DsNodes(i)(jloc) = dsj;
              }
	  }

	// storing matrix dfj on the boundary
	int nb_faces_loc = mesh.Element(Ne).GetNbBoundary();
	MatDfjNodes(Ne).Reallocate(nb_faces_loc);
	MatDfjNodes(Ne)(NumLoc).Reallocate(r+1);
	for (int j = 0; j <= r; j++)
	  {
	    int node = mesh.GetNodalNumber(Ne, NumLoc, j);
	    MatDfjNodes(Ne)(NumLoc)(j) = Matrices.GetPointNodal(node);
	  }	
	
	// and on the opposite side too
	int num_edge = mesh.Element(Ne).numBoundary(NumLoc);
	int Ne2 = mesh.Boundary(num_edge).numElement(0);
	if (Ne2 == Ne)
	  Ne2 = mesh.Boundary(num_edge).numElement(1);
	
	int NumLoc2 = mesh.Element(Ne2).GetPositionBoundary(num_edge);
	nb_faces_loc = mesh.Element(Ne2).GetNbBoundary();
        mesh.GetVerticesElement(Ne2, S) ;
        mesh.FjElemNodal(S, Points, mesh, Ne2) ;
        mesh.DFjElemNodal(S, Points, Matrices, mesh, Ne2) ;
	
	MatDfjNodes(Ne2).Reallocate(nb_faces_loc);
	MatDfjNodes(Ne2)(NumLoc2).Reallocate(r+1);
	for (int j = 0; j <= r; j++)
	  {
	    int node = mesh.GetNodalNumber(Ne2, NumLoc2, j);
	    MatDfjNodes(Ne2)(NumLoc2)(j) = Matrices.GetPointNodal(node);
	  }
      }
    
    // adding vertices to the 2-D mesh
    mesh.ResizeVertices(nb_vertices);
    nb_vertices = nb_old_vertices;
    for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
      if (PointToDuplicate(i))
        {
          // vertices are shifted with respect to normale
          mesh.Vertex(nb_vertices) = mesh_surf.Vertex(i);
          mesh.Vertex(NumVertices(i)) -= Decalage*delta*Normales(i); 
          mesh.Vertex(nb_vertices) += Decalage*delta*Normales(i); 
          nb_vertices++;
        }
        
    // copying faces
    int nb_edges_ref = mesh_surf.GetNbBoundaryRef();
    // and duplicating faces for volume mesh
    int nb_old_edge2d = mesh.GetNbBoundaryRef(); 
    int nb_edge2d = nb_old_edge2d;
    mesh.ResizeBoundariesRef(nb_old_edge2d + nb_edges_ref);
    NumFaceMinus.Reallocate(nb_edges_ref);
    NumFacePlus.Reallocate(nb_edges_ref);
    
    // changing boundary conditions for gamma+ and gamma-
    mesh.GetNewReference(RefGammaMoins, RefGammaPlus);    
    
    if (exact_solution)
      {
        mesh.SetBoundaryCondition(RefGammaMoins, BoundaryConditionEnum::LINE_INSIDE) ;
        mesh.SetBoundaryCondition(RefGammaPlus, BoundaryConditionEnum::LINE_INSIDE) ;
      }
    else
      {
        mesh.SetBoundaryCondition(RefGammaMoins, BoundaryConditionEnum::LINE_TRANSMISSION) ;
        mesh.SetBoundaryCondition(RefGammaPlus, BoundaryConditionEnum::LINE_TRANSMISSION) ;
      }
    
    // creating the new edges
    for (int i = 0; i < nb_edges_ref; i++)
      {
        int nf = NumEdge(i);
        NumFaceMinus(i) = nf;
        NumFacePlus(i) = nb_edge2d;
        int n1 = IndexVert(mesh.BoundaryRef(nf).numVertex(0));
        int n2 = IndexVert(mesh.BoundaryRef(nf).numVertex(1));
        
        // On change la reference de l ancienne face
        mesh.BoundaryRef(nf).SetReference(RefGammaMoins) ;
        
        // on cree la nouvelle face sur Gamma Plus
        mesh.BoundaryRef(nb_edge2d).Init(n1, n2, RefGammaPlus);
        
        // points intermediaires
        for (int k = 0; k < r-1; k++)
          {
            R2 ptC = mesh.GetPointInsideEdge(nf, k);
            R2 ptM = ptC, ptP = ptC;
            Add(-Decalage*delta, NormaleEdge(i, k), ptM);
            Add(Decalage*delta, NormaleEdge(i, k), ptP);
            
            mesh.SetPointInsideEdge(nf, k, ptM);
            mesh.SetPointInsideEdge(nb_edge2d, k, ptP);	    
          }
        
        nb_edge2d++;
      }
    
    // we are changing vertex numbers of edges which are adjacent to Gamma
    for (int i = 0; i < nb_old_edge2d; i++)
      {
        bool edge_intersecting = false;
        bool edge_not_on_plane = false;
        int ref = mesh.BoundaryRef(i).GetReference();
        
        int n1 = mesh.BoundaryRef(i).numVertex(0);
        int n2 = mesh.BoundaryRef(i).numVertex(1);
        if (IndexVert(n1) >= nb_old_vertices)
          edge_intersecting = true;
        else
          edge_not_on_plane = true;
        
        if (IndexVert(n2) >= nb_old_vertices)
          edge_intersecting = true;
        else
          edge_not_on_plane = true;
        
        int ne = mesh.BoundaryRef(i).numElement(0);
        if ((edge_intersecting) && (edge_not_on_plane))
          if (mesh.Element(ne).GetReference() == RefOmegaPlus)
            mesh.BoundaryRef(i).Init(IndexVert(n1), IndexVert(n2), ref);
      }
    
    mesh.SetCurveType(RefGammaMoins, mesh.CURVE_FILE) ;
    mesh.SetCurveType(RefGammaPlus, mesh.CURVE_FILE) ;
      
    // modification of elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int nb_vert = mesh.Element(i).GetNbVertices();
        int ref = mesh.Element(i).GetReference() ;
        IVect num(nb_vert); bool elt_near_plate = false;
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = mesh.Element(i).numVertex(j);            
            if (IndexVert(num(j)) >= nb_old_vertices)
              elt_near_plate = true;
          }
        
        if (elt_near_plate)
          {
            if (ref == RefOmegaPlus) 
              {
                int ref = mesh.Element(i).GetReference();
                for (int j = 0; j < nb_vert; j++)
                  num(j) = IndexVert(num(j));
                
                mesh.Element(i).Init(num, ref);
              }
          }
      }
    
    // creating elements inside the layer if exact solution is required
    if (exact_solution)
      {
        int nb_elt = mesh.GetNbElt();
        mesh.ResizeElements(nb_elt+nb_edges_ref);
        for (int i = 0; i < nb_edges_ref; i++)
          {
            int nf = NumEdge(i);
            IVect num(4);
            num(0) = mesh.BoundaryRef(nf).numVertex(0);
            num(1) = mesh.BoundaryRef(nf).numVertex(1);
            num(2) = IndexVert(num(1));
            num(3) = IndexVert(num(0));
            
            mesh.Element(nb_elt+i).Init(num, ref_omega_layer);
          }
        
        mesh.ReorientElements();
	
	// for exact solution, we clear transmission_references
	// since no transmission conditions are applied
	// but the exact solution is computed with a thin layer
	transmission_references.Clear();
      }

    // updating faces/edges
    mesh.SortBoundariesRef();
    mesh.FindConnectivity();
    mesh.AddBoundaryEdges();
    mesh.FindConnectivity();
    mesh.ProjectPointsOnCurves();
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif

    if (nb_proc == 1)
      if (var_problem.print_level >= 6)
        mesh.Write("domain.mesh");
    
    if (var_problem.print_level >= 2)
      cout << "Mesh successfully parted " << endl;
    
    
    /*VectReal_wp step_subdiv(5); step_subdiv.Fill(); Mlt(0.25, step_subdiv);
    mesh.SubdivideMesh(step_subdiv);
    mesh.Write("domain_split.mesh");
    
    mesh.ClearElements();
    mesh.Write("faces.mesh");
    exit(0);*/
  }


  //! initialisation of transmission conditions  
  template<>
  void VarTransmission_Base<Dimension2>::TreatTransmission(const IVect& Epart)
  {
    if (exact_solution)
      transmission_references.Clear();

    if (transmission_references.GetM() <= 0)
      return;

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    // implementation similar to 3-D
    const Mesh<Dimension2>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    
    int nb_edges_ref = mesh_surf.GetNbBoundaryRef();
    NodleSurf.Reallocate(nb_edges_ref);
    int order = mesh_num.GetOrder();
    int nb_vertices = mesh_surf.GetNbVertices();
    int nb_edges = mesh_surf.GetNbBoundaryRef();
    int nb_dof_vertex = mesh_num.number_map.GetNbDofVertex(order);
    int nb_dof_edge = mesh_num.number_map.GetNbDofEdge(order);
    
    // retrieving dofs on vertices
    // here we construct an array OffsetDofVertexNumber
    // containing the offset for each vertex
    // the dof number of the j-th dof of vertex i is equal to OffsetDofVertexNumber(i)+j
    // the dof numbers are sorted such that an increasing sequence in surface dofs
    // correspond to an increasing sequence in volume dofs
    IVect OffsetDofVertexNumber(nb_vertices);    
    IVect offset_dof(nb_vertices);
    IVect permut(nb_vertices);
    int offset = 0, off_prec = 0; 
    if (nb_dof_vertex > 0)
      {
        for (int i = 0; i < nb_vertices; i++)
          offset_dof(i) = mesh_num.OffsetDofVertexNumber(OldNumVertexSurf(i));
        
        permut.Fill();
        Sort(offset_dof, permut);
        
        off_prec = offset_dof(0);
        offset_dof(0) = offset;
        for (int i = 1; i < nb_vertices; i++)
          {            
            if (offset_dof(i) != off_prec)
              offset += nb_dof_vertex;
                
            off_prec = offset_dof(i);
            offset_dof(i) = offset;
          }
        
        offset += nb_dof_vertex;
        
        for (int i = 0; i < nb_vertices; i++)
          OffsetDofVertexNumber(permut(i)) = offset_dof(i);
      }
    
    // retrieving dofs on edges
    // here we construct an array OffsetDofEdgeNumber
    // containing the offset for each edge
    // the dof number of the j-th dof of edge i is equal to OffsetDofEdgexNumber(i)+j
    // the dof numbers are sorted such that an increasing sequence in surface dofs
    // correspond to an increasing sequence in volume dofs
    IVect OffsetDofEdgeNumber(nb_edges);    
    IVect nb_dof(nb_edges);
    offset_dof.Reallocate(nb_edges);
    permut.Reallocate(nb_edges+1);
    for (int i = 0; i < nb_edges; i++)
      {
        offset_dof(i) = mesh_num.OffsetDofEdgeNumber(NumFaceMinus(i));
        int re = mesh_num.GetOrderEdge(NumFaceMinus(i));
        nb_dof(i) = mesh_num.number_map.GetNbDofEdge(re);
      }
    
    permut.Fill();
    Sort(offset_dof, nb_dof, permut);
    
    off_prec = offset_dof(0);
    offset_dof(0) = offset;
    for (int i = 1; i < nb_edges; i++)
      {                    
        if (offset_dof(i) != off_prec)
          offset += nb_dof(i);
        
        off_prec = offset_dof(i);
        offset_dof(i) = offset;
      }
    
    offset += nb_dof(nb_edges-1);
    
    for (int i = 0; i < nb_edges; i++)
      OffsetDofEdgeNumber(permut(i)) = offset_dof(i);
    
    int nodl_surf = offset;
    DdlVolMinus.Reallocate(nodl_surf);
    DdlVolPlus.Reallocate(nodl_surf);
    
    // filling NodleSurf
    for (int i = 0; i < nb_edges_ref; i++)
      {
        NodleSurf(i).Reallocate(2*nb_dof_vertex + nb_dof_edge);
        
        int n1 = mesh_surf.BoundaryRef(i).numVertex(0);
        int n2 = mesh_surf.BoundaryRef(i).numVertex(1);
        for (int k = 0; k < nb_dof_vertex; k++)
          NodleSurf(i)(k) = OffsetDofVertexNumber(n1) + k;
        
        for (int k = 0; k < nb_dof_edge; k++)
          NodleSurf(i)(nb_dof_vertex + k) = OffsetDofEdgeNumber(i) + k;
        
        for (int k = 0; k < nb_dof_vertex; k++)
          NodleSurf(i)(nb_dof_vertex + nb_dof_edge + k) = OffsetDofVertexNumber(n2) + k;
      }
    
    // first vertices
    for (int i = 0; i < OldNumVertexSurf.GetM(); i++)
      {
        int nv_minus = OldNumVertexSurf(i);
        int nv_plus = IndexVert(nv_minus);
        for (int k = 0; k < nb_dof_vertex; k++)
          {
            DdlVolMinus(OffsetDofVertexNumber(i) + k)
              = mesh_num.OffsetDofVertexNumber(nv_minus) + k;
            
            DdlVolPlus(OffsetDofVertexNumber(i) + k)
              = mesh_num.OffsetDofVertexNumber(nv_plus) + k;
          }
      }
    
    // then edges
    for (int i = 0; i < nb_edges_ref; i++)
      {
        int ne_minus = NumFaceMinus(i);
        int ne_plus = NumFacePlus(i);
        int offset_minus = mesh_num.OffsetDofEdgeNumber(ne_minus);
        int offset_plus = mesh_num.OffsetDofEdgeNumber(ne_plus);
        for (int k = 0; k < nb_dof_edge; k++)
          {
            DdlVolMinus(OffsetDofEdgeNumber(i) + k) = offset_minus + k;
            DdlVolPlus(OffsetDofEdgeNumber(i) + k) = offset_plus + k;
          }
      }
    
    nb_ddl_volume_on_surface = DdlVolMinus.GetM();
    mesh_surf.Clear();
    
    if (nb_proc > 1)
      {
        int nb_edges_ref = NumFaceMinus.GetM();
        NodleVolMinus.Reallocate(nb_edges_ref);
        NodleVolPlus.Reallocate(nb_edges_ref);
        ProcVolMinus.Reallocate(nb_edges_ref);
        ProcVolPlus.Reallocate(nb_edges_ref);
        for (int i = 0; i < nb_edges_ref; i++)
          {
            int nb_dof_edge = NodleSurf(i).GetM();
            NodleVolMinus(i).Reallocate(nb_dof_edge);
            NodleVolPlus(i).Reallocate(nb_dof_edge);
            int ne_minus = mesh.BoundaryRef(NumFaceMinus(i)).numElement(0);
            int ne_plus = mesh.BoundaryRef(NumFacePlus(i)).numElement(0);
            ProcVolMinus(i) = Epart(ne_minus);
            ProcVolPlus(i) = Epart(ne_plus);
            for (int j = 0; j < nb_dof_edge; j++)
              {
                int num_dof = NodleSurf(i)(j);
                NodleVolMinus(i)(j) = DdlVolMinus(num_dof);
                NodleVolPlus(i)(j) = DdlVolPlus(num_dof);
	      }
	  }

      }    
  }
  

#ifdef MONTJOIE_WITH_THREE_DIM  
  //! separation of mesh into two parts  
  template<>
  void VarTransmission_Base<Dimension3>::PartMeshTransmission()
  {
    if (transmission_references.GetM() <= 0)
      return;

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    Mesh<Dimension3>& mesh = var_problem.mesh;
    
    bool store_normale_gamma = true;
    
    int ref_plate = 1;
    IVect ref_cond(mesh.GetNbReferences()+1); ref_cond.Fill(0);
    for (int i = 0; i < transmission_references.GetM(); i++)
      ref_cond(transmission_references(i)) = ref_plate;
        
    if (gamma_is_plate)
      {
        // reorienting faces so that 2-D jacobians will be positive without need of 
        // reorienting 2-D element. Otherwise, we would have to handle different orientations
        // between 2-D element and 3-D face
        mesh.ReorientFaces(ref_plate, ref_cond, R3(0, 0, 1));
        mesh.FindConnectivity();
      }
    
    // we are extracting the surface mesh
    IVect NumFace, NumVertices, NumElement, NumLocalFace;
    // creation du maillage mesh_surf
    // NumFace donne le numero globale de la face i de mesh_surf dans le maillage volumique
    // idem pour NumVertices
    // NumElement : element volumique auquel la face appartient (un element sur deux au hasard ))
    // NumLocalFace numero local de la face dans l element 
    mesh.GetBoundaryMesh(ref_plate, mesh_surf, NumFace, NumVertices, 
                         NumElement, NumLocalFace, ref_cond);
  
    OldNumFaceSurf = NumFace;
    OldNumVertexSurf = NumVertices;
    mesh_surf.FindConnectivity();
    //mesh_surf.Write("surf.mesh");
    
    if (gamma_is_plate)
      {
        // changing order of approximation if necessary
        int order = mesh_plate.GetGeometryOrder();
        mesh_surf.SetGeometryOrder(order);

        // copying vertices
        mesh_plate.Clear();
        mesh_plate.ReallocateVertices(mesh_surf.GetNbVertices());
        for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
          mesh_plate.Vertex(i).Init(mesh_surf.Vertex(i)(0), mesh_surf.Vertex(i)(1));
        
        // copying edges on the boundary
        int nb_edges_ref = 0;
        for (int i = 0; i < mesh_surf.GetNbEdgesRef(); i++)
          if (mesh_surf.EdgeRef(i).GetNbFaces() == 1)
            nb_edges_ref++;
        
        mesh_plate.ReallocateBoundariesRef(nb_edges_ref);
        nb_edges_ref = 0;
        R2 ptA; R3 pt3d;
        for (int i = 0; i < mesh_surf.GetNbEdgesRef(); i++)
          if (mesh_surf.EdgeRef(i).GetNbFaces() == 1)
            {
              mesh_plate.BoundaryRef(nb_edges_ref).
                Init(mesh_surf.EdgeRef(i).numVertex(0),
                     mesh_surf.EdgeRef(i).numVertex(1), ref_plate);
              
              for (int k = 0; k < order-1; k++)
                {
                  pt3d = mesh_surf.GetPointInsideEdge(i, k);
                  ptA(0) = pt3d(0); ptA(1) = pt3d(1);
                  mesh_plate.SetPointInsideEdge(nb_edges_ref, k, ptA);
                }
              
              nb_edges_ref++;
            }               
      }
    
    Vector<bool> FaceToDuplicate(mesh.GetNbBoundary());
    FaceToDuplicate.Fill(false);
    for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
      FaceToDuplicate(NumFace(i)) = true;
    
    // detection of vertices to duplicate
    Vector<bool> PointToDuplicate(mesh_surf.GetNbVertices());
    PointToDuplicate.Fill(true);
    //DISP(duplicate_extremity);
    if (!duplicate_extremity)
      {
        for (int i = 0; i < mesh_surf.GetNbEdges(); i++)
          if (mesh_surf.GetEdge(i).GetNbFaces() == 1)
            {
              int n1 = mesh_surf.GetEdge(i).numVertex(0);
              PointToDuplicate(n1) = false;
              int n2 = mesh_surf.GetEdge(i).numVertex(1);
              PointToDuplicate(n2) = false;
            }
      }
    
    IndexVert.Reallocate(mesh.GetNbVertices()); IndexVert.Fill();
    int nb_old_vertices = mesh.GetNbVertices();
    int nb_vertices = nb_old_vertices;
    for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
      if (PointToDuplicate(i))
        {
          IndexVert(NumVertices(i)) = nb_vertices;
          nb_vertices++;
        }

    // references of domains must be different
    if (!gamma_is_plate)
      {
        if (RefOmegaMoins == RefOmegaPlus)
          {
            cout << "Use different references for Omega- and Omega+ " << endl;
            abort();
          }
      }
    
    int r = mesh_surf.GetGeometryOrder();
    int nb_vert_surf = mesh_surf.GetNbVertices();
    int nb_face_surf = mesh_surf.GetNbFaces();
    // allocating arrays that will contain normale associated with vertices,
    // nodal points of edges, and nodal points of faces
    VectR3 Normales(nb_vert_surf) ;
    Matrix<R3> NormaleEdge; Vector<VectR3> NormaleFace;
    if (r > 1)
      {
        NormaleEdge.Reallocate(mesh.GetNbEdgesRef(), r-1);
        NormaleFace.Reallocate(nb_face_surf);
      }
    
    // computation of normales for each vertex, edge and face
    // normales are outward to domain Omega-
    VectR3 S ; SetPoints<Dimension3> Points ; SetMatrices<Dimension3> Matrices ;
    Matrix3_3 mat_dfj, dfjm1; R3 normale_fj, vec_unit; Real_wp dsj(0);
    Matrix<int> RotationTri, RotationQuad;
    mesh.GetGeometryFaceRotation(RotationTri, RotationQuad);
    if (store_normale_gamma)
      {
        NormaleNodes.Reallocate(nb_face_surf);
        DsNodes.Reallocate(nb_face_surf);
        PosNodes.Reallocate(nb_face_surf);
      }
    
    MatDfjNodes.Reallocate(mesh.GetNbElt());    
    for (int i = 0; i < nb_face_surf; i++)
      {
        int Ne = NumElement(i) ;
        int NumLoc =  NumLocalFace(i) ;
        mesh.GetVerticesElement(Ne,S) ;
        mesh.FjElemNodal(S, Points,mesh, Ne) ;
        mesh.DFjElemNodal(S, Points,Matrices, mesh, Ne) ;
        int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
        int nb_nodes = nb_vert, nb_inside = 0;
        if (nb_vert == 3)
          {
            nb_nodes = (r+1)*(r+2)/2;
            nb_inside = (r-1)*(r-2)/2;
          }
        else
          {
            nb_nodes = (r+1)*(r+1);
            nb_inside = (r-1)*(r-1);
          }

        if (nb_inside > 0)
          NormaleFace(i).Reallocate(nb_inside);

        if (store_normale_gamma)
          {
            NormaleNodes(i).Reallocate(nb_nodes);
            PosNodes(i).Reallocate(nb_nodes);
            DsNodes(i).Reallocate(nb_nodes);        
          }
        
        int rot = mesh.Element(Ne).GetOrientationFace(NumLoc);
        for (int j = 0; j < nb_nodes; j++)
          {
            int jrot = j;
            if (nb_vert == 3)
              jrot = RotationTri(rot, j);
            else
              jrot = RotationQuad(rot, j);
            
            int node = mesh.GetNodalNumber(Ne, NumLoc, jrot);
            mat_dfj = Matrices.GetPointNodal(node);
            GetInverse(mat_dfj, dfjm1);
            mesh.GetNormale(dfjm1, normale_fj, dsj, Ne, NumLoc);
            
            // normale pour une sphere
            if (gamma_is_sphere)
              {
                vec_unit = Points.GetPointNodal(node);
                Mlt(1.0/Norm2(vec_unit), vec_unit);
                if (DotProd(vec_unit, normale_fj) < 0)
                  Mlt(-1.0, vec_unit);
                
                normale_fj = vec_unit;
              }

            if (RefOmegaMoins > RefOmegaPlus)
              Mlt(-1.0, normale_fj);
            
            if (j < nb_vert)
              Normales(mesh_surf.BoundaryRef(i).numVertex(j)) = normale_fj;
            else if (j < nb_vert*r)
              {
                int na = (j-nb_vert)/(r-1);
                int k = (j-nb_vert)%(r-1);
                int nf = NumFace(i);
                int num_edge = mesh.BoundaryRef(nf).numEdge(na);
                if (num_edge < mesh.GetNbEdgesRef())
                  {                
                    if (mesh.BoundaryRef(nf).GetOrientationEdge(na))         
                      NormaleEdge(num_edge, k) = normale_fj;
                    else
                      NormaleEdge(num_edge, r-2-k) = normale_fj;
                  }
              }
            else
              {
                NormaleFace(i)(j-nb_vert*r) = normale_fj;
              }
            
            if (store_normale_gamma)
              {
                NormaleNodes(i)(j) = normale_fj;
                DsNodes(i)(j) = dsj;
                PosNodes(i)(j) = Points.GetPointNodal(node);
              }
          }

	// storing matrix dfj on the boundary
	int nb_faces_loc = mesh.Element(Ne).GetNbBoundary();
	MatDfjNodes(Ne).Reallocate(nb_faces_loc);
	MatDfjNodes(Ne)(NumLoc).Reallocate(nb_nodes);
	for (int j = 0; j < nb_nodes; j++)
	  {
	    int node = mesh.GetNodalNumber(Ne, NumLoc, j);
	    MatDfjNodes(Ne)(NumLoc)(j) = Matrices.GetPointNodal(node);
	  }	
	
	// and on the opposite side too
	int num_face = mesh.Element(Ne).numBoundary(NumLoc);
	int Ne2 = mesh.Boundary(num_face).numElement(0);
	if (Ne2 == Ne)
	  Ne2 = mesh.Boundary(num_face).numElement(1);
	
	int NumLoc2 = mesh.Element(Ne2).GetPositionBoundary(num_face);
	nb_faces_loc = mesh.Element(Ne2).GetNbBoundary();
        mesh.GetVerticesElement(Ne2, S) ;
        mesh.FjElemNodal(S, Points, mesh, Ne2) ;
        mesh.DFjElemNodal(S, Points, Matrices, mesh, Ne2) ;
	
	MatDfjNodes(Ne2).Reallocate(nb_faces_loc);
	MatDfjNodes(Ne2)(NumLoc2).Reallocate(nb_nodes);
	for (int j = 0; j < nb_nodes; j++)
	  {
	    int node = mesh.GetNodalNumber(Ne2, NumLoc2, j);
	    MatDfjNodes(Ne2)(NumLoc2)(j) = Matrices.GetPointNodal(node);
	  }
        
      }
    
    // adding vertices to the 3-D mesh
    mesh.ResizeVertices(nb_vertices);
    nb_vertices = nb_old_vertices;
    for (int i = 0; i < mesh_surf.GetNbVertices(); i++)
      if (PointToDuplicate(i))
        {
          // vertices are shifted with respect to normale
          mesh.Vertex(nb_vertices) = mesh_surf.Vertex(i);
          mesh.Vertex(NumVertices(i)) -= Decalage*delta*Normales(i); 
          mesh.Vertex(nb_vertices) += Decalage*delta*Normales(i); 
          nb_vertices++;
        }
        
    // copying faces
    int nb_faces_ref = mesh_surf.GetNbBoundaryRef();
    // and duplicating faces for volume mesh
    int nb_old_face3d = mesh.GetNbBoundaryRef(); 
    int nb_face3d = nb_old_face3d;
    mesh.ResizeBoundariesRef(nb_old_face3d + nb_faces_ref);
    NumFaceMinus.Reallocate(nb_faces_ref);
    NumFacePlus.Reallocate(nb_faces_ref);
    
    // changing boundary conditions for gamma+ and gamma-
    mesh.GetNewReference(RefGammaMoins, RefGammaPlus);    
    if (gamma_is_plate)
      {
        //RefGammaMoins = 1;
        //RefGammaPlus = 1;
        mesh.SetBoundaryCondition(RefGammaMoins, BoundaryConditionEnum::LINE_NEUMANN) ;
        mesh.SetBoundaryCondition(RefGammaPlus, BoundaryConditionEnum::LINE_NEUMANN) ;
      }
    else
      {
        mesh.SetBoundaryCondition(RefGammaMoins, BoundaryConditionEnum::LINE_TRANSMISSION) ;
        
        if (exact_solution)
          {
            mesh.SetBoundaryCondition(RefGammaMoins, BoundaryConditionEnum::LINE_INSIDE) ;
            mesh.SetBoundaryCondition(RefGammaPlus, BoundaryConditionEnum::LINE_INSIDE) ;
          }
        else
          mesh.SetBoundaryCondition(RefGammaPlus, BoundaryConditionEnum::LINE_TRANSMISSION) ;
      }
    
    // creating the new faces
    int ref_among_other = -1;
    //DISP(nb_faces_ref);
    if (gamma_is_plate)
      mesh_plate.ReallocateElements(nb_faces_ref);
    
    for (int i = 0; i < nb_faces_ref; i++)
      {
        int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
        int ref = mesh_surf.BoundaryRef(i).GetReference();
        ref_among_other = ref;
        IVect num(nb_vert);
        IVect num_surf(nb_vert);
        
        // and 3-D face
        int nf = NumFace(i);
        NumFaceMinus(i) = nf;
        NumFacePlus(i) = nb_face3d;
        for (int j = 0; j < nb_vert; j++)
          {
            num_surf(j) = mesh_surf.BoundaryRef(i).numVertex(j);
            num(j) = IndexVert(mesh.BoundaryRef(nf).numVertex(j));
          }
        
        if (gamma_is_plate)
          {
            // initialization of 2-D element
            mesh_plate.Element(i).Init(num_surf, ref);
            
            // On change la reference de l ancienne face
            mesh.BoundaryRef(nf).SetReference(RefGammaMoins) ;
            
            // on cree la nouvelle face sur Gamma Plus
            mesh.BoundaryRef(nb_face3d).Init(num, RefGammaPlus);
          }
        else
          {
            mesh.BoundaryRef(nf).SetReference(RefGammaMoins) ;
            mesh.BoundaryRef(nb_face3d).Init(num, RefGammaPlus);
          }
        
        nb_face3d++;
      }

    // we find the equation of the plate z = cte
    Vector<Real_wp> coef_plane;    
    Real_wp z_cte = 0;
    if (gamma_is_plate)
      {
        mesh_surf.FindParametersPlane(ref_among_other, coef_plane);
        if ((abs(coef_plane(0)) > epsilon_machine) || (abs(coef_plane(1)) > epsilon_machine))
          {
            cout << "We handle only plates located on plane z = cte "<< endl;
            cout << "the equation of your plate is equal to "<< endl;
            cout << coef_plane(0) << " x + " << coef_plane(1) << " y + " << coef_plane(2)
                 << " z + " << coef_plane(3) << " = 0 " << endl;
            abort();
          }
        else
          {
            z_cte = -coef_plane(3)/coef_plane(2);
          }
      }
    
    // creating new edges
    int nb_edges_ref = mesh.GetNbEdgesRef();
    int nb_old_edges_ref = nb_edges_ref;
    for (int i = 0; i < nb_old_edges_ref; i++)
      {
        bool edge_on_face = true;
        int n1 = mesh.EdgeRef(i).numVertex(0);
        int n2 = mesh.EdgeRef(i).numVertex(1);
        if (IndexVert(n1) < nb_old_vertices)
          edge_on_face = false;
        
        if (IndexVert(n2) < nb_old_vertices)
          edge_on_face = false;
        
        if (edge_on_face)
          nb_edges_ref++;
        
      }
    
    if (nb_edges_ref > nb_old_edges_ref)
      {
        mesh.ResizeEdgesRef(nb_edges_ref);
        int nb_edges_ref = nb_old_edges_ref;
        for (int i = 0; i < nb_old_edges_ref; i++)
          {
            bool edge_on_face = true;
            bool edge_intersecting = false;
            int n1 = mesh.EdgeRef(i).numVertex(0);
            int n2 = mesh.EdgeRef(i).numVertex(1);
            int ref = mesh.EdgeRef(i).GetReference();
            if (IndexVert(n1) >= nb_old_vertices)
              edge_intersecting = true;
            else
              edge_on_face = false;
            
            if (IndexVert(n2) >= nb_old_vertices)
              edge_intersecting = true;
            else
              edge_on_face = false;
            
            if (edge_on_face)
              {
                if (gamma_is_plate)
                  {
                    mesh.EdgeRef(nb_edges_ref).Init(IndexVert(n1), IndexVert(n2), ref);
                    mesh.EdgeRef(i).SetReference(ref);
                  }
                else
                  {
                    mesh.EdgeRef(nb_edges_ref).Init(IndexVert(n1), IndexVert(n2), RefGammaPlus);
                    mesh.EdgeRef(i).SetReference(RefGammaMoins);
                  }
                
                for (int k = 0; k < r-1; k++)
                  {
                    R3 ptC = mesh.GetPointInsideEdge(i, k);
                    R3 ptM = ptC, ptP = ptC;
                    Add(-Decalage*delta, NormaleEdge(i, k), ptM);
                    Add(Decalage*delta, NormaleEdge(i, k), ptP);
                    
                    mesh.SetPointInsideEdge(i, k, ptM);
                    mesh.SetPointInsideEdge(nb_edges_ref, k, ptP);
                  }
                nb_edges_ref++;
              }
            else if (edge_intersecting)
              {
                int num_elem = mesh.EdgeRef(i).numElement(0);
                if (mesh.Element(num_elem).GetReference() == RefOmegaPlus)
                  mesh.EdgeRef(i).Init(IndexVert(n1), IndexVert(n2), 0);
                else
                  mesh.EdgeRef(i).SetReference(0);
                
                mesh.SetStraightEdge(i);
                
              }
            
          }
        
      }
    
    // shifting internal nodes (for curved meshes)
    for (int i = 0; i < nb_faces_ref; i++)
      {
        int nm = NumFaceMinus(i);
        int np = NumFacePlus(i);
        
        int ref = mesh_surf.BoundaryRef(i).GetReference();
        if (mesh.GetCurveType(ref) > 0)
          {
            int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
            int nb_inside = 0;
            if (nb_vert == 3)
              nb_inside = (r-1)*(r-2)/2;
            else
              nb_inside = (r-1)*(r-1);
            
            if (nb_inside > 0)
              {
                VectR3 points_m(nb_inside), points_p(nb_inside);
                for (int j = 0; j < nb_inside; j++)
                  {
                    points_p(j) = mesh.GetPointInsideFace(nm, j);
                    points_m(j) = points_p(j);
                    Add(-Decalage*delta, NormaleFace(i)(j), points_m(j));
                    Add(Decalage*delta, NormaleFace(i)(j), points_p(j));
                  }
                
                mesh.SetPointInsideFace(nm, points_m);
                mesh.SetPointInsideFace(np, points_p);
              }
          }
        else
          {
            mesh.SetStraightFace(nm);
            mesh.SetStraightFace(np);
          }
        
      }
    
    // we are changing vertex numbers of faces which are adjacent to Gamma
    R3 center;
    for (int i = 0; i < nb_old_face3d; i++)
      {
        bool face_intersecting = false;
        bool face_not_on_plane = false;
        int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
        int ref = mesh.BoundaryRef(i).GetReference();
        IVect num(nb_vert);
        center.Zero();
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = mesh.BoundaryRef(i).numVertex(j);            
            center += mesh.Vertex(num(j));
            if (IndexVert(num(j)) >= nb_old_vertices)
              face_intersecting = true;
            else
              face_not_on_plane = true;

          }

        Mlt(1.0/Real_wp(nb_vert), center);
        
        if (FaceToDuplicate(i))
          face_not_on_plane = false;
        
        int ne = mesh.BoundaryRef(i).numElement(0);
        if ((face_intersecting) && (face_not_on_plane))
          if ( (gamma_is_plate && (center(2) > z_cte))
               || ( (!gamma_is_plate) && (mesh.Element(ne).GetReference() == RefOmegaPlus) ) )
            {
              for (int j = 0; j < nb_vert; j++)
                num(j) = IndexVert(num(j));
              
              mesh.BoundaryRef(i).Init(num, ref);
              mesh.SetStraightFace(i);
            }
      }
    
    if (!gamma_is_plate)
      {
        mesh.SetCurveType(RefGammaMoins, mesh.CURVE_FILE) ;
        mesh.SetCurveType(RefGammaPlus, mesh.CURVE_FILE) ;
      }
    
    // modification of elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int nb_vert = mesh.Element(i).GetNbVertices();
        int ref = mesh.Element(i).GetReference() ;
        IVect num(nb_vert); bool elt_near_plate = false;
        center.Zero();
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = mesh.Element(i).numVertex(j);            
            center += mesh.Vertex(num(j));
            if (IndexVert(num(j)) >= nb_old_vertices)
              elt_near_plate = true;
          }
        
        Mlt(1.0/Real_wp(nb_vert), center);
        
        if (elt_near_plate)
          {
            if ((gamma_is_plate && (center(2) > z_cte))
                || ((!gamma_is_plate) && (ref == RefOmegaPlus)) )
              {
                int ref = mesh.Element(i).GetReference();
                for (int j = 0; j < nb_vert; j++)
                  num(j) = IndexVert(num(j));
                
                mesh.Element(i).Init(num, ref);
              }
          }
      }
    
    // creating elements inside the layer if exact solution is required
    if (exact_solution)
      {
        int nb_elt = mesh.GetNbElt();
        mesh.ResizeElements(nb_elt+nb_faces_ref);
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
            int nf = NumFace(i);
            IVect num(2*nb_vert);
            for (int j = 0; j < nb_vert; j++)
              {
                num(j) = mesh.BoundaryRef(nf).numVertex(j);
                num(nb_vert+j) = IndexVert(num(j));
              }

            mesh.Element(nb_elt+i).Init(num, ref_omega_layer);
          }
        
        mesh.ReorientElements();

	// for exact solution, we clear transmission_references
	// since no transmission conditions are applied
	// but the exact solution is computed with a thin layer
	transmission_references.Clear();
      }
    
    if (!gamma_is_plate)
      {
        mesh_surf.SortBoundariesRef();
        mesh.SortBoundariesRef();
        mesh_surf.FindConnectivity();        
      }
    
    // updating faces/edges
    mesh.FindConnectivity();
    mesh.AddBoundaryFaces();
    mesh.ProjectPointsOnCurves();

    if (gamma_is_plate)
      {
        mesh_plate.ReorientElements();
        mesh_plate.FindConnectivity();
        //mesh_plate.Write("plaque.mesh");
      }
    
    if (nb_proc == 1)
      if (var_problem.print_level >= 6)
        mesh.Write("domain.mesh");
    
    /*{
      VectReal_wp step_subdiv(5); step_subdiv.Fill(); Mlt(0.25, step_subdiv);
      mesh.SubdivideMesh(step_subdiv);
      mesh.Write("domain_split.mesh");
      
      mesh.ClearElements();
      mesh.Write("faces.mesh");
      exit(0);
    }*/
  }


  //! initialisation of transmission conditions  
  template<>
  void VarTransmission_Base<Dimension3>::TreatTransmission(const IVect& Epart)
  {
    if (exact_solution)
      transmission_references.Clear();
    
    if (transmission_references.GetM() <= 0)
      return;
    
    Mesh<Dimension3>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension3>& mesh_num = var_problem.GetMeshNumbering(0);

#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif
        
    if (gamma_is_plate)
      {
        int nb_faces = NumFaceMinus.GetM();
        // order for acoustic
        int ra = mesh_num.GetOrder();
        // order for mecanic
        //int rm = mesh_plate.GetOrder();
        NodleAcous.Reallocate(nb_faces);
        NodleMeca.Reallocate(nb_faces);
        int nb_vertices_vol = mesh.GetNbVertices();
        //int nb_edges_vol = mesh.GetNbEdges();
        
        int nb_vertices = mesh_plate.GetNbVertices();
        int nb_edges = mesh_plate.GetNbEdges();
        int nodl_acous = nb_vertices + nb_edges*(ra-1);
        JacobianAcous.Reallocate(nb_faces);
        // computing dofs for acoustic
        for (int i = 0; i < nb_faces; i++)
          {
            int nb_vert = mesh_plate.Element(i).GetNbVertices();
            int nb_dof  = 0;
            if (nb_vert == 3)
              nb_dof = (ra+1)*(ra+2)/2;
            else
              nb_dof = (ra+1)*(ra+1);
            
            NodleAcous(i).Reallocate(nb_dof);
            nodl_acous += nb_dof - nb_vert*ra;
          }
        
        DdlVolMinus.Reallocate(nodl_acous); DdlVolMinus.Fill(-1);
        DdlVolPlus.Reallocate(nodl_acous); DdlVolPlus.Fill(-1);
        nodl_acous = nb_vertices + nb_edges*(ra-1);
        // int nodl_vol = nb_vertices_vol + (ra-1)*nb_edges_vol;
        SetPoints<Dimension2> pts_surf; VectR2 s2; R2 pt2d;
        Matrix2_2 mat_dfj; SetMatrices<Dimension2> mat_surf;
        for (int i = 0; i < nb_faces; i++)
          {        
            if (nb_proc == 1)
              NodleMeca(i) = mesh_plate_num.Element(i).GetNodle();
            
            int nb_vert = mesh_plate.Element(i).GetNbVertices();
            //int ref = mesh_plate.Element(i).GetReference();
            int nb_dof = NodleAcous(i).GetM();
            // vertices
            for (int j = 0; j < nb_vert; j++)
              {
                NodleAcous(i)(j) = mesh_plate.Element(i).numVertex(j);
                DdlVolMinus(NodleAcous(i)(j)) = mesh.BoundaryRef(NumFaceMinus(i)).numVertex(j);
                DdlVolPlus(NodleAcous(i)(j)) = mesh.BoundaryRef(NumFacePlus(i)).numVertex(j);
              }
            
            // edges
            for (int j = 0; j < nb_vert; j++)
              {
                int ne = mesh_plate.Element(i).numEdge(j);
                int nem = mesh.BoundaryRef(NumFaceMinus(i)).numEdge(j);
                int nep = mesh.BoundaryRef(NumFacePlus(i)).numEdge(j);
                bool way_surf = mesh_plate.Element(i).GetOrientationEdge(j);
                int off_acous = nb_vert + (ra-1)*j;
                if (way_surf)
                  for (int k = 0; k < ra-1; k++)
                    NodleAcous(i)(off_acous + k) = nb_vertices + ne*(ra-1) + k;
                else
                  for (int k = 0; k < ra-1; k++)
                    NodleAcous(i)(off_acous + k) = nb_vertices + ne*(ra-1) + ra-2-k;
                
                for (int k = 0; k < ra-1; k++)
                  DdlVolMinus(nb_vertices + ne*(ra-1) + k) = nb_vertices_vol + nem*(ra-1) + k;
                
                for (int k = 0; k < ra-1; k++)
                  DdlVolPlus(nb_vertices + ne*(ra-1) + k) = nb_vertices_vol + nep*(ra-1) + k;
              }
            
            // then inside of faces
            for (int j = nb_vert*ra; j < nb_dof; j++)
              {
                int j2 = j - ra*nb_vert;
                DdlVolMinus(nodl_acous) = mesh_num.OffsetDofFaceNumber(NumFaceMinus(i)) + j2;
                DdlVolPlus(nodl_acous) = mesh_num.OffsetDofFaceNumber(NumFacePlus(i)) + j2;
                NodleAcous(i)(j) = nodl_acous++;
              }
            
            // computation of jacobian
            int Nquad = (ra+1)*(ra+1);
            JacobianAcous(i).Reallocate(Nquad);
            mesh_plate.GetVerticesElement(i, s2);
            
            const VectReal_wp& FluxWeight = mesh_num.number_map
	      .GetFluxWeight(ra, mesh.BoundaryRef(NumFacePlus(i)));
            
	    mesh_plate.FjElemNodal(s2, pts_surf, mesh_plate, i);
            for (int j = 0; j < Nquad; j++)
              {
                R2 pt_quad = mesh_num.number_map
		  .GetQuadraturePoint(ra, j, mesh.BoundaryRef(NumFacePlus(i)));
                
		mesh_plate.DFj(s2, pts_surf, pt_quad,
                               mat_dfj, mesh_plate, i);
                
                Real_wp poids = 2.0*FluxWeight(j);
                JacobianAcous(i)(j) = Det(mat_dfj)*poids;
              }
          }
        
        //DISP(DdlVolMinus); DISP(DdlVolPlus);
        
        // retrieving dofs with Dirichlet condition
        DirichletAcous.Reallocate((ra+1)*mesh_plate.GetNbBoundaryRef());
        int nb_dir = 0;
        for (int i = 0; i < mesh_plate.GetNbBoundaryRef(); i++)
          {
            int ne = i;
            DirichletAcous(nb_dir++) = mesh_plate.BoundaryRef(i).numVertex(0);
            DirichletAcous(nb_dir++) = mesh_plate.BoundaryRef(i).numVertex(1);
            for (int k = 0; k < ra - 1; k++)
              DirichletAcous(nb_dir++) = nb_vertices + ne*(ra-1) + k;
          }
        
        Assemble(nb_dir, DirichletAcous); DirichletAcous.Resize(nb_dir);
        
        NodleSurf = NodleAcous;
      }
    else
      {        
        int nb_faces_ref = mesh_surf.GetNbBoundaryRef();
        NodleSurf.Reallocate(nb_faces_ref);
        int order = mesh_num.GetOrder();
        int nb_vertices = mesh_surf.GetNbVertices();
        int nb_edges = mesh_surf.GetNbEdges();
        int nb_dof_vertex = mesh_num.number_map.GetNbDofVertex(order);
        int nb_dof_edge = mesh_num.number_map.GetNbDofEdge(order);
        
        // we fill NumEdgeMinus and NumEdgePlus (edge numbers across the interface)
        IVect NumEdgeMinus(nb_edges), NumEdgePlus(nb_edges);
        NumEdgeMinus.Fill(-1); NumEdgePlus.Fill(-1);
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int nm = NumFaceMinus(i);
            int np = NumFacePlus(i);
            int nb_vert = mesh.BoundaryRef(nm).GetNbVertices();
            for (int j = 0; j < nb_vert; j++)
              {
                int ne = mesh_surf.BoundaryRef(i).numEdge(j);
                int ne_minus = mesh.BoundaryRef(nm).numEdge(j);
                int ne_plus = mesh.BoundaryRef(np).numEdge(j);
                NumEdgeMinus(ne) = ne_minus;
                NumEdgePlus(ne) = ne_plus;                
              }
            //DISP(i); DISP(NumEdgeMinus); DISP(NumEdgePlus);
          }
        //DISP(NumEdgeMinus); DISP(NumEdgePlus);
        
        // retrieving dofs on vertices
        IVect OffsetDofVertexNumber(nb_vertices);    
        IVect offset_dof(nb_vertices);
        IVect permut(nb_vertices);
        int offset = 0, off_prec = 0; 
        if (nb_dof_vertex > 0)
          {
            for (int i = 0; i < nb_vertices; i++)
              offset_dof(i) = mesh_num.OffsetDofVertexNumber(OldNumVertexSurf(i));
            
            permut.Fill();
            Sort(offset_dof, permut);
            
            off_prec = offset_dof(0);
            offset_dof(0) = offset;
            for (int i = 1; i < nb_vertices; i++)
              {            
                if (offset_dof(i) != off_prec)
                  offset += nb_dof_vertex;
                
                off_prec = offset_dof(i);
                offset_dof(i) = offset;
              }
            
            offset += nb_dof_vertex;
            
            for (int i = 0; i < nb_vertices; i++)
              OffsetDofVertexNumber(permut(i)) = offset_dof(i);
          }
        
        // dofs on edges
        IVect OffsetDofEdgeNumber(nb_edges);    
        IVect nb_dof(nb_edges);
        offset_dof.Reallocate(nb_edges);
        permut.Reallocate(nb_edges+1);
        for (int i = 0; i < nb_edges; i++)
          {
            offset_dof(i) = mesh_num.OffsetDofEdgeNumber(NumEdgeMinus(i));
            int re = mesh_num.GetOrderEdge(NumEdgeMinus(i));
            nb_dof(i) = mesh_num.number_map.GetNbDofEdge(re);
          }
        
        permut.Fill();
        Sort(offset_dof, nb_dof, permut);
        
        off_prec = offset_dof(0);
        offset_dof(0) = offset;
        for (int i = 1; i < nb_edges; i++)
          {                    
            if (offset_dof(i) != off_prec)
              offset += nb_dof(i);
            
            off_prec = offset_dof(i);
            offset_dof(i) = offset;
          }
        
        offset += nb_dof(nb_edges-1);
        
        for (int i = 0; i < nb_edges; i++)
          OffsetDofEdgeNumber(permut(i)) = offset_dof(i);
        
        // dofs on faces
        IVect OffsetDofFaceNumber(nb_faces_ref);    
        nb_dof.Reallocate(nb_faces_ref);
        offset_dof.Reallocate(nb_faces_ref);
        permut.Reallocate(nb_faces_ref);
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int nf = OldNumFaceSurf(i);
            offset_dof(i) = mesh_num.OffsetDofFaceNumber(nf);
            int rf = mesh_num.GetOrderFace(nf);
            nb_dof(i) = mesh_num.number_map.GetNbDofElement(rf, mesh.BoundaryRef(nf));
          }
        
        permut.Fill();
        Sort(offset_dof, nb_dof, permut);
        
        off_prec = offset_dof(0);
        offset_dof(0) = offset;
        for (int i = 1; i < nb_faces_ref; i++)
          {                    
            if (offset_dof(i) != off_prec)
              offset += nb_dof(i);
            
            off_prec = offset_dof(i);
            offset_dof(i) = offset;
          }
        
        offset += nb_dof(nb_faces_ref-1);
        
        for (int i = 0; i < nb_faces_ref; i++)
          OffsetDofFaceNumber(permut(i)) = offset_dof(i);
        
        int nodl_surf = offset;
        DdlVolMinus.Reallocate(nodl_surf);
        DdlVolPlus.Reallocate(nodl_surf);
        // filling NodleSurf
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int nb_vert = mesh_surf.BoundaryRef(i).GetNbVertices();
            int nb_dof_inside = mesh_num.number_map
	      .GetNbDofElement(order, mesh.BoundaryRef(OldNumFaceSurf(i)));
            
            int offset_loc = (nb_dof_vertex + nb_dof_edge)*nb_vert;
            NodleSurf(i).Reallocate(offset_loc + nb_dof_inside);
            for (int j = 0; j < nb_vert; j++)
              {
                int nv = mesh_surf.BoundaryRef(i).numVertex(j);
                for (int k = 0; k < nb_dof_vertex; k++)
                  NodleSurf(i)(j*nb_dof_vertex+k) = OffsetDofVertexNumber(nv) + k;
                
                int ne = mesh_surf.BoundaryRef(i).numEdge(j);
                if (mesh_surf.BoundaryRef(i).GetOrientationEdge(j))
                  for (int k = 0; k < nb_dof_edge; k++)
                    NodleSurf(i)(nb_dof_vertex*nb_vert + j*nb_dof_edge + k)
                      = OffsetDofEdgeNumber(ne) + k;
                else
                  for (int k = 0; k < nb_dof_edge; k++)
                    NodleSurf(i)(nb_dof_vertex*nb_vert + j*nb_dof_edge + k)
                      = OffsetDofEdgeNumber(ne) + nb_dof_edge-1-k;
              }
            
            for (int k = 0; k < nb_dof_inside; k++)
              NodleSurf(i)(offset_loc + k) = OffsetDofFaceNumber(i) + k;
            
          }
        
        // first vertices
        for (int i = 0; i < OldNumVertexSurf.GetM(); i++)
          {
            int nv_minus = OldNumVertexSurf(i);
            int nv_plus = IndexVert(nv_minus);
            for (int k = 0; k < nb_dof_vertex; k++)
              {
                DdlVolMinus(OffsetDofVertexNumber(i) + k)
                  = mesh_num.OffsetDofVertexNumber(nv_minus) + k;
                DdlVolPlus(OffsetDofVertexNumber(i) + k)
                  = mesh_num.OffsetDofVertexNumber(nv_plus) + k;
              }
          }
        
        // then edges
        for (int i = 0; i < nb_edges; i++)
          {
            int ne_minus = NumEdgeMinus(i);
            int ne_plus = NumEdgePlus(i);
            int offset_minus = mesh_num.OffsetDofEdgeNumber(ne_minus);
            int offset_plus = mesh_num.OffsetDofEdgeNumber(ne_plus);
            
            for (int k = 0; k < nb_dof_edge; k++)
              {
                DdlVolMinus(OffsetDofEdgeNumber(i) + k) = offset_minus + k;
                DdlVolPlus(OffsetDofEdgeNumber(i) + k) = offset_plus + k;
              }
          }
        
        // then faces
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int nf_minus = NumFaceMinus(i);
            int nf_plus = NumFacePlus(i);
            int nb_dof_inside = mesh_num.number_map.
              GetNbDofElement(order, mesh.BoundaryRef(OldNumFaceSurf(i)));
            
            int vol_minus = mesh_num.OffsetDofFaceNumber(nf_minus);
            int vol_plus = mesh_num.OffsetDofFaceNumber(nf_plus);
            for (int k = 0; k < nb_dof_inside; k++)
              {
                DdlVolMinus(OffsetDofFaceNumber(i) + k) = vol_minus + k;
                DdlVolPlus(OffsetDofFaceNumber(i) + k) = vol_plus + k;
              }
            
          }
      }

    
    /*var_problem.mesh_num.CheckMeshNumber();

    // checking correspondence of dofs
    VectR3 CoorDofs(mesh_num.GetNbDof());
    VectR3 pts_loc;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
      mesh.GetDofsElement(i, pts_loc, var.GetReferenceElement(i));
        for (int j = 0; j < mesh_num.Element(i).GetNbDof(); j++)
          CoorDofs(mesh_num.Element(i).GetNumberDof(j)) = pts_loc(j);
      }
    
    for (int i = 0; i < this->DdlVolMinus.GetM(); i++)
      if (CoorDofs(this->DdlVolMinus(i)) != CoorDofs(this->DdlVolPlus(i)))
        {
          cout << "Dofs are not matching" << endl;
          abort();
        }

    string name_file = "CoorDofsN" + to_str(nb_proc) + ".dat";
    ofstream file_out(name_file.data());
    for (int i = 0; i < CoorDofs.GetM(); i++)
      {
        PrintNoBrackets(file_out, CoorDofs(i));
        file_out << '\n';
      }
    
    file_out.close();
   
    for (int i = 0; i < NumFaceMinus.GetM(); i++)
      {
        int ne_minus = mesh.BoundaryRef(NumFaceMinus(i)).numElement(0);
        int ne_plus = mesh.BoundaryRef(NumFacePlus(i)).numElement(0);
        int nb_dof_face = NodleSurf(i).GetM();
        for (int j = 0; j < nb_dof_face; j++)
          {
            int dof_plus = DdlVolPlus(NodleSurf(i)(j));
            int dof_minus = DdlVolMinus(NodleSurf(i)(j));
            bool pres_dof = false;
            for (int k = 0; k < mesh_num.Element(ne_minus).GetNbDof(); k++)
              if (mesh_num.Element(ne_minus).GetNumberDof(k) == dof_minus)
                pres_dof = true;
            
            if (!pres_dof)
              {
                cout << "Problem on numbering" << endl;
                abort();
              }
            
            pres_dof = false;
            for (int k = 0; k < mesh_num.Element(ne_plus).GetNbDof(); k++)
              if (mesh_num.Element(ne_plus).GetNumberDof(k) == dof_plus)
                pres_dof = true;
            
            if (!pres_dof)
              {
                cout << "Problem on numbering" << endl;
                abort();
              }
          }    
      }
    */ 
    nb_ddl_volume_on_surface = DdlVolMinus.GetM();
    mesh_surf.Clear();
    
    //DISP(DdlVolMinus); DISP(DdlVolPlus);
    
    //mesh.Write("domain.mesh");
    if (nb_proc > 1)
      {
        int nb_faces_ref = NumFaceMinus.GetM();
        NodleVolMinus.Reallocate(nb_faces_ref);
        NodleVolPlus.Reallocate(nb_faces_ref);
        ProcVolMinus.Reallocate(nb_faces_ref);
        ProcVolPlus.Reallocate(nb_faces_ref);
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int nb_dof_face = NodleSurf(i).GetM();
            NodleVolMinus(i).Reallocate(nb_dof_face);
            NodleVolPlus(i).Reallocate(nb_dof_face);
            int ne_minus = mesh.BoundaryRef(NumFaceMinus(i)).numElement(0);
            int ne_plus = mesh.BoundaryRef(NumFacePlus(i)).numElement(0);
            ProcVolMinus(i) = Epart(ne_minus);
            ProcVolPlus(i) = Epart(ne_plus);
            for (int j = 0; j < nb_dof_face; j++)
              {
                int num_dof = NodleSurf(i)(j);
                NodleVolMinus(i)(j) = DdlVolMinus(num_dof);
                NodleVolPlus(i)(j) = DdlVolPlus(num_dof);
                /*if (NodleVolPlus(i)(j) == 109)
                  {
                    DISP(i); DISP(j); DISP(ne_minus); DISP(ne_plus);
                    DISP(Epart(ne_minus)); DISP(Epart(ne_plus));
                    DISP(ne_minus); DISP(ne_plus); DISP(NodleVolPlus(i));
                    DISP(mesh_num.Element(ne_plus).GetNodle());
                    } */
              }
          }
      }
  }
  
#endif
   
  //! arrays are distributed between processors in order to ensure parallel execution
  /*!
    // number of faces involved for the processor proc
    int nb_surfaces_involved(0);
    // number of dofs on the surface (for the current processor)
    int nb_loc_ddl_on_surface(0);
    // array containing all what is needed for processor proc
    // (this is the array that will be sent to the processor)
    IVect InfoSurfInvolved, InfoSurfMatching;
   */
  template<class Dimension>
  void VarTransmission_Base<Dimension>
  ::SendTransmissionDofs(const IVect& Epart,
                         int& nb_surfaces_involved, int& nb_loc_ddl_on_surface,
                         IVect& InfoSurfInvolved, IVect& InfoSurfMatching)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_problem.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif

    // if gamma_is_plate is true, the distribution is performed in VibroAcousticSession
    if (gamma_is_plate)
      return;
    
    if (exact_solution)
      transmission_references.Clear();
    
    if (transmission_references.GetM() <= 0)
      return;
    
    if (nb_proc == 1)
      return;
    
#ifdef SELDON_WITH_MPI
    if (rank_proc == 0)
      {

        // in order to sort local surface dofs in the same order as the global dofs
        Vector<bool> SurfDofUsed(nb_ddl_volume_on_surface);
        IVect InfoSurf;
        // list of global surface dofs for each processor (sorted list)
        Vector<IVect> ListSurfDofs(nb_proc);
        for (int proc = 0; proc < nb_proc; proc++)
          {
            // counting size of arrays
            SurfDofUsed.Fill(false);
            int nb_surf = 0, offset = 0;
            for (int i = 0; i < NodleVolMinus.GetM(); i++)
              {
                int p_plus = ProcVolPlus(i);
                int p_minus = ProcVolMinus(i);
                if ((p_plus == proc) || (p_minus == proc))
                  {
                    nb_surf++;
                    offset += 5 + 3*NodleVolMinus(i).GetM();
                    
                    // marking global surface dofs that are used by the processor proc
                    for (int j = 0; j < NodleSurf(i).GetM(); j++)
                      SurfDofUsed(NodleSurf(i)(j)) = true;
                  }
              }
            
            // counting dofs on the surface used by the processor proc
            int nb_ddl_surf = 0;
            for (int i = 0; i < nb_ddl_volume_on_surface; i++)
              if (SurfDofUsed(i))
                nb_ddl_surf++;
            
            // filling array to send to the processor proc
            InfoSurf.Reallocate(offset + nb_ddl_surf);
            offset = 0;
            for (int i = 0; i < NodleVolMinus.GetM(); i++)
              {
                int p_plus = ProcVolPlus(i);
                int p_minus = ProcVolMinus(i);
                // if one element across the face belongs to the processor proc
                if ((p_plus == proc) || (p_minus == proc))
                  {
                    // we send to the processor proc
                    // informations about this face 
                    InfoSurf(offset++) = NumFaceMinus(i); // global face number on Gamma-
                    InfoSurf(offset++) = NumFacePlus(i); // global face number on Gamma+
                    InfoSurf(offset++) = p_minus; // processor on Gamma-
                    InfoSurf(offset++) = p_plus; // processor on Gamma+
                    InfoSurf(offset++) = NodleVolMinus(i).GetM(); // number of dofs on the face
                    
                    // global volume dofs on Gamma- side
                    for (int j = 0; j < NodleVolMinus(i).GetM(); j++)
                      InfoSurf(offset++) = NodleVolMinus(i)(j);
                    
                    // global volume dofs on Gamma+ side
                    for (int j = 0; j < NodleVolPlus(i).GetM(); j++)
                      InfoSurf(offset++) = NodleVolPlus(i)(j);
                    
                    // global surface dofs
                    for (int j = 0; j < NodleSurf(i).GetM(); j++)
                      InfoSurf(offset++) = NodleSurf(i)(j);
                  }
              }
            
            // finally, list of global surface dofs
            ListSurfDofs(proc).Reallocate(nb_ddl_surf);
            nb_ddl_surf = 0;
            for (int i = 0; i < nb_ddl_volume_on_surface; i++)
              if (SurfDofUsed(i))
                {
                  InfoSurf(offset++) = i;
                  ListSurfDofs(proc)(nb_ddl_surf) = i;
                  nb_ddl_surf++;
                }
            
            // sending nb_surf and size of array InfoSurf
            IVect num(6);
            num(0) = nb_surf; num(1) = offset;
            num(2) = nb_ddl_volume_on_surface; num(3) = nb_ddl_surf;
            num(4) = RefGammaMoins; num(5) = RefGammaPlus;
            if (proc == 0)
              {
                nb_surfaces_involved = nb_surf;
                InfoSurfInvolved = InfoSurf;
                nb_loc_ddl_on_surface = nb_ddl_surf;
              }
            else
              {
                MPI_Send(num.GetData(), num.GetM(), MPI_INTEGER, proc, 42, comm);
                MPI_Send(InfoSurf.GetData(), InfoSurf.GetM(), MPI_INTEGER, proc, 43, comm);
              }
          }

        if (var_problem.print_level >= 2)
          cout << "Starting to compute MatchingDofOrig_Subdomain "
               << "for transmission surface ... " << endl;
        
        // expensive loop to construct MatchingDofOrig_Subdomain
        InfoSurf.Reallocate((nb_ddl_volume_on_surface+1)*nb_proc);
        for (int proc = 0; proc < nb_proc; proc++)
          {
            // in InfoSurf, we store for each processor the number of dofs interacting
            // and the list of dofs numbers interacting with the other processor
            int offset = 0;
            for (int proc2 = 0; proc2 < nb_proc; proc2++)
              if (proc != proc2)
                {
                  int nb_dof_interac = 0;
                  int offset_begin = offset++;
                  
                  // the two lists are sorted, we find all the dofs (in a single loop)
                  // such that ListSurfDofs(proc)(k) = ListSurfDofs(proc2)(k2)
                  int k = 0;
                  for (int j = 0; j < ListSurfDofs(proc).GetM(); j++)
                    {
                      int dof = ListSurfDofs(proc)(j);
                      while ( (k < ListSurfDofs(proc2).GetM())
                              && (ListSurfDofs(proc2)(k) < dof))
                        k++;
                      
                      if (k < ListSurfDofs(proc2).GetM())
                        if (dof == ListSurfDofs(proc2)(k))
                          {
                            InfoSurf(offset++) = j;
                            nb_dof_interac++;
                          }
                    }
                  
                  InfoSurf(offset_begin) = nb_dof_interac;
                }
            
            // sending InfoSurf to the processor proc
            if (proc == 0)
              InfoSurfMatching = InfoSurf;
            else
              {
                //DISP(proc); DISP(offset); DISP(InfoSurf.GetM());
                MPI_Send(&offset, 1, MPI_INTEGER, proc, 40, comm);
                MPI_Send(InfoSurf.GetData(), offset, MPI_INTEGER, proc, 41, comm);
              }
          }

        if (var_problem.print_level >= 2)
          cout << "MatchingDofOrig_Subdomain computed for transmission surface " << endl;
      }
    else
      {
        MPI_Status status;
        // receiving the size of array InfoSurfInvolved, number of dofs (local, global)
        // reference for gamma- and gamma+
        IVect num(6);
        MPI_Recv(num.GetData(), num.GetM(), MPI_INTEGER, 0, 42, comm, &status);
        nb_surfaces_involved = num(0);
        int size_array = num(1);
        nb_ddl_volume_on_surface = num(2);
        nb_loc_ddl_on_surface = num(3);
        RefGammaMoins = num(4);
        RefGammaPlus = num(5);
        
        // receiving array InfoSurfInvolved
        InfoSurfInvolved.Reallocate(size_array);
        MPI_Recv(InfoSurfInvolved.GetData(), size_array, MPI_INTEGER, 0, 43, comm, &status);
        
        // receiving array InfoSurfMatching
        MPI_Recv(&size_array, 1, MPI_INTEGER, 0, 40, comm, &status);
        //DISP(size_array);
        
        InfoSurfMatching.Reallocate(size_array);
        MPI_Recv(InfoSurfMatching.GetData(), size_array, MPI_INTEGER, 0, 41, comm, &status);
      }    

    // exchanging jacobian matrices on the boundary
    Vector<int64_t> xtmp;
    if (rank_proc == 0)
      {
        // loop over processors
        for (int rank = 0; rank < nb_proc; rank++)
          {
            // counting the number of elements to send
            int nb_elt = 0, nb_nodes = 0;
            for (int i = 0; i < MatDfjNodes.GetM(); i++)
              if ((Epart(i) == rank) && (MatDfjNodes(i).GetM() > 0))
                {
                  nb_elt++;
                  nb_nodes += MatDfjNodes(i).GetM();
                }
            
            if (rank == 0)
              {
                nb_elt_local = nb_elt;
                nb_nodes_local = nb_nodes;
              }
            else
              {
                MPI_Send(&nb_elt, 1, MPI_INTEGER, rank, 92, comm);
                MPI_Send(&nb_nodes, 1, MPI_INTEGER, rank, 93, comm);
              }
            
            // loop over elements
            IVect nb_pts_quad(nb_elt + nb_nodes + 1);
            nb_pts_quad.Fill(0); nb_elt = 0;
            int nb = 0, nb_total_pts = 0;
            for (int i = 0; i < MatDfjNodes.GetM(); i++)
              if ((Epart(i) == rank) && (MatDfjNodes(i).GetM() > 0))
                {
                  nb_pts_quad(nb++) = i;
                  for (int j = 0; j < MatDfjNodes(i).GetM(); j++)
                    {
                      nb_pts_quad(nb++) = MatDfjNodes(i)(j).GetM();
                      nb_total_pts += MatDfjNodes(i)(j).GetM();
                    }
                }
            
            nb_total_pts *= Dimension::dim_N*Dimension::dim_N;
            nb_pts_quad(nb++) = nb_total_pts;
            
            if (rank == 0)
              nb_pts_quad_local = nb_pts_quad;
            else
              MPI_Send(nb_pts_quad.GetData(), nb_pts_quad.GetM(),
                       MPI_INTEGER, rank, 94, comm);
            
            VectReal_wp coef_dfj_quad(nb_total_pts);
            nb = 0;
            for (int i = 0; i < MatDfjNodes.GetM(); i++)
              if ((Epart(i) == rank) && (MatDfjNodes(i).GetM() > 0))
                {
                  for (int j = 0; j < MatDfjNodes(i).GetM(); j++)
                    for (int k = 0; k < MatDfjNodes(i)(j).GetM(); k++)
                      {
                        for (int p = 0; p < Dimension::dim_N; p++)
                          for (int q = 0; q < Dimension::dim_N; q++)
                            coef_dfj_quad(nb++) = MatDfjNodes(i)(j)(k)(p, q);
                      }
                }
            
            if (rank == 0)
              coef_dfj_quad_local = coef_dfj_quad;
            else
              MpiSsend(comm, coef_dfj_quad, xtmp, coef_dfj_quad.GetM(), rank, 95);
          }
      }
    else
      {
        MPI_Status status;
        MPI_Recv(&nb_elt_local, 1, MPI_INTEGER, 0, 92, comm, &status);
        MPI_Recv(&nb_nodes_local, 1, MPI_INTEGER, 0, 93, comm, &status);
        
        nb_pts_quad_local.Reallocate(nb_elt_local+nb_nodes_local+1);
        MPI_Recv(nb_pts_quad_local.GetData(), nb_pts_quad_local.GetM(),
                 MPI_INTEGER, 0, 94, comm, &status);
        
        int N = nb_pts_quad_local(nb_elt_local+nb_nodes_local);
        coef_dfj_quad_local.Reallocate(N);
        MpiRecv(comm, coef_dfj_quad_local, xtmp, N, 0, 95, status);
      }
    
    MatDfjNodes.Clear();
#endif
    
  }


  //! arrays are distributed between processors in order to ensure parallel execution
  /*!
    // number of faces involved for the processor proc
    int nb_surfaces_involved(0);
    // number of dofs on the surface (for the current processor)
    int nb_loc_ddl_on_surface(0);
    // array containing all what is needed for processor proc
    // (this is the array that will be sent to the processor)
    IVect InfoSurfInvolved, InfoSurfMatching;
   */
  template<class Dimension>
  void VarTransmission_Base<Dimension>
  ::DistributeTransmissionDofs(int nb_surfaces_involved, int nb_loc_ddl_on_surface,
                               IVect& InfoSurfInvolved, IVect& InfoSurfMatching)
  {    
#ifdef SELDON_WITH_MPI    
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    if (gamma_is_plate)
      return;
    
    if (transmission_references.GetM() <= 0)
      return;
    
    if (nb_proc == 1)
      return;

#ifdef SELDON_WITH_MPI    
    //var_problem.mesh.Write("subdomain"+to_str(rank_proc)+".mesh");
    // clearing arrays 
    NumFaceMinus.Clear(); NumFacePlus.Clear();
    NodleVolMinus.Clear(); NodleVolPlus.Clear();
    NodleSurf.Clear(); DdlSurfLambda.Clear();
    DdlVolPlus.Clear(); DdlVolMinus.Clear();
    ProcVolPlus.Clear(); ProcVolMinus.Clear();
    MatchingNumber_Subdomain.Clear();
    MatchingDofOrig_Subdomain.Clear();
    OldNumFaceSurf.Clear(); OldNumVertexSurf.Clear(); IndexVert.Clear();

    // retrieving for each face : processor number on gamma- and gamma+
    // and degrees of freedom on gamma- and gamma+
    IVect NumGlobalFaceMinus(nb_surfaces_involved);
    IVect NumGlobalFacePlus(nb_surfaces_involved);
    IVect ProcSurfMinus, ProcSurfPlus;
    ProcSurfMinus.Reallocate(nb_surfaces_involved);
    ProcSurfPlus.Reallocate(nb_surfaces_involved);
    
    NodleSurf.Reallocate(nb_surfaces_involved);
    NodleVolMinus.Reallocate(nb_surfaces_involved);
    NodleVolPlus.Reallocate(nb_surfaces_involved);
    
    int offset = 0;
    for (int i = 0; i < nb_surfaces_involved; i++)
      {
        // global numbers for faces on Gamma- and Gamma+
        NumGlobalFaceMinus(i) = InfoSurfInvolved(offset++);
        NumGlobalFacePlus(i) = InfoSurfInvolved(offset++);
        // processors on Gamma- and Gamma+ 
        ProcSurfMinus(i) = InfoSurfInvolved(offset++);
        ProcSurfPlus(i) = InfoSurfInvolved(offset++);
        int nb_dof = InfoSurfInvolved(offset++);
        // dof numbers on Gamma- and Gamma+ (for dofs associated with E)
        NodleVolMinus(i).Reallocate(nb_dof);
        NodleVolPlus(i).Reallocate(nb_dof);
        NodleSurf(i).Reallocate(nb_dof);
        for (int j = 0; j < nb_dof; j++)
          NodleVolMinus(i)(j) = InfoSurfInvolved(offset++);
        
        for (int j = 0; j < nb_dof; j++)
          NodleVolPlus(i)(j) = InfoSurfInvolved(offset++);
        
        // dof numbers for lambda
        for (int j = 0; j < nb_dof; j++)
          NodleSurf(i)(j) = InfoSurfInvolved(offset++);
      }
    
    // global dofs numbers for lambda
    DdlSurfLambda.Reallocate(nb_loc_ddl_on_surface);
    for (int i = 0; i < nb_loc_ddl_on_surface; i++)
      DdlSurfLambda(i) = InfoSurfInvolved(offset++);    

    InfoSurfInvolved.Clear();
    
    // global to local numbering for volume dofs
    IVect Glob_to_local(var_problem.GetNbGlobalMeshDof());
    Glob_to_local.Fill(-1);
    //DISP(var.mesh_num.GlobDofNumber_Subdomain);
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    for (int i = 0; i < mesh_num.GetNbDof(); i++)
      Glob_to_local(mesh_num.GlobDofNumber_Subdomain(i)) = i;
    
    // global to local numbering for surface dofs
    IVect Glob_surf_to_local(nb_ddl_volume_on_surface);
    Glob_surf_to_local.Fill(-1);
    for (int i = 0; i < nb_loc_ddl_on_surface; i++)
      Glob_surf_to_local(DdlSurfLambda(i)) = i;
    
    // constructing DdlVolMinus, DdlVolPlus, DdlSurfLambda
    DdlVolMinus.Reallocate(nb_loc_ddl_on_surface);
    DdlVolPlus.Reallocate(nb_loc_ddl_on_surface);
    ProcVolMinus.Reallocate(nb_loc_ddl_on_surface);
    ProcVolPlus.Reallocate(nb_loc_ddl_on_surface);
    DdlVolMinus.Fill(-1); DdlVolPlus.Fill(-1);
    ProcVolMinus.Fill(-1); ProcVolPlus.Fill(-1);
    
    
    for (int i = 0; i < NodleVolMinus.GetM(); i++)
      for (int j = 0; j < NodleVolMinus(i).GetM(); j++)
        {
          int n = NodleSurf(i)(j);
          int dof = Glob_surf_to_local(n);
          int proc = ProcSurfMinus(i);
          if (ProcVolMinus(dof) != rank_proc)
            {
              ProcVolMinus(dof) = proc;
              if (proc == rank_proc)
                DdlVolMinus(dof) = Glob_to_local(NodleVolMinus(i)(j));
              else
                DdlVolMinus(dof) = NodleVolMinus(i)(j);
            }

          proc = ProcSurfPlus(i);
          if (ProcVolPlus(dof) != rank_proc)
            {
              ProcVolPlus(dof) = proc;
              if (proc == rank_proc)
                DdlVolPlus(dof) = Glob_to_local(NodleVolPlus(i)(j));
              else
                DdlVolPlus(dof) = NodleVolPlus(i)(j);
            }
        }
    
    for (int i = 0; i < nb_loc_ddl_on_surface; i++)
      {
        if ((this->DdlVolPlus(i) == -1) || (this->DdlVolMinus(i) == -1))
          {
            cout << "Invalid dof" << endl;
            abort();
          }
      }
    
    // constructing MatchingDofOrig_Surface
    int nb_dom = 0; offset = 0;
    for (int proc = 0; proc < nb_proc; proc++)
      if (proc != rank_proc)
        {
          if (InfoSurfMatching(offset) > 0)
            nb_dom++;
          
          offset += InfoSurfMatching(offset)+1;
        }
    
    MatchingNumber_Subdomain.Reallocate(nb_dom);
    MatchingDofOrig_Subdomain.Reallocate(nb_dom);
    nb_dom = 0; offset = 0;
    for (int proc = 0; proc < nb_proc; proc++)
      if (proc != rank_proc)
        {
          int nb_dof = InfoSurfMatching(offset); offset++;
          if (nb_dof > 0)
            {
              MatchingNumber_Subdomain(nb_dom) = proc;
              MatchingDofOrig_Subdomain(nb_dom).Reallocate(nb_dof);
              for (int j = 0; j < nb_dof; j++)
                MatchingDofOrig_Subdomain(nb_dom)(j) = InfoSurfMatching(offset++);
              
              nb_dom++;
            }
        }

    InfoSurfMatching.Clear();    

    // checking correspondence of dofs
    /*
      typename Dimension::VectR_N CoorDofs(var.mesh_num.GetNbDof());
    typename Dimension::VectR_N pts_loc;
    for (int i = 0; i < var.mesh.GetNbElt(); i++)
      {
      var.mesh.GetDofsElement(i, pts_loc, var.GetReferenceElement(i));
        
        for (int j = 0; j < var.mesh_num.Element(i).GetNbDof(); j++)
          CoorDofs(var.mesh_num.Element(i).GetNumberDof(j)) = pts_loc(j);
      }
    
    string name_file = "CoorDofs" + to_str(rank_proc) 
    + "N" + to_str(nb_proc) + ".dat";
    
    ofstream file_out(name_file.data());
    for (int i = 0; i < CoorDofs.GetM(); i++)
      {
        PrintNoBrackets(file_out, CoorDofs(i));
        file_out << '\n';
      }
    
    file_out.close();
    */

    // filling MatDfjNodes
    MatDfjNodes.Reallocate(var_problem.mesh.GetNbElt());
    int nb = 0, ind = 0, num_elem = 0;
    for (int i = 0; i < nb_elt_local; i++)
      {
        int num_elem_glob = nb_pts_quad_local(nb++);
        while (num_elem_glob > var_problem.mesh.GlobElementNumber_Subdomain(num_elem))
          num_elem++;
        
        if (num_elem_glob == var_problem.mesh.GlobElementNumber_Subdomain(num_elem))
          {
            int nb_faces_loc = var_problem.mesh.Element(num_elem).GetNbBoundary();
            MatDfjNodes(num_elem).Reallocate(nb_faces_loc);
            for (int j = 0; j < nb_faces_loc; j++)
              {
                int nb_quad = nb_pts_quad_local(nb++);
                MatDfjNodes(num_elem)(j).Reallocate(nb_quad);
                for (int k = 0; k < nb_quad; k++)
                  {
                    for (int p = 0; p < Dimension::dim_N; p++)
                      for (int q = 0; q < Dimension::dim_N; q++)
                        MatDfjNodes(num_elem)(j)(k)(p, q)
                          = coef_dfj_quad_local(ind++);
                  }
              }
          }
        else
          {
            cout << "Unknown element" << endl;
            abort();
          }
      }

#endif
      
  }


  template<class Dimension>
  void VarTransmission_Base<Dimension>::Clear()
  {
    transmission_references.Clear();
    
    NumFaceMinus.Clear();
    NumFacePlus.Clear();
    DdlVolMinus.Clear();
    DdlVolPlus.Clear();
    
    NodleSurf.Clear(); 
    NodleAcous.Clear();
    NodleMeca.Clear();
    NodleVolMinus.Clear();
    NodleVolPlus.Clear();
    
    ProcVolMinus.Clear();
    ProcVolPlus.Clear();
    JacobianAcous.Clear();
    DirichletAcous.Clear();
    
    OldNumFaceSurf.Clear();
    OldNumVertexSurf.Clear();
    IndexVert.Clear();
  }


  template<class Dimension>
  void VarTransmission_Base<Dimension>
  ::CopyArray(IVect& ddl_vol_minus, IVect& ddl_vol_plus, 
	      Vector<IVect>& NodleAcous_, Vector<IVect>& NodleMeca_,
	      Vector<VectReal_wp>& JacobianAcous_, IVect& DirichletAcous_)
  {
    ddl_vol_minus = DdlVolMinus;
    ddl_vol_plus = DdlVolPlus;
    
    //NumFaceMinus = var_volume.var_harmonic.var_transmission.NumFaceMinus;
    //NumFacePlus = var_volume.var_harmonic.var_transmission.NumFacePlus;
    NumFaceMinus.Clear();
    NumFacePlus.Clear();
    
    NodleAcous_ = NodleAcous;
    NodleMeca_ = NodleMeca;
    
    NodleAcous.Clear();
    NodleMeca.Clear();
    
    JacobianAcous_ = JacobianAcous;
    DirichletAcous_ = DirichletAcous;
    JacobianAcous.Clear();
    DirichletAcous.Clear();    
  }
  
} // end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_CXX
#endif

  
