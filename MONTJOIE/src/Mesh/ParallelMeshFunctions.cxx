#ifndef MONTJOIE_FILE_PARALLEL_MESH_FUNCTIONS_CXX

namespace Montjoie
{
  
  //! sending a mesh to a processor
  /*!
    \param[in] glob_mesh global mesh (a part of this mesh is actually sent)
    \param[in] glob_mesh_num global numberings
    \param[in] NumElement list of element numbers associated with the sub-mesh to send
    \param[in] Epart Epart(i) is the processor number associated with element i
    \param[in] NumLoc NumLoc(i) is the local element number of global element i
    \param[inout] mesh sub-mesh to send
    \param[inout] proc mesh will be send to this processor
    \param[inout] comm MPI communicator
    \param[in] tag MPI communications will use tag+something as MPI tag
    \param[out] param other parameters that are sent
  */
  template<>
  void SendMeshToProcessor(const Mesh<Dimension2>& glob_mesh,
                           const Vector<MeshNumbering<Dimension2>* >& glob_mesh_num,
			   const IVect& NumElement, const IVect& Epart, const IVect& NumLoc,
                           Mesh<Dimension2>& mesh, int proc, const MPI_Comm& comm, ParamParallelMesh<Dimension2>& Param, int tag)
  {
    IVect& num = Param.num; VectReal_wp& all_param = Param.all_param;
    VectR2& PointsEdgeRef = Param.PointsEdgeRef; IVect& NumBoundary = Param.NumBoundary;
    IVect& NumDofNeighbor = Param.NumDofNeighbor; IVect& NumElem = Param.NumElem;
    IVect& NumDof = Param.NumDof; Vector<IVect>& MinimalProc = Param.MinimalProc;
    VectR2& TranslatPeriodicEdge = Param.TranslatPeriodicEdge;
    IVect& OffsetDofV = Param.OffsetDofV; IVect& NumDofPML = Param.NumDofPML;
    
    int rank; MPI_Comm_rank(comm, &rank);
    
    // the array all_param is constructed and will contain
    // all the parameters (such as parameters of the curves) of the mesh
    // that we want to send to other processors
    int nb_ref = mesh.GetNbReferences();
    int nb_param = 5*nb_ref+3;
    for (int i = 0; i < glob_mesh.GetNbPmlAreas(); i++)
      nb_param += glob_mesh.GetPmlArea(i).GetNbParameters();
    
    VectReal_wp param;
    for (int i = 1; i <= nb_ref; i++)
      {
        mesh.GetCurveParameter(i, param);
        nb_param += param.GetM();
      }

    nb_param += 3*glob_mesh.GetNbPeriodicReferences();
    
    Vector<int64_t> all_param_tmp;
    IVect ref_neighbor = mesh.GetOriginalNeighborReference();
    all_param.Reallocate(nb_param); all_param.Fill(0);    
    nb_param = 5*nb_ref+1;
    all_param(0) = nb_param;
    for (int i = 1; i <= mesh.GetNbReferences(); i++)
      {
        mesh.GetCurveParameter(i, param);
        for (int j = 0; j < param.GetM(); j++)
          all_param(nb_param+j) = param(j);
        
        nb_param += param.GetM();
        all_param(i) = nb_param;
        all_param(nb_ref+i) = mesh.GetBodyNumber(i);
        all_param(2*nb_ref+i) = mesh.GetBoundaryCondition(i);
        all_param(3*nb_ref+i) = mesh.GetCurveType(i);
        all_param(4*nb_ref+i) = ref_neighbor(i);
      }    
    
    all_param(nb_param) = glob_mesh.GetPeriodicAlpha(); nb_param++;
    all_param(nb_param) = glob_mesh.GetNbPeriodicReferences(); nb_param++;
    for (int i = 0; i < glob_mesh.GetNbPeriodicReferences(); i++)
      {
	all_param(nb_param) = glob_mesh.GetPeriodicReference(i, 0);
	all_param(nb_param+1) = glob_mesh.GetPeriodicReference(i, 1);
	all_param(nb_param+2) = glob_mesh.GetPeriodicityTypeReference(i);
	nb_param += 3;
      }
    
    for (int i = 0; i < glob_mesh.GetNbPmlAreas(); i++)
      glob_mesh.GetPmlArea(i).FillParameters(all_param, nb_param);
    
    // sending number of vertices, referenced edges and elements
    int nb_vertices = mesh.GetNbVertices();
    int nb_edges_ref = mesh.GetNbBoundaryRef();
    int nb_elt = mesh.GetNbElt();
    int nb_int_elt = 0;
    for (int i = 0; i < nb_elt; i++)
      nb_int_elt += mesh.Element(i).GetNbVertices();

    // counting dofs
    int nb_mesh_num = glob_mesh_num.GetM();
    int nb_matching_proc = 0, nb_matching_dofs = 0;
    
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    Vector<bool> ProcMatching(nb_proc);
    ProcMatching.Fill(false);
    
    Vector<Vector<IVect> > MatchingDofs(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      MatchingDofs(i).Reallocate(nb_proc);
    
    for (int n = 0; n < glob_mesh_num.GetM(); n++)
      if (glob_mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        {
          // returning to positive numbers
          for (int i = 0; i < MinimalProc(n).GetM(); i++)
            MinimalProc(n)(i) = abs(MinimalProc(n)(i));
          
          // we put opposite signs in MinimalProc for dofs associated with the processor proc            
          for (int i = 0; i < nb_elt; i++)
            {
              int ne = NumElement(i);
              for (int j = 0; j < glob_mesh_num(n)->Element(ne).GetNbDof(); j++)
                {
                  int num_dof = glob_mesh_num(n)->Element(ne).GetNumberDof(j);
                  if ((num_dof >= 0) && (MinimalProc(n)(num_dof) > 0))
                    MinimalProc(n)(num_dof) = -MinimalProc(n)(num_dof);
                }
            }
          
          // the number of dofs shared with other processors is estimated
          IVect nb_dofs_proc(nb_proc);
          nb_dofs_proc.Zero();
          for (int i = 0; i < glob_mesh.GetNbElt(); i++)
            for (int j = 0; j < glob_mesh_num(n)->Element(i).GetNbDof(); j++)
              {
                int num_dof = glob_mesh_num(n)->Element(i).GetNumberDof(j);
                if ((num_dof >= 0) && (MinimalProc(n)(num_dof) < 0))
                  nb_dofs_proc(Epart(i))++;
              }
          
          for (int p = 0; p < MatchingDofs(n).GetM(); p++)
            if (p != proc)
              MatchingDofs(n)(p).Reallocate(nb_dofs_proc(p));
          
          nb_dofs_proc.Zero();
          for (int i = 0; i < glob_mesh.GetNbElt(); i++)
            {
              int p = Epart(i);
              if (p != proc)
                for (int j = 0; j < glob_mesh_num(n)->Element(i).GetNbDof(); j++)
                  {
                    int num_dof = glob_mesh_num(n)->Element(i).GetNumberDof(j);
                    if ((num_dof >= 0) && (MinimalProc(n)(num_dof) < 0))
                      {
                        MatchingDofs(n)(p)(nb_dofs_proc(p)) = num_dof;
                        nb_dofs_proc(p)++;
                      }
                  }
            }
          
          // duplicates are removed and numbers sorted
          for (int p = 0; p < MatchingDofs(n).GetM(); p++)
            if (MatchingDofs(n)(p).GetM() > 0)
              {
                if (!ProcMatching(p))
                  {
                    ProcMatching(p) = true;
                    nb_matching_proc++;
                  }
                
                RemoveDuplicate(MatchingDofs(n)(p));
                nb_matching_dofs += MatchingDofs(n)(p).GetM();
              }        
        }
    
    num.Reallocate(20 + 2*nb_mesh_num);
    num.Fill(0);
    num(0) = nb_vertices;
    num(1) = nb_edges_ref;
    num(2) = nb_elt;
    num(3) = nb_ref;
    num(4) = nb_param;
    num(5) = nb_int_elt;
    num(6) = mesh.GetGeometryOrder();
    num(7) = nb_mesh_num;
    num(8) = nb_matching_proc;
    num(9) = nb_matching_dofs;
    num(11) = OffsetDofV.GetM();
    if (glob_mesh_num.GetM() > 0)
      {
        num(12) = glob_mesh_num(0)->compute_dof_pml;
        num(13) = glob_mesh_num(0)->GetVariableOrder();
        num(17) = glob_mesh_num(0)->GetFormulationForPeriodicCondition();
      }
    
    num(14) = glob_mesh.GetNbPmlAreas();
    num(15) = glob_mesh.IsOnlyQuadrilateral();
    num(16) = glob_mesh.GetNbBoundaryRef();

    for (int n = 0; n < nb_mesh_num; n++)
      {
        num(20 + n) = glob_mesh_num(n)->GetOrder();
        num(20 + n + nb_mesh_num) = glob_mesh_num(n)->number_map.FormulationDG();
      }
    
    if (rank != proc)
      {
        int nb_param_int = num.GetM();
        MPI_Send(&nb_param_int, 1, MPI_INTEGER, proc, tag, comm); tag++;
        MPI_Send(num.GetData(), nb_param_int, MPI_INTEGER, proc, tag, comm); tag++;
      }
    
    // sending parameters of curves
    if (rank != proc)
      {
        MpiSsend(comm, all_param, all_param_tmp, nb_param, proc, tag);
        tag++;
      }
    
    // sending vertices and internal points of edges
    int r = mesh.GetGeometryOrder();
    PointsEdgeRef.Reallocate(nb_vertices + (r-1)*nb_edges_ref);
    for (int i = 0; i < nb_vertices; i++)
      PointsEdgeRef(i) = mesh.Vertex(i);
    
    int nb = nb_vertices;
    if (r > 1)
      {
        for (int i = 0; i < nb_edges_ref; i++)
          for (int k = 0; k < r-1; k++)
            PointsEdgeRef(nb++) = mesh.GetPointInsideEdge(i, k);
      }
    
    if (rank != proc)
      {
        MpiSsend(comm, reinterpret_cast<Real_wp*>(PointsEdgeRef.GetData()),
                 all_param_tmp, 2*nb, proc, tag);
        tag++;
      }    
    
    // sending edges
    int nb_param_per_edge = 14 + 2*nb_mesh_num;
    NumBoundary.Reallocate(nb_param_per_edge*nb_edges_ref);
    NumBoundary.Fill(-1);
    nb = 0; int nb_dof_neighbor = 0;
    for (int i = 0; i < nb_edges_ref; i++)
      {
        NumBoundary(nb) = mesh.BoundaryRef(i).numVertex(0);
        NumBoundary(nb+1) = mesh.BoundaryRef(i).numVertex(1);
        NumBoundary(nb+2) = mesh.BoundaryRef(i).GetReference();
        int ne = mesh.BoundaryRef(i).numElement(0);
        int num_loc = mesh.Element(ne).GetPositionBoundary(i);
        int num_elem = NumElement(ne);
        int iglob = glob_mesh.Element(num_elem).numEdge(num_loc);
        if (glob_mesh_num.GetM() > 0)
          NumBoundary(nb+3) = glob_mesh_num(0)->GetOrderQuadrature(iglob);
        
	NumBoundary(nb+9) = iglob;
        if (glob_mesh_num.GetM() > 0)
          {
            NumBoundary(nb+12) = glob_mesh_num(0)->GetPeriodicBoundary(iglob);
            NumBoundary(nb+13) = glob_mesh_num(0)->GetPeriodicityTypeForBoundary(iglob);
          }
        
        int ref = mesh.BoundaryRef(i).GetReference();
        
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            // referenced edge at the interface between subdomains            
            int num_elt2 = glob_mesh.Boundary(iglob).numElement(0);
            if (num_elt2 == num_elem)
              num_elt2 = glob_mesh.Boundary(iglob).numElement(1);
            
            //int num_loc2 = glob_mesh.Element(num_elt2).numEdge(num_loc);
            int num_loc2 = -1, rot = -1;
            if (glob_mesh_num.GetM() > 0)
              glob_mesh_num(0)->GetBoundaryRotation(iglob, num_elem, num_elt2, num_loc, num_loc2, rot);
            
            NumBoundary(nb+4) = num_elt2;
            NumBoundary(nb+5) = Epart(num_elt2); 
            NumBoundary(nb+6) = NumLoc(num_elt2);
            NumBoundary(nb+7) = num_loc2;
            NumBoundary(nb+8) = rot;
            for (int n = 0; n < glob_mesh_num.GetM(); n++)
              {
                NumBoundary(nb+14+n) = glob_mesh_num(n)->GetOrderElement(num_elt2);
                NumBoundary(nb+14+n+nb_mesh_num) = glob_mesh_num(n)->GetNbLocalDof(num_elt2);
                nb_dof_neighbor += glob_mesh_num(n)->GetNbLocalDof(num_elt2);
              }
            
            NumBoundary(nb+10) = glob_mesh.GetTypeElement(num_elt2);            
            NumBoundary(nb+11) = glob_mesh.Element(num_elt2).GetReference();            
          }
        
        nb += nb_param_per_edge;
      }
    
    if (rank != proc)
      {
        MPI_Send(NumBoundary.GetData(), nb, MPI_INTEGER, proc, tag, comm);    
        tag++;
      }

    // and periodic edges
    if ((glob_mesh.GetNbPeriodicReferences() > 0) && (glob_mesh_num.GetM() > 0))
      {
	TranslatPeriodicEdge.Reallocate(nb_edges_ref);
	for (int i = 0; i < nb_edges_ref; i++)
	  {
            int num = NumBoundary(nb_param_per_edge*i+12);
            if (num >= 0)
              TranslatPeriodicEdge(i) = glob_mesh_num(0)->GetTranslationPeriodicBoundary(num);
            else
              TranslatPeriodicEdge(i).Zero();
          }
	
	if (rank != proc)
	  {
	    MpiSsend(comm, reinterpret_cast<Real_wp*>(TranslatPeriodicEdge.GetData()), all_param_tmp,
		     2*nb_edges_ref, proc, tag);
	    tag++;
	  }
      }
    
    // dofs associated with neighboring elements
    if (glob_mesh_num.GetM() > 0)
      {
        if (glob_mesh_num(0)->compute_dof_pml)
          NumDofNeighbor.Reallocate(2*nb_dof_neighbor);
        else
          NumDofNeighbor.Reallocate(nb_dof_neighbor);
    
        nb = 0;
        for (int i = 0; i < nb_edges_ref; i++)
          {
            int ref = mesh.BoundaryRef(i).GetReference();
            if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
              for (int n = 0; n < glob_mesh_num.GetM(); n++)
                {
                  int num_elt2 = NumBoundary(nb_param_per_edge*i + 4);
                  IVect Nodle = glob_mesh_num(n)->Element(num_elt2).GetNodle();
                  if (glob_mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
                    {
                      Nodle.Resize(glob_mesh_num(n)->GetNbLocalDof(num_elt2));
                      for (int j = 1; j < Nodle.GetM(); j++)
                        Nodle(j) = Nodle(0) + j;
                    }
                  
                  for (int j = 0; j < Nodle.GetM(); j++)
                    NumDofNeighbor(nb+j) = Nodle(j);
                
                  nb += Nodle.GetM();
                  
                  if (glob_mesh_num(n)->compute_dof_pml)
                    {
                      for (int j = 0; j < Nodle.GetM(); j++)
                        {
                          if (Nodle(j) >= 0)
                            NumDofNeighbor(nb+j) = glob_mesh_num(n)->GetDofPML(Nodle(j));
                          else
                            NumDofNeighbor(nb+j) = -1;
                        }
                      
                      nb += Nodle.GetM();
                    }
                }
          }
        
        if (rank != proc)
          {
            MPI_Send(NumDofNeighbor.GetData(), nb, MPI_INTEGER, proc, tag, comm);    
            tag++;
          }
      }
    
    // and elements
    NumElem.Reallocate((5 + 3*nb_mesh_num)*nb_elt + (3+nb_mesh_num)*nb_int_elt);
    nb = 0; int nb_dof = 0;
    for (int i = 0; i < nb_elt; i++)
      {
        NumElem(nb) = mesh.Element(i).GetNbVertices();
        NumElem(nb+1) = mesh.Element(i).GetReference();
        NumElem(nb+2) = mesh.Element(i).GetNumberPML();
        NumElem(nb+3) = mesh.Element(i).GetTypePML();
        int num_elem = NumElement(i);
        NumElem(nb+4) = NumElement(i);
        for (int n = 0; n < glob_mesh_num.GetM(); n++)
          {
            NumElem(nb+5+n) = glob_mesh_num(n)->GetOrderElement(num_elem);
            NumElem(nb+5+n+nb_mesh_num) = glob_mesh_num(n)->GetOrderInside(num_elem);   
            NumElem(nb+5+n+2*nb_mesh_num) = glob_mesh_num(n)->GetNbLocalDof(num_elem);
            if (glob_mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
              nb_dof += NumElem(nb+5+n+2*nb_mesh_num);
          }
        
        nb += 5 + 3*nb_mesh_num;
      }
    
    for (int i = 0; i < nb_elt; i++)
      {
        int ne = NumElement(i);
        int nv = mesh.Element(i).GetNbVertices();
        for (int j = 0; j < nv; j++)
          {
            NumElem(nb+j) = mesh.Element(i).numVertex(j); 
            NumElem(nb+nv+j) = glob_mesh.Element(ne).numVertex(j);
            
            int num_edge = glob_mesh.Element(ne).numEdge(j);
            NumElem(nb+2*nv+j) = num_edge;
            for (int n = 0; n < glob_mesh_num.GetM(); n++)
              NumElem(nb+(3+n)*nv+j) = glob_mesh_num(n)->GetOrderEdge(num_edge);
          }
        
        nb += (3+nb_mesh_num)*nv;
      }
    
    if (rank != proc)
      {
        MPI_Send(NumElem.GetData(), nb, MPI_INTEGER, proc, tag, comm);    
        tag++;
      }
        
    // then dofs
    int nb_mesh_num_discont = 0, nb_mesh_num_cont = 0;
    for (int n = 0; n < glob_mesh_num.GetM(); n++)
      {
        if (glob_mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
          nb_mesh_num_discont++;
        else
          nb_mesh_num_cont++;
      }
    
    if (glob_mesh_num.GetM() > 0)
      {
        NumDof.Reallocate(nb_dof + 2*nb_mesh_num_cont*nb_matching_proc + nb_matching_dofs + nb_elt*nb_mesh_num_discont);
        NumDofPML.Reallocate(nb_dof + nb_elt*nb_mesh_num_discont);
      }
    
    nb = 0 ; int nb_pml = 0;
    for (int n = 0; n < glob_mesh_num.GetM(); n++)
      {  
        if (glob_mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          {
            int nb_old = nb;
            for (int i = 0; i < nb_elt; i++)
              {
                int ne = NumElement(i);
                IVect Nodle = glob_mesh_num(n)->Element(ne).GetNodle();
                for (int j = 0; j < Nodle.GetM(); j++)
                  NumDof(nb+j) = Nodle(j);
                
                nb += Nodle.GetM();
              }
            
            for (int p = 0; p < MatchingDofs(n).GetM(); p++)
              if (ProcMatching(p))
                {
                  if (MatchingDofs(n)(p).GetM() > 0)
                    {
                      NumDof(nb) = p;
                      NumDof(nb+1) = MatchingDofs(n)(p).GetM();
                      for (int j = 0; j < MatchingDofs(n)(p).GetM(); j++)
                        NumDof(nb+2+j) = MatchingDofs(n)(p)(j);
                      
                      nb += 2+MatchingDofs(n)(p).GetM();
                    }
                  else
                    {
                      NumDof(nb) = p;
                      NumDof(nb+1) = 0;
                      nb += 2;
                    }
                }

            if (glob_mesh_num(n)->compute_dof_pml)
              {
                nb = nb_old;
                // we send dofs for PML if needed
                for (int i = 0; i < nb_elt; i++)
                  {
                    int ne = NumElement(i);
                    for (int j = 0; j < glob_mesh_num(n)->Element(ne).GetNbDof(); j++)
                      {
                        if (NumDof(nb+j) >= 0)
                          NumDofPML(nb_pml+j) = glob_mesh_num(n)->GetDofPML(NumDof(nb+j));
                        else
                          NumDofPML(nb_pml+j) = -1;
                      }
                    
                    nb_pml += glob_mesh_num(n)->Element(ne).GetNbDof();
                    nb += glob_mesh_num(n)->Element(ne).GetNbDof();
                  }
                
                for (int p = 0; p < MatchingDofs(n).GetM(); p++)
                  if (MatchingDofs(n)(p).GetM() > 0)
                    nb += 2+MatchingDofs(n)(p).GetM();
              }
          }
        else
          {
            int nb_old = nb;
            for (int i = 0; i < nb_elt; i++)
              {
                int ne = NumElement(i);
                NumDof(nb + i) = glob_mesh_num(n)->Element(ne).GetNumberDof(0);
              }
            
            nb += nb_elt;
            
            if (glob_mesh_num(n)->compute_dof_pml)
              {
                nb = nb_old;
                for (int i = 0; i < nb_elt; i++)
                  {
                    if (NumDof(nb + i) >= 0)
                      NumDofPML(nb_pml + i) = glob_mesh_num(n)->GetDofPML(NumDof(nb+i));
                    else
                      NumDofPML(nb_pml + i) = -1;
                  }
                
                nb_pml += nb_elt;
                nb += nb_elt;
              }
          }
      }
    
    if (rank != proc)
      {
        MPI_Send(NumDof.GetData(), nb, MPI_INTEGER, proc, tag, comm);    
        tag++;
      }
    
    if (rank != proc)
      {
        MPI_Send(NumDofPML.GetData(), nb_pml, MPI_INTEGER, proc, tag, comm);    
        tag++;
      }
    
    // and offsets for vectorial dofs
    if (OffsetDofV.GetM() > 0)
      {
        if (rank != proc)
          {
            MPI_Send(OffsetDofV.GetData(), OffsetDofV.GetM(), MPI_INTEGER, proc, tag, comm);
            tag++;
          }
      }
  }


  template<>
  void RecvMeshFromProcessor(Mesh<Dimension2>& mesh, Vector<MeshNumbering<Dimension2>* >& mesh_num,
			     Vector<IVect>& ConnecEdge, IVect& MatchingProc, Vector<IVect>& MatchingDofs,
                             ParamParallelMesh<Dimension2>& Param, int proc, MPI_Comm& comm, int tag)
  {
    IVect& nodl_mesh = Param.nodl_mesh;
    IVect& num = Param.num; VectReal_wp& all_param = Param.all_param;
    VectR2& PointsEdgeRef = Param.PointsEdgeRef; IVect& NumBoundary = Param.NumBoundary;
    IVect& NumDofNeighbor = Param.NumDofNeighbor; IVect& NumElem = Param.NumElem;
    IVect& NumDof = Param.NumDof; IVect& OffsetDofV = Param.OffsetDofV;
    VectR2& TranslatPeriodicEdge = Param.TranslatPeriodicEdge;
    IVect& NumDofPML = Param.NumDofPML;
    Vector<IVect>& offset_matching = Param.offset_matching;
    
    MPI_Status status;
    int rank; MPI_Comm_rank(comm, &rank);
    
    if (mesh.print_level >= 7)
      cout << rank << " Receiving numbers of vertices, edges, etc " << endl;
    
    // receiving number of vertices, edges and elements
    if (rank != proc)
      {
        int nb_param_int;
        MPI_Recv(&nb_param_int, 1, MPI_INTEGER, proc, tag, comm, &status); tag++;

        num.Reallocate(nb_param_int);        
        MPI_Recv(num.GetData(), nb_param_int, MPI_INTEGER, proc, tag, comm, &status);    
        tag++;
      }
    
    Vector<int64_t> all_param_tmp;
    int nb_vertices = num(0);
    int nb_edges_ref = num(1);
    int nb_elt = num(2);
    int nb_ref = num(3);
    int nb_param = num(4);
    int nb_int_elt = num(5);
    bool only_quad = num(15);
    mesh.SetGeometryOrder(num(6), only_quad);
    int nb_mesh_num = num(7);
    int nb_matching_proc = num(8);
    int nb_matching_dofs = num(9);
    int size_offset_v = num(11);
    int nb_pml_areas = num(14);
    bool variable_order = false;
    if (num(13) != MeshNumbering<Dimension2>::CONSTANT_ORDER)
      variable_order = true;

    int nb_glob_edges_ref = num(16);
    if (nb_mesh_num != mesh_num.GetM())
      {
        cout << "Number of mesh numbering does not coincide " << nb_mesh_num << " different from " << mesh_num.GetM() << endl;
        abort();
      }

    for (int n = 0; n < nb_mesh_num; n++)
      {
        mesh_num(n)->compute_dof_pml = num(12);
        mesh_num(n)->SetFormulationForPeriodicCondition(num(17));
        mesh_num(n)->SetOrder(num(20+n));
        mesh_num(n)->number_map.SetFormulationDG(num(20+n+nb_mesh_num));
      }
    
    if (nb_ref > mesh.GetNbReferences())
      mesh.ResizeNbReferences(nb_ref);
    
    if (mesh.print_level >= 7)
      cout << rank << " Receiving curves " << endl;
    
    // parameters of curves
    if (rank != proc)
      {
        all_param.Reallocate(nb_param);
        MpiRecv(comm, all_param, all_param_tmp, nb_param, proc, tag, status);
        tag++;
      }

    if (mesh.print_level >= 7)
      cout << rank << " Receiving curves completed" << endl;
    
    VectReal_wp param;
    IVect ref_neighbor(nb_ref+1); ref_neighbor.Fill(-3);
    nb_param = 5*nb_ref + 1;
    for (int i = 1; i <= nb_ref; i++)
      {
        int p0 = toInteger(all_param(i-1));
        int p1 = toInteger(all_param(i));
        
        mesh.SetBodyNumber(i, toInteger(all_param(nb_ref+i)));
        mesh.SetBoundaryCondition(i, toInteger(all_param(2*nb_ref+i)));
        mesh.SetCurveType(i, toInteger(all_param(3*nb_ref+i)));
        if (p1 > p0)
          {
            param.Reallocate(p1-p0);
            for (int p = p0; p < p1; p++)
              param(p-p0) = all_param(p);
            
            mesh.SetCurveParameter(i, param);
            nb_param += param.GetM();
          }        
        
        ref_neighbor(i) = toInteger(all_param(4*nb_ref+i));
      }
    
    Real_wp alpha = all_param(nb_param); nb_param++;
    mesh.SetPeriodicAlpha(alpha);

    mesh.ClearPeriodicCondition();
    int nb_periodic_ref = toInteger(all_param(nb_param)); nb_param++;
    for (int i = 0; i < nb_periodic_ref; i++)
      {
	mesh.AddPeriodicCondition(TinyVector<int, 2>(toInteger(all_param(nb_param)),
						     toInteger(all_param(nb_param+1))),
				  toInteger(all_param(nb_param+2)));

	nb_param += 3;
      }
    
    mesh.ReallocatePmlAreas(nb_pml_areas);
    for (int i = 0; i < nb_pml_areas; i++)
      mesh.GetPmlArea(i).SetRegion(all_param, nb_param);
    
    mesh.SetOriginalNeighborReference(ref_neighbor);
    
    // receiving vertices and internal points of edges
    mesh.ReallocateBoundariesRef(nb_edges_ref);
    mesh.ReallocateVertices(nb_vertices);
    int r = mesh.GetGeometryOrder();
    if (mesh.print_level >= 7)
      cout << rank << " Receiving internal points of the mesh " << endl;
    
    if (rank != proc)
      {
        PointsEdgeRef.Reallocate(nb_vertices + (r-1)*nb_edges_ref);
        MpiRecv(comm, reinterpret_cast<Real_wp*>(PointsEdgeRef.GetData()), all_param_tmp,
                2*PointsEdgeRef.GetM(), proc, tag, status);
        tag++;
      }
    
    for (int i = 0; i < nb_vertices; i++)
      mesh.Vertex(i) = PointsEdgeRef(i);
    
    int nb = nb_vertices;
    if (r > 1)
      {
        for (int i = 0; i < nb_edges_ref; i++)
          for (int k = 0; k < r-1; k++)
            mesh.SetPointInsideEdge(i, k, PointsEdgeRef(nb++));
      }
    
    // receiving edge numbers
    PointsEdgeRef.Clear();
    ConnecEdge.Reallocate(nb_edges_ref);
    if (mesh.print_level >= 7)
      cout << rank << " Receiving referenced edges " << endl;

    int nb_param_per_edge = 14 + 2*nb_mesh_num;
    if (rank != proc)
      {
        NumBoundary.Reallocate(nb_param_per_edge*nb_edges_ref);
        MPI_Recv(NumBoundary.GetData(), NumBoundary.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;
      }
    
    IVect OrderQuad(nb_edges_ref);
    nb = 0; int nb_dof_neighbor = 0;
    for (int n = 0; n < nb_mesh_num; n++)
      mesh_num(n)->InitPeriodicBoundary();

    IVect Global_to_local_edge(nb_glob_edges_ref);
    Global_to_local_edge.Fill(-1);
    for (int i = 0; i < nb_edges_ref; i++)
      Global_to_local_edge(NumBoundary(nb_param_per_edge*i + 9)) = i;
    
    for (int i = 0; i < nb_edges_ref; i++)
      {
        int n1 = NumBoundary(nb);
        int n2 = NumBoundary(nb+1);
        int ref = NumBoundary(nb+2);
        OrderQuad(i) = NumBoundary(nb+3);
        mesh.BoundaryRef(i).Init(n1, n2, ref);
	int nf_opp = NumBoundary(nb+12);
	int type_per = NumBoundary(nb+13);
	
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            int num_elt2 = NumBoundary(nb+4);
            int proc2 = NumBoundary(nb+5);
            int ne_loc2 = NumBoundary(nb+6);
            int num_loc2 = NumBoundary(nb+7);
            int rot = NumBoundary(nb+8);
            int type_elt = NumBoundary(nb+10);
            int iglob = NumBoundary(nb+9);
            int ref2 = NumBoundary(nb+11);

            int size_connec = 10 + nb_mesh_num;
            for (int n = 0; n < nb_mesh_num; n++)
              {
                //int re = NumBoundary(nb+14+n);
                int nb_dof = NumBoundary(nb+14+n+nb_mesh_num);                
                nb_dof_neighbor += nb_dof;
                size_connec += nb_dof + 3;
              }

            ConnecEdge(i).Reallocate(size_connec);
            ConnecEdge(i)(0) = iglob;
            ConnecEdge(i)(1) = num_elt2;
            ConnecEdge(i)(2) = proc2;
            ConnecEdge(i)(3) = ne_loc2;
            ConnecEdge(i)(4) = num_loc2;
            ConnecEdge(i)(5) = rot;
            ConnecEdge(i)(6) = type_elt;
            ConnecEdge(i)(7) = ref2;
            ConnecEdge(i)(8) = 0;
            ConnecEdge(i)(9) = 0;

            int offset_connec = 10+nb_mesh_num;
            for (int n = 0; n < nb_mesh_num; n++)
              {
                int re = NumBoundary(nb+14+n);
                int nb_dof = NumBoundary(nb+14+n+nb_mesh_num);
                ConnecEdge(i)(10+n) = offset_connec;
                ConnecEdge(i)(offset_connec) = re;
                ConnecEdge(i)(offset_connec+1) = nb_dof;
                ConnecEdge(i)(offset_connec+2) = 0;
                offset_connec += nb_dof + 3;
              }
          }
        else
	  {
	    ConnecEdge(i).Clear();
	    
	    // we find the local number for nf_opp
	    if (nf_opp >= 0)
	      nf_opp = Global_to_local_edge(nf_opp);
	  }
	
	if (nf_opp != -1)
          for (int n = 0; n < nb_mesh_num; n++)
            {
              mesh_num(n)->SetPeriodicBoundary(i, nf_opp);
              mesh_num(n)->SetPeriodicityTypeForBoundary(i, type_per);
            }
        
        nb += nb_param_per_edge;
      }

    Global_to_local_edge.Clear();
    
    // retrieving periodic edges
    if (mesh.print_level >= 7)
      cout << rank << " Receiving periodic edge translations " << endl;

    if (mesh.GetNbPeriodicReferences() > 0)
      {
	if (rank != proc)
	  {
	    TranslatPeriodicEdge.Reallocate(nb_edges_ref);
	    MpiRecv(comm, reinterpret_cast<Real_wp*>(TranslatPeriodicEdge.GetData()),
		    all_param_tmp, 2*nb_edges_ref, proc, tag, status);
	    
	    tag++;
	  }
	
	for (int n = 0; n < nb_mesh_num; n++)
          mesh_num(n)->GetTranslationPeriodicBoundary() = TranslatPeriodicEdge;
      }
    
    if (mesh.print_level >= 7)
      cout << rank << " Receiving dof number of neighboring elements " << endl;
    
    // receiving dof numbers of neighboring elements
    if (rank != proc)
      {
        if (mesh_num(0)->compute_dof_pml)
          NumDofNeighbor.Reallocate(2*nb_dof_neighbor);
        else
          NumDofNeighbor.Reallocate(nb_dof_neighbor);
        
        MPI_Recv(NumDofNeighbor.GetData(),
                 NumDofNeighbor.GetM(), MPI_INTEGER, proc, tag, 
                 comm, &status);
        tag++;
      }
    
    nb = 0;
    for (int i = 0; i < nb_edges_ref; i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            Vector<bool> element_in_pml(nb_mesh_num);
            element_in_pml.Fill(false); int nb_dof_all_pml = 0;
            for (int n = 0; n < nb_mesh_num; n++)
              if (mesh_num(n)->compute_dof_pml)
                {
                  int p = ConnecEdge(i)(10+n);
                  int nb_dof = ConnecEdge(i)(p+1);
                  for (int j = 0; j < nb_dof; j++)
                    if (NumDofNeighbor(nb + nb_dof + j) >= 0)
                      element_in_pml(n) = true;

                  if (element_in_pml(n))
                    nb_dof_all_pml += nb_dof;
                }
            
            if (nb_dof_all_pml <= 0)
              {
                for (int n = 0; n < nb_mesh_num; n++)
                  {
                    int p = ConnecEdge(i)(10+n);
                    int nb_dof = ConnecEdge(i)(p+1);                    
                    for (int j = 0; j < nb_dof; j++)
                      ConnecEdge(i)(p+3+j) =  NumDofNeighbor(nb+j);
            
                    nb += nb_dof;
                    if (mesh_num(n)->compute_dof_pml)
                      nb += nb_dof;
                  }
              }
            else
              {
                Vector<int> old_connec = ConnecEdge(i);
                ConnecEdge(i).Resize(old_connec.GetM() + nb_dof_all_pml);
                int offset_connec = 10 + nb_mesh_num, old_off = offset_connec;
                for (int n = 0; n < nb_mesh_num; n++)
                  {
                    int nb_dof = old_connec(old_off+1);
                    ConnecEdge(i)(10+n) = offset_connec;
                    ConnecEdge(i)(offset_connec) = old_connec(old_off);
                    ConnecEdge(i)(offset_connec+1) = nb_dof;
                    for (int j = 0; j < nb_dof; j++)
                      ConnecEdge(i)(offset_connec+3+j) =  NumDofNeighbor(nb+j);

                    nb += nb_dof;
                    if (element_in_pml(n))
                      {
                        ConnecEdge(i)(offset_connec+2) = nb_dof;
                        offset_connec += nb_dof;
                        for (int j = 0; j < nb_dof; j++)
                          ConnecEdge(i)(offset_connec+3+j) =  NumDofNeighbor(nb+j);
                      }
                    else
                      ConnecEdge(i)(offset_connec+2) = 0;

                    if (mesh_num(n)->compute_dof_pml)
                      nb += nb_dof;

                    offset_connec += nb_dof+3;
                    old_off += nb_dof+3;
                  }
              }
          }
      }
    
    // receiving elements
    if (mesh.print_level >= 7)
      cout << rank << " Receiving elements " << endl;
    
    mesh.ReallocateElements(nb_elt);
    if (rank != proc)
      {
        NumElem.Reallocate((5 + 3*nb_mesh_num)*nb_elt + (3+nb_mesh_num)*nb_int_elt);
        MPI_Recv(NumElem.GetData(), NumElem.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;
      }
    
    nb = 0;
    Vector<IVect> OrderElt(nb_mesh_num), OrderInside(nb_mesh_num), NbDofElt(nb_mesh_num);
    for (int n = 0; n < nb_mesh_num; n++)
      {
        OrderElt(n).Reallocate(nb_elt); OrderInside(n).Reallocate(nb_elt);
        NbDofElt(n).Reallocate(nb_elt);
      }
    
    int offset = (5+3*nb_mesh_num)*nb_elt, nb_dof = 0;
    mesh.GlobElementNumber_Subdomain.Reallocate(nb_elt);
    mesh.GlobVertexNumber_Subdomain.Reallocate(nb_vertices);
    for (int i = 0; i < nb_elt; i++)
      {
        int nb_vert = NumElem(nb);
        int ref = NumElem(nb+1);
        int pml_num = NumElem(nb+2);
        int pml_type = NumElem(nb+3);
        num.Reallocate(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = NumElem(offset+j);
            mesh.GlobVertexNumber_Subdomain(num(j)) = NumElem(offset + nb_vert + j);
          }
	
        mesh.Element(i).Init(num, ref);
	mesh.Element(i).SetPML(pml_num, pml_type);

        mesh.GlobElementNumber_Subdomain(i) = NumElem(nb+4);

        for (int n = 0; n < nb_mesh_num; n++)
          {
            OrderElt(n)(i) = NumElem(nb+5+n);
            OrderInside(n)(i) = NumElem(nb+5+n+nb_mesh_num);
            NbDofElt(n)(i) = NumElem(nb+5+n+2*nb_mesh_num);
            if (mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
              nb_dof += NbDofElt(n)(i);
          }
        
        nb += 5 + 3*nb_mesh_num;
        offset += (3+nb_mesh_num)*nb_vert;
      }
    
    // updating connectivity
    mesh.ReorientElements();
    mesh.FindConnectivity();
    mesh.ProjectPointsOnCurves();
    
    // taking into account periodicity
    if (nb_mesh_num > 0)
      for (int i = 0; i < nb_edges_ref; i++)
        {
          int ref = NumBoundary(nb_param_per_edge*i + 2);
          
          if (mesh.GetBoundaryCondition(ref) != BoundaryConditionEnum::LINE_NEIGHBOR)
            if (mesh_num(0)->GetPeriodicBoundary(i) >= 0)
              {
                int i2 = mesh_num(0)->GetPeriodicBoundary(i);
                int ne2 = mesh.BoundaryRef(i2).numElement(0);
                mesh.BoundaryRef(i).AddElement(ne2);
                mesh.Boundary(i).AddElement(ne2);
              }	  
        }

    if (variable_order)
      {
        for (int n = 0; n < nb_mesh_num; n++)
          {
            mesh_num(n)->SetVariableOrder(mesh_num(n)->USER_ORDER);
        
            for (int i = 0; i < nb_elt; i++)
              {
                mesh_num(n)->SetOrderElement(i, OrderElt(n)(i));
                mesh_num(n)->SetOrderInside(i, OrderInside(n)(i));
              }
          }

        nb = (5+3*nb_mesh_num)*nb_elt;
        for (int i = 0; i < nb_elt; i++)
          {
            int nb_vert = mesh.Element(i).GetNbEdges();            
            for (int j = 0; j < nb_vert; j++)
              {
                int ne = mesh.Element(i).numEdge(j);
                for (int n = 0; n < nb_mesh_num; n++)
                  mesh_num(n)->SetOrderEdge(ne, NumElem(nb + (3+n)*nb_vert + j));
              }
            
            nb += (3+nb_mesh_num)*nb_vert;
          }
        
        // order for referenced boundaries
        for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
          for (int n = 0; n < nb_mesh_num; n++)
            mesh_num(n)->SetOrderQuadrature(i, OrderQuad(i));
        
        // order for internal boundaries
        for (int n = 0; n < nb_mesh_num; n++)
          for (int i = mesh.GetNbBoundaryRef(); i < mesh.GetNbBoundary(); i++)
            {
              int rf = 0;
              for (int k = 0; k < mesh.Boundary(i).GetNbElements(); k++)
                rf = max(rf, mesh_num(n)->GetOrderElement(mesh.Boundary(i).numElement(k)));
              
              mesh_num(n)->SetOrderQuadrature(i, rf);
            }        
      }
    
    nb = (5+nb_mesh_num)*nb_elt;
    mesh.GlobEdgeNumber_Subdomain.Reallocate(mesh.GetNbEdges());
    mesh.GlobEdgeNumber_Subdomain.Fill(-1);
    for (int i = 0; i < nb_elt; i++)
      {
        int nb_vert = mesh.Element(i).GetNbEdges();
        for (int j = 0; j < nb_vert; j++)
          {
            int ne = mesh.Element(i).numEdge(j);
            mesh.GlobEdgeNumber_Subdomain(ne) = NumElem(nb + 2*nb_vert + j);
          }
        
        nb += (3+nb_mesh_num)*nb_vert;
      }
    
    // retrieving dof numbers
    if (mesh.print_level >= 7)
      cout << rank << " Receiving dof numbers " << endl;
    
    int nb_mesh_num_discont = 0, nb_mesh_num_cont = 0;
    for (int n = 0; n < mesh_num.GetM(); n++)
      {
        if (mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
          nb_mesh_num_discont++;
        else
          nb_mesh_num_cont++;
      }
    
    if (mesh_num.GetM() > 0)
      {
        NumDof.Reallocate(nb_dof + 2*nb_mesh_num_cont*nb_matching_proc + nb_matching_dofs + nb_elt*nb_mesh_num_discont);
        NumDofPML.Reallocate(nb_dof + nb_elt*nb_mesh_num_discont);
      }
    
    if (rank != proc)
      {
        MPI_Recv(NumDof.GetData(), NumDof.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;

        MPI_Recv(NumDofPML.GetData(), NumDofPML.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;
      }

    MatchingDofs.Reallocate(nb_matching_proc);
    MatchingProc.Reallocate(nb_matching_proc);
    // counting size of arrays and allocating them
    nb = 0; IVect size_matching(nb_matching_proc);
    size_matching.Reallocate(nb_matching_proc);
    size_matching.Zero();
    for (int n = 0; n < nb_mesh_num; n++)
      {
        if (mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          {
            for (int i = 0; i < nb_elt; i++)
              nb += NbDofElt(n)(i);
            
            for (int p = 0; p < nb_matching_proc; p++)
              {
                int n_dof = NumDof(nb+1);
                MatchingProc(p) = NumDof(nb);
                size_matching(p) += n_dof;
                nb += 2 + n_dof;
              }
          }
        else
          nb += nb_elt;
      }
    
    offset_matching.Reallocate(nb_matching_proc);
    for (int p = 0; p < nb_matching_proc; p++)
      {
        MatchingDofs(p).Reallocate(size_matching(p));
        MatchingDofs(p).Fill(-1);
        
        offset_matching(p).Reallocate(nb_mesh_num+1);
        offset_matching(p).Zero();
      }
    
    size_matching.Zero();
    nb = 0; int nb_pml = 0;
    for (int n = 0; n < nb_mesh_num; n++)
      {
        if (mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
          {
            nb_dof = 0;
            for (int i = 0; i < nb_elt; i++)
              nb_dof += NbDofElt(n)(i);
            
            mesh_num(n)->GlobDofNumber_Subdomain.Reallocate(nb_dof);
            mesh_num(n)->GlobDofNumber_Subdomain.Fill(-1);
            int offset = 0;
            mesh_num(n)->ReallocateElements(nb_elt);
            for (int i = 0; i < nb_elt; i++)
              {
                int off_g = NumDof(nb + i);
                mesh_num(n)->Element(i).ReallocateDof(1);
                mesh_num(n)->Element(i).SetNumberDof(0, offset);
                for (int j = 0; j < NbDofElt(n)(i); j++)
                  mesh_num(n)->GlobDofNumber_Subdomain(offset+j) = off_g + j;
                
                offset += NbDofElt(n)(i);
              }            
            
            mesh_num(n)->SetNbDof(offset);
            nb += nb_elt;
            
            if (mesh_num(n)->compute_dof_pml)
              {
                IVect IndexDof(nb_elt);
                IndexDof.Fill(-1);
                
                nb_dof = 0;
                for (int i = 0; i < nb_elt; i++)
                  if (NumDofPML(nb_pml + i) >= 0)
                    {
                      IndexDof(i) = nb_dof;
                      nb_dof += NbDofElt(n)(i);
                    }
                
                mesh_num(n)->GlobDofPML_Subdomain.Reallocate(nb_dof);
                mesh_num(n)->GlobDofPML_Subdomain.Fill(-1);
                mesh_num(n)->ReallocateDofPML(nb_dof);
                for (int i = 0; i < nb_elt; i++)
                  if (NumDofPML(nb_pml+i) >= 0)
                    {
                      int npml = NumDofPML(nb_pml+i);
                      int off = mesh_num(n)->Element(i).GetNumberDof(0);
                      for (int j = 0; j < NbDofElt(n)(i); j++)
                        {
                          mesh_num(n)->GlobDofPML_Subdomain(IndexDof(i) + j) = npml + j;
                          mesh_num(n)->SetDofPML(off+j, IndexDof(i)+j);
                        }
                    }

                nb_pml += nb_elt;
              }
	  }
        else
          {
            IVect IndexDof(nodl_mesh(n));
            IndexDof.Fill(-1);
            
            mesh_num(n)->ReallocateElements(nb_elt);
            for (int i = 0; i < nb_elt; i++)
              {
                mesh_num(n)->Element(i).ReallocateDof(NbDofElt(n)(i));            
                for (int j = 0; j < NbDofElt(n)(i); j++)
                  {
                    int n0 = NumDof(nb+j);
                    mesh_num(n)->Element(i).SetNumberDof(j, n0);
                    if (n0 >= 0)
                      IndexDof(n0) = 1;
                  }
                
                nb += NbDofElt(n)(i);
              }
            
            for (int p = 0; p < nb_matching_proc; p++)
              {
                int n_dof = NumDof(nb+1);
                for (int j = 0; j < n_dof; j++)
                  MatchingDofs(p)(size_matching(p) + j) = NumDof(nb+2+j);   

                size_matching(p) += n_dof;
                offset_matching(p)(n+1) = size_matching(p);
                nb += 2 + n_dof;
              }
            
            int nb_dof_true = 0;
            for (int i = 0; i < IndexDof.GetM(); i++)
              if (IndexDof(i) > 0)
                IndexDof(i) = nb_dof_true++;
            
            mesh_num(n)->SetNbDof(nb_dof_true);
            mesh_num(n)->GlobDofNumber_Subdomain.Reallocate(nb_dof_true);
            mesh_num(n)->GlobDofNumber_Subdomain.Fill(-1);
            for (int i = 0; i < nb_elt; i++)
              for (int j = 0; j < NbDofElt(n)(i); j++)
                {
                  int n0 = mesh_num(n)->Element(i).GetNumberDof(j);
                  if (n0 >= 0)
                    {
                      mesh_num(n)->GlobDofNumber_Subdomain(IndexDof(n0)) = n0;
                      mesh_num(n)->Element(i).SetNumberDof(j, IndexDof(n0));
                    }
                }

            // PML dofs
            int nb_pml_old = nb_pml;
            if (mesh_num(n)->compute_dof_pml)
              {
                IndexDof.Fill(-1);                
                for (int i = 0; i < nb_elt; i++)
                  {
                    for (int j = 0; j < NbDofElt(n)(i); j++)
                      {
                        int n0 = NumDofPML(nb_pml + j);
                        if (n0 >= 0)
                          IndexDof(n0) = 1;
                      }
                    
                    nb_pml += NbDofElt(n)(i);
                  }
                
                nb_dof = 0;
                for (int i = 0; i < IndexDof.GetM(); i++)
                  if (IndexDof(i) > 0)
                    IndexDof(i) = nb_dof++;
                
                mesh_num(n)->GlobDofPML_Subdomain.Reallocate(nb_dof);
                mesh_num(n)->GlobDofPML_Subdomain.Fill(-1);
                mesh_num(n)->ReallocateDofPML(nb_dof);
                nb_pml = nb_pml_old;
                for (int i = 0; i < nb_elt; i++)
                  {                
                    for (int j = 0; j < NbDofElt(n)(i); j++)
                      {
                        int n0 = mesh_num(n)->Element(i).GetNumberDof(j);
                        int npml = NumDofPML(nb_pml+j);
                        if (npml >= 0)
                          {
                            mesh_num(n)->GlobDofPML_Subdomain(IndexDof(npml)) = npml;
                            mesh_num(n)->SetDofPML(n0, IndexDof(npml));
                          }
                      }                
                    nb_pml += NbDofElt(n)(i);
                  }
              }
          }
      }
    
    if (mesh.print_level >= 7)
      cout << rank << " Receiving offsets for vectorial dofs " << endl;

    // receiving offsets for vectorial dofs
    if (size_offset_v > 0)
      {
        if (rank != proc)
          {
            OffsetDofV.Reallocate(nb_elt);
            MPI_Recv(OffsetDofV.GetData(), nb_elt, MPI_INTEGER, proc, tag,
                     comm, &status);
            tag++;
          }
      }
  }
  
   
#ifdef MONTJOIE_WITH_THREE_DIM
  //! sending a mesh to a processor
  template<>
  void SendMeshToProcessor(const Mesh<Dimension3>& glob_mesh,
                           const Vector<MeshNumbering<Dimension3>* >& glob_mesh_num,
			   const IVect& NumElement, const IVect& Epart, const IVect& NumLoc,
                           Mesh<Dimension3>& mesh, int proc, const MPI_Comm& comm,
                           ParamParallelMesh<Dimension3>& Param, int tag)
  {
    IVect& num = Param.num; VectReal_wp& all_param = Param.all_param;
    VectR3& PointsEdgeRef = Param.PointsEdgeRef; IVect& NumBoundary = Param.NumBoundary;
    IVect& NumDofNeighbor = Param.NumDofNeighbor; IVect& NumElem = Param.NumElem;
    IVect& NumDof = Param.NumDof; Vector<IVect>& MinimalProc = Param.MinimalProc;
    VectR3& TranslatPeriodicFace = Param.TranslatPeriodicEdge;
    IVect& OffsetDofV = Param.OffsetDofV; IVect& NumDofPML = Param.NumDofPML;
    
    int rank; MPI_Comm_rank(comm, &rank);
    
    int nb_ref = mesh.GetNbReferences();
    int nb_param = 5*nb_ref+3;
    for (int i = 0; i < glob_mesh.GetNbPmlAreas(); i++)
      nb_param += glob_mesh.GetPmlArea(i).GetNbParameters();
    
    VectReal_wp param;
    for (int i = 1; i <= nb_ref; i++)
      {
        mesh.GetCurveParameter(i, param);
        nb_param += param.GetM();
      }

    nb_param += 3*glob_mesh.GetNbPeriodicReferences();
	
    Vector<int64_t> all_param_tmp;    
    IVect ref_neighbor = mesh.GetOriginalNeighborReference();
    all_param.Reallocate(nb_param); all_param.Fill(0);    
    nb_param = 5*nb_ref+1;
    all_param(0) = nb_param;
    for (int i = 1; i <= mesh.GetNbReferences(); i++)
      {
        mesh.GetCurveParameter(i, param);
        for (int j = 0; j < param.GetM(); j++)
          all_param(nb_param+j) = param(j);
        
        nb_param += param.GetM();
        all_param(i) = nb_param;
        all_param(nb_ref+i) = mesh.GetBodyNumber(i);
        all_param(2*nb_ref+i) = mesh.GetBoundaryCondition(i);
        all_param(3*nb_ref+i) = mesh.GetCurveType(i);
        all_param(4*nb_ref+i) = ref_neighbor(i);
      }    
    
    all_param(nb_param) = glob_mesh.GetPeriodicAlpha(); nb_param++;
    all_param(nb_param) = glob_mesh.GetNbPeriodicReferences(); nb_param++;
    for (int i = 0; i < glob_mesh.GetNbPeriodicReferences(); i++)
      {
	all_param(nb_param) = glob_mesh.GetPeriodicReference(i, 0);
	all_param(nb_param+1) = glob_mesh.GetPeriodicReference(i, 1);
	all_param(nb_param+2) = glob_mesh.GetPeriodicityTypeReference(i);
	nb_param += 3;
      }

    for (int i = 0; i < glob_mesh.GetNbPmlAreas(); i++)
      glob_mesh.GetPmlArea(i).FillParameters(all_param, nb_param);

    // sending number of vertices, referenced edges and elements
    int nb_vertices = mesh.GetNbVertices();
    int nb_edges = mesh.GetNbEdges();
    int nb_edges_ref = mesh.GetNbEdgesRef();
    int nb_faces_ref = mesh.GetNbBoundaryRef();
    int nb_tri_ref = 0, nb_quad_ref = 0;
    for (int i = 0; i < nb_faces_ref; i++)
      {
	if (mesh.BoundaryRef(i).GetNbVertices() == 3)
	  nb_tri_ref++;
	else
	  nb_quad_ref++;
      }
    
    int nb_int_face = 3*nb_tri_ref + 4*nb_quad_ref;
    
    int nb_elt = mesh.GetNbElt();
    int nb_intE_elt = 0, nb_intV_elt = 0, nb_intF_elt = 0;
    for (int i = 0; i < nb_elt; i++)
      {
        nb_intV_elt += mesh.Element(i).GetNbVertices();
        nb_intE_elt += mesh.Element(i).GetNbEdges();
        nb_intF_elt += mesh.Element(i).GetNbFaces();
      }
        
    // counting dofs
    int nb_mesh_num = glob_mesh_num.GetM();
    int nb_matching_proc = 0, nb_matching_dofs = 0;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    Vector<bool> ProcMatching(nb_proc);
    ProcMatching.Fill(false);
    
    Vector<Vector<IVect> > MatchingDofs(nb_mesh_num);
    for (int i = 0; i < nb_mesh_num; i++)
      MatchingDofs(i).Reallocate(nb_proc);

    for (int n = 0; n < glob_mesh_num.GetM(); n++)
      if (glob_mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
        {
          // returning to positive numbers
          for (int i = 0; i < MinimalProc(n).GetM(); i++)
            MinimalProc(n)(i) = abs(MinimalProc(n)(i));
          
          // we put opposite signs in MinimalProc for dofs associated with the processor proc            
          for (int i = 0; i < nb_elt; i++)
            {
              int ne = NumElement(i);
              for (int j = 0; j < glob_mesh_num(n)->Element(ne).GetNbDof(); j++)
                {
                  int num_dof = glob_mesh_num(n)->Element(ne).GetNumberDof(j);
                  if ((num_dof >= 0) && (MinimalProc(n)(num_dof) > 0))
                    MinimalProc(n)(num_dof) = -MinimalProc(n)(num_dof);
                }
            }
          
          // the number of dofs shared with other processors is estimated
          IVect nb_dofs_proc(nb_proc);
          nb_dofs_proc.Zero();
          for (int i = 0; i < glob_mesh.GetNbElt(); i++)
            for (int j = 0; j < glob_mesh_num(n)->Element(i).GetNbDof(); j++)
              {
                int num_dof = glob_mesh_num(n)->Element(i).GetNumberDof(j);
                if ((num_dof >= 0) && (MinimalProc(n)(num_dof) < 0))
                  nb_dofs_proc(Epart(i))++;
              }
          
          for (int p = 0; p < MatchingDofs(n).GetM(); p++)
            if (p != proc)
              MatchingDofs(n)(p).Reallocate(nb_dofs_proc(p));
          
          nb_dofs_proc.Zero();
          for (int i = 0; i < glob_mesh.GetNbElt(); i++)
            {
              int p = Epart(i);
              if (p != proc)
                for (int j = 0; j < glob_mesh_num(n)->Element(i).GetNbDof(); j++)
                  {
                    int num_dof = glob_mesh_num(n)->Element(i).GetNumberDof(j);
                    if ((num_dof >= 0) && (MinimalProc(n)(num_dof) < 0))
                      {
                        MatchingDofs(n)(p)(nb_dofs_proc(p)) = num_dof;
                        nb_dofs_proc(p)++;
                      }
                  }
            }
          
          // duplicates are removed and numbers sorted
          for (int p = 0; p < MatchingDofs(n).GetM(); p++)
            if (MatchingDofs(n)(p).GetM() > 0)
              {
                if (!ProcMatching(p))
                  {
                    ProcMatching(p) = true;
                    nb_matching_proc++;
                  }
                
                RemoveDuplicate(MatchingDofs(n)(p));
                nb_matching_dofs += MatchingDofs(n)(p).GetM();
              }          
        }
    
    num.Reallocate(30 + 2*nb_mesh_num);
    num.Fill(0);
    
    num(0) = nb_vertices;
    num(1) = nb_edges_ref;
    num(2) = nb_faces_ref;
    num(3) = nb_tri_ref;
    num(4) = nb_quad_ref;
    
    num(5) = nb_elt;
    num(6) = nb_intV_elt;    
    num(7) = nb_intE_elt;    
    num(8) = nb_intF_elt; 
    num(9) = nb_ref;
    num(10) = nb_param;
    
    num(11) = mesh.GetGeometryOrder();
    num(12) = nb_mesh_num;
    num(13) = nb_matching_proc;
    num(14) = nb_matching_dofs;
    num(15) = 0;
    num(16) = glob_mesh.GetNbBoundaryRef();
    num(17) = nb_edges;
    num(19) = OffsetDofV.GetM();
    if (glob_mesh_num.GetM() > 0)
      {
        num(20) = glob_mesh_num(0)->compute_dof_pml;
        num(21) = glob_mesh_num(0)->GetVariableOrder();
        num(24) = glob_mesh_num(0)->GetFormulationForPeriodicCondition();
      }
    
    num(22) = glob_mesh.GetNbPmlAreas();
    num(23) = glob_mesh.IsOnlyHexahedral();

    for (int n = 0; n < nb_mesh_num; n++)
      {
        num(30 + n) = glob_mesh_num(n)->GetOrder();
        num(30 + n + nb_mesh_num) = glob_mesh_num(n)->number_map.FormulationDG();
      }

    if (rank != proc)
      {
        int nb_param_int = num.GetM();
        MPI_Send(&nb_param_int, 1, MPI_INTEGER, proc, tag, comm); tag++;
        MPI_Send(num.GetData(), nb_param_int, MPI_INTEGER, proc, tag, comm); tag++;
      }
    
    // sending parameters of curves
    if (rank != proc)
      {
        MpiSsend(comm, all_param, all_param_tmp, nb_param, proc, tag);
        tag++;
      }
    
    // sending vertices and internal points of edges and faces
    int r = mesh.GetGeometryOrder();
    PointsEdgeRef.Reallocate(nb_vertices + (r-1)*nb_edges_ref
                             + (r-1)*(r-1)*nb_quad_ref + (r-1)*(r-2)/2*nb_tri_ref);
    
    for (int i = 0; i < nb_vertices; i++)
      PointsEdgeRef(i) = mesh.Vertex(i);
    
    int nb = nb_vertices;
    if (r > 1)
      {
        for (int i = 0; i < nb_edges_ref; i++)
          for (int k = 0; k < r-1; k++)
            PointsEdgeRef(nb++) = mesh.GetPointInsideEdge(i, k);
	
	int nb_nodes = 0, nb_nodes_tri = (r-1)*(r-2)/2, nb_nodes_quad = (r-1)*(r-1);
	for (int i = 0; i < nb_faces_ref; i++)
	  {
	    if (mesh.BoundaryRef(i).GetNbVertices() == 3)
	      nb_nodes = nb_nodes_tri;
	    else
	      nb_nodes = nb_nodes_quad;
	    
	    for (int k = 0; k < nb_nodes; k++)
	      PointsEdgeRef(nb++) = mesh.GetPointInsideFace(i, k);
	  }	    
      }        
    
    if (rank != proc)
      {
        MpiSsend(comm, reinterpret_cast<Real_wp*>(PointsEdgeRef.GetData()),
                 all_param_tmp, 3*nb, proc, tag);
        tag++;
      }
    
    // faces
    int nb_param_per_face = 14 + 2*nb_mesh_num;
    NumBoundary.Reallocate(nb_param_per_face*nb_faces_ref + nb_int_face + 3*nb_edges_ref);
    NumBoundary.Fill(-1);
    nb = 0;
    int nb_dof_neighbor = 0;
    for (int i = 0; i < nb_faces_ref; i++)
      {
	NumBoundary(nb) = mesh.BoundaryRef(i).GetNbVertices();
	NumBoundary(nb+1) = mesh.BoundaryRef(i).GetReference();

        int ne = mesh.BoundaryRef(i).numElement(0);
        int num_face = i;
        int num_loc = -1;
        int num_elem = -1;
        if (ne >= 0)
          {
            num_loc = mesh.Element(ne).GetPositionBoundary(num_face);
            num_elem = NumElement(ne);
          }
        
        int iglob = 0;
        if (num_elem >= 0)
          iglob = glob_mesh.Element(num_elem).numBoundary(num_loc);

        if (glob_mesh_num.GetM() > 0)
          NumBoundary(nb+2) = glob_mesh_num(0)->GetOrderQuadrature(iglob);
        
	NumBoundary(nb+9) = iglob;
        if (glob_mesh_num.GetM() > 0)
          {
            NumBoundary(nb+12) = glob_mesh_num(0)->GetPeriodicBoundary(iglob);
            NumBoundary(nb+13) = glob_mesh_num(0)->GetPeriodicityTypeForBoundary(iglob);
          }
	
        int ref = mesh.BoundaryRef(i).GetReference();
        
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            // referenced edge at the interface between subdomains
            int num_elt2 = glob_mesh.Boundary(iglob).numElement(0);            
            if (num_elt2 == num_elem)
              num_elt2 = glob_mesh.Boundary(iglob).numElement(1);
            
            //int num_loc2 = glob_mesh.Element(num_elt2).numEdge(num_loc);
            int num_loc2 = -1, rot = -1;
            if (glob_mesh_num.GetM() > 0)
              glob_mesh_num(0)->GetBoundaryRotation(iglob, num_elem, num_elt2, num_loc, num_loc2, rot);
            
            NumBoundary(nb+4) = num_elt2;
            NumBoundary(nb+5) = Epart(num_elt2);
            NumBoundary(nb+6) = NumLoc(num_elt2);
            NumBoundary(nb+7) = num_loc2;
            NumBoundary(nb+8) = rot;
            for (int n = 0; n < glob_mesh_num.GetM(); n++)
              {                
                NumBoundary(nb+14+n) = glob_mesh_num(n)->GetOrderElement(num_elt2);
                NumBoundary(nb+14+nb_mesh_num+n) = glob_mesh_num(n)->GetNbLocalDof(num_elt2);
                nb_dof_neighbor += glob_mesh_num(n)->GetNbLocalDof(num_elt2);
              }
            
            NumBoundary(nb+10) = glob_mesh.GetTypeElement(num_elt2);
            NumBoundary(nb+11) = glob_mesh.Element(num_elt2).GetReference();
          }
        
        nb += nb_param_per_face;
      }
    
    for (int i = 0; i < nb_faces_ref; i++)
      {
        for (int j = 0; j < mesh.BoundaryRef(i).GetNbVertices(); j++)
          NumBoundary(nb+j) = mesh.BoundaryRef(i).numVertex(j); 
	
        nb += mesh.BoundaryRef(i).GetNbVertices();
      }

    for (int i = 0; i < nb_edges_ref; i++)
      {
        NumBoundary(nb) = mesh.EdgeRef(i).numVertex(0);
        NumBoundary(nb+1) = mesh.EdgeRef(i).numVertex(1);
        NumBoundary(nb+2) = mesh.EdgeRef(i).GetReference();
	nb += 3;
      }
    
    if (rank != proc)
      {
        MPI_Send(NumBoundary.GetData(), nb, MPI_INTEGER, proc, tag, comm);
        tag++;
      }

    // and periodic faces
    if ((glob_mesh.GetNbPeriodicReferences() > 0) && (glob_mesh_num.GetM() > 0))
      {
	TranslatPeriodicFace.Reallocate(nb_faces_ref);
	for (int i = 0; i < nb_faces_ref; i++)
	  {
            int num = NumBoundary(nb_param_per_face*i+12);
            if (num >= 0)
              TranslatPeriodicFace(i)
                = glob_mesh_num(0)->GetTranslationPeriodicBoundary(num);
            else
              TranslatPeriodicFace(i).Zero();
          }
        
	if (rank != proc)
	  {
	    MpiSsend(comm, reinterpret_cast<Real_wp*>(TranslatPeriodicFace.GetData()), all_param_tmp,
		     3*nb_faces_ref, proc, tag);
	    tag++;
	  }
      }

    // dofs associated with neighboring elements
    if (glob_mesh_num.GetM() > 0)
      {
        if (glob_mesh_num(0)->compute_dof_pml)
          NumDofNeighbor.Reallocate(2*nb_dof_neighbor);
        else
          NumDofNeighbor.Reallocate(nb_dof_neighbor);
        
        nb = 0;
        for (int i = 0; i < nb_faces_ref; i++)
          {
            int ref = mesh.BoundaryRef(i).GetReference();
            if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
              for (int n = 0; n < glob_mesh_num.GetM(); n++)
                {
                  int num_elt2 = NumBoundary(nb_param_per_face*i + 4);
                  IVect Nodle = glob_mesh_num(n)->Element(num_elt2).GetNodle();
                  if (glob_mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
                    {
                      Nodle.Resize(glob_mesh_num(n)->GetNbLocalDof(num_elt2));
                      for (int j = 1; j < Nodle.GetM(); j++)
                        Nodle(j) = Nodle(0) + j;
                    }
                  
                  for (int j = 0; j < Nodle.GetM(); j++)
                    NumDofNeighbor(nb+j) = Nodle(j);
            
                  nb += Nodle.GetM();
                  
                  if (glob_mesh_num(n)->compute_dof_pml)
                    {
                      for (int j = 0; j < Nodle.GetM(); j++)
                        {
                          if (Nodle(j) >= 0)
                            NumDofNeighbor(nb+j) = glob_mesh_num(n)->GetDofPML(Nodle(j));
                          else
                            NumDofNeighbor(nb+j) = -1;
                        }
                      
                      nb += Nodle.GetM();
                    }
                }
          }
        
        if (rank != proc)
          {
            MPI_Send(NumDofNeighbor.GetData(), nb, MPI_INTEGER, proc, tag, comm);
            tag++;
          }
      }
    
    // and elements
    NumElem.Reallocate((5 + 3*nb_mesh_num)*nb_elt + 2*nb_intV_elt + (1+nb_mesh_num)*nb_intE_elt + (1+nb_mesh_num)*nb_intF_elt);
    nb = 0; int nb_dof = 0;
    for (int i = 0; i < nb_elt; i++)
      {
        NumElem(nb) = mesh.Element(i).GetNbVertices();
        NumElem(nb+1) = mesh.Element(i).GetReference();
        NumElem(nb+2) = mesh.Element(i).GetNumberPML();
        NumElem(nb+3) = mesh.Element(i).GetTypePML();
        int ne = NumElement(i);
        NumElem(nb+4) = ne;
        for (int n = 0; n < glob_mesh_num.GetM(); n++)
          {
            NumElem(nb+5+n) = glob_mesh_num(n)->GetOrderElement(ne);
            NumElem(nb+5+n+nb_mesh_num) = glob_mesh_num(n)->GetOrderInside(ne);   
            NumElem(nb+5+n+2*nb_mesh_num) = glob_mesh_num(n)->GetNbLocalDof(ne);
            if (glob_mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
              nb_dof += NumElem(nb+5+n+2*nb_mesh_num);
          }
        
        nb += 5 + 3*nb_mesh_num;
      }
    
    if (mesh.print_level >= 7)
      cout << " we send " << nb_dof << " degrees of freedom " << endl;
    
    for (int i = 0; i < nb_elt; i++)
      {
        int ne = NumElement(i);
        int nv = mesh.Element(i).GetNbVertices();
        for (int j = 0; j < nv; j++)
          {
            NumElem(nb+j) = mesh.Element(i).numVertex(j); 
            NumElem(nb+nv+j) = glob_mesh.Element(ne).numVertex(j);
          }
        
        nb += 2*nv;
        
        nv = mesh.Element(i).GetNbEdges();
        for (int j = 0; j < nv; j++)
          {
            int num_edge = glob_mesh.Element(ne).numEdge(j);
            NumElem(nb+j) = num_edge;
            for (int n = 0; n < glob_mesh_num.GetM(); n++)
              NumElem(nb+(1+n)*nv+j) = glob_mesh_num(n)->GetOrderEdge(num_edge);
          }
        
        nb += (1+nb_mesh_num)*nv;
        
        nv = mesh.Element(i).GetNbFaces();
        for (int j = 0; j < nv; j++)
          {
            int num_face = glob_mesh.Element(ne).numFace(j);
            NumElem(nb+j) = num_face;
            for (int n = 0; n < glob_mesh_num.GetM(); n++)
              NumElem(nb+(1+n)*nv+j) = glob_mesh_num(n)->GetOrderFace(num_face);
          }
        
        nb += (1+nb_mesh_num)*nv;
      }
    
    if (rank != proc)
      {
        MPI_Send(NumElem.GetData(), nb, MPI_INTEGER, proc, tag, comm);
        tag++;
      }
    
    // then dofs
    int nb_mesh_num_discont = 0, nb_mesh_num_cont = 0;
    for (int n = 0; n < glob_mesh_num.GetM(); n++)
      {
        if (glob_mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
          nb_mesh_num_discont++;
        else
          nb_mesh_num_cont++;
      }
    
    if (glob_mesh_num.GetM() > 0)
      {
        NumDof.Reallocate(nb_dof + 2*nb_mesh_num_cont*nb_matching_proc + nb_matching_dofs + nb_elt*nb_mesh_num_discont);
        NumDofPML.Reallocate(nb_dof + nb_elt*nb_mesh_num_discont);
      }
    
    nb = 0 ; int nb_pml = 0;
    for (int n = 0; n < glob_mesh_num.GetM(); n++)
      {  
        if (glob_mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          {
            int nb_old = nb;
            for (int i = 0; i < nb_elt; i++)
              {
                int ne = NumElement(i);
                IVect Nodle = glob_mesh_num(n)->Element(ne).GetNodle();
                for (int j = 0; j < Nodle.GetM(); j++)
                  NumDof(nb+j) = Nodle(j);
                
                nb += Nodle.GetM();
              }
            
            for (int p = 0; p < MatchingDofs(n).GetM(); p++)
              if (ProcMatching(p))
                {
                  if (MatchingDofs(n)(p).GetM() > 0)
                    {
                      NumDof(nb) = p;
                      NumDof(nb+1) = MatchingDofs(n)(p).GetM();
                      for (int j = 0; j < MatchingDofs(n)(p).GetM(); j++)
                        NumDof(nb+2+j) = MatchingDofs(n)(p)(j);
                      
                      nb += 2+MatchingDofs(n)(p).GetM();
                    }
                  else
                    {
                      NumDof(nb) = p;
                      NumDof(nb+1) = 0;
                      nb += 2;
                    }
                }

            if (glob_mesh_num(n)->compute_dof_pml)
              {
                nb = nb_old;
                // we send dofs for PML if needed
                for (int i = 0; i < nb_elt; i++)
                  {
                    int ne = NumElement(i);
                    for (int j = 0; j < glob_mesh_num(n)->Element(ne).GetNbDof(); j++)
                      {
                        if (NumDof(nb+j) >= 0)
                          NumDofPML(nb_pml+j) = glob_mesh_num(n)->GetDofPML(NumDof(nb+j));
                        else
                          NumDofPML(nb_pml+j) = -1;
                      }
                    
                    nb_pml += glob_mesh_num(n)->Element(ne).GetNbDof();
                    nb += glob_mesh_num(n)->Element(ne).GetNbDof();
                  }
                
                for (int p = 0; p < MatchingDofs(n).GetM(); p++)
                  if (MatchingDofs(n)(p).GetM() > 0)
                    nb += 2+MatchingDofs(n)(p).GetM();
              }
          }
        else
          {
            int nb_old = nb;
            for (int i = 0; i < nb_elt; i++)
              {
                int ne = NumElement(i);
                NumDof(nb + i) = glob_mesh_num(n)->Element(ne).GetNumberDof(0);
              }
            
            nb += nb_elt;
            
            if (glob_mesh_num(n)->compute_dof_pml)
              {
                nb = nb_old;
                for (int i = 0; i < nb_elt; i++)
                  {
                    if (NumDof(nb + i) >= 0)
                      NumDofPML(nb_pml + i) = glob_mesh_num(n)->GetDofPML(NumDof(nb+i));
                    else
                      NumDofPML(nb_pml + i) = -1;
                  }
                
                nb_pml += nb_elt;
                nb += nb_elt;
              }
          }
      }
    
    if (rank != proc)
      {
        MPI_Send(NumDof.GetData(), nb, MPI_INTEGER, proc, tag, comm);    
        tag++;
      }
    
    if (rank != proc)
      {
        MPI_Send(NumDofPML.GetData(), nb_pml, MPI_INTEGER, proc, tag, comm);    
        tag++;
      }
    
    // and offsets for vectorial dofs
    if (OffsetDofV.GetM() > 0)
      {
        if (rank != proc)
          {
            MPI_Send(OffsetDofV.GetData(), OffsetDofV.GetM(), MPI_INTEGER, proc, tag, comm);
            tag++;
          }
      }
  }


  template<>
  void RecvMeshFromProcessor(Mesh<Dimension3>& mesh, Vector<MeshNumbering<Dimension3>* >& mesh_num,
			     Vector<IVect>& ConnecFace, IVect& MatchingProc, Vector<IVect>& MatchingDofs,
                             ParamParallelMesh<Dimension3>& Param, int proc, MPI_Comm& comm, int tag)
  {
    IVect& nodl_mesh = Param.nodl_mesh;
    IVect& num = Param.num; VectReal_wp& all_param = Param.all_param;
    VectR3& PointsEdgeRef = Param.PointsEdgeRef; IVect& NumBoundary = Param.NumBoundary;
    IVect& NumDofNeighbor = Param.NumDofNeighbor;
    IVect& NumElem = Param.NumElem; IVect& NumDof = Param.NumDof;
    VectR3& TranslatPeriodicFace = Param.TranslatPeriodicEdge;
    IVect& OffsetDofV = Param.OffsetDofV; IVect& NumDofPML = Param.NumDofPML;
    Vector<IVect>& offset_matching = Param.offset_matching;
    
    MPI_Status status;
    int rank; MPI_Comm_rank(comm, &rank);

    if (mesh.print_level >= 7)
      cout << rank << " Receiving numbers of vertices, edges, etc " << endl;
    
    // receiving number of vertices, edges and elements
    if (rank != proc)
      {
        int nb_param_int;
        MPI_Recv(&nb_param_int, 1, MPI_INTEGER, proc, tag, comm, &status); tag++;        
        num.Reallocate(nb_param_int);
        
        MPI_Recv(num.GetData(), nb_param_int, MPI_INTEGER, proc, tag,
                 comm, &status);    
        tag++;
      }
    
    Vector<int64_t> all_param_tmp;    
    int nb_vertices = num(0);
    int nb_edges_ref = num(1);
    int nb_faces_ref = num(2);
    int nb_tri_ref = num(3);
    int nb_quad_ref = num(4);
    int nb_int_face = 3*nb_tri_ref + 4*nb_quad_ref;

    int nb_elt = num(5);
    int nb_intV_elt = num(6);
    int nb_intE_elt = num(7);
    int nb_intF_elt = num(8);
    int nb_ref = num(9);
    int nb_param = num(10);
    
    int order_geom = num(11);
    int nb_mesh_num = num(12);
    int nb_matching_proc = num(13);
    int nb_matching_dofs = num(14);
    int nb_glob_faces_ref = num(16);
    //int nb_edges = num(17);
    int size_offset_v = num(19);    
    bool variable_order = false;
    if (num(21) != MeshNumbering<Dimension3>::CONSTANT_ORDER)
      variable_order = true;

    int nb_pml_areas = num(22);
    bool only_quad = num(23);
    mesh.SetGeometryOrder(order_geom, only_quad);

    if (nb_mesh_num != mesh_num.GetM())
      {
        cout << "Number of mesh numbering does not coincide " << nb_mesh_num << " different from " << mesh_num.GetM() << endl;
        abort();
      }
    
    for (int n = 0; n < nb_mesh_num; n++)
      {
        mesh_num(n)->compute_dof_pml = num(20);
        mesh_num(n)->SetFormulationForPeriodicCondition(num(24));
        mesh_num(n)->SetOrder(num(30+n));
        mesh_num(n)->number_map.SetFormulationDG(num(30+n+nb_mesh_num));
      }
    
    if (nb_ref > mesh.GetNbReferences())
      mesh.ResizeNbReferences(nb_ref);
    
    if (mesh.print_level >= 7)
      cout << rank << " Receiving curves " << endl;
    
    // parameters of curves
    if (rank != proc)
      {
        all_param.Reallocate(nb_param);
        MpiRecv(comm, all_param, all_param_tmp, nb_param, proc, tag, status);
        tag++;
      }

    if (mesh.print_level >= 7)
      cout << rank << " Receiving curves completed" << endl;
    
    VectReal_wp param;
    IVect ref_neighbor(nb_ref+1); ref_neighbor.Fill(-3);
    nb_param = 5*nb_ref + 1;
    for (int i = 1; i <= nb_ref; i++)
      {
        int p0 = toInteger(all_param(i-1));
        int p1 = toInteger(all_param(i));
        
        mesh.SetBodyNumber(i, toInteger(all_param(nb_ref+i)));
        mesh.SetBoundaryCondition(i, toInteger(all_param(2*nb_ref+i)));
        mesh.SetCurveType(i, toInteger(all_param(3*nb_ref+i)));
        if (p1 > p0)
          {
            param.Reallocate(p1-p0);
            for (int p = p0; p < p1; p++)
              param(p-p0) = all_param(p);
            
            mesh.SetCurveParameter(i, param);
          }
        
        ref_neighbor(i) = toInteger(all_param(4*nb_ref+i));
      }

    Real_wp alpha = all_param(nb_param); nb_param++;
    mesh.SetPeriodicAlpha(alpha);

    mesh.ClearPeriodicCondition();
    int nb_periodic_ref = toInteger(all_param(nb_param)); nb_param++;
    for (int i = 0; i < nb_periodic_ref; i++)
      {
	mesh.AddPeriodicCondition(TinyVector<int, 2>(toInteger(all_param(nb_param)),
                                                     toInteger(all_param(nb_param+1))),
				  toInteger(all_param(nb_param+2)));
        
	nb_param += 3;
      }
    
    mesh.ReallocatePmlAreas(nb_pml_areas);
    for (int i = 0; i < nb_pml_areas; i++)
      mesh.GetPmlArea(i).SetRegion(all_param, nb_param);
    
    mesh.SetOriginalNeighborReference(ref_neighbor);
    
    // receiving vertices and internal points of edges
    mesh.ReallocateEdgesRef(nb_edges_ref);
    mesh.ReallocateBoundariesRef(nb_faces_ref);
    mesh.ReallocateVertices(nb_vertices);
    int r = mesh.GetGeometryOrder();
    int nb_nodes_tri = (r-1)*(r-2)/2, nb_nodes_quad = (r-1)*(r-1);
    if (mesh.print_level >= 7)
      cout << rank << " Receiving internal points of the mesh " << endl;
    
    if (rank != proc)
      {
        PointsEdgeRef.Reallocate(nb_vertices + (r-1)*nb_edges_ref
                                 + (r-1)*(r-1)*nb_quad_ref + (r-1)*(r-2)/2*nb_tri_ref);
        
        MpiRecv(comm, reinterpret_cast<Real_wp*>(PointsEdgeRef.GetData()), all_param_tmp,
                3*PointsEdgeRef.GetM(), proc, tag, status);
        tag++;
      }
    
    for (int i = 0; i < nb_vertices; i++)
      mesh.Vertex(i) = PointsEdgeRef(i);
    
    int nb = nb_vertices;
    if (r > 1)
      {
        for (int i = 0; i < nb_edges_ref; i++)
          for (int k = 0; k < r-1; k++)
            mesh.SetPointInsideEdge(i, k, PointsEdgeRef(nb++));
      }
    
    int nbf = nb;
    
    // receiving connectivity of faces
    if (mesh.print_level >= 7)
      cout << rank << " Receiving referenced faces " << endl;

    int nb_param_per_face = 14 + 2*nb_mesh_num;
    ConnecFace.Reallocate(nb_faces_ref);
    if (rank != proc)
      {
        NumBoundary.Reallocate(nb_param_per_face*nb_faces_ref + nb_int_face + 3*nb_edges_ref);
        MPI_Recv(NumBoundary.GetData(), NumBoundary.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;
      }
    
    nb = 0;
    int offset = nb + nb_param_per_face*nb_faces_ref, nb_dof_neighbor = 0;
    for (int n = 0; n < nb_mesh_num; n++)
      mesh_num(n)->InitPeriodicBoundary();

    IVect Global_to_local_face(nb_glob_faces_ref);
    Global_to_local_face.Fill(-1);
    for (int i = 0; i < nb_faces_ref; i++)
      Global_to_local_face(NumBoundary(nb_param_per_face*i + 9)) = i;
    
    IVect OrderQuad(nb_faces_ref); VectR3 Pts;
    for (int i = 0; i < nb_faces_ref; i++)
      {
	int nb_vert = NumBoundary(nb);
	int ref = NumBoundary(nb+1);
	int rf = NumBoundary(nb+2);
	num.Reallocate(nb_vert);
	for (int j = 0; j < nb_vert; j++)
          num(j) = NumBoundary(offset+j);
	
       	mesh.BoundaryRef(i).Init(num, ref);
	OrderQuad(i) = rf;
	int nb_nodes = nb_nodes_quad;
	if (nb_vert == 3)
	  nb_nodes = nb_nodes_tri;
	
	Pts.Reallocate(nb_nodes);
	for (int k = 0; k < nb_nodes; k++)
	  Pts(k) = PointsEdgeRef(nbf + k);
	
	mesh.SetPointInsideFace(i, Pts);

	int nf_opp = NumBoundary(nb+12);
	int type_per = NumBoundary(nb+13);
	
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            int num_elt2 = NumBoundary(nb+4);
            int proc2 = NumBoundary(nb+5);
            int ne_loc2 = NumBoundary(nb+6);
            int num_loc2 = NumBoundary(nb+7);
            int rot = NumBoundary(nb+8);
            int type_elt = NumBoundary(nb+10);
            int iglob = NumBoundary(nb+9);
            int ref2 = NumBoundary(nb+11);

            int size_connec = 10 + nb_mesh_num;
            for (int n = 0; n < nb_mesh_num; n++)
              {
                //int re = NumBoundary(nb+14+n);
                int nb_dof = NumBoundary(nb+14+n+nb_mesh_num);                
                nb_dof_neighbor += nb_dof;
                size_connec += nb_dof + 3;
              }

            ConnecFace(i).Reallocate(size_connec);
            ConnecFace(i)(0) = iglob;
            ConnecFace(i)(1) = num_elt2;
            ConnecFace(i)(2) = proc2;
            ConnecFace(i)(3) = ne_loc2;
            ConnecFace(i)(4) = num_loc2;
            ConnecFace(i)(5) = rot;
            ConnecFace(i)(6) = type_elt;
            ConnecFace(i)(7) = ref2;
            ConnecFace(i)(8) = 0;
	    ConnecFace(i)(9) = 0;

            int offset_connec = 10+nb_mesh_num;
            for (int n = 0; n < nb_mesh_num; n++)
              {
                int re = NumBoundary(nb+14+n);
                int nb_dof = NumBoundary(nb+14+n+nb_mesh_num);
                ConnecFace(i)(10+n) = offset_connec;
                ConnecFace(i)(offset_connec) = re;
                ConnecFace(i)(offset_connec+1) = nb_dof;
                ConnecFace(i)(offset_connec+2) = 0;
                offset_connec += nb_dof + 3;
              }
          }
        else
	  {
	    ConnecFace(i).Clear();

	    // we find the local number for nf_opp
	    if (nf_opp >= 0)
	      nf_opp = Global_to_local_face(nf_opp);
	  }
	
	if (nf_opp != -1)
          for (int n = 0; n < nb_mesh_num; n++)
            {
              mesh_num(n)->SetPeriodicBoundary(i, nf_opp);
              mesh_num(n)->SetPeriodicityTypeForBoundary(i, type_per);
            }
        
	nbf += nb_nodes;
	nb += nb_param_per_face;
	offset += nb_vert;
      }

    Global_to_local_face.Clear();
    
    nb = offset;
    for (int i = 0; i < nb_edges_ref; i++)
      {
	mesh.EdgeRef(i).Init(NumBoundary(nb), NumBoundary(nb+1), NumBoundary(nb+2));
	nb += 3;
      }
    
    PointsEdgeRef.Clear();

    // retrieving periodic faces
    if (mesh.print_level >= 7)
      cout << rank << " Receiving periodic face translations " << endl;

    if (mesh.GetNbPeriodicReferences() > 0)
      {
	if (rank != proc)
	  {
	    TranslatPeriodicFace.Reallocate(nb_faces_ref);
	    MpiRecv(comm, reinterpret_cast<Real_wp*>(TranslatPeriodicFace.GetData()),
		    all_param_tmp, 3*nb_faces_ref, proc, tag, status);
	    
	    tag++;
	  }

	for (int n = 0; n < nb_mesh_num; n++)
          mesh_num(n)->GetTranslationPeriodicBoundary() = TranslatPeriodicFace;
      }
    
    // receiving dof numbers of neighboring elements
    if (mesh.print_level >= 7)
      cout << rank << " Receiving dof number of neighboring elements " << endl;
    
    if (rank != proc)
      {
        if (mesh_num(0)->compute_dof_pml)
          NumDofNeighbor.Reallocate(2*nb_dof_neighbor);
        else
          NumDofNeighbor.Reallocate(nb_dof_neighbor);
        
        MPI_Recv(NumDofNeighbor.GetData(),
                 NumDofNeighbor.GetM(), MPI_INTEGER, proc, tag,
                 comm, &status);
        tag++;
      }
    
    nb = 0;
    for (int i = 0; i < nb_faces_ref; i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            Vector<bool> element_in_pml(nb_mesh_num);
            element_in_pml.Fill(false); int nb_dof_all_pml = 0;
            for (int n = 0; n < nb_mesh_num; n++)
              if (mesh_num(n)->compute_dof_pml)
                {
                  int p = ConnecFace(i)(10+n);
                  int nb_dof = ConnecFace(i)(p+1);
                  for (int j = 0; j < nb_dof; j++)
                    if (NumDofNeighbor(nb + nb_dof + j) >= 0)
                      element_in_pml(n) = true;
                  
                  if (element_in_pml(n))
                    nb_dof_all_pml += nb_dof;
                }
            
            if (nb_dof_all_pml <= 0)
              {
                for (int n = 0; n < nb_mesh_num; n++)
                  {
                    int p = ConnecFace(i)(10+n);
                    int nb_dof = ConnecFace(i)(p+1);                    
                    for (int j = 0; j < nb_dof; j++)
                      ConnecFace(i)(p+3+j) =  NumDofNeighbor(nb+j);
            
                    nb += nb_dof;
                    if (mesh_num(n)->compute_dof_pml)
                      nb += nb_dof;
                  }
              }
            else
              {
                Vector<int> old_connec = ConnecFace(i);
                ConnecFace(i).Resize(old_connec.GetM() + nb_dof_all_pml);
                int offset_connec = 10 + nb_mesh_num, old_off = offset_connec;
                for (int n = 0; n < nb_mesh_num; n++)
                  {
                    int nb_dof = old_connec(old_off+1);
                    ConnecFace(i)(10+n) = offset_connec;
                    ConnecFace(i)(offset_connec) = old_connec(old_off);
                    ConnecFace(i)(offset_connec+1) = nb_dof;
                    for (int j = 0; j < nb_dof; j++)
                      ConnecFace(i)(offset_connec+3+j) =  NumDofNeighbor(nb+j);

                    nb += nb_dof;
                    if (element_in_pml(n))
                      {
                        ConnecFace(i)(offset_connec+2) = nb_dof;
                        offset_connec += nb_dof;
                        for (int j = 0; j < nb_dof; j++)
                          ConnecFace(i)(offset_connec+3+j) =  NumDofNeighbor(nb+j);
                      }
                    else
                      ConnecFace(i)(offset_connec+2) = 0;

                    if (mesh_num(n)->compute_dof_pml)
                      nb += nb_dof;
                    
                    offset_connec += nb_dof+3;
                    old_off += nb_dof+3;
                  }
              }
          }
      }
    
    // elements  
    if (mesh.print_level >= 7)
      cout << rank << " Receiving elements " << endl;
    
    mesh.ReallocateElements(nb_elt);
    if (rank != proc)
      {
        NumElem.Reallocate((5+3*nb_mesh_num)*nb_elt + 2*nb_intV_elt + (1+nb_mesh_num)*nb_intE_elt + (1+nb_mesh_num)*nb_intF_elt);
        MPI_Recv(NumElem.GetData(), NumElem.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;
      }
    
    nb = 0;
    Vector<IVect> OrderElt(nb_mesh_num), OrderInside(nb_mesh_num), NbDofElt(nb_mesh_num);
    for (int n = 0; n < nb_mesh_num; n++)
      {
        OrderElt(n).Reallocate(nb_elt); OrderInside(n).Reallocate(nb_elt);
        NbDofElt(n).Reallocate(nb_elt);
      }

    offset = (5+3*nb_mesh_num)*nb_elt; int nb_dof = 0;
    mesh.GlobElementNumber_Subdomain.Reallocate(nb_elt);
    mesh.GlobVertexNumber_Subdomain.Reallocate(nb_vertices);
    for (int i = 0; i < nb_elt; i++)
      {
        int nb_vert = NumElem(nb);
        int ref = NumElem(nb+1);
        int pml_num = NumElem(nb+2);
        int pml_type = NumElem(nb+3);
        num.Reallocate(nb_vert);
        for (int j = 0; j < nb_vert; j++)
          {
            num(j) = NumElem(offset+j);
            mesh.GlobVertexNumber_Subdomain(num(j)) = NumElem(offset + nb_vert + j);
          }
        
        mesh.Element(i).Init(num, ref);
	mesh.Element(i).SetPML(pml_num, pml_type);

        mesh.GlobElementNumber_Subdomain(i) = NumElem(nb+4);
        
        for (int n = 0; n < nb_mesh_num; n++)
          {
            OrderElt(n)(i) = NumElem(nb+5+n);
            OrderInside(n)(i) = NumElem(nb+5+n+nb_mesh_num);
            NbDofElt(n)(i) = NumElem(nb+5+n+2*nb_mesh_num);
            if (mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
              nb_dof += NbDofElt(n)(i);
          }
        
        nb += 5 + 3*nb_mesh_num;
        offset += 2*nb_vert + (1+nb_mesh_num)*mesh.Element(i).GetNbEdges() + (1+nb_mesh_num)*mesh.Element(i).GetNbFaces();
      }
    
    // updating connectivity
    mesh.ReorientElements();
    mesh.FindConnectivity();
    mesh.ProjectPointsOnCurves();

    // taking into account periodicity
    for (int i = 0; i < nb_faces_ref; i++)
      {
	int ref = NumBoundary(nb_param_per_face*i + 1);

	if (mesh.GetBoundaryCondition(ref) != BoundaryConditionEnum::LINE_NEIGHBOR)
	  if (mesh_num(0)->GetPeriodicBoundary(i) >= 0)
	    {
	      int i2 = mesh_num(0)->GetPeriodicBoundary(i);
	      int ne2 = mesh.BoundaryRef(i2).numElement(0);
	      mesh.BoundaryRef(i).AddElement(ne2);
	      mesh.Boundary(i).AddElement(ne2);
	    }
      }
    
    if (variable_order)
      {
        for (int n = 0; n < nb_mesh_num; n++)
          {
            mesh_num(n)->SetVariableOrder(mesh_num(n)->USER_ORDER);
        
            for (int i = 0; i < nb_elt; i++)
              {
                mesh_num(n)->SetOrderElement(i, OrderElt(n)(i));
                mesh_num(n)->SetOrderInside(i, OrderInside(n)(i));
              }
          }
        
        nb = (5+3*nb_mesh_num)*nb_elt;
        for (int i = 0; i < nb_elt; i++)
          {
            nb += 2*mesh.Element(i).GetNbVertices();
            int nv = mesh.Element(i).GetNbEdges();
            for (int j = 0; j < nv; j++)
              {
                int ne = mesh.Element(i).numEdge(j);
                for (int n = 0; n < nb_mesh_num; n++)
                  mesh_num(n)->SetOrderEdge(ne, NumElem(nb + (1+n)*nv + j));
              }
            
            nb += (1+nb_mesh_num)*nv;
            nv = mesh.Element(i).GetNbFaces();
            for (int j = 0; j < nv; j++)
              {
                int ne = mesh.Element(i).numFace(j);
                for (int n = 0; n < nb_mesh_num; n++)
                  mesh_num(n)->SetOrderFace(ne, NumElem(nb + (1+n)*nv + j));
              }
            
            nb += (1+nb_mesh_num)*nv;
          }
        
        for (int i = 0; i < nb_faces_ref; i++)
          for (int n = 0; n < nb_mesh_num; n++)
            mesh_num(n)->SetOrderQuadrature(i, OrderQuad(i));
        
        // order for intern boundary
        for (int n = 0; n < nb_mesh_num; n++)
          for (int i = nb_faces_ref; i < mesh.GetNbBoundary(); i++)
            {
              int rf = 0;
              for (int k = 0; k < mesh.Boundary(i).GetNbElements(); k++)
                {
                  if (mesh.Boundary(i).numElement(k) >= 0)
                    rf = max(rf, OrderElt(n)(mesh.Boundary(i).numElement(k)));
                }
              
              mesh_num(n)->SetOrderQuadrature(i, rf);
            }
      }
    
    nb = (5+3*nb_mesh_num)*nb_elt;
    mesh.GlobEdgeNumber_Subdomain.Reallocate(mesh.GetNbEdges());
    mesh.GlobFaceNumber_Subdomain.Reallocate(mesh.GetNbFaces());
    mesh.GlobEdgeNumber_Subdomain.Fill(-1);
    mesh.GlobFaceNumber_Subdomain.Fill(-1);
    for (int i = 0; i < nb_elt; i++)
      {
        nb += 2*mesh.Element(i).GetNbVertices();
        int nv = mesh.Element(i).GetNbEdges();
        for (int j = 0; j < nv; j++)
          {
            int ne = mesh.Element(i).numEdge(j);
            mesh.GlobEdgeNumber_Subdomain(ne) = NumElem(nb + j);
          }
        
        nb += (1+nb_mesh_num)*nv;
        nv = mesh.Element(i).GetNbFaces();
        for (int j = 0; j < nv; j++)
          {
            int ne = mesh.Element(i).numFace(j);
            mesh.GlobFaceNumber_Subdomain(ne) = NumElem(nb + j);
          }
        
        nb += (1+nb_mesh_num)*nv;
      }
    
    // retrieving dof numbers
    if (mesh.print_level >= 7)
      cout << rank << " Receiving dof numbers " << endl;
    
    int nb_mesh_num_discont = 0, nb_mesh_num_cont = 0;
    for (int n = 0; n < mesh_num.GetM(); n++)
      {
        if (mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
          nb_mesh_num_discont++;
        else
          nb_mesh_num_cont++;
      }
    
    if (mesh_num.GetM() > 0)
      {
        NumDof.Reallocate(nb_dof + 2*nb_mesh_num_cont*nb_matching_proc + nb_matching_dofs + nb_elt*nb_mesh_num_discont);
        NumDofPML.Reallocate(nb_dof + nb_elt*nb_mesh_num_discont);
      }
    
    if (rank != proc)
      {
        MPI_Recv(NumDof.GetData(), NumDof.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;

        MPI_Recv(NumDofPML.GetData(), NumDofPML.GetM(),
                 MPI_INTEGER, proc, tag, comm, &status);
        
        tag++;
      }

    MatchingDofs.Reallocate(nb_matching_proc);
    MatchingProc.Reallocate(nb_matching_proc);
    // counting size of arrays and allocating them
    nb = 0; IVect size_matching(nb_matching_proc);
    size_matching.Reallocate(nb_matching_proc);
    size_matching.Zero();
    for (int n = 0; n < nb_mesh_num; n++)
      {
        if (mesh_num(n)->number_map.FormulationDG() != ElementReference_Base::DISCONTINUOUS)
          {
            for (int i = 0; i < nb_elt; i++)
              nb += NbDofElt(n)(i);
            
            for (int p = 0; p < nb_matching_proc; p++)
              {
                int n_dof = NumDof(nb+1);
                MatchingProc(p) = NumDof(nb);
                size_matching(p) += n_dof;
                nb += 2 + n_dof;
              }
          }
        else
          nb += nb_elt;
      }
    
    offset_matching.Reallocate(nb_matching_proc);
    for (int p = 0; p < nb_matching_proc; p++)
      {
        MatchingDofs(p).Reallocate(size_matching(p));
        MatchingDofs(p).Fill(-1);
        
        offset_matching(p).Reallocate(nb_mesh_num+1);
        offset_matching(p).Zero();
      }
    
    size_matching.Zero();
    nb = 0; int nb_pml = 0;
    for (int n = 0; n < nb_mesh_num; n++)
      {
        if (mesh_num(n)->number_map.FormulationDG() == ElementReference_Base::DISCONTINUOUS)
          {
            nb_dof = 0;
            for (int i = 0; i < nb_elt; i++)
              nb_dof += NbDofElt(n)(i);
            
            mesh_num(n)->GlobDofNumber_Subdomain.Reallocate(nb_dof);
            mesh_num(n)->GlobDofNumber_Subdomain.Fill(-1);
            int offset = 0;
            mesh_num(n)->ReallocateElements(nb_elt);
            for (int i = 0; i < nb_elt; i++)
              {
                int off_g = NumDof(nb + i);
                mesh_num(n)->Element(i).ReallocateDof(1);
                mesh_num(n)->Element(i).SetNumberDof(0, offset);
                for (int j = 0; j < NbDofElt(n)(i); j++)
                  mesh_num(n)->GlobDofNumber_Subdomain(offset+j) = off_g + j;
                
                offset += NbDofElt(n)(i);
              }            
            
            mesh_num(n)->SetNbDof(offset);
            nb += nb_elt;
            
            if (mesh_num(n)->compute_dof_pml)
              {
                IVect IndexDof(nb_elt);
                IndexDof.Fill(-1);
                
                nb_dof = 0;
                for (int i = 0; i < nb_elt; i++)
                  if (NumDofPML(nb_pml + i) >= 0)
                    {
                      IndexDof(i) = nb_dof;
                      nb_dof += NbDofElt(n)(i);
                    }
                
                mesh_num(n)->GlobDofPML_Subdomain.Reallocate(nb_dof);
                mesh_num(n)->GlobDofPML_Subdomain.Fill(-1);
                mesh_num(n)->ReallocateDofPML(nb_dof);
                for (int i = 0; i < nb_elt; i++)
                  if (NumDofPML(nb_pml+i) >= 0)
                    {
                      int npml = NumDofPML(nb_pml+i);
                      int off = mesh_num(n)->Element(i).GetNumberDof(0);
                      for (int j = 0; j < NbDofElt(n)(i); j++)
                        {
                          mesh_num(n)->GlobDofPML_Subdomain(IndexDof(i) + j) = npml + j;
                          mesh_num(n)->SetDofPML(off+j, IndexDof(i)+j);
                        }
                    }

                nb_pml += nb_elt;
              }
	  }
        else
          {
            IVect IndexDof(nodl_mesh(n));
            IndexDof.Fill(-1);
            
            mesh_num(n)->ReallocateElements(nb_elt);
            for (int i = 0; i < nb_elt; i++)
              {
                mesh_num(n)->Element(i).ReallocateDof(NbDofElt(n)(i));            
                for (int j = 0; j < NbDofElt(n)(i); j++)
                  {
                    int n0 = NumDof(nb+j);
                    mesh_num(n)->Element(i).SetNumberDof(j, n0);
                    if (n0 >= 0)
                      IndexDof(n0) = 1;
                  }
                
                nb += NbDofElt(n)(i);
              }
            
            for (int p = 0; p < nb_matching_proc; p++)
              {
                int n_dof = NumDof(nb+1);
                for (int j = 0; j < n_dof; j++)
                  MatchingDofs(p)(size_matching(p) + j) = NumDof(nb+2+j);   

                size_matching(p) += n_dof;
                offset_matching(p)(n+1) = size_matching(p);
                nb += 2 + n_dof;
              }
            
            int nb_dof_true = 0;
            for (int i = 0; i < IndexDof.GetM(); i++)
              if (IndexDof(i) > 0)
                IndexDof(i) = nb_dof_true++;
            
            mesh_num(n)->SetNbDof(nb_dof_true);
            mesh_num(n)->GlobDofNumber_Subdomain.Reallocate(nb_dof_true);
            mesh_num(n)->GlobDofNumber_Subdomain.Fill(-1);
            for (int i = 0; i < nb_elt; i++)
              for (int j = 0; j < NbDofElt(n)(i); j++)
                {
                  int n0 = mesh_num(n)->Element(i).GetNumberDof(j);
                  if (n0 >= 0)
                    {
                      mesh_num(n)->GlobDofNumber_Subdomain(IndexDof(n0)) = n0;
                      mesh_num(n)->Element(i).SetNumberDof(j, IndexDof(n0));
                    }
                }

            // PML dofs
            int nb_pml_old = nb_pml;
            if (mesh_num(n)->compute_dof_pml)
              {
                IndexDof.Fill(-1);                
                for (int i = 0; i < nb_elt; i++)
                  {
                    for (int j = 0; j < NbDofElt(n)(i); j++)
                      {
                        int n0 = NumDofPML(nb_pml + j);
                        if (n0 >= 0)
                          IndexDof(n0) = 1;
                      }
                    
                    nb_pml += NbDofElt(n)(i);
                  }
                
                nb_dof = 0;
                for (int i = 0; i < IndexDof.GetM(); i++)
                  if (IndexDof(i) > 0)
                    IndexDof(i) = nb_dof++;
                
                mesh_num(n)->GlobDofPML_Subdomain.Reallocate(nb_dof);
                mesh_num(n)->GlobDofPML_Subdomain.Fill(-1);
                mesh_num(n)->ReallocateDofPML(nb_dof);
                nb_pml = nb_pml_old;
                for (int i = 0; i < nb_elt; i++)
                  {                
                    for (int j = 0; j < NbDofElt(n)(i); j++)
                      {
                        int n0 = mesh_num(n)->Element(i).GetNumberDof(j);
                        int npml = NumDofPML(nb_pml+j);
                        if (npml >= 0)
                          {
                            mesh_num(n)->GlobDofPML_Subdomain(IndexDof(npml)) = npml;
                            mesh_num(n)->SetDofPML(n0, IndexDof(npml));
                          }
                      }                
                    nb_pml += NbDofElt(n)(i);
                  }
              }
          }
      }
    
    // receiving offsets for vectorial dofs
    if (size_offset_v > 0)
      {
        if (rank != proc)
          {
            OffsetDofV.Reallocate(size_offset_v);
            MPI_Recv(OffsetDofV.GetData(), OffsetDofV.GetM(), MPI_INTEGER, proc, tag,
                     comm, &status);
            tag++;
          }
      }
    
    if (mesh.print_level >= 7)
      cout << rank << " Receiving completed " << endl;    
  }
#endif

}

#define MONTJOIE_FILE_PARALLEL_MESH_FUNCTIONS_CXX
#endif
