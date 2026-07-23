#ifndef MONTJOIE_FILE_EVALUATION_CFL_WAVE_CXX

namespace Montjoie
{

  //! sets nominal time step and stability function 
  void EigenvalueComparisonStabFunc::SetFunction(const Real_wp& dt, const UnivariatePolynomial<Real_wp>& P)
  {
    dt_nom = dt;
    Pol = P;
  }


  //! compares two eigenvalues Lr + 1j*Li and Lr2 + 1j*Li2
  int EigenvalueComparisonStabFunc::CompareEigenvalue(const Real_wp& Lr, const Real_wp& Li, const Real_wp& Lr2, const Real_wp& Li2)
  {
    // we compute z = \lambda \Delta t
    Complex_wp z = dt_nom*Complex_wp(Lr, Li);
    Complex_wp z2 = dt_nom*Complex_wp(Lr2, Li2);    

    // the modulus of R(z) is computed
    Real_wp R = abs(Pol.Evaluate(z));
    Real_wp R2 = abs(Pol.Evaluate(z2));
    
    // the eigenvalues are classified by increasing modulus of R(z)
    if (abs(R-R2) <= 1e-6)
      return 0;
    else if (R > R2)
      return -1;
    
    return 1;
  }
  
  
  //! finds the time step dt such that the eigenvalue Lr + 1j*Li belongs to the stability domain
  Real_wp EigenvalueComparisonStabFunc::FindTimeStep(const Real_wp& Lr, const Real_wp& Li) const
  {
    // bisection method to find dt such that |R(\lambda dt)| = 1
    Real_wp dt0 = 0.0, dt1 = dt_nom, dt = dt0;

    /* int nb_points = 10000;
    ofstream file_out("module.dat"); file_out.precision(15);
    for (int i = 0; i < nb_points; i++)
      {
        dt = dt0 + Real_wp(i+1)*(dt1-dt0) / nb_points;
        Complex_wp z = dt*Complex_wp(Lr, Li);
        Real_wp R = abs(Pol.Evaluate(z));
        file_out << dt << " " << R << '\n';
      }
    
      file_out.close(); */
    
    while (abs(dt0 - dt1) > 1e-10*dt_nom)
      {
        dt = (dt0 + dt1)/2;
        Complex_wp z = dt*Complex_wp(Lr, Li);
        Real_wp R = abs(Pol.Evaluate(z));
        if (R > 1)
          dt1 = dt;
        else
          dt0 = dt;                
      }
    
    return dt0;
  }
  
  
  /****************************
   * EigenProblemTimeMontjoie *
   ****************************/


  //! default constructor
  void EigenProblemTimeMontjoie::InitDefaultValues()
  {
    this->n_ = var_time.GetNumberOfUnknowns();
    this->SetDiagonalMass();
    this->print_level = var_problem.print_level;
    this->complex_system = false;

#ifdef SELDON_WITH_MPI
    this->SetCommunicator(var_comm.comm_group_mode);
#endif
    
    this->Init(this->n_);
    this->stopping_criterion = var_time.default_threshold_cfl;
  }
  
  
  //! Cholesky factorisation of mass matrix
  void EigenProblemTimeMontjoie::FactorizeCholeskyMass()
  {
    // already done by var_time.ComputeMassMatrix()
  }
  
  
  //! computation of mass matrix
  void EigenProblemTimeMontjoie::ComputeDiagonalMass()
  {
    Vector<Real_wp>& D = this->sqrt_diagonal_mass;
    // for first-order formulation, the mass matrix is integrated to the stiffness matrix
    if (var_problem.FirstOrderFormulation())
      {
	D.Reallocate(this->n_);
	D.Fill(Real_wp(1));
      }
    else
      {
	D = var_time.Glob_mat_Dh->GetDiagonalDh();
        if (D.GetM() == 0)
          {
            const VectReal_wp& invD = var_time.Glob_mat_Dh->GetInverseDiagonalDh();
            if (invD.GetM() == 0)
              {
                cout << "Diagonal not found" << endl;
                abort();
              }
            
            D.Reallocate(invD.GetM());
            for (int i = 0; i < D.GetM(); i++)
              D(i) = Real_wp(1) / invD(i);
          }
      }
    
    // D.WriteText("Dh.dat");
  }
  
  
  //! multiplication by stiffness matrix
  void EigenProblemTimeMontjoie
  ::MltStiffness(const Real_wp&, const Real_wp&,
		 const Vector<Real_wp>& X, Vector<Real_wp>& Y)
  {
    abort();
  }
  
  
  //! matrix vector product with stiffness matrix
  void EigenProblemTimeMontjoie
  ::MltStiffness(const Vector<Real_wp>& X, Vector<Real_wp>& Y)
  {
    // Dirichlet values are stored
    VectReal_wp& Xt = const_cast<Vector<Real_wp>& >(X);
    VectReal_wp Xdir(var_boundary.GetNbDirichletDof());
    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      {
        int num_dof = var_boundary.GetDirichletDofNumber(i);
        Xdir(i) = X(num_dof);
        SetComplexZero(Xt(num_dof));
      }
    
    if (var_problem.FirstOrderFormulation())
      var_time.EvaluateFunction(Real_wp(0), X, Y, true, false);
    else
      var_time.ApplyOperatorKh(-1.0, 0, X, 0.0, Y);

    // Dirichlet components are set to zero
    for (int i = 0; i < var_boundary.GetNbDirichletDof(); i++)
      {
        int num_dof = var_boundary.GetDirichletDofNumber(i);
        Xt(num_dof) = Xdir(i);
        SetComplexZero(Y(num_dof));
      }
  }

  void EigenProblemTimeMontjoie
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Real_wp>& X, Vector<Real_wp>& Y)
  {
    cout << "Not implemented" << endl;
    abort();
  }

  void EigenProblemTimeMontjoie
  ::MltStiffness(const Real_wp&, const Real_wp&,
                 const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
  {
    cout << "Not implemented" << endl;
    abort();
  }

  void EigenProblemTimeMontjoie
  ::MltStiffness(const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
  {
    cout << "Not implemented" << endl;
    abort();
  }

  void EigenProblemTimeMontjoie
  ::MltStiffness(const SeldonTranspose& trans, const Vector<Complex_wp>& X, Vector<Complex_wp>& Y)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  //! returns true if the eigenvalue problem is symmetric (mass and stiffness matrix)
  bool EigenProblemTimeMontjoie
  ::IsSymmetricProblem() const 
  {
    if (var_problem.FirstOrderFormulation())
      return false;
    
    return true;
  }


  //! returns true if the eigenvalue problem is symmetric (mass and stiffness matrix)
  bool EigenProblemTimeMontjoie
  ::IsHermitianProblem() const 
  {
    if (var_problem.FirstOrderFormulation())
      return false;
    
    return true;
  }

  
  //! debug function
  void EigenProblemTimeMontjoie::CheckStiffnessMatrix()
  {
    int nodl = this->n_;
    Matrix<Real_wp, General, ArrayRowSparse> Ah_full(nodl, nodl);
    Vector<Real_wp> Ones(nodl), Ah_Ones(nodl);
    for (int i = 0; i < nodl; i++)
      {
	Ones.Fill(0); Ones(i) = Real_wp(1);
	cout << "colonne " << i << endl;
	MltStiffness(Ones, Ah_Ones);
	for (int j = 0; j < nodl; j++)
	  if (Ah_Ones(j) != Real_wp(0))
	    Ah_full.AddInteraction(j, i, Ah_Ones(j));
      }
    
    Ah_full.WriteText("KhTest.dat");
  }


  /*******************
   * MeshGeometryCFL *
   *******************/
  
  template<class Dimension>
  void MeshGeometryCFL<Dimension>
  ::Init(const DistributedProblem<Dimension>& var_problem)
  {
    const Mesh<Dimension>& mesh = var_problem.mesh;
    IndexVertex.Reallocate(mesh.GetNbVertices());
    ElementOnSubdomain.Reallocate(mesh.GetNbElt());
    IndexVertex.Fill(-1);
    ElementOnSubdomain.Fill(false);
    
    IndexNeighbor.Reallocate(mesh.GetNbBoundaryRef());
    IndexNeighbor.Fill(-1);

#ifdef SELDON_WITH_MPI   
    int nb_proc = var_problem.GetNbProcPerMode();
    if (nb_proc <= 1)
      return;

    // loop over referenced faces to find elements to send to other processors
    Vector<int> NbVertexSendPerProc(nb_proc), NbVertexRecvPerProc(nb_proc);
    NbVertexSendPerProc.Zero(); NbVertexRecvPerProc.Zero();
    int nb_vertices_faces = 0, nb_faces_neigh = 0;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            // element to send
            int p = var_problem.GetProcessorNeighboringFace(i);
            int ne = mesh.BoundaryRef(i).numElement(0);
            int nb_vert = mesh.Element(ne).GetNbVertices();
            NbVertexSendPerProc(p) += nb_vert;

            // element to receive
            int type_elt = var_problem.GetTypeEltNeighboringFace(i);
            nb_vert = 1 + Dimension::dim_N + type_elt;
            if (type_elt == 3)
              nb_vert = 8;
            
            NbVertexRecvPerProc(p) += nb_vert;
            nb_vertices_faces += nb_vert;
            nb_faces_neigh++;
          }
      }
    
    // we allocate vertices numbers and coordinates (in arrays NumVertexXXX and VertexXXX)
    Vector<Vector<int> > NumVertexSend(nb_proc), NumVertexRecv(nb_proc);
    Vector<typename Dimension::VectR_N> VertexSend(nb_proc), VertexRecv(nb_proc); 
    for (int p = 0; p < nb_proc; p++)
      {
        NumVertexSend(p).Reallocate(NbVertexSendPerProc(p));
        NumVertexRecv(p).Reallocate(NbVertexRecvPerProc(p));
        VertexSend(p).Reallocate(NbVertexSendPerProc(p));
        VertexRecv(p).Reallocate(NbVertexRecvPerProc(p));
      }

    // filling vertices to be sent
    NbVertexSendPerProc.Zero();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            // element to send
            int p = var_problem.GetProcessorNeighboringFace(i);
            int ne = mesh.BoundaryRef(i).numElement(0);
            int nb_vert = mesh.Element(ne).GetNbVertices();
            int offset = NbVertexSendPerProc(p);
            for (int j = 0; j < nb_vert; j++)
              {
                int nv = mesh.Element(ne).numVertex(j);
                NumVertexSend(p)(offset+j) = mesh.GlobVertexNumber_Subdomain(nv);
                VertexSend(p)(offset+j) = mesh.Vertex(nv);
              }
            
            NbVertexSendPerProc(p) += nb_vert;
          }
      }

    // MPI communications
    Vector<int64_t> xtmp;
    Vector<Vector<int64_t> > xtmp_send(nb_proc), xtmp_recv(nb_proc);
    Vector<MPI_Request> request_send(2*nb_proc), request_recv(2*nb_proc);
    const MPI_Comm& comm = var_problem.comm_group_mode;
    for (int p = 0; p < nb_proc; p++)
      if (NbVertexSendPerProc(p) != 0)
        {
          request_send(p) = MpiIsend(comm, NumVertexSend(p), xtmp, NbVertexSendPerProc(p), p, 131);
          request_send(nb_proc + p) = MpiIsend(comm, &VertexSend(p)(0)(0), xtmp_send(p),
                                               Dimension::dim_N*NbVertexSendPerProc(p), p, 132);
        }
    
    for (int p = 0; p < nb_proc; p++)
      if (NbVertexRecvPerProc(p) != 0)
        {
          request_recv(p) = MpiIrecv(comm, NumVertexRecv(p), xtmp, NbVertexRecvPerProc(p), p, 131);
          request_recv(nb_proc + p) = MpiIrecv(comm, &VertexRecv(p)(0)(0), xtmp_recv(p),
                                               Dimension::dim_N*NbVertexRecvPerProc(p), p, 132);
        }
    
    MPI_Status status;
    for (int p = 0; p < nb_proc; p++)
      if (NbVertexSendPerProc(p) != 0)
        {
          MPI_Wait(&request_send(p), &status);
          MPI_Wait(&request_send(nb_proc+p), &status);
        }

    for (int p = 0; p < nb_proc; p++)
      if (NbVertexRecvPerProc(p) != 0)
        {
          MPI_Wait(&request_recv(p), &status);
          MPI_Wait(&request_recv(nb_proc+p), &status);
        }
    
    for (int p = 0; p < nb_proc; p++)
      if (NbVertexRecvPerProc(p) != 0)
        MpiCompleteIrecv(&VertexRecv(p)(0)(0), xtmp_recv(p), Dimension::dim_N*NbVertexRecvPerProc(p));
    
    xtmp_send.Clear(); xtmp_recv.Clear();
    VertexSend.Clear(); NumVertexSend.Clear();

    // then filling VertexFace and NumVertexFace
    Vector<int> NumVertexFace(nb_vertices_faces), NumElement(nb_faces_neigh);
    typename Dimension::VectR_N VertexFace(nb_vertices_faces);
    NbVertexRecvPerProc.Zero(); nb_vertices_faces = 0; nb_faces_neigh = 0;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            int p = var_problem.GetProcessorNeighboringFace(i);
            int type_elt = var_problem.GetTypeEltNeighboringFace(i);
            int nb_vert = 1 + Dimension::dim_N + type_elt;
            if (type_elt == 3)
              nb_vert = 8;
            
            int offset = NbVertexRecvPerProc(p);
            for (int j = 0; j < nb_vert; j++)
              {
                NumVertexFace(nb_vertices_faces + j) = NumVertexRecv(p)(offset+j);
                VertexFace(nb_vertices_faces + j) = VertexRecv(p)(offset+j);
              }
            
            NumElement(nb_faces_neigh) = var_problem.GetElementNumberNeighboringFace(i);
            NbVertexRecvPerProc(p) += nb_vert;
            nb_vertices_faces += nb_vert;
            nb_faces_neigh++;
          }
      }
    
    NumVertexRecv.Clear(); VertexRecv.Clear();
    
    // checking that global numbers are sorted
    for (int i = 0; i < mesh.GetNbVertices()-1; i++)
      if (mesh.GlobVertexNumber_Subdomain(i+1) < mesh.GlobVertexNumber_Subdomain(i))
        {
          cout << "invalid numbers " << endl;
          abort();
        }

    for (int i = 0; i < mesh.GetNbElt()-1; i++)
      if (mesh.GlobElementNumber_Subdomain(i+1) < mesh.GlobElementNumber_Subdomain(i))
        {
          cout << "invalid numbers " << endl;
          abort();
        }
    
    // NumElement is sorted
    Vector<int> permut_elt(NumElement.GetM());
    permut_elt.Fill();
    Sort(NumElement, permut_elt);
    
    Vector<int> NumLocalElt(nb_faces_neigh);
    NumLocalElt.Fill(-1);

    int nb_elt = 0, pos0 = 0;
    while (pos0 < NumElement.GetM())
      {
        int pos1 = pos0 + 1;
        while ((pos1 < NumElement.GetM()) && (NumElement(pos1) == NumElement(pos0)))
          pos1++;
        
        for (int k = pos0; k < pos1; k++)
          NumLocalElt(permut_elt(k)) = nb_elt;
        
        nb_elt++;
        pos0 = pos1;
      }
    
    // NumVertexFace is sorted
    Vector<int> permut_face(NumVertexFace.GetM());
    permut_face.Fill();
    Sort(NumVertexFace, permut_face);
    
    // creating PointsNeighbor
    PointsNeighbor.Reallocate(NumVertexFace.GetM());
    Vector<int> NumLocalVertex(NumVertexFace.GetM());
    NumLocalVertex.Fill(-1);
    
    int nb_vertex = 0;
    pos0 = 0; int pos_loc = 0;
    while (pos0 < NumVertexFace.GetM())
      {
        int pos1 = pos0 + 1;
        while ((pos1 < NumVertexFace.GetM()) && (NumVertexFace(pos1) == NumVertexFace(pos0)))
          pos1++;
        
        int nv = NumVertexFace(pos0);
        while ((pos_loc < mesh.GetNbVertices()) && (mesh.GlobVertexNumber_Subdomain(pos_loc) < nv))
          pos_loc++;
        
        bool dist_point = true;
        if ((pos_loc < mesh.GetNbVertices()) && (mesh.GlobVertexNumber_Subdomain(pos_loc) == nv))
          dist_point = false;
        
        if (dist_point)
          {
            PointsNeighbor(nb_vertex) = VertexFace(permut_face(pos0));
            for (int k = pos0; k < pos1; k++)
              NumLocalVertex(permut_face(k)) = -nb_vertex-1;
            
            nb_vertex++;
          }
        else
          {
            for (int k = pos0; k < pos1; k++)
              NumLocalVertex(permut_face(k)) = pos_loc;
          }

        pos0 = pos1;
      }

    PointsNeighbor.Resize(nb_vertex);

    // creating IndexNeighbor and ElementNeighbor
    ReferenceNeighbor.Reallocate(nb_elt);
    ElementNeighbor.Reallocate(nb_elt);
    Vector<bool> ElementTreated(nb_elt);
    ElementTreated.Fill(false);
    nb_vertices_faces = 0; nb_faces_neigh = 0;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int ref = mesh.BoundaryRef(i).GetReference();
        if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR)
          {
            int type_elt = var_problem.GetTypeEltNeighboringFace(i);
            int nb_vert = 1 + Dimension::dim_N + type_elt;
            if (type_elt == 3)
              nb_vert = 8;
            
            int ne = NumLocalElt(nb_faces_neigh);
            IndexNeighbor(i) = ne;
            if (!ElementTreated(ne))
              {
                ReferenceNeighbor(ne) = var_problem.GetRefDomainNeighboringFace(i);
                ElementNeighbor(ne).Reallocate(nb_vert);
                for (int j = 0; j < nb_vert; j++)
                  ElementNeighbor(ne)(j) = NumLocalVertex(nb_vertices_faces+j);
              }

            nb_vertices_faces += nb_vert;
            nb_faces_neigh++;
          }
      }

    NeighIndexVertex.Reallocate(PointsNeighbor.GetM());
    NeighElementOnSubdomain.Reallocate(ElementNeighbor.GetM());
    NeighIndexVertex.Fill(-1);
    NeighElementOnSubdomain.Fill(false);

#endif
  }

  //! extracting a small mesh from the global mesh so that CFL will be computed on this mesh
  template<class Dimension>
  void MeshGeometryCFL<Dimension>
  ::ExtractSubMesh(const Mesh<Dimension>& mesh, const MeshNumbering<Dimension>& mesh_num,
                   Mesh<Dimension>& sub_mesh, MeshNumbering<Dimension>& sub_mesh_num,
                   const IVect& num, const IVect& num_vertex,
                   const IVect& num_neigh, const IVect& num_vertex_neigh)
  {
    int nb_vertices_subdomain = num_vertex.GetM();
    int nb_elt_subdomain = num.GetM();
    
    // no curved elements
    sub_mesh.SetGeometryOrder(1);
    sub_mesh.ReallocateVertices(nb_vertices_subdomain + num_vertex_neigh.GetM());
    sub_mesh.ReallocateElements(nb_elt_subdomain + num_neigh.GetM());
    for (int i = 0; i < num_vertex.GetM(); i++)
      sub_mesh.Vertex(i) = mesh.Vertex(num_vertex(i));

    for (int i = 0; i < num_vertex_neigh.GetM(); i++)
      sub_mesh.Vertex(nb_vertices_subdomain + i) = PointsNeighbor(num_vertex_neigh(i));
    
    for (int i = 0; i < num.GetM(); i++)
      {
        int n0 = num(i);
        int nb_vert = mesh.Element(n0).GetNbVertices();
        IVect numv(nb_vert);
        int ref = mesh.Element(n0).GetReference();
        for (int k = 0; k < nb_vert; k++)
          numv(k) = IndexVertex(mesh.Element(n0).numVertex(k));
        
        sub_mesh.Element(i).Init(numv, ref);
      }

    for (int i = 0; i < num_neigh.GetM(); i++)
      {
        int n0 = num_neigh(i);
        int nb_vert = ElementNeighbor(n0).GetM();
        IVect numv(nb_vert);
        int ref = ReferenceNeighbor(n0);
        for (int k = 0; k < nb_vert; k++)
          {
            int nv = ElementNeighbor(n0)(k);
            if (nv >= 0)
              numv(k) = IndexVertex(nv);
            else
              numv(k) = nb_vertices_subdomain + NeighIndexVertex(-nv-1);
          }
        
        sub_mesh.Element(num.GetM() + i).Init(numv, ref);
      }
    
    // changing boundaries
    sub_mesh.ReorientElements();
    sub_mesh.FindConnectivity();
    sub_mesh.AddBoundaryEdges();
    sub_mesh.FindConnectivity();
    
    // Neumann boundary condition for internal boundaries
    //int cond_intern = BoundaryConditionEnum::LINE_ABSORBING;
    int cond_intern = BoundaryConditionEnum::LINE_NEUMANN;
    for (int ref = 1; ref <= sub_mesh.GetNbReferences(); ref++)
      if ((mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_INSIDE) ||
          (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR))
        sub_mesh.SetBoundaryCondition(ref, cond_intern);
    
    // filling sub_mesh_num
    int order = mesh_num.GetOrderElement(num(0));
    sub_mesh_num.SetOrder(order);
    bool order_variable = false;
    for (int k = 0; k < nb_elt_subdomain; k++)
      {
        int r2 = mesh_num.GetOrderElement(num(k));
        sub_mesh_num.SetOrderElement(k, r2);
        if (r2 != order)
          order_variable = true;
      }
    
    if (order_variable)
      {
        for (int j = 0; j < sub_mesh.GetNbBoundary(); j++)
          {
            int n0 = sub_mesh.Boundary(j).numElement(0);
            sub_mesh_num.SetOrderQuadrature(j, sub_mesh_num.GetOrderElement(n0));
          }
        
        for (int k = 0; k < num.GetM(); k++)
          for (int j = 0; j < mesh.Element(num(k)).GetNbBoundary(); j++)
            {
              int rf = mesh_num.GetOrderQuadrature(mesh.Element(num(k)).numBoundary(j));
              int ne = sub_mesh.Element(k).numBoundary(j);
              sub_mesh_num.SetOrderQuadrature(ne, rf);
            }
      }
    else
      sub_mesh_num.SetVariableOrder(mesh_num.CONSTANT_ORDER);
  }        
  
  
  //! creates a small mesh containing element i and its neighbors
  template<class Dimension>
  void MeshGeometryCFL<Dimension>
  ::CreateSmallMesh(int i, const Mesh<Dimension>& mesh, const MeshNumbering<Dimension>& mesh_num,
                    Mesh<Dimension>& small_mesh, MeshNumbering<Dimension>& small_mesh_num)
  {
    // first we place the element i in this small mesh :
    IVect num(1); num(0) = i;
    int nb_elt_subdomain = 1;
    int nb_vertices_subdomain = mesh.Element(i).GetNbVertices();
    IVect num_vertex(nb_vertices_subdomain);
    ElementOnSubdomain(i) = true;
    for (int k = 0; k < nb_vertices_subdomain; k++)
      {
        IndexVertex(mesh.Element(i).numVertex(k)) = k;
        num_vertex(k) = mesh.Element(i).numVertex(k);
      }
    
    // then we take surrounding elements            
    int nb_vert_neigh = 0;
    IVect num_neigh, num_vertex_neigh;
    for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
      {
        int nf = mesh.Element(i).numBoundary(j);
        if (mesh.Boundary(nf).GetNbElements() == 2)
          {
            int n0 = mesh.Boundary(nf).numElement(0);
            if (n0 == i)
              n0 = mesh.Boundary(nf).numElement(1);
            
            num.PushBack(n0);
            ElementOnSubdomain(n0) = true;
            nb_elt_subdomain++;
            int nb_vert = mesh.Element(n0).GetNbVertices();
            for (int k = 0; k < nb_vert; k++)
              {
                int nv = mesh.Element(n0).numVertex(k);
                if (IndexVertex(nv) == -1)
                  {
                    IndexVertex(nv) = nb_vertices_subdomain;
                    nb_vertices_subdomain++;
                    num_vertex.PushBack(nv);
                  }
              }
          }
        else if ((nf < mesh.GetNbBoundaryRef())
                 && (IndexNeighbor(nf) >= 0))
          {
            // element on another processor
            int n0 = IndexNeighbor(nf);
            NeighElementOnSubdomain(n0) = true;
            num_neigh.PushBack(n0);
            int nb_vert = ElementNeighbor(n0).GetM();
            for (int k = 0; k < nb_vert; k++)
              {
                int nv = ElementNeighbor(n0)(k);
                if (nv >= 0)
                  {
                    if (IndexVertex(nv) == -1)
                      {
                        IndexVertex(nv) = nb_vertices_subdomain;
                        nb_vertices_subdomain++;
                        num_vertex.PushBack(nv);
                      }
                  }
                else
                  {
                    // vertex on another processor
                    if (NeighIndexVertex(-nv-1) == -1)
                      {
                        NeighIndexVertex(-nv-1) = nb_vert_neigh++;
                        num_vertex_neigh.PushBack(-nv-1);
                      }
                  }
              }
          }
      }
    
    // creation of the small mesh
    ExtractSubMesh(mesh, mesh_num, small_mesh, small_mesh_num,
                   num, num_vertex, num_neigh, num_vertex_neigh);
    
    // IndexVertex and ElementOnSubdomain are placed in the initial configuation
    for (int k = 0; k < nb_elt_subdomain; k++)
      ElementOnSubdomain(num(k)) = false; 
    
    for (int k = 0; k < nb_vertices_subdomain; k++)
      IndexVertex(num_vertex(k)) = -1;                 

    for (int k = 0; k < num_neigh.GetM(); k++)
      NeighElementOnSubdomain(num_neigh(k)) = false;
    
    for (int k = 0; k < nb_vert_neigh; k++)
      NeighIndexVertex(num_vertex_neigh(k)) = -1;
  }
  

  /*************
   * Functions *
   *************/
  

  //! evaluation of CFL 
  Real_wp EvaluateCFL_Generic(EigenProblemTimeMontjoie& pb_eig, const UnivariatePolynomial<Real_wp>& P,
                              const Real_wp& dt, bool first_order_pb)
  {
    EigenvalueComparisonStabFunc class_comp;
    class_comp.SetFunction(dt, P);
    
    pb_eig.SetNbAskedEigenvalues(10);
    pb_eig.SetNbArnoldiVectors(30);
    pb_eig.SetNbAdditionalEigenvalues(10);
    pb_eig.SetComputationalMode(pb_eig.REGULAR_MODE);

    int type_sort = EigenProblem_Base<Real_wp>::SORTED_MODULUS;
    int type_solver = TypeEigenvalueSolver::ARPACK;
    if (first_order_pb)
      {
        type_sort = EigenProblem_Base<Real_wp>::SORTED_USER;
        pb_eig.SetUserComparisonClass(&class_comp);
        type_solver = TypeEigenvalueSolver::SLEPC;
      }
    
    pb_eig.SetTypeSpectrum(pb_eig.LARGE_EIGENVALUES, 0.0, type_sort);
    
    Vector<Real_wp> eigen_values, eigen_imag;
    Matrix<Real_wp, General, ColMajor> eigen_vec;
    GetEigenvaluesEigenvectors(pb_eig, eigen_values, eigen_imag, eigen_vec, type_solver);
    //DISP(eigen_values); DISP(eigen_imag);

    Real_wp cfl = 0;
    if (first_order_pb)
      {
        // time step is searched by bisection
        cfl = class_comp.FindTimeStep(eigen_values(0), eigen_imag(0));
      }
    else
      {
        // maximal eigenvalue
        for (int i = 0; i < pb_eig.GetNbAskedEigenvalues(); i++)
          cfl = max(cfl, sqrt(eigen_values(i)*eigen_values(i) + eigen_imag(i)*eigen_imag(i)));
        
	cfl = P(0)/sqrt(cfl);
      }
    
    return cfl;
  }
  
  
  //! evaluating CFL for each element of the mesh
  /*!
    For each element of the mesh, we extract the element and its neighbours
    and compute the CFL of this small mesh
   */
  template<class Dimension>
  void ComputeLocalTimeStep(VarInstationary_Dim<Dimension>& var_time,
                            VectReal_wp& dt)
  {
    DistributedProblem<Dimension>& var_problem = var_time.GetDistributedProblem();
    
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    const Mesh<Dimension>& mesh = var_problem.mesh;
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    dt.Reallocate(mesh.GetNbElt());
    dt.Zero();
    
    VarInstationary_Dim<Dimension>* var;
    var = var_time.GetNewHyperbolicProblem();
    //DISP(var_time.GetTimeStep());

    DistributedProblem<Dimension>& var_p = var->GetDistributedProblem();
    var_p.CopyInputData(var_problem);
    
#ifdef SELDON_WITH_MPI
    var_p.comm_group_mode = MPI_COMM_SELF;
#endif
    
    var_p.GetMeshNumbering(0).number_map = mesh_num.number_map;
        
    var->SetInputData(string("OrderTimeScheme"), var_time.param_time_scheme);
    var->SetTimeStep(var_time.GetTimeStep());
    
    var_p.GetOutputProblem().grid_to_be_computed = false;
    var_p.mesh.print_level = -1;
    var_p.print_level = -1;
    var->ForceComputationTimeStep(false);
    var->default_threshold_cfl = var_time.default_threshold_cfl;
    
    Real_wp dt_min = 1e30, dt_max = 0;
    
    if (var_time.GetFileStoringLocalDt().size() > 2)
      {
        // case where time steps are stored in a file
        VectReal_wp glob_dt;
        glob_dt.ReadText(var_time.GetFileStoringLocalDt());
        
        if (nb_proc == 1)
          {
            dt = glob_dt;
            if (dt.GetM() != mesh.GetNbElt())
              {
                cout << "Wrong file containing local time steps ? " << endl;
                DISP(dt.GetM());
                DISP(mesh.GetNbElt());
                abort();
              }
          }
        else
          {
            for (int i = 0; i < mesh.GetNbElt(); i++)
              dt(i) = glob_dt(mesh.GlobElementNumber_Subdomain(i));
          }
        
        for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            dt_min = min(dt_min, dt(i));
            dt_max = max(dt_max, dt(i));
          }
      }
    else if ((var_time.GetCFL_ComputationalMethod() == var_time.CFL_MESH_SIZE_ELT)
	     || (var_time.GetCFL_ComputationalMethod() == var_time.CFL_MESH_SIZE_FACE)
	     || (var_time.GetCFL_ComputationalMethod() == var_time.CFL_MESH_SIZE_VERTEX))
      {
	// case where the local time step is computed by the geometrical size of elements
	
	// first we compute the lengths of all the edges
	VectReal_wp length_edge(mesh.GetNbEdges());
	for (int i = 0; i < mesh.GetNbEdges(); i++)
	  {
	    int n0 = mesh.GetEdge(i).numVertex(0), n1 = mesh.GetEdge(i).numVertex(1);
	    length_edge(i) = mesh.Vertex(n0).Distance(mesh.Vertex(n1));
	  }

	// then for each element, we take the minimal length
	VectReal_wp size_elt(mesh.GetNbElt());
	size_elt.Zero();
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	      {
		int ne = mesh.Element(i).numEdge(j);
		if (j == 0)
		  size_elt(i) = length_edge(ne);
		else
		  size_elt(i) = min(size_elt(i), length_edge(ne));		
	      }
	  }

	// then depending on the strategy, we take into account neighboring elements
	if (var_time.GetCFL_ComputationalMethod() == var_time.CFL_MESH_SIZE_FACE)
	  {
	    // for each face, we set h = min(two elements adjacent)
	    VectReal_wp length_face(mesh.GetNbBoundary());
	    for (int i = 0; i < length_face.GetM(); i++)
	      {
		int ne = mesh.Boundary(i).numElement(0);
		length_face(i) = size_elt(ne);
		if (mesh.Boundary(i).GetNbElements() > 1)
		  {
		    ne = mesh.Boundary(i).numElement(1);
		    length_face(i) = min(length_face(i), size_elt(ne));
		  }
	      }

	    // then each element has the minimal length of its faces
	    for (int i = 0; i < mesh.GetNbElt(); i++)
	      {
		for (int j = 0; j < mesh.Element(i).GetNbBoundary(); j++)
		  size_elt(i) = min(size_elt(i), length_face(mesh.Element(i).numBoundary(j)));
	      }
	  }
	else if (var_time.GetCFL_ComputationalMethod() == var_time.CFL_MESH_SIZE_VERTEX)
	  {
	    // similar process with vertices
	    VectReal_wp length_vertex(mesh.GetNbVertices());
	    length_vertex.Fill(1e300);
	    for (int i = 0; i < mesh.GetNbElt(); i++)
	      {
		for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
		  {
		    int nv = mesh.Element(i).numVertex(j);
		    length_vertex(nv) = min(length_vertex(nv), size_elt(i));
		  }
	      }
	    
	    for (int i = 0; i < mesh.GetNbElt(); i++)
	      {
		for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
		  size_elt(i) = min(size_elt(i), length_vertex(mesh.Element(i).numVertex(j)));
	      }
	  }

	// we set dt = h (Geometrical point of view)
	dt = size_elt;

	for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            dt_min = min(dt_min, dt(i));
            dt_max = max(dt_max, dt(i));
          }
	
	if (nb_proc > 1)
          {
            cout << "not implemented " << endl;
            abort();
          }
        else
          {
            dt.WriteText("local_dt2D.don");            
          }
      }
    else
      {
        MeshGeometryCFL<Dimension> var_geom;
        var_geom.Init(var_problem);
        
        string suffix = "_P" + to_str(var_problem.GetRankProcMode()) + ".mesh";
	// case where the local time step is evaluated by computing eigenvalues associated with small meshes
        int old_percent = 0, new_percent = 0;
        for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            // displays a progress bar
            new_percent = toInteger(round(Real_wp(i)/(mesh.GetNbElt()-1)*80));
            if (rank_proc == 0)
              if (var_problem.print_level >= 2)
                for (int percent = old_percent; percent < new_percent; percent++)
                  { cout<<"#"; cout.flush(); }
            
            old_percent = new_percent;
            
            // small mesh with surrounding elements
            var_p.mesh.Clear();
            
            var_geom.CreateSmallMesh(i, mesh, mesh_num, var_p.mesh, var_p.GetMeshNumbering(0));
            	    
	    var_p.CopyFiniteElement(var_problem);
	    
            // var.var_harmonic.exit_if_no_boundary_condition = false;
	    //var_p.mesh.Write("sub"+to_str(i)+suffix);
            var_p.GetMeshNumbering(0).NumberMesh();
	    var_p.ComputeNumberOfDofs();
            
            dt(i) = var->EvaluateCFL();
            
            if (dt(i) < dt_min)
              dt_min = dt(i);
            
            if (dt(i) > dt_max)
              dt_max = dt(i);
	    
            /*DistributedMatrix<Real_wp, General, ArrayRowSparse> Kh, Mh;
	    GlobalGenericMatrix<Real_wp> nat_mat;

	    nat_mat.SetCoefMass(1.0);
	    nat_mat.SetCoefDamping(0.0);
	    nat_mat.SetCoefStiffness(0.0);
	    var_p.GetComputationProblem().AddMatrixWithBC(Mh, nat_mat);
	    
	    nat_mat.SetCoefMass(0.0);
	    nat_mat.SetCoefDamping(0.0);
	    nat_mat.SetCoefStiffness(1.0);
	    var_p.GetComputationProblem().AddMatrixWithBC(Kh, nat_mat);

	    if (rank_proc == 0)
              {
                Mh.WriteText("Mh.dat");
                Kh.WriteText("Kh.dat");
              }
            */
            
	    //int test_input; cout << "waiting" << endl; cin >> test_input;
          }
        
        if (rank_proc == 0)
          if (var_problem.print_level >= 2)
            cout << endl;
        
        if (nb_proc > 1)
          {
#ifdef SELDON_WITH_MPI
            Vector<int64_t> xtmp;
            if (rank_proc == 0)
              {
                Vector<int> num_glob, nb_elt(nb_proc);
                int nb_global_elt = dt.GetM();
                Vector<Real_wp> global_dt, local_dt; MPI_Status status;
                for (int i = 1; i < nb_proc; i++)
                  {
                    MPI_Recv(&nb_elt(i), 1, MPI_INTEGER, i, 132, var_problem.comm_group_mode, &status);
                    nb_global_elt += nb_elt(i);
                  }
                
                global_dt.Reallocate(nb_global_elt);
                global_dt.Zero();
                for (int j = 0; j < mesh.GetNbElt(); j++)
                  global_dt(mesh.GlobElementNumber_Subdomain(j)) = dt(j);
                
                for (int i = 1; i < nb_proc; i++)
                  {
                    num_glob.Reallocate(nb_elt(i));
                    local_dt.Reallocate(nb_elt(i));
                    MPI_Recv(num_glob.GetData(), nb_elt(i), MPI_INTEGER, i, 133, var_problem.comm_group_mode, &status);
                    MpiRecv(var_problem.comm_group_mode, local_dt, xtmp, nb_elt(i), i, 134, status);
                    
                    for (int j = 0; j < nb_elt(i); j++)
                      global_dt(num_glob(j)) = local_dt(j);
                  }

                global_dt.WriteText("local_dt.don");
              }
            else
              {
                int nb_elt = mesh.GlobElementNumber_Subdomain.GetM();
                MPI_Send(&nb_elt, 1, MPI_INTEGER, 0, 132, var_problem.comm_group_mode);
                MPI_Send(mesh.GlobElementNumber_Subdomain.GetData(), mesh.GetNbElt(),
                          MPI_INTEGER, 0, 133, var_problem.comm_group_mode);
                
                MpiSsend(var_problem.comm_group_mode, dt, xtmp, dt.GetM(), 0, 134);
              }
#endif
          }
        else
          {
            dt.WriteText("local_dt.don");            
          }
      }

 #ifdef SELDON_WITH_MPI
    Vector<int64_t> xtmp;
    Real_wp dt_tmp(dt_min);
    MpiReduce(var_problem.comm_group_mode, &dt_tmp, xtmp, &dt_min, 1, MPI_MIN, 0);
    dt_tmp = dt_max;
    MpiReduce(var_problem.comm_group_mode, &dt_tmp, xtmp, &dt_max, 1, MPI_MAX, 0);
#endif
    
    if (rank_proc == 0)
      {
        cout << "Maximal time step = " << dt_max << endl;
        cout << "Minimal time step = " << dt_min << endl;
        cout << "Ratio between time steps = " << dt_max/dt_min << endl;
      }

    delete var;
  }

  
  Real_wp EvaluateCFL_Unsteady(VarComputationProblem& var, VarInstationary_Base& var_time,
			       EigenProblemTimeMontjoie& pb_eig, const UnivariatePolynomial<Real_wp>& P, bool init_computation)
  {
    if (init_computation)
      {
        var.ComputeMassMatrix();
        var_time.ComputeMassMatrix();
	var_time.ComputeStiffnessMatrix();
        var_time.InitTimeIterations();
      }
    
    VarProblem_Base& var_base = dynamic_cast<VarProblem_Base&>(var);
    bool first_order = var_base.FirstOrderFormulation();
    return EvaluateCFL_Generic(pb_eig, P, var_time.GetTimeStep(), first_order);
  }
  
}

#define MONTJOIE_FILE_EVALUATION_CFL_WAVE_CXX
#endif

