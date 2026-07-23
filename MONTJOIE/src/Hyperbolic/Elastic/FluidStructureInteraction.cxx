#ifndef MONTJOIE_FILE_FLUID_STRUCTURE_INTERACTION_CXX

namespace Montjoie
{
  
  /**********************************
   * FluidStructureInteraction_Base *
   **********************************/

  
  //! constructor
  template<class Dimension>
  FluidStructureInteraction_Base<Dimension>::
  FluidStructureInteraction_Base(VarProblem<Dimension>& var_f,
				 VarProblem<Dimension>& var_s) 
    : var_fluid(var_f), var_solid(var_s)
  {
    nodl = 0;
  }
  

  //! construction of the mesh, and separation of the mesh into fluid part and solid part
  template<class Dimension>
  void FluidStructureInteraction_Base<Dimension>
  ::ComputeMeshAndFiniteElement(const string& name_element, bool first_order_formulation)
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif
    
    // the lowest order is used for the geometry
    Mesh<Dimension>& mesh = var_fluid.mesh;
    int rf = mesh.GetGeometryOrder();
    int rs = var_solid.mesh.GetGeometryOrder();
    int r = min(rf, rs);
    mesh.SetGeometryOrder(r);

    // enforcing first-order formulation if needed
    if (first_order_formulation)
      {
	var_fluid.SetFirstOrderFormulation(true);
        var_solid.SetFirstOrderFormulation(true);
      }
    
    if (rank_proc == 0)
      {
	// we construct the general mesh on the root processor
	mesh.ConstructMesh(0, var_fluid.mesh_data(0));
	
	int ref_max = 1;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  ref_max = max(ref_max, mesh.Element(i).GetReference());
	
	// is_fluid(i) will be true if the reference i is associated with the fluid
	Vector<bool> is_fluid(ref_max+1);
	is_fluid.Fill(false);
	for (int i = 0; i < ref_fluid.GetM(); i++)
	  if (ref_fluid(i) <= ref_max)
	    is_fluid(ref_fluid(i)) = true;
	
	// we add referenced faces for faces at the interface (if they don't exist)
	int nb_new_faces = 0, nb_face_on_interface = 0;
	for (int i = 0; i < mesh.GetNbBoundary(); i++)
	  if (mesh.Boundary(i).GetNbElements() == 2)
	    {
	      int n1 = mesh.Boundary(i).numElement(0);
	      int n2 = mesh.Boundary(i).numElement(1);
	      int ref1 = mesh.Element(n1).GetReference();
	      int ref2 = mesh.Element(n2).GetReference();
	      if (is_fluid(ref1) != is_fluid(ref2))
		{
		  nb_face_on_interface++;
		  if (i >= mesh.GetNbBoundaryRef())
		    nb_new_faces++;
		}
	    }
	
	int nb_old = mesh.GetNbBoundaryRef();
	int new_ref = -1;
	if (nb_new_faces > 0)
	  {
	    new_ref = mesh.GetNewReference();
	    mesh.ResizeBoundariesRef(nb_old + nb_new_faces);
	  }
	
	Vector<bool> is_interface(mesh.GetNbReferences()+1);
	is_interface.Fill(false);
	if (new_ref != -1)
	  is_interface(new_ref) = true;
	
	// we retrieve for each boundary on the interface
	// the element number and local position
	// for the fluid and solid
	this->NumElementFluid.Reallocate(nb_face_on_interface);
	this->NumElementSolid.Reallocate(nb_face_on_interface);
	this->NumLocFluid.Reallocate(nb_face_on_interface);
	this->NumLocSolid.Reallocate(nb_face_on_interface);
	nb_face_on_interface = 0;
	nb_new_faces = 0;
	for (int i = 0; i < mesh.GetNbBoundary(); i++)
	  if (mesh.Boundary(i).GetNbElements() == 2)
	    {
	      int n1 = mesh.Boundary(i).numElement(0);
	      int n2 = mesh.Boundary(i).numElement(1);
	      int ref1 = mesh.Element(n1).GetReference();
	      int ref2 = mesh.Element(n2).GetReference();
	      if (is_fluid(ref1) != is_fluid(ref2))
		{
		  int pos1 = mesh.Element(n1).GetPositionBoundary(i);
		  int pos2 = mesh.Element(n2).GetPositionBoundary(i);
		  if (is_fluid(ref1))
		    {
		      this->NumElementFluid(nb_face_on_interface) = n1;
		      this->NumLocFluid(nb_face_on_interface) = pos1;
		      this->NumElementSolid(nb_face_on_interface) = n2;
		      this->NumLocSolid(nb_face_on_interface) = pos2;
		    }
		  else
		    {
		      this->NumElementFluid(nb_face_on_interface) = n2;
		      this->NumLocFluid(nb_face_on_interface) = pos2;
		      this->NumElementSolid(nb_face_on_interface) = n1;
		      this->NumLocSolid(nb_face_on_interface) = pos1;
		    }
		  
		  nb_face_on_interface++;
		  if (i >= nb_old)
		    {
		      mesh.Boundary(i).SetReference(new_ref);
		      mesh.BoundaryRef(nb_old + nb_new_faces) = mesh.Boundary(i);
		      nb_new_faces++;
		    }
		  else
		    {
		      is_interface(mesh.BoundaryRef(i).GetReference()) = true;
		    }
		}
	    }    
	
	int nb_ref = 0;
	for (int i = 1; i < is_interface.GetM(); i++)
	  if (is_interface(i))
	    nb_ref++;
	
	// constructing ref_interface
	ref_interface.Reallocate(nb_ref);
	nb_ref = 0;
	for (int i = 1; i < is_interface.GetM(); i++)
	  if (is_interface(i))
	    ref_interface(nb_ref++) = i;
	
	// separating the mesh into fluid part and solid part
	IVect solid_boundary_condition = var_solid.mesh.GetBoundaryCondition();
	var_solid.mesh = mesh;
	for (int ref = 1; ref < solid_boundary_condition.GetM(); ref++)
	  var_solid.mesh.SetBoundaryCondition(ref, solid_boundary_condition(ref));
	
	for (int i = 1; i < is_interface.GetM(); i++)
	  if (is_interface(i))
	    {
	      var_solid.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_TRANSMISSION);
	      var_fluid.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_TRANSMISSION);
	    }
	
	int nb_elt = 0;
	IVect IndexElt(mesh.GetNbElt()); IndexElt.Fill(-1);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = mesh.Element(i).GetReference();
	    if (is_fluid(ref))
	      {
		mesh.Element(nb_elt) = mesh.Element(i);
		IndexElt(i) = nb_elt;
		nb_elt++;
	      }
	  }
	
	mesh.ResizeElements(nb_elt);
	
	// changing element numbers 
	for (int i = 0; i < nb_face_on_interface; i++)
	  this->NumElementFluid(i) = IndexElt(this->NumElementFluid(i));
	
	nb_elt = 0;
	IndexElt.Reallocate(var_solid.mesh.GetNbElt());
	IndexElt.Fill(-1);
	for (int i = 0; i < var_solid.mesh.GetNbElt(); i++)
	  {
	    int ref = var_solid.mesh.Element(i).GetReference();
	    if (!is_fluid(ref))
	      {
		var_solid.mesh.Element(nb_elt) = var_solid.mesh.Element(i);
		IndexElt(i) = nb_elt;
		nb_elt++;
	      }
	  }
	
	var_solid.mesh.ResizeElements(nb_elt);
	
	// changing element numbers     
	for (int i = 0; i < nb_face_on_interface; i++)
	  this->NumElementSolid(i) = IndexElt(this->NumElementSolid(i));
	
	var_fluid.mesh.ForceCoherenceMesh();
	var_solid.mesh.ForceCoherenceMesh();
	var_fluid.mesh.ReorientElements();
	var_solid.mesh.ReorientElements();
	var_fluid.mesh.FindConnectivity();
	var_solid.mesh.FindConnectivity();
	
	// writing the two meshes if needed
	if (var_fluid.print_level >= 4)
	  {
	    var_fluid.mesh.Write("fluid.mesh");
	    var_solid.mesh.Write("solid.mesh");
	  }
	
	var_fluid.mesh_data.Clear();
	var_solid.mesh_data.Clear();
        
	// updating orders of quadrature
        MeshNumbering<Dimension>& mesh_num_fluid = var_fluid.GetMeshNumbering(0);
        MeshNumbering<Dimension>& mesh_num_solid = var_solid.GetMeshNumbering(0);
	if (mesh_num_fluid.GetOrder() != mesh_num_solid.GetOrder())
	  {
	    mesh_num_fluid.SetVariableOrder(MeshNumbering<Dimension2>::USER_ORDER);
	    mesh_num_solid.SetVariableOrder(MeshNumbering<Dimension2>::USER_ORDER);
	    
	    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
	      {
		int ne_s = this->NumElementSolid(i);
		int num_face_s = var_solid.mesh.Element(ne_s).numBoundary(this->NumLocSolid(i));
		int rs = mesh_num_solid.GetOrderQuadrature(num_face_s);
		
		int ne_f = this->NumElementFluid(i);
		int num_face_f = var_fluid.mesh.Element(ne_f).numBoundary(this->NumLocFluid(i));
		int rf = mesh_num_fluid.GetOrderQuadrature(num_face_f);
		
		int r = max(rf, rs);
		mesh_num_fluid.SetOrderQuadrature(num_face_f, r);
		mesh_num_solid.SetOrderQuadrature(num_face_s, r);	    
	      }	
	    
	    //var_fluid.UpdateInterpolationElement();
	    //var_solid.UpdateInterpolationElement();
	  }
      }

#ifdef SELDON_WITH_MPI
    int nb_elt_fluid_global = var_fluid.mesh.GetNbElt();
    int nb_elt_solid_global = var_solid.mesh.GetNbElt();
#endif
    
    // now the fluid and solid meshes are split into the different processors
    var_fluid.ComputeMeshAndFiniteElement(name_element);
    var_solid.ComputeMeshAndFiniteElement(name_element);

    if (nb_proc == 1)
      return;
    
#ifdef SELDON_WITH_MPI

    // for parallel execution, we need to construct additional arrays
    //var_fluid.mesh.Write("fluid"+to_str(rank_proc)+".mesh");
    //var_solid.mesh.Write("solid"+to_str(rank_proc)+".mesh");
    MeshNumbering<Dimension>& mesh_num_fluid = var_fluid.GetMeshNumbering(0);
    MeshNumbering<Dimension>& mesh_num_solid = var_solid.GetMeshNumbering(0);

    MPI_Status status;
    Vector<int64_t> xtmp;
    if (rank_proc == 0)
      {	
	// we recover epart and local elements numbers for the fluid mesh and solid mesh
	Vector<int> EpartFluid(nb_elt_fluid_global), NumLocalFluid(nb_elt_fluid_global);
	Vector<int> EpartSolid(nb_elt_solid_global), NumLocalSolid(nb_elt_solid_global);
	EpartFluid.Zero(); NumLocalFluid.Fill(-1);
	EpartSolid.Zero(); NumLocalSolid.Fill(-1);
	Vector<int> num(2);	
	for (int i = 0; i < nb_proc; i++)
	  {
	    IVect NumEltFluid, NumEltSolid;
	    if (i == 0)
	      {
		NumEltFluid = var_fluid.mesh.GlobElementNumber_Subdomain;
		NumEltSolid = var_solid.mesh.GlobElementNumber_Subdomain;
	      }
	    else
	      {
		MpiRecv(comm, num, xtmp, 2, i, 134, status);
		NumEltFluid.Reallocate(num(0));
		NumEltSolid.Reallocate(num(1));	    
		MpiRecv(comm, NumEltFluid, xtmp, num(0), i, 135, status);
		MpiRecv(comm, NumEltSolid, xtmp, num(1), i, 136, status);
	      }
	    
	    for (int j = 0; j < NumEltFluid.GetM(); j++)
	      {
		EpartFluid(NumEltFluid(j)) = i;
		NumLocalFluid(NumEltFluid(j)) = j;
	      }
	    
	    for (int j = 0; j < NumEltSolid.GetM(); j++)
	      {
		EpartSolid(NumEltSolid(j)) = i;
		NumLocalSolid(NumEltSolid(j)) = j;
	      }
	  }

	Vector<IVect> NumElt_Fluid(nb_proc), NumElt_Solid(nb_proc);
	Vector<IVect> NumLoc_Fluid(nb_proc), NumLoc_Solid(nb_proc);
	Vector<IVect> ProcElt_Fluid(nb_proc), ProcElt_Solid(nb_proc);
	// then loop over elements of the interface
	for (int i = 0; i < this->NumElementFluid.GetM(); i++)
	  {
	    int ne_f = this->NumElementFluid(i);
	    int ne_s = this->NumElementSolid(i);
	    int pf = EpartFluid(ne_f);
	    int ps = EpartSolid(ne_s);
	    NumElt_Fluid(pf).PushBack(NumLocalFluid(ne_f));
	    NumLoc_Fluid(pf).PushBack(this->NumLocFluid(i));
	    ProcElt_Fluid(pf).PushBack(ps);
	    
	    NumElt_Solid(ps).PushBack(NumLocalSolid(ne_s));
	    NumLoc_Solid(ps).PushBack(this->NumLocSolid(i));
	    ProcElt_Solid(ps).PushBack(pf);
	  }
	
	for (int i = 0; i < nb_proc; i++)
	  {
	    if (i == 0)
	      {
		this->NumElementFluid = NumElt_Fluid(0);
		this->NumLocFluid = NumLoc_Fluid(0);
		this->ProcElementFluid = ProcElt_Fluid(0);
		this->NumElementSolid = NumElt_Solid(0);
		this->NumLocSolid = NumLoc_Solid(0);
		this->ProcElementSolid = ProcElt_Solid(0);
	      }
	    else
	      {
		num(0) = NumElt_Fluid(i).GetM();
		num(1) = NumElt_Solid(i).GetM();
		MpiSsend(comm, num, xtmp, 2, i, 137);

		if (num(0) > 0)
		  {
		    MpiSsend(comm, NumElt_Fluid(i), xtmp, num(0), i, 138);
		    MpiSsend(comm, NumLoc_Fluid(i), xtmp, num(0), i, 139);
		    MpiSsend(comm, ProcElt_Fluid(i), xtmp, num(0), i, 140);
		  }
		
		if (num(1) > 0)
		  {
		    MpiSsend(comm, NumElt_Solid(i), xtmp, num(1), i, 141);
		    MpiSsend(comm, NumLoc_Solid(i), xtmp, num(1), i, 142);
		    MpiSsend(comm, ProcElt_Solid(i), xtmp, num(1), i, 143);
		  }
	      }
	  }
	
      }
    else
      {        
	Vector<int> num(2);
	num(0) = var_fluid.mesh.GlobElementNumber_Subdomain.GetM();
	num(1) = var_solid.mesh.GlobElementNumber_Subdomain.GetM();
	
	MpiSsend(comm, num, xtmp, 2, 0, 134);
	MpiSsend(comm, var_fluid.mesh.GlobElementNumber_Subdomain, xtmp, num(0), 0, 135);	
	MpiSsend(comm, var_solid.mesh.GlobElementNumber_Subdomain, xtmp, num(1), 0, 136);
	
	MpiRecv(comm, num, xtmp, 2, 0, 137, status);
	this->NumElementFluid.Reallocate(num(0));
	this->NumLocFluid.Reallocate(num(0));
	this->ProcElementFluid.Reallocate(num(0));

	this->NumElementSolid.Reallocate(num(1));
	this->NumLocSolid.Reallocate(num(1));
	this->ProcElementSolid.Reallocate(num(1));
	
	if (num(0) > 0)
	  {
	    MpiRecv(comm, this->NumElementFluid, xtmp, num(0), 0, 138, status);
	    MpiRecv(comm, this->NumLocFluid, xtmp, num(0), 0, 139, status);
	    MpiRecv(comm, this->ProcElementFluid, xtmp, num(0), 0, 140, status);
	  }
	
	if (num(1) > 0)
	  {
	    MpiRecv(comm, this->NumElementSolid, xtmp, num(1), 0, 141, status);
	    MpiRecv(comm, this->NumLocSolid, xtmp, num(1), 0, 142, status);
	    MpiRecv(comm, this->ProcElementSolid, xtmp, num(1), 0, 143, status);
	  }
      }

    Vector<IVect> EntierToSendFluid(nb_proc), EntierToSendSolid(nb_proc);
    Vector<IVect> EntierToRecvFluid(nb_proc), EntierToRecvSolid(nb_proc);
    IVect NbEntierFluid(nb_proc), NbEntierSolid(nb_proc);
    NbEntierFluid.Zero(); NbEntierSolid.Zero();
    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
      {
	int ne_f = this->NumElementFluid(i);
	int p = this->ProcElementFluid(i);
	if (p != rank_proc)
	  NbEntierFluid(p) += 8 + mesh_num_fluid.GetNbLocalDof(ne_f);
      }

    for (int i = 0; i < this->NumElementSolid.GetM(); i++)
      {
	int ne_s = this->NumElementSolid(i);
	int p = this->ProcElementSolid(i);
	if (p != rank_proc)
	  NbEntierSolid(p) += 8 + mesh_num_solid.GetNbLocalDof(ne_s);
      }
    
    for (int i = 0; i < nb_proc; i++)
      {
	EntierToSendFluid(i).Reallocate(NbEntierFluid(i));
	EntierToSendSolid(i).Reallocate(NbEntierSolid(i));
      }
    
    NbEntierFluid.Zero(); NbEntierSolid.Zero();
    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
      {
	int ne_f = this->NumElementFluid(i);
	int p = this->ProcElementFluid(i);
	int num_loc = this->NumLocFluid(i);
	if (p != rank_proc)
	  {
	    int& nb = NbEntierFluid(p);
	    int rot = var_fluid.mesh.Element(ne_f).GetOrientationBoundary(num_loc);
	    IVect Nodle = var_fluid.GetDofNumberOnElement(ne_f);
	    EntierToSendFluid(p)(nb++) = Nodle.GetM();
	    EntierToSendFluid(p)(nb++) = var_fluid.mesh.GlobElementNumber_Subdomain(ne_f);
	    EntierToSendFluid(p)(nb++) = ne_f;
	    EntierToSendFluid(p)(nb++) = num_loc;
	    EntierToSendFluid(p)(nb++) = rot;
	    EntierToSendFluid(p)(nb++) = mesh_num_fluid.GetOrderElement(ne_f);
	    EntierToSendFluid(p)(nb++) = var_fluid.mesh.GetTypeElement(ne_f);
	    EntierToSendFluid(p)(nb++) = Nodle.GetM();
	    for (int j = 0; j < Nodle.GetM(); j++)
	      EntierToSendFluid(p)(nb++) = mesh_num_fluid.GlobDofNumber_Subdomain(Nodle(j));
	  }
	else
	  {
	    EntierToRecvSolid(p).PushBack(ne_f);
	    EntierToRecvSolid(p).PushBack(num_loc);
	  }
      }

    for (int i = 0; i < this->NumElementSolid.GetM(); i++)
      {
	int ne_s = this->NumElementSolid(i);
	int p = this->ProcElementSolid(i);
	int num_loc = this->NumLocSolid(i);
	if (p != rank_proc)
	  {
	    int& nb = NbEntierSolid(p);
	    int rot = var_solid.mesh.Element(ne_s).GetOrientationBoundary(num_loc);
	    IVect Nodle = var_solid.GetDofNumberOnElement(ne_s);
	    EntierToSendSolid(p)(nb++) = Nodle.GetM();
	    EntierToSendSolid(p)(nb++) = var_solid.mesh.GlobElementNumber_Subdomain(ne_s);
	    EntierToSendSolid(p)(nb++) = ne_s;
	    EntierToSendSolid(p)(nb++) = num_loc;
	    EntierToSendSolid(p)(nb++) = rot;
	    EntierToSendSolid(p)(nb++) = mesh_num_solid.GetOrderElement(ne_s);
	    EntierToSendSolid(p)(nb++) = var_solid.mesh.GetTypeElement(ne_s);
	    EntierToSendSolid(p)(nb++) = Nodle.GetM();
	    for (int j = 0; j < Nodle.GetM(); j++)
	      EntierToSendSolid(p)(nb++) = mesh_num_solid.GlobDofNumber_Subdomain(Nodle(j));
	  }
	else
	  {
	    EntierToRecvFluid(p).PushBack(ne_s);
	    EntierToRecvFluid(p).PushBack(num_loc);
	  }
      }

    Vector<int> NbEntierRecvSolid(nb_proc), NbEntierRecvFluid(nb_proc);
    NbEntierRecvFluid.Zero(); NbEntierRecvSolid.Zero();    
    Vector<MPI_Request> request_send(2*nb_proc), request_recv(2*nb_proc);
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbEntierFluid(i) > 0)
	  request_send(2*i) = MpiIsend(comm, &NbEntierFluid(i), xtmp, 1, i, 144);

	if (NbEntierSolid(i) > 0)
	  request_send(2*i+1) = MpiIsend(comm, &NbEntierSolid(i), xtmp, 1, i, 145);
	
	if (NbEntierSolid(i) > 0)
	  request_recv(2*i) = MpiIrecv(comm, &NbEntierRecvSolid(i), xtmp, 1, i, 144);

	if (NbEntierFluid(i) > 0)
	  request_recv(2*i+1) = MpiIrecv(comm, &NbEntierRecvFluid(i), xtmp, 1, i, 145);
      }
    
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbEntierFluid(i) > 0)
	  MPI_Wait(&request_send(2*i), &status);

	if (NbEntierSolid(i) > 0)
	  MPI_Wait(&request_send(2*i+1), &status);

	if (NbEntierSolid(i) > 0)
	  MPI_Wait(&request_recv(2*i), &status);
        
	if (NbEntierFluid(i) > 0)
	  MPI_Wait(&request_recv(2*i+1), &status);
      }
    
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbEntierFluid(i) > 0)
	  request_send(2*i) = MpiIsend(comm, EntierToSendFluid(i), xtmp, NbEntierFluid(i), i, 146);

	if (NbEntierSolid(i) > 0)
	  request_send(2*i+1) = MpiIsend(comm, EntierToSendSolid(i), xtmp, NbEntierSolid(i), i, 147);

	if (NbEntierSolid(i) > 0)
	  {
	    EntierToRecvSolid(i).Reallocate(NbEntierRecvSolid(i));
	    request_recv(2*i) = MpiIrecv(comm, EntierToRecvSolid(i), xtmp, NbEntierRecvSolid(i), i, 146);
	  }

	if (NbEntierFluid(i) > 0)
	  {
	    EntierToRecvFluid(i).Reallocate(NbEntierRecvFluid(i));
	    request_recv(2*i+1) = MpiIrecv(comm, EntierToRecvFluid(i), xtmp, NbEntierRecvFluid(i), i, 147);
	  }
      }

    for (int i = 0; i < nb_proc; i++)
      {
	if (NbEntierFluid(i) > 0)
	  MPI_Wait(&request_send(2*i), &status);

	if (NbEntierSolid(i) > 0)
	  MPI_Wait(&request_send(2*i+1), &status);

	if (NbEntierSolid(i) > 0)
	  MPI_Wait(&request_recv(2*i), &status);
        
	if (NbEntierFluid(i) > 0)
	  MPI_Wait(&request_recv(2*i+1), &status);
      }
    
    NbEntierRecvSolid.Zero();
    NbEntierRecvFluid.Zero();
    
    ConnecSolid.Reallocate(NumElementSolid.GetM());
    for (int i = 0; i < NumElementSolid.GetM(); i++)
      {
	int p = ProcElementSolid(i);
	int& nb = NbEntierRecvSolid(p);
	if (p != rank_proc)
	  {
	    int nb_dof = EntierToRecvSolid(p)(nb);
	    ConnecSolid(i).Reallocate(nb_dof + 8);
	    for (int j = 0; j < ConnecSolid(i).GetM(); j++)
	      ConnecSolid(i)(j) = EntierToRecvSolid(p)(nb++);	
	  }
	else
	  {
	    ConnecSolid(i).Reallocate(2);
	    ConnecSolid(i)(0) = EntierToRecvSolid(p)(nb++);
	    ConnecSolid(i)(1) = EntierToRecvSolid(p)(nb++);
	  }
      }

    ConnecFluid.Reallocate(NumElementFluid.GetM());
    for (int i = 0; i < NumElementFluid.GetM(); i++)
      {
	int p = ProcElementFluid(i);
	int& nb = NbEntierRecvFluid(p);
	if (p != rank_proc)
	  {
	    int nb_dof = EntierToRecvFluid(p)(nb);
	    ConnecFluid(i).Reallocate(nb_dof + 8);
	    for (int j = 0; j < ConnecFluid(i).GetM(); j++)
	      ConnecFluid(i)(j) = EntierToRecvFluid(p)(nb++);	
	  }
	else
	  {
	    ConnecFluid(i).Reallocate(2);
	    ConnecFluid(i)(0) = EntierToRecvFluid(p)(nb++);
	    ConnecFluid(i)(1) = EntierToRecvFluid(p)(nb++);
	  }
      }

    
    /* DISP(NumElementFluid); DISP(NumLocFluid); DISP(ProcElementFluid);
    DISP(NumElementSolid); DISP(NumLocSolid); DISP(ProcElementSolid);
    for (int i = 0; i < NumElementSolid.GetM(); i++)
      {
	DISP(i); DISP(ConnecSolid(i));
      }

    for (int i = 0; i < NumElementFluid.GetM(); i++)
      {
	DISP(i); DISP(ConnecFluid(i));
      }

    abort(); */
#endif
    
  }
  
  
  //! Retrieves informations about fluid element 
  /*!
    L'element fluide recherche est en face de l'element solide i
    Si le calcul est en parallele, cet element peut etre place sur un processeur different
    proc_fluid : processor where the fluid element is
    rot2 : orientation of the fluid element
    NodleAcous : degrees of freedom of the fluid element
    num_loc_acous : local position of the edge/face inside the fluid element
   */
  template<class Dimension>
  const ElementReference<Dimension, 1>&
  FluidStructureInteraction_Base<Dimension>::
  GetFluidData(int i, int& proc_fluid,
	       int& rot2, IVect& NodleAcous, int& num_loc_acous) const
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    if (nb_proc == 1)
      {
	proc_fluid = rank_proc;
	int ne_f = this->NumElementFluid(i);
	num_loc_acous = this->NumLocFluid(i);
	rot2 = var_fluid.mesh.Element(ne_f).GetOrientationBoundary(this->NumLocFluid(i));
	NodleAcous = var_fluid.GetDofNumberOnElement(ne_f);
	return var_fluid.GetReferenceElementH1(ne_f);
      }
    else
      {
	proc_fluid = this->ProcElementSolid(i);
	if (proc_fluid == rank_proc)
	  {
	    int ne_f = this->ConnecSolid(i)(0);
	    num_loc_acous = this->ConnecSolid(i)(1);
	    rot2 = var_fluid.mesh.Element(ne_f).GetOrientationBoundary(num_loc_acous);
	    NodleAcous = var_fluid.GetDofNumberOnElement(ne_f);
	    return var_fluid.GetReferenceElementH1(ne_f);
	  }
	else
	  {
	    num_loc_acous = this->ConnecSolid(i)(3);
	    rot2 = this->ConnecSolid(i)(4);
	    NodleAcous.Reallocate(this->ConnecSolid(i)(0));
	    for (int j = 0; j < NodleAcous.GetM(); j++)
	      NodleAcous(j) = this->ConnecSolid(i)(8+j);
	  }
      }
    
    return var_fluid.GetReferenceElementH1(0);
  }
  
  
  //! returns the orientation of the fluid element located in front of solid element i
  template<class Dimension>
  int FluidStructureInteraction_Base<Dimension>::GetFluidRotation(int i, int& proc_fluid) const
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    proc_fluid = rank_proc;
    if (nb_proc == 1)
      {
	int ne_f = this->NumElementFluid(i);
	return var_fluid.mesh.Element(ne_f).GetOrientationBoundary(this->NumLocFluid(i));
      }

    proc_fluid = this->ProcElementSolid(i);
    if (proc_fluid == rank_proc)
      {
	int ne_f = this->ConnecSolid(i)(0);
	int num_loc_acous = this->ConnecSolid(i)(1);
	return var_fluid.mesh.Element(ne_f).GetOrientationBoundary(num_loc_acous);
      }
    
    return this->ConnecSolid(i)(4);
  }
  

  //! initialization of physical values (rho_f, c^2, rho, tensor C)
  template<class Dimension>
  void FluidStructureInteraction_Base<Dimension>::InitIndices(int nb)
  {
    var_fluid.InitIndices(nb);
    var_solid.InitIndices(nb);
  }


  //! Constructs arrays needed for computation of finite element matrix
  template<class Dimension>
  void FluidStructureInteraction_Base<Dimension>
  ::ConstructAll(const string& input_file, const string& name_element, bool first)
  {
    // maximal number of different physical domains
    this->InitIndices(PhysicalConstant::nb_max_indices);
    var_fluid.SetTypeEquation("none");
    var_solid.SetTypeEquation("none");
    
    // data file is read
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(input_file, lines_data_file, var_fluid.comm_group_mode);
#else
    ReadLinesFile(input_file, lines_data_file);
#endif

    ReadInputFile(lines_data_file, *this);
    
    // mesh and finite element is constructed
    this->ComputeMeshAndFiniteElement(name_element, first);

    // other initializations
    var_fluid.PerformOtherInitializations();
    var_solid.PerformOtherInitializations();

    // mass and stiffness matrices are computed
    var_fluid.ComputeMassMatrix();    
    var_solid.ComputeMassMatrix();
    
    var_fluid.GetBoundaryConditionProblem().ComputeQuasiPeriodicPhase();
    var_solid.GetBoundaryConditionProblem().ComputeQuasiPeriodicPhase();
  }
  

  /*********************************
   * FluidStructureInteraction_Elas *
   *********************************/


  //! constructor
  template<class Complexe, class Dimension>
  FluidStructureInteraction_Elas<Complexe, Dimension>::
  FluidStructureInteraction_Elas(VarHarmonic_Base<Complexe, Dimension>& var_f,
				 VarHarmonic_Base<Complexe, Dimension>& var_s) 
    : FluidStructureInteraction_Base<Dimension>(var_f, var_s),
      var_fluid(var_f), var_solid(var_s)
  {
  }
  
  
  template<class Complexe, class Dimension>
  void FluidStructureInteraction_Elas<Complexe, Dimension>
  ::ConstructAll(const string& input_file, const string& name_element, bool first)
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    FluidStructureInteraction_Base<Dimension>::ConstructAll(input_file, name_element, first);
    const MeshNumbering<Dimension>& mesh_num_fluid = var_fluid.GetMeshNumbering(0);
    const MeshNumbering<Dimension>& mesh_num_solid = var_solid.GetMeshNumbering(0);
    
    NbPointsQuadFluid.Reallocate(nb_proc);
    NbPointsQuadSolid.Reallocate(nb_proc);
    NbPointsQuadSolid.Zero(); NbPointsQuadFluid.Zero();
    for (int i = 0; i < this->NumElementSolid.GetM(); i++)
      {
	int ne_s = this->NumElementSolid(i);
	int num_face = var_solid.mesh.Element(ne_s).numBoundary(this->NumLocSolid(i));
	int nb_pts = mesh_num_solid.GetNbPointsQuadratureBoundary(num_face);
	int proc = rank_proc;
	if (nb_proc > 1)
	  proc = this->ProcElementSolid(i);
	
	NbPointsQuadSolid(proc) += nb_pts;
      }

    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
      {
	int ne_f = this->NumElementFluid(i);
	int num_face = var_fluid.mesh.Element(ne_f).numBoundary(this->NumLocFluid(i));
	int nb_pts = mesh_num_fluid.GetNbPointsQuadratureBoundary(num_face);
	int proc = rank_proc;
	if (nb_proc > 1)
	  proc = this->ProcElementFluid(i);
	
	NbPointsQuadFluid(proc) += nb_pts;
      }
  }
  

  template<class Complexe, class Dimension>
  void FluidStructureInteraction_Elas<Complexe, Dimension>::
  SendUnFluid(Vector<VectComplexe>& UnQuad)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif
    
    if (nb_proc == 1)
      return;

#ifdef SELDON_WITH_MPI
    Vector<VectComplexe> UnFluid(UnQuad);

    Vector<int64_t> xtmp;
    Vector<MPI_Request> request_send(nb_proc), request_recv(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbPointsQuadFluid(i) > 0)
	  request_send(i) = MpiIsend(comm, UnFluid(i), xtmp, UnFluid(i).GetM(), i, 152);
	
	UnQuad(i).Reallocate(NbPointsQuadSolid(i));
	if (NbPointsQuadSolid(i) > 0)
	  request_recv(i) = MpiIrecv(comm, UnQuad(i), xtmp, UnQuad(i).GetM(), i, 152);	
      }
    
    MPI_Status status;
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbPointsQuadFluid(i) > 0)
	  MPI_Wait(&request_send(i), &status);

	if (NbPointsQuadSolid(i) > 0)
	  MPI_Wait(&request_recv(i), &status);
      }    
#endif

  }


  template<class Complexe, class Dimension>
  void FluidStructureInteraction_Elas<Complexe, Dimension>::
  SendUnSolid(Vector<VectComplexe>& UnQuad)
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
#else
    int nb_proc(1);
#endif

    if (nb_proc == 1)
      return;

#ifdef SELDON_WITH_MPI
    Vector<VectComplexe> UnSolid(UnQuad);
    
    Vector<int64_t> xtmp;
    Vector<MPI_Request> request_send(nb_proc), request_recv(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbPointsQuadSolid(i) > 0)
	  request_send(i) = MpiIsend(comm, UnSolid(i), xtmp, UnSolid(i).GetM(), i, 151);
	
	UnQuad(i).Reallocate(NbPointsQuadFluid(i));
	if (NbPointsQuadFluid(i) > 0)
	  request_recv(i) = MpiIrecv(comm, UnQuad(i), xtmp, UnQuad(i).GetM(), i, 151);	
      }
    
    MPI_Status status;
    for (int i = 0; i < nb_proc; i++)
      {
	if (NbPointsQuadSolid(i) > 0)
	  MPI_Wait(&request_send(i), &status);

	if (NbPointsQuadFluid(i) > 0)
	  MPI_Wait(&request_recv(i), &status);
      }    
#endif

  }
  

  //! computing ProdPn = beta*ProdPn + alpha*B*Un
  template<class Complexe, class Dimension>
  void FluidStructureInteraction_Elas<Complexe, Dimension>::
  ApplyCouplingOperatorB(const Complexe& alpha, const VectComplexe& Un,
                         const Complexe& beta, VectComplexe& ProdPn)
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    Complexe zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    if (beta == zero)
      ProdPn.Fill(0);
    else if (beta != one)
      Mlt(beta, ProdPn);
    
    const Mesh<Dimension>& mesh = var_fluid.mesh;
    const MeshNumbering<Dimension>& mesh_num = var_fluid.GetMeshNumbering(0);
    const MeshNumbering<Dimension>& mesh_num_solid = var_solid.GetMeshNumbering(0);
    
    const Mesh<Dimension>& mesh_elas = var_solid.mesh;
    
    TinyVector<VectComplexe, Dimension::dim_N> UnSurf;
    
    TinyVector<Complexe, Dimension::dim_N> vec_u;
    TinyVector<VectComplexe, Dimension::dim_N> UnLoc;
    TinyVector<VectComplexe, 1> PnLoc;

    // variables storing u \cdot n for parallel execution
    Vector<VectComplexe> UnQuad(nb_proc);
    Vector<int> NumQuadP(nb_proc);
    NumQuadP.Zero();
    for (int i = 0; i < nb_proc; i++)
      UnQuad(i).Reallocate(NbPointsQuadSolid(i));
    
    int proc;
    for (int i = 0; i < this->NumElementSolid.GetM(); i++)
      {
	int ne_s = this->NumElementSolid(i);
	const ElementReference<Dimension, 1>& Fb = var_solid.GetReferenceElementH1(ne_s);
	
	var_solid.GetLocalUnknownVector(Un, ne_s, UnLoc);
	
	int num_face = mesh_elas.Element(ne_s).numBoundary(this->NumLocSolid(i));
	int r = mesh_num_solid.GetOrderQuadrature(num_face);
	int nb_pts = mesh_num_solid.GetNbPointsQuadratureBoundary(num_face);

	num_face = mesh_elas.Element(ne_s).numBoundary(this->NumLocSolid(i));
	for (int k = 0; k < Dimension::dim_N; k++)
	  {
	    UnSurf(k).Reallocate(nb_pts);
	    UnSurf(k).Zero();
	    Fb.ApplyShTranspose(this->NumLocSolid(i), UnLoc(k), UnSurf(k), r);
	  }

	int rot1 = mesh_elas.Element(ne_s).GetOrientationBoundary(this->NumLocSolid(i));
	int rot2 = this->GetFluidRotation(i, proc);
	int nv = mesh_elas.Boundary(num_face).GetNbVertices();
	int rot = mesh_num.GetRotationFace(rot1, rot2, nv);
	const VectReal_wp& PoidsFlux = mesh_num_solid.number_map.
	  GetFluxWeight(r, mesh_elas.Boundary(num_face));
	
	const Matrix<int>& FacesQuadRotation = mesh_num_solid.number_map.
	  GetRotationQuadraturePoints(r, mesh_elas.Boundary(num_face));

	for (int k = 0; k < nb_pts; k++)
	  {
	    Real_wp poids = 2.0*var_solid.Glob_dsj(num_face)(k)*PoidsFlux(k);
	    R_N normale = var_solid.Glob_normale(num_face)(k);
	    CopyVector(UnSurf, k, vec_u);
	    UnSurf(0)(k) = poids*DotProd(vec_u, normale);
	  }
	
	for (int k = 0; k < nb_pts; k++)
	  UnQuad(proc)(NumQuadP(proc) + FacesQuadRotation(rot, k)) = UnSurf(0)(k);
	
	NumQuadP(proc) += nb_pts;
      }

    this->SendUnSolid(UnQuad);
    
    NumQuadP.Zero();
    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
      {
	int ne_f = this->NumElementFluid(i);
	const ElementReference<Dimension, 1>& Fb_acous = var_fluid.GetReferenceElementH1(ne_f);
	int num_face = mesh.Element(ne_f).numBoundary(this->NumLocFluid(i));
	int r = mesh_num.GetOrderQuadrature(num_face);
	int nb_pts = mesh_num.GetNbPointsQuadratureBoundary(num_face);
	int proc = rank_proc;
	if (nb_proc > 1)
	  proc = this->ProcElementFluid(i);
	
	UnSurf(0).Reallocate(nb_pts);
	for (int j = 0; j < nb_pts; j++)
	  UnSurf(0)(j) = UnQuad(proc)(NumQuadP(proc) + j);
	
	NumQuadP(proc) += nb_pts;
	PnLoc(0).Reallocate(Fb_acous.GetNbDof());
	PnLoc(0).Zero();
	
	Fb_acous.ApplySh(alpha, this->NumLocFluid(i), UnSurf(0), PnLoc(0), r);
    
	var_fluid.AddLocalUnknownVector(one, PnLoc, ne_f, ProdPn);	
      }
  }
  
  
  //! computing ProdUn = beta*ProdUn + alpha*B^T Pn
  template<class Complexe, class Dimension>
  void FluidStructureInteraction_Elas<Complexe, Dimension>::
  ApplyCouplingTransposeB(const Complexe& alpha, const VectComplexe& Pn,
                          const Complexe& beta, VectComplexe& ProdUn)
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int nb_proc(1), rank_proc(0);
#endif

    Complexe zero, one;
    SetComplexZero(zero); SetComplexOne(one);
    if (beta == zero)
      ProdUn.Zero();
    else if (beta != one)
      Mlt(beta, ProdUn);
    
    const Mesh<Dimension>& mesh_acous = var_fluid.mesh;
    const Mesh<Dimension>& mesh = var_solid.mesh;
    const MeshNumbering<Dimension>& mesh_num = var_solid.GetMeshNumbering(0);
    const MeshNumbering<Dimension>& mesh_num_fluid = var_fluid.GetMeshNumbering(0);

    TinyVector<Complexe, Dimension::dim_N> vec_u;
    TinyVector<VectComplexe, Dimension::dim_N> PnSurf;
    VectComplexe PnRot;    
    TinyVector<VectComplexe, Dimension::dim_N> UnLoc;
    TinyVector<VectComplexe, 1> PnLoc;

    // variables storing p for parallel execution
    Vector<VectComplexe> PnQuad(nb_proc);
    Vector<int> NumQuadP(nb_proc);
    NumQuadP.Zero();
    for (int i = 0; i < nb_proc; i++)
      PnQuad(i).Reallocate(NbPointsQuadFluid(i));

    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
      {
	int ne_f = this->NumElementFluid(i);
	const ElementReference<Dimension, 1>& Fb = var_fluid.GetReferenceElementH1(ne_f);

	var_fluid.GetLocalUnknownVector(Pn, ne_f, PnLoc);
	
	int num_face = mesh_acous.Element(ne_f).numBoundary(this->NumLocFluid(i));
	int r = mesh_num_fluid.GetOrderQuadrature(num_face);
	int nb_pts = mesh_num_fluid.GetNbPointsQuadratureBoundary(num_face);
		
	PnRot.Reallocate(nb_pts);
	PnRot.Zero();
	Fb.ApplyShTranspose(this->NumLocFluid(i), PnLoc(0), PnRot, r);
	
	int proc = rank_proc;
	if (nb_proc > 1)
	  proc = this->ProcElementFluid(i);
	
	for (int j = 0; j < nb_pts; j++)
	  PnQuad(proc)(NumQuadP(proc) + j) = PnRot(j);
	
	NumQuadP(proc) += nb_pts;
      }

    this->SendUnFluid(PnQuad);
    
    NumQuadP.Zero();
    for (int i = 0; i < this->NumElementSolid.GetM(); i++)
      {
	int ne_s = this->NumElementSolid(i);
	const ElementReference<Dimension, 1>& Fb_elas = var_solid.GetReferenceElementH1(ne_s);
	int num_face = mesh.Element(ne_s).numBoundary(this->NumLocSolid(i));
    
	int r = mesh_num.GetOrderQuadrature(num_face);
	int nb_pts = mesh_num.GetNbPointsQuadratureBoundary(num_face);
    
	int rot1 = mesh.Element(ne_s).GetOrientationBoundary(this->NumLocSolid(i));
	int proc;
	int nv = mesh.Boundary(num_face).GetNbVertices();
	int rot2 = this->GetFluidRotation(i, proc);
	int rot = mesh_num.GetRotationFace(rot1, rot2, nv);
	const VectReal_wp& PoidsFlux = mesh_num.number_map.
	  GetFluxWeight(r, mesh.Boundary(num_face));
	
	const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	  GetRotationQuadraturePoints(r, mesh.Boundary(num_face));
	
	for (int k = 0; k < Dimension::dim_N; k++)
	  {
	    PnSurf(k).Reallocate(nb_pts);
	    PnSurf(k).Zero();
	  }
	
	for (int k = 0; k < nb_pts; k++)
	  PnSurf(0)(k) = PnQuad(proc)(NumQuadP(proc) + FacesQuadRotation(rot, k));

	NumQuadP(proc) += nb_pts;
	for (int k = 0; k < nb_pts; k++)
	  {
	    Real_wp poids = 2.0*var_solid.Glob_dsj(num_face)(k)*PoidsFlux(k);
	    R_N normale = var_solid.Glob_normale(num_face)(k);
	    vec_u = poids*normale*PnSurf(0)(k);
	    CopyVector(vec_u, k, PnSurf);
	  }
	
	for (int k = 0; k < Dimension::dim_N; k++)
	  {
	    UnLoc(k).Reallocate(Fb_elas.GetNbDof());
	    UnLoc(k).Zero();
	    Fb_elas.ApplySh(alpha, this->NumLocSolid(i), PnSurf(k), UnLoc(k), r);
	  }
	
	var_solid.AddLocalUnknownVector(one, UnLoc, ne_s, ProdUn);
      }
  }
  

  //! adding coupling operator alpha B to matrix A
  template<class Complexe, class Dimension>
  void FluidStructureInteraction_Elas<Complexe, Dimension>::
  AddOperatorB(const Complexe& alpha, VirtualMatrix<Complexe>& mat_sp)
  {
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = var_fluid.comm_group_mode;
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
#else
    int rank_proc(0);
#endif

    Complexe one; SetComplexOne(one);
    const Mesh<Dimension>& mesh_elas = var_solid.mesh;
    const MeshNumbering<Dimension>& mesh_elas_num = var_solid.GetMeshNumbering(0);

    int Nsol = var_solid.GetNbDof();
    int Nsol_glob = var_solid.GetNbGlobalDof();
    int Nscal = mesh_elas_num.GetNbDof();
    // int Nscal_glob = var_solid.GetNbGlobalMeshDof();
    
    for (int i = 0; i < this->NumElementSolid.GetM(); i++)
      {
	int ne_s = this->NumElementSolid(i);
	const ElementReference<Dimension, 1>& Fb_elas = var_solid.GetReferenceElementH1(ne_s);
	int num_face = mesh_elas.Element(ne_s).numBoundary(this->NumLocSolid(i));
	
	int r = mesh_elas_num.GetOrderQuadrature(num_face);
	int nb_pts = mesh_elas_num.GetNbPointsQuadratureBoundary(num_face);
	
	int rot1 = mesh_elas.Element(ne_s).GetOrientationBoundary(this->NumLocSolid(i));
	int nv = mesh_elas.Boundary(num_face).GetNbVertices();
	
	int rot2 = -1; IVect NodleAcous; int num_loc_acous = -1;
	int proc_fluid = rank_proc;
	const ElementReference<Dimension, 1>& Fb_acous
	  = this->GetFluidData(i, proc_fluid,
			       rot2, NodleAcous, num_loc_acous);
	
	int rot = mesh_elas_num.GetRotationFace(rot1, rot2, nv);
	const VectReal_wp& PoidsFlux = mesh_elas_num.number_map.
	  GetFluxWeight(r, mesh_elas.Boundary(num_face));
	
	const Matrix<int>& FacesQuadRotation = mesh_elas_num.number_map.
	  GetRotationQuadraturePoints(r, mesh_elas.Boundary(num_face));
	
	VectReal_wp phi(Fb_elas.GetNbDof()), val_phi(nb_pts);
	TinyVector<VectComplexe, Dimension::dim_N> feval;
	TinyVector<Complexe, Dimension::dim_N> vec_u;
	VectComplexe contrib(Fb_acous.GetNbDof());
	contrib.Zero();
	phi.Zero(); val_phi.Zero();
	
	for (int k = 0; k < Dimension::dim_N; k++)
	  {
	    feval(k).Reallocate(nb_pts);
	    feval(k).Zero();
	  }
	
	VectComplexe val_interac(Fb_acous.GetNbDof());
	int c_interac;
	IVect col_interac(Fb_acous.GetNbDof());
	val_interac.Zero(); col_interac.Fill(-1);
	
	//DISP(i); DISP(ne_s); DISP(ne_f);     DISP(num_face);
	// loop over dofs for elastic function test
	for (int iloc = 0; iloc < Fb_elas.GetNbDof(); iloc++)
	  {
	    phi(iloc) = 1.0;
	    int num_dof = mesh_elas_num.Element(ne_s).GetNumberDof(iloc);
	    
	    Fb_elas.ApplyShTranspose(this->NumLocSolid(i), phi, val_phi, r);
	    bool null_fct = true;
	    for (int k = 0; k < nb_pts; k++)
	      if (abs(val_phi(k)) > epsilon_machine)
		null_fct = false;
	    
	    if (!null_fct)
	      {
		for (int k = 0; k < nb_pts; k++)
		  {
		    Real_wp poids = 2.0*var_solid.Glob_dsj(num_face)(k)*PoidsFlux(k)*val_phi(k);
		    vec_u = alpha*var_solid.Glob_normale(num_face)(k);
		    int krot = FacesQuadRotation(rot, k);
		    vec_u *= poids;
		    CopyVector(vec_u, krot, feval);
		  }
		
		for (int k = 0; k < Dimension::dim_N; k++)
		  {
		    contrib.Zero();		
		    Fb_acous.ApplySh(one, num_loc_acous, feval(k), contrib, r);
		    c_interac = 0;
		    for (int j = 0; j < Fb_acous.GetNbDof(); j++)
		      if (abs(contrib(j)) > var_solid.GetThresholdMatrix())
			{
			  col_interac(c_interac) = Nsol + NodleAcous(j);
			  if (proc_fluid != rank_proc)
			    col_interac(c_interac) = Nsol_glob + NodleAcous(j);
			  
			  val_interac(c_interac) = contrib(j);
			  if (proc_fluid == rank_proc)
			    mat_sp.AddInteraction(col_interac(c_interac), k*Nscal + num_dof,
						  val_interac(c_interac));
			  else
			    {
			      mat_sp.AddRowDistantInteraction(col_interac(c_interac),
							      k*Nscal + num_dof, proc_fluid,
							      val_interac(c_interac));
			      
			      mat_sp.AddDistantInteraction(k*Nscal + num_dof, col_interac(c_interac),
							   proc_fluid, val_interac(c_interac));
			    }
			  
			  c_interac++;
			}
		    
		    if (proc_fluid == rank_proc)
		      mat_sp.AddInteractionRow(k*Nscal + num_dof, c_interac, col_interac, val_interac);
		  }
	      }
	    
	    phi(iloc) = Real_wp(0);
	  }
      }
  }
  
  
  /*************************************
   * FluidStructureInteractionHarmonic *
   *************************************/


  //! Default constructor
  template<class TypeEqFluid, class TypeEqSolid>
  FluidStructureInteractionHarmonic<TypeEqFluid, TypeEqSolid>::FluidStructureInteractionHarmonic()
    : FluidStructureInteraction_Elas<Complexe, Dimension>(var_fluid, var_solid)
  {
  }
  

  //! keywords associated with fluid structure interaction
  template<class TypeEqFluid, class TypeEqSolid>
  void FluidStructureInteractionHarmonic<TypeEqFluid, TypeEqSolid>::
  SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (!keyword.compare("ReferenceFluid"))
      {
        // providing all the references for the fluid
        this->ref_fluid.Reallocate(param.GetM());
        for (int i = 0; i < this->ref_fluid.GetM(); i++)
          this->ref_fluid(i) = to_num<int>(param(i));
      }
    else if (keyword.find("Solid") == 0)
      {
        // for all keywords beginning with Solid
        // we call the method SetInputData of var_solid
        string keyword2 = keyword.substr(5);
        var_solid.SetInputData(keyword2, param);
      }
    else if (keyword.find("Fluid") == 0)
      {
        // for all keywords beginning with Fluid
        // we call the method SetInputData of var_fluid
        string keyword2 = keyword.substr(5);
        var_fluid.SetInputData(keyword2, param);
      }
    else
      {
        var_fluid.SetInputData(keyword, param);
        var_solid.SetInputData(keyword, param);
      }        
  }


  template<class TypeEqFluid, class TypeEqSolid>
  void FluidStructureInteractionHarmonic<TypeEqFluid, TypeEqSolid>
  ::RunAll(const string& input_file, const string& name_element)
  {
    FluidStructureInteraction_Elas<Complexe, Dimension>::ConstructAll(input_file, name_element, false);

    DistributedMatrix<Complexe, Symmetric, ArrayRowSymSparse> mat_sp;
    
    int Np = var_fluid.GetNbDof();
    int Nu = var_solid.GetNbDof();
    int N = Np + Nu;
    
    mat_sp.Reallocate(N, N);
    
    GlobalGenericMatrix<Complexe> nat_mat;

    var_solid.AddMatrixWithBC(mat_sp, nat_mat, 0, 0);
    var_fluid.AddMatrixWithBC(mat_sp, nat_mat, Nu, Nu);
        
    Complexe coef;
    SetComplexReal(var_fluid.GetOmega(), coef);
    this->AddOperatorB(coef, mat_sp);

#ifdef SELDON_WITH_MPI
    Vector<int> GlobalNum, OverlapNum, OverlapProc;
    Vector<int> MatchingProc;
    Vector<IVect> MatchingDofNumber;
    int nb_proc; MPI_Comm_size(var_fluid.comm_group_mode, &nb_proc);

    const MeshNumbering<Dimension>& mesh_num_fluid = var_fluid.GetMeshNumbering(0);
    const MeshNumbering<Dimension>& mesh_num_solid = var_solid.GetMeshNumbering(0);        
    // arrays needed for the distributed matrix are constructed
    if (nb_proc > 1)
      {
	MPI_Comm& comm = var_solid.comm_group_mode;
	int Np_glob = var_fluid.GetNbGlobalDof();
	int Nu_glob = var_solid.GetNbGlobalDof();
        
        IVect& MatchingNumber_Solid = var_solid.GetProcMatchingNeighbor();
        Vector<IVect>& MatchingDofOrig_Solid = var_solid.GetOriginalMatchingDofNeighbor();
        IVect& MatchingNumber_Fluid = var_fluid.GetProcMatchingNeighbor();
        Vector<IVect>& MatchingDofOrig_Fluid = var_fluid.GetOriginalMatchingDofNeighbor();

	GlobalNum.Reallocate(Np+Nu);
	for (int i = 0; i < Nu; i++)
	  GlobalNum(i) = mesh_num_solid.GlobDofNumber_Subdomain(i);
	
	for (int i = 0; i < Np; i++)
	  GlobalNum(Nu + i) = Nu_glob + mesh_num_fluid.GlobDofNumber_Subdomain(i);
	
	OverlapNum.Reallocate(var_fluid.GetNbOverlappedDof() + var_solid.GetNbOverlappedDof());
	for (int i = 0; i < var_solid.GetNbOverlappedDof(); i++)
	  OverlapNum(i) = var_solid.GetOverlappedDofNumber(i);
	
	int offset = var_solid.GetNbOverlappedDof();
	for (int i = 0; i < var_fluid.GetNbOverlappedDof(); i++)
	  OverlapNum(offset+i) = Nu + var_fluid.GetOverlappedDofNumber(i);
	
	OverlapProc.Reallocate(OverlapNum.GetM());
	for (int i = 0; i < var_solid.GetNbOverlappedDof(); i++)
	  OverlapProc(i) = var_solid.GetOverlappedProcNumber(i);
	
	for (int i = 0; i < var_fluid.GetNbOverlappedDof(); i++)
	  OverlapProc(offset+i) = var_fluid.GetOverlappedProcNumber(i);
	
	int Nvol = mesh_num_solid.GetNbDof();
	int i0 = 0;	

	for (int i = 0; i < MatchingNumber_Solid.GetM(); i++)
	  {
	    int proc = MatchingNumber_Solid(i);
	    while ((i0 < MatchingNumber_Fluid.GetM())
		   && (MatchingNumber_Fluid(i0) < proc))
	      {
		int p = MatchingNumber_Fluid(i0);
		MatchingProc.PushBack(p);
		IVect DofNumber(MatchingDofOrig_Fluid(i0));
		for (int j = 0; j < DofNumber.GetM(); j++)
		  DofNumber(j) += Nu;
		
		MatchingDofNumber.PushBack(DofNumber);
		i0++;
	      }
	    
	    offset = MatchingDofOrig_Solid(i).GetM();
	    IVect DofNumber(Dimension::dim_N*offset);
	    for (int j = 0; j < offset; j++)
	      {
		DofNumber(j) = MatchingDofOrig_Solid(i)(j);
		DofNumber(offset+j) = Nvol + DofNumber(j);
		if (Dimension::dim_N == 3)
		  DofNumber(2*offset+j) = 2*Nvol + DofNumber(j);
	      }

	    MatchingProc.PushBack(proc);
	    
	    if ((i0 < MatchingNumber_Fluid.GetM())
		   && (MatchingNumber_Fluid(i0) == proc))
	      {
		offset = DofNumber.GetM();
		DofNumber.Resize(DofNumber.GetM() + MatchingDofOrig_Fluid(i0).GetM());
		for (int j = 0; j < MatchingDofOrig_Fluid(i0).GetM(); j++)
		  DofNumber(offset + j) = Nu + MatchingDofOrig_Fluid(i0)(j);
		
		i0++;
	      }	    
	    
	    MatchingDofNumber.PushBack(DofNumber);
	  }

	while (i0 < MatchingNumber_Fluid.GetM())
	  {
	    int p = MatchingNumber_Fluid(i0);
	    MatchingProc.PushBack(p);
	    IVect DofNumber(MatchingDofOrig_Fluid(i0));
	    for (int j = 0; j < DofNumber.GetM(); j++)
	      DofNumber(j) += Nu;
	    
	    MatchingDofNumber.PushBack(DofNumber);
	    i0++;
	  }
	
	mat_sp.Init(Np_glob+Nu_glob, &GlobalNum, &OverlapNum, &OverlapProc,
		    Np+Nu, 1, &MatchingProc, &MatchingDofNumber, comm);
      }
#endif
    
    All_MatrixLU<Complexe> mat_lu;
    
    mat_lu.ShowMessages();
    mat_lu.Factorize(mat_sp);
    
    Vector<Complexe> rhs_elas(Nu);
    Vector<Complexe> rhs_acous(Np);  
    Vector<Complexe> rhs(N);
        
    var_fluid.ComputeRightHandSide(rhs_acous);
    var_solid.ComputeRightHandSide(rhs_elas);
    
    for (int i = 0; i < Nu; i++)
      rhs(i) = rhs_elas(i);
    
    for (int i = 0; i < Np; i++)
      rhs(Nu+i) = rhs_acous(i);
    
    mat_lu.Solve(rhs);
    
    for (int i = 0; i < Nu; i++)
      rhs_elas(i) = rhs(i);
    
    for (int i = 0; i < Np; i++)
      rhs_acous(i) = rhs(Nu+i);
    
    var_fluid.WriteDatas(rhs_acous);
    var_solid.WriteDatas(rhs_elas); 
  }
  
    
  /*****************************
   * FluidStructureInteraction *
   *****************************/


  //! keywords associated with fluid structure interaction
  template<class TypeEqFluid, class TypeEqSolid>
  void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::
  SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (!keyword.compare("ReferenceFluid"))
      {
        // providing all the references for the fluid
        this->ref_fluid.Reallocate(param.GetM());
        for (int i = 0; i < this->ref_fluid.GetM(); i++)
          this->ref_fluid(i) = to_num<int>(param(i));
      }
    else if (keyword.find("Solid") == 0)
      {
        // for all keywords beginning with Solid
        // we call the method SetInputData of var_solid
        string keyword2 = keyword.substr(5);
        var_solid.SetInputData(keyword2, param);
      }
    else if (keyword.find("Fluid") == 0)
      {
        // for all keywords beginning with Fluid
        // we call the method SetInputData of var_fluid
        string keyword2 = keyword.substr(5);
        var_fluid.SetInputData(keyword2, param);
      }
    else
      {
        var_fluid.SetInputData(keyword, param);
        var_solid.SetInputData(keyword, param);
      }        
  }
  
  
  //! construction of mesh, finite elements, matrices so that time-scheme can be launched
  template<class TypeEqFluid, class TypeEqSolid>
  void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>
  ::ConstructAll(const string& input_file, const string& name_element)
  {
    FluidStructureInteraction_Elas<Real_wp, Dimension>::ConstructAll(input_file, name_element);
    
    // right hand sides are computed
    var_fluid.ComputeRightHandSide();
    var_solid.ComputeRightHandSide();
    
    var_fluid.ComputeMassMatrix();
    var_solid.ComputeMassMatrix();
    
    var_fluid.ComputeStiffnessMatrix();
    var_solid.ComputeStiffnessMatrix();
    
    // initialization before time iterations
    var_fluid.CheckIdentityJacobianInPML();
    var_solid.CheckIdentityJacobianInPML();
    var_fluid.InitTimeIterations();
    var_solid.InitTimeIterations();

    int Np = var_fluid.GetNbScalarUnknowns();
    int Nv = var_fluid.GetNbVectorialUnknowns();
    int Nu = var_solid.GetNbScalarUnknowns();
    int Ns = var_solid.GetNbVectorialUnknowns();
    
    // total number of unknowns
    this->nodl = Np + Ns + Nv + Nu;
    
  }
  
  
  //! computation of Y = g(t, X) (evolution system can be written under the form dX/dt = g(t, X))
  template<class TypeEqFluid, class TypeEqSolid>
  void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
			     const VectReal_wp& X, VectReal_wp& Y, bool invert, bool source)
  {
    int Np = var_fluid.GetNbScalarUnknowns();
    int Nv = var_fluid.GetNbVectorialUnknowns();
    int Nu = var_solid.GetNbScalarUnknowns();
    int Ns = var_solid.GetNbVectorialUnknowns();
    
    VectReal_wp Pn, Vn, Un, Sn;
    VectReal_wp ProdPn, ProdVn, ProdUn, ProdSn;
    
    Real_wp* x_ptr = X.GetData();
    Real_wp* y_ptr = Y.GetData();
    Y.Fill(0);
    
    Pn.SetData(Np, x_ptr); x_ptr += Np; ProdPn.SetData(Np, y_ptr); y_ptr += Np;
    Vn.SetData(Nv, x_ptr); x_ptr += Nv; ProdVn.SetData(Nv, y_ptr); y_ptr += Nv;
    Un.SetData(Nu, x_ptr); x_ptr += Nu; ProdUn.SetData(Nu, y_ptr); y_ptr += Nu;
    Sn.SetData(Ns, x_ptr); x_ptr += Ns; ProdSn.SetData(Ns, y_ptr); y_ptr += Ns;
    
    var_fluid.AddScalarSourceAtTime(1.0, tn, nb_deriv, ProdPn);
    this->ApplyStiffnessPressure(1.0, Vn, 1.0, ProdPn);
    this->ApplyCouplingOperatorB(1.0, Un, 1.0, ProdPn);
    var_fluid.ApplyOperatorSh(-1.0, tn, Pn, 1.0, ProdPn);
    if (invert)
      this->SolveMassPressure(ProdPn);

    this->ApplyStiffnessDisplacementFluid(-1.0, Pn, 1.0, ProdVn);
    if (invert)
      this->SolveMassDisplacementFluid(ProdVn);

    if (source)
      var_solid.AddScalarSourceAtTime(1.0, tn, nb_deriv, ProdUn);
    
    this->ApplyStiffnessDisplacementSolid(-1.0, Sn, 1.0, ProdUn);
    this->ApplyCouplingTransposeB(-1.0, Pn, 1.0, ProdUn);
    var_solid.ApplyOperatorSh(-1.0, tn, Un, 1.0, ProdUn);
    if (invert)
      this->SolveMassDisplacementSolid(ProdUn);
    
    this->ApplyStiffnessStress(1.0, Un, 1.0, ProdSn);
    if (invert)
      this->SolveMassStress(ProdSn);
    
    var_fluid.var_harmonic.ImposeNullDirichletCondition(ProdPn);
    var_solid.var_harmonic.ImposeNullDirichletCondition(ProdUn);
    
    Pn.Nullify(); Un.Nullify(); Vn.Nullify(); Sn.Nullify();
    ProdPn.Nullify(); ProdUn.Nullify(); ProdVn.Nullify(); ProdSn.Nullify();
  }
  
  
  template<class TypeEqFluid, class TypeEqSolid>
  void FluidStructureInteraction<TypeEqFluid, TypeEqSolid>::
  EvaluateFunction(const Real_wp& tn, const VectReal_wp& X, VectReal_wp& Y,
                   bool invert_mass, bool source)
  {
    this->EvaluateDerivativeFunction(tn, 0, X, Y, invert_mass, source);
  }


  /********************************
   * FluidStructureLeapFrogScheme *
   ********************************/
  
  
  //! null initial conditions for leap-frog scheme
  template<class GenericPb>
  void FluidStructureLeapFrogScheme
  ::SetInitialCondition(const Real_wp& t0, const Real_wp& dt_, GenericPb& var)
  {
    dt = dt_;
    
    int Np = var.var_fluid.GetNbScalarUnknowns();
    int Nv = var.var_fluid.GetNbVectorialUnknowns();
    int Nu = var.var_solid.GetNbScalarUnknowns();
    int Ns = var.var_solid.GetNbVectorialUnknowns();

    Un_half.Reallocate(Nu); ProdUn.Reallocate(Nu);
    Vn_half.Reallocate(Nv); ProdVn.Reallocate(Nv);
    Pn.Reallocate(Np); ProdPn.Reallocate(Np);
    Sn.Reallocate(Ns); ProdSn.Reallocate(Ns);
    
    Un_half.Fill(0); ProdUn.Fill(0);
    Vn_half.Fill(0); ProdVn.Fill(0);
    Pn.Fill(0); ProdPn.Fill(0);
    Sn.Fill(0); ProdSn.Fill(0);
  }
  
  
  //! basic leap-frog scheme
  template<class GenericPb>
  void FluidStructureLeapFrogScheme::Advance(const Real_wp& t, int n, GenericPb& var)
  {
    // U^n+1/2 = U^n-1/2 + dt M_u^-1 (F_u^n - R_u Sigma^n - B* P^n)
    var.ApplyStiffnessDisplacementSolid(-1.0, Sn, 0.0, ProdUn);    
    var.ApplyCouplingTransposeB(-1.0, Pn, 1.0, ProdUn);
    var.var_solid.AddScalarSourceAtTime(1.0, t, 0, ProdUn);
    var.var_solid.ApplyOperatorDhMinusdtSh(1.0, t, Un_half, dt, ProdUn);
     var.var_solid.SolveOperatorDhPlusdtSh(ProdUn);
    Copy(ProdUn, Un_half);
    var.var_solid.var_harmonic.ImposeNullDirichletCondition(Un_half);
    
    // Sigma^n+1 = Sigma^n + dt M_sigma^-1 R_sigma U^n+1/2
    var.ApplyStiffnessStress(1.0, Un_half, 0.0, ProdSn);
    var.SolveMassStress(ProdSn);
    Add(dt, ProdSn, Sn);
    
    // P^n+1 = P^n + dt M_p^-1 ( F_p^n+1/2 + R_p V^n+1/2 + B U^n+1/2)
    var.ApplyStiffnessPressure(1.0, Vn_half, 0.0, ProdPn);
    var.ApplyCouplingOperatorB(1.0, Un_half, 1.0, ProdPn);
    var.var_fluid.AddScalarSourceAtTime(1.0, t+0.5*dt, 0, ProdPn);
    var.var_fluid.ApplyOperatorDhMinusdtSh(1.0, t+0.5*dt, Pn, dt, ProdPn);
    var.var_fluid.SolveOperatorDhPlusdtSh(ProdPn);
    Copy(ProdPn, Pn);
    var.var_fluid.var_harmonic.ImposeNullDirichletCondition(Pn);
    
    // V^n+3/2 = V^n+1/2 - dt M_v^-1 R_v P^n+1
    var.ApplyStiffnessDisplacementFluid(-1.0, Pn, 0.0, ProdVn);
    var.SolveMassDisplacementFluid(ProdVn);
    Add(dt, ProdVn, Vn_half);    
  }


  /*******************************
   * AcousticAcousticInteraction *
   *******************************/


  //! keywords associated with fluid structure interaction
  template<class Dimension>
  void AcousticAcousticInteraction<Dimension>::
  SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (!keyword.compare("ReferenceFluid"))
      {
        // providing all the references for the fluid
        this->ref_fluid.Reallocate(param.GetM());
        for (int i = 0; i < this->ref_fluid.GetM(); i++)
          this->ref_fluid(i) = to_num<int>(param(i));
      }
    else if (keyword.find("Solid") == 0)
      {
        // for all keywords beginning with Solid
        // we call the method SetInputData of var_solid
        string keyword2 = keyword.substr(5);
        var_solid.SetInputData(keyword2, param);
      }
    else if (keyword.find("Fluid") == 0)
      {
        // for all keywords beginning with Fluid
        // we call the method SetInputData of var_fluid
        string keyword2 = keyword.substr(5);
        var_fluid.SetInputData(keyword2, param);
      }
    else
      {
        var_fluid.SetInputData(keyword, param);
        var_solid.SetInputData(keyword, param);
      }        
  }
  
  
  //! construction of mesh, finite elements, matrices so that time-scheme can be launched
  template<class Dimension>
  void AcousticAcousticInteraction<Dimension>
  ::ConstructAll(const string& input_file, const string& name_element)
  {
    // maximal number of different physical domains
    this->InitIndices(PhysicalConstant::nb_max_indices);
    
    // data file is read
    Vector<string> lines_data_file;
    
#ifdef SELDON_WITH_MPI
    ReadLinesFile(input_file, lines_data_file, var_fluid.var_harmonic.comm_group_mode);
#else
    ReadLinesFile(input_file, lines_data_file);
#endif

    ReadInputFile(lines_data_file, *this);
    
    // mesh and finite element is constructed
    this->ComputeMeshAndFiniteElement(name_element, false);

    // other initializations
    var_fluid.var_harmonic.PerformOtherInitializations();
    var_solid.var_harmonic.PerformOtherInitializations();
    
    // right hand sides are computed
    var_fluid.ComputeRightHandSide();
    var_solid.ComputeRightHandSide();
    
    // mass and stiffness matrices are computed
    var_fluid.var_harmonic.ComputeMassMatrix();
    var_solid.var_harmonic.ComputeMassMatrix();
    
    var_fluid.var_harmonic.ComputeQuasiPeriodicPhase();
    var_solid.var_harmonic.ComputeQuasiPeriodicPhase();
    
    var_fluid.ComputeMassMatrix();
    var_solid.ComputeMassMatrix();

    var_fluid.ComputeStiffnessMatrix();
    var_solid.ComputeStiffnessMatrix();
    
    // initialization before time iterations
    var_fluid.CheckIdentityJacobianInPML();
    var_solid.CheckIdentityJacobianInPML();
    var_fluid.InitTimeIterations();
    var_solid.InitTimeIterations();
    
    var_fluid.var_harmonic.FindDofsOnReference(this->ref_interface, DdlLambda);
    DISP(DdlLambda);
    
    // total number of unknowns
    this->nodl = var_fluid.var_harmonic.GetNbDof() + var_solid.var_harmonic.GetNbDof() + DdlLambda.GetM();
    
  }

  
  //! adding coupling operator alpha B to matrix A
  template<class Dimension>
  void AcousticAcousticInteraction<Dimension>
  ::AddOperatorB(const Real_wp& alpha, VirtualMatrix<Real_wp>& mat_sp)
  {
    const Mesh<Dimension>& mesh_elas = var_solid.var_harmonic.mesh;
    const MeshNumbering<Dimension>& mesh_elas_num = var_solid.var_harmonic.GetMeshNumbering(0);

    const Mesh<Dimension>& mesh_acous = var_fluid.var_harmonic.mesh;
    const MeshNumbering<Dimension>& mesh_acous_num = var_fluid.var_harmonic.GetMeshNumbering(0);

    EllipticProblem<LaplaceEquation<Dimension> > & var_acous = var_fluid.var_harmonic;    
    EllipticProblem<LaplaceEquation<Dimension> > & var_elas = var_solid.var_harmonic;

    Vector<int> IndexLambda(var_acous.GetNbDof());
    IndexLambda.Fill(-1);
    for (int i = 0; i < DdlLambda.GetM(); i++)
      IndexLambda(DdlLambda(i)) = i;
    
    for (int i = 0; i < this->NumElementFluid.GetM(); i++)
      {
	int ne_f = this->NumElementFluid(i);
	const ElementReference<Dimension, 1>& Fb_acous = var_acous.GetReferenceElementH1(ne_f);
	
	int ne_s = this->NumElementSolid(i);
	const ElementReference<Dimension, 1>& Fb_elas = var_elas.GetReferenceElementH1(ne_s);
	
	int Nsol = var_elas.GetNbDof();
	int offset_lambda = Nsol + var_acous.GetNbDof();
	
	int num_face = mesh_elas.Element(ne_s).numBoundary(this->NumLocSolid(i));
	int r = mesh_acous_num.GetOrderQuadrature(num_face);
	int nb_pts = mesh_acous_num.GetNbPointsQuadratureBoundary(num_face);
	
	int rot1 = mesh_elas.Element(ne_s).GetOrientationBoundary(this->NumLocSolid(i));
	int rot2 = mesh_acous.Element(ne_f).GetOrientationBoundary(this->NumLocFluid(i));
	int nv = mesh_elas.Boundary(num_face).GetNbVertices();
	int rot = mesh_acous_num.GetRotationFace(rot1, rot2, nv);
	const VectReal_wp& PoidsFlux = mesh_elas_num.number_map.
	  GetFluxWeight(r, mesh_elas.Boundary(num_face));
	
	const Matrix<int>& FacesQuadRotation = mesh_elas_num.number_map.
	  GetRotationQuadraturePoints(r, mesh_elas.Boundary(num_face));
	
	VectReal_wp phi(Fb_elas.GetNbDof()), val_phi(nb_pts), feval(nb_pts);
	VectReal_wp contrib(Fb_acous.GetNbDof());
	contrib.Zero(); feval.Zero();
	phi.Zero(); val_phi.Zero();
	
	VectReal_wp val_interac(Fb_acous.GetNbDof());
	int c_interac;
	IVect col_interac(Fb_acous.GetNbDof());
	val_interac.Zero(); col_interac.Fill(-1);
	
	//DISP(i); DISP(ne_s); DISP(ne_f);     DISP(num_face);
	// loop over dofs for elastic function test
	for (int iloc = 0; iloc < Fb_elas.GetNbDof(); iloc++)
	  {
	    phi(iloc) = 1.0;
	    int num_dof = mesh_elas_num.Element(ne_s).GetNumberDof(iloc);
	    
	    Fb_elas.ApplyShTranspose(this->NumLocSolid(i), phi, val_phi, r);
	    bool null_fct = true;
	    for (int k = 0; k < nb_pts; k++)
	      if (abs(val_phi(k)) > epsilon_machine)
		null_fct = false;
	    
	    if (!null_fct)
	      {
		for (int k = 0; k < nb_pts; k++)
		  {
		    Real_wp poids = 2.0*var_elas.Glob_dsj(num_face)(k)*PoidsFlux(k)*val_phi(k)*alpha;
		    int krot = FacesQuadRotation(rot, k);
		    feval(krot) = poids;
		  }
		
		contrib.Zero();		
		Fb_acous.ApplySh(1.0, this->NumLocFluid(i), feval, contrib, r);
		c_interac = 0;
		for (int j = 0; j < Fb_acous.GetNbDof(); j++)
		  if (abs(contrib(j)) > var_elas.GetThresholdMatrix())
		    {
		      int nlambda = IndexLambda(mesh_acous_num.Element(ne_f).GetNumberDof(j));
		      if (nlambda < 0)
			{
			  cout << "Impossible case" << endl;
			  abort();
			}
		      
		      col_interac(c_interac) = offset_lambda + nlambda;
		      val_interac(c_interac) = contrib(j);
		      mat_sp.AddInteraction(col_interac(c_interac), num_dof, val_interac(c_interac));
		      c_interac++;
		    }
		
		mat_sp.AddInteractionRow(num_dof, c_interac, col_interac, val_interac);
	      }
	    
	    phi(iloc) = 0.0;
	  }

	num_face = mesh_acous.Element(ne_f).numBoundary(this->NumLocFluid(i));

	// loop over dofs for acoustic function test
	phi.Reallocate(Fb_acous.GetNbDof());
	phi.Zero(); val_phi.Zero();
	for (int iloc = 0; iloc < Fb_acous.GetNbDof(); iloc++)
	  {
	    phi(iloc) = 1.0;
	    int num_dof = mesh_acous_num.Element(ne_f).GetNumberDof(iloc);
	    
	    Fb_acous.ApplyShTranspose(this->NumLocFluid(i), phi, val_phi, r);
	    bool null_fct = true;
	    for (int k = 0; k < nb_pts; k++)
	      if (abs(val_phi(k)) > epsilon_machine)
		null_fct = false;
	    
	    if (!null_fct)
	      {
		for (int k = 0; k < nb_pts; k++)
		  {
		    Real_wp poids = -2.0*var_acous.Glob_dsj(num_face)(k)*PoidsFlux(k)*val_phi(k)*alpha;
		    feval(k) = poids;
		  }
		
		contrib.Zero();		
		Fb_acous.ApplySh(1.0, this->NumLocFluid(i), feval, contrib, r);
		c_interac = 0;
		for (int j = 0; j < Fb_acous.GetNbDof(); j++)
		  if (abs(contrib(j)) > var_elas.GetThresholdMatrix())
		    {
		      int nlambda = IndexLambda(mesh_acous_num.Element(ne_f).GetNumberDof(j));
		      if (nlambda < 0)
			{
			  cout << "Impossible case" << endl;
			  abort();
			}
		      
		      col_interac(c_interac) = offset_lambda + nlambda;
		      val_interac(c_interac) = contrib(j);
		      mat_sp.AddInteraction(col_interac(c_interac), Nsol + num_dof, val_interac(c_interac));
		      c_interac++;
		    }

		mat_sp.AddInteractionRow(Nsol + num_dof, c_interac, col_interac, val_interac);
	      }
	    
	    phi(iloc) = 0.0;
	  }	
      }
  }
  
  
  //! Computes mass and stiffness matrices (without damping)
  template<class Dimension>
  void AcousticAcousticInteraction<Dimension>
  ::ComputeMatrices(DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& M,
		    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>& K)
  {
    cout << "Number of dofs in fluid domain = " << var_fluid.var_harmonic.GetNbDof() << endl;
    cout << "Number of dofs in solid domain = " << var_solid.var_harmonic.GetNbDof() << endl;
    cout << "Number of dofs for Lambda = " << DdlLambda.GetM() << endl;
    
    M.Reallocate(this->nodl, this->nodl);
    K.Reallocate(this->nodl, this->nodl);
    
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(0);
    nat_mat.SetCoefStiffness(1);
    nat_mat.SetCoefDamping(0);
    
    int Ns = var_solid.var_harmonic.GetNbDof();

    var_solid.var_harmonic.SetCoefficientDirichlet(Real_wp(0));
    var_fluid.var_harmonic.SetCoefficientDirichlet(Real_wp(0));
    var_solid.var_harmonic.AddMatrixWithBC(K, nat_mat, 0, 0);
    var_fluid.var_harmonic.AddMatrixWithBC(K, nat_mat, Ns, Ns);

    AddOperatorB(Real_wp(1), K);
    
    nat_mat.SetCoefMass(1);
    nat_mat.SetCoefStiffness(0);
    nat_mat.SetCoefDamping(0);

    var_solid.var_harmonic.SetCoefficientDirichlet(Real_wp(1));
    var_fluid.var_harmonic.SetCoefficientDirichlet(Real_wp(1));
    var_solid.var_harmonic.AddMatrixWithBC(M, nat_mat, 0, 0);
    var_fluid.var_harmonic.AddMatrixWithBC(M, nat_mat, Ns, Ns);
    
  }

}

#define MONTJOIE_FILE_FLUID_STRUCTURE_INTERACTION_CXX
#endif

