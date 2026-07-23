#ifndef FILE_MODEL_THIN_SLOT_HELMHOLTZ_CXX

namespace Montjoie
{

#ifdef MONTJOIE_WITH_THIN_SLOT_MODEL
#ifdef MONTJOIE_WITH_TWO_DIM  

  int SlotModelParameters<Dimension2>::type_model(0);
  
  //! computation of \f$\int_{t_0}^{t_1} \varphi_k(x) dt \f$
  /*!
    \param[in] t0 
    \param[in] t1
    \param[in] ds length of the edge, the integral is multiplied by this length
    \param[in] k basis function number
  */
  template<>
  Real_wp VarHelmholtz_Base<Complex_wp, Dimension2>
  ::FirstModeSlot2D(const Real_wp & t0, const Real_wp& t1,
                    const Real_wp& ds, int k, const ElementReference<Dimension2, 1>& Fb) const
  {
    Real_wp dt = t1 - t0;
    Real_wp FM = 0;
    
    const ElementReference<Dimension1, 1>& edge = Fb.GetSurfaceFiniteElement(0);
    
    for (int m = 0; m <= edge.GetOrder(); m++)
      FM += edge.Weights(m)*edge.GetValuePhi1D(k, edge.Points(m)*(t1-t0) + t0);
    
    FM *= ds*abs(dt);
    
    return FM;
  }


  //! localization of thin slots in the mesh 
  template<>
  void VarHelmholtz_Base<Complex_wp, Dimension2>::LocalizeSlots()
  {
    // DISP(this->list_slots);
    // localization of slots
    R2 ptA, ptB, ptS, vec_u, vec_v; Real_wp local_pos, prod_vec, prod_scal;
    int nb_slots = this->list_slots.GetM();
    Vector<IVect> Slot_Extremity(nb_slots);
    Vector<IVect> Slot_NumEdge(nb_slots); Vector<VectReal_wp> Slot_LocalPos(nb_slots);
    for (int i = 0; i < var_problem.mesh.GetNbBoundaryRef(); i++)
      {
	// loop on the two points of the edge
	// if these points own to the slots of the mesh
	for (int l = 0; l < 2; l++)
	  {
	    local_pos = l;
	    int nv = var_problem.mesh.BoundaryRef(i).numVertex(l);
	    ptA = var_problem.mesh.Vertex(nv);
	    for (int nslot = 0; nslot < nb_slots; nslot++)
	      for (int k = 0; k < 2; k++)
		{
		  if ((this->list_slots(nslot).ptC(k) != ptA)
                      &&(this->list_slots(nslot).ptD(k) != ptA))
		    {
		      vec_u = ptA - this->list_slots(nslot).ptC(k);
		      vec_v = ptA - this->list_slots(nslot).ptD(k);
		      prod_vec = vec_u(0)*vec_v(1) - vec_u(1)*vec_v(0);
		      prod_scal = DotProd(vec_u,vec_v);
		      if ((abs(prod_vec) < R2::threshold)&&(prod_scal <= R2::threshold))
			{
			  // the point is inside the slot, we add it
			  Slot_NumEdge(nslot).PushBack(i);
			  Slot_LocalPos(nslot).PushBack(local_pos);
			  Slot_Extremity(nslot).PushBack(k);
			}
		    }
		}
	  }
	
	// now we test if the points of the slot are inside the edge
	int nv = var_problem.mesh.BoundaryRef(i).numVertex(0);
	ptA = var_problem.mesh.Vertex(nv);
	nv = var_problem.mesh.BoundaryRef(i).numVertex(1);
	ptB = var_problem.mesh.Vertex(nv);
	Real_wp invAB = 1.0/ptA.Distance(ptB);
	for (int nslot = 0; nslot < nb_slots; nslot++)
	  for (int k = 0; k < 2; k++)
	    for (int l = 0; l < 2; l++)
	      {
		if (l == 0)
		  ptS = this->list_slots(nslot).ptC(k);
		else
		  ptS = this->list_slots(nslot).ptD(k);
		
		vec_u = ptS - ptA;
		vec_v = ptS - ptB;
		prod_vec = vec_u(0)*vec_v(1) - vec_u(1)*vec_v(0);
		prod_scal = DotProd(vec_u,vec_v);
		if ((abs(prod_vec) < R2::threshold)&&(prod_scal <= R2::threshold))
		  {
		    // the point of the slot is inside the edge
		    // finding local position in the edge
		    local_pos = Norm2(vec_u)*invAB;
		    // adding the edge to the slot
		    Slot_NumEdge(nslot).PushBack(i);
		    Slot_LocalPos(nslot).PushBack(local_pos);
		    Slot_Extremity(nslot).PushBack(k);
		  }
	      }
      }
    
    // now sorting the edges of each slot
    for (int nslot = 0; nslot < nb_slots; nslot++)
      Sort(Slot_NumEdge(nslot).GetM(), Slot_NumEdge(nslot),
           Slot_LocalPos(nslot), Slot_Extremity(nslot));
    
    Vector<VectReal_wp> Slot_SecondPos(nb_slots);
    // eliminating edges with only one point 
    // (edges that are at the extremities of the slot, but no overlapping)
    for (int nslot = 0; nslot < nb_slots; nslot++)
      {
	int iedge = 0, nb_edges = 0;
	Slot_SecondPos(nslot).Reallocate(Slot_NumEdge(nslot).GetM());
	while (iedge < (Slot_NumEdge(nslot).GetM()-1))
	  {
	    if (Slot_NumEdge(nslot)(iedge) == Slot_NumEdge(nslot)(iedge+1))
	      {
		// case of an edge overlapping the slot
		Slot_NumEdge(nslot)(nb_edges) = Slot_NumEdge(nslot)(iedge);
		Slot_Extremity(nslot)(nb_edges) = Slot_Extremity(nslot)(iedge);
		Slot_LocalPos(nslot)(nb_edges) = Slot_LocalPos(nslot)(iedge);
		Slot_SecondPos(nslot)(nb_edges) = Slot_LocalPos(nslot)(iedge+1);
		iedge += 2; nb_edges++;
	      }
	    else
	      iedge++;
	  }
	
	// resizing arrays to the true number of edges belonging to the slot
	Slot_NumEdge(nslot).Resize(nb_edges);
	Slot_Extremity(nslot).Resize(nb_edges);
	Slot_LocalPos(nslot).Resize(nb_edges);
	Slot_SecondPos(nslot).Resize(nb_edges);
	
	// copy of the array in list_slots
	this->list_slots(nslot).Slot_NumEdge = Slot_NumEdge(nslot);
	this->list_slots(nslot).Slot_Extremity = Slot_Extremity(nslot);
	this->list_slots(nslot).Slot_LocalPos = Slot_LocalPos(nslot);
	this->list_slots(nslot).Slot_SecondPos = Slot_SecondPos(nslot);
	
	// DISP(nslot);
	// DISP(this->list_slots(nslot).Slot_NumEdge); DISP(this->list_slots(nslot).Slot_Extremity);
	// DISP(this->list_slots(nslot).Slot_LocalPos); DISP(this->list_slots(nslot).Slot_SecondPos);
      }
  }
  

  //! For 1D-2D model of a slot, some interactions are added on the matrix
  /*!
    \param[in] mat2D matrix to which interactions are added
  */
  template<>
  void VarHelmholtz_Base<Complex_wp, Dimension2>::AddContributions1D_Slot(VirtualMatrix<Complex_wp>& mat2D)
  {    
    Real_wp FM1, FM2; R2 ptA, ptB;

    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    // loop on all the slots
    for (int nslot = 0; nslot < this->list_slots.GetM(); nslot++)
      {
	// offset for the 1-D dofs
	int offset = this->list_slots(nslot).offset_nodl;
	Matrix<Complex_wp, Symmetric, ArrayRowSymSparse>& mat1D = this->list_slots(nslot).mat_direct1D;
	int dof1D_ext(0), num_elt1D(0);
	Real_wp epsilon = this->list_slots(nslot).epsilon;
	
        int order = mesh_num.GetOrder();
        int nb_elt_int = this->list_slots(nslot).Slot_NumEdge.GetM()+2;
	IVect col_interac(nb_elt_int*order+2);
	Vector<Complex_wp> val_interac(nb_elt_int*order+2);
	int c_interac = 0;
	
	// loop on all the edges belonging to the slot
	for (int iedge = 0; iedge < this->list_slots(nslot).Slot_NumEdge.GetM(); iedge++)
	  {
	    if (this->list_slots(nslot).Slot_Extremity(iedge) == 0)
	      {
		// edge on the left side of the slot
		num_elt1D = 0;
		dof1D_ext = this->list_slots(nslot).var1D.mesh.GetNumberDof(num_elt1D, 0);
	      }
	    else
	      {
		// edge on the right side of the slot
		num_elt1D = this->list_slots(nslot).var1D.mesh.GetNbElt() - 1;
		dof1D_ext = this->list_slots(nslot).var1D.mesh.GetNumberDof(num_elt1D, order);
	      }
	    
	    int num_ref = this->list_slots(nslot).Slot_NumEdge(iedge);
	    ptA = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref).numVertex(0));
	    ptB = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref).numVertex(1));
	    Real_wp dsj = ptA.Distance(ptB);
	    Real_wp t0 = this->list_slots(nslot).Slot_LocalPos(iedge);
	    Real_wp t1 = this->list_slots(nslot).Slot_SecondPos(iedge);
	    
	    // getting informations about the edge
            int num_edge = num_ref;
            int iquad = var_problem.mesh.Boundary(num_edge).numElement(0);
            int num_loc = var_problem.mesh.Element(iquad).GetPositionBoundary(num_edge);            
            order = mesh_num.GetOrderElement(iquad);
            const ElementReference<Dimension2, 1>& Fb = var_problem.GetReferenceElementH1(iquad);
            int nb_dof_edge = order+1;
            if (var_problem.mesh.Element(iquad).numVertex(num_loc)
                != var_problem.mesh.BoundaryRef(num_ref).numVertex(0))
              {
                Real_wp ti = 1.0-t1;
                t1 = 1.0-t0;
                t0 = ti;
              }
            
            IVect num_dof(nb_dof_edge);            
	    // first, we add the 2-D <-> 2-D contributions	    
	    for(int k = 0; k < nb_dof_edge; k++)
	      {  
		col_interac.Fill(0);
		val_interac.Fill(0);
		c_interac = 0;
		
		FM1 = this->FirstModeSlot2D(t0, t1, dsj, k, Fb);
		
		// we get global number dof
		int num_dof_loc = Fb.GetLocalNumber(num_loc, k);
		num_dof(k) = mesh_num.Element(iquad).GetNumberDof(num_dof_loc);
		
		for (int iedge2 = 0; iedge2 < this->list_slots(nslot).Slot_NumEdge.GetM(); iedge2++)
		  if (this->list_slots(nslot).Slot_Extremity(iedge2)
                      == this->list_slots(nslot).Slot_Extremity(iedge))
		    {
		      int num_ref2 = this->list_slots(nslot).Slot_NumEdge(iedge2);
		      ptA = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref2).numVertex(0));
		      ptB = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref2).numVertex(1));
		      Real_wp dsj2 = ptA.Distance(ptB);
		      Real_wp t0b = this->list_slots(nslot).Slot_LocalPos(iedge2);
		      Real_wp t1b = this->list_slots(nslot).Slot_SecondPos(iedge2);
		      
		      // getting informations about the edge
                      int num_edge2 = num_ref2;
                      int iquad2 = var_problem.mesh.Boundary(num_edge2).numElement(0);
                      int num_loc2 = var_problem.mesh.Element(iquad2).GetPositionBoundary(num_edge2);
                      int order2 = mesh_num.GetOrderElement(iquad2);
                      int nb_dof_edge2 = order2 + 1;
		      const ElementReference<Dimension2, 1>& Fb2 = var_problem.GetReferenceElementH1(iquad2);
                      if (var_problem.mesh.Element(iquad2).numVertex(num_loc2)
                          != var_problem.mesh.BoundaryRef(num_ref2).numVertex(0))
                        {
                          Real_wp ti = 1.0-t1b;
                          t1b = 1.0-t0b;
                          t0b = ti;
                        }

                      for(int l = 0; l < nb_dof_edge2; l++)
			{
			  num_dof_loc = Fb2.GetLocalNumber(num_loc2, l);
			  col_interac(c_interac)
                            = mesh_num.Element(iquad2).GetNumberDof(num_dof_loc);
			  
                          FM2 = this->FirstModeSlot2D(t0b, t1b, dsj2, l, Fb2);
			  to_complex(FM1*FM2/epsilon*mat1D(dof1D_ext,dof1D_ext),
                                     val_interac(c_interac));
			  //DISP(FM1);DISP(FM2);DISP(mat1D(0,0));
			  c_interac++;
			}
		    }
		
		// 2-D 1-D contributions
		for(int q = 0; q < nb_dof_edge; q++)
		  {
		    int num1D = this->list_slots(nslot).var1D.mesh.GetNumberDof(num_elt1D, q);
		    if (num1D != dof1D_ext)
		      {
			col_interac(c_interac) = num1D + offset - 1;
			to_complex(FM1*mat1D(num1D, dof1D_ext), val_interac(c_interac));
			c_interac++;
		      }
		  }
		
		// adding interactions to the 2-D matrix
		mat2D.AddInteractionRow(num_dof(k), c_interac , col_interac, val_interac);
	      }
	    
	    // loop on 1-D dofs
	    num_dof.Fill(0);
	    for(int p = 0; p < nb_dof_edge; p++)
	      {
		int num1D = this->list_slots(nslot).var1D.mesh.GetNumberDof(num_elt1D, p);
		num_dof(p) = offset + num1D-1;
		if (num1D != dof1D_ext)
		  {
		    c_interac = 0;
		    col_interac.Fill(0);
		    val_interac.Fill(0);
		    
		    // interactions with 2-D
		    for(int l = 0; l < nb_dof_edge; l++)
		      {
			int num_dof_loc = Fb.GetLocalNumber(num_loc, l);
			col_interac(c_interac)
                          = mesh_num.Element(iquad).GetNumberDof(num_dof_loc);
			
                        FM1 = this->FirstModeSlot2D(t0, t1, dsj, l, Fb);
			to_complex(FM1*mat1D(dof1D_ext, num1D), val_interac(c_interac));
			c_interac++;
		      }
		    
		    mat2D.AddInteractionRow(num_dof(p), c_interac , col_interac, val_interac);
		  }
	      }
	  }
	
	// adding the 1-D matrix
	// interactions with 1-D
	for (int p = 1; p < mat1D.GetM()-1; p++)
	  {
	    c_interac = 0;
	    int num1D = offset + p -1;
	    for(int q = 0; q < mat1D.GetRowSize(p); q++)
	      {
		int num_dof_loc = mat1D.Index(p, q);
		if ((num_dof_loc != 0)&&(num_dof_loc != mat1D.GetM()-1))
		  {
		    col_interac(c_interac) = offset + num_dof_loc-1;
		    to_complex(epsilon*mat1D.Value(p, q), val_interac(c_interac));
		    c_interac++;
		  }
	      }
            
	    mat2D.AddInteractionRow(num1D, c_interac , col_interac, val_interac);
	  }
      }
  }

  
  //! modification of the matrix in the case of approximate thin-slot models
  template<>
  void VarHelmholtz_Base<Complex_wp, Dimension2>::TreatThinSlotCondition(VirtualMatrix<Complex_wp>& mat_sp)
  {
    //cout << "coucou slot " << endl;
    if (this->list_slots.GetM() <= 0)
      return;
        
    // localization of the slots
    LocalizeSlots();
    
    if (list_slots(0).type_model == list_slots(0).MODEL_MESH1D)
      {
	// case of coupling 1-D mesh with 2-D mesh
	this->AddContributions1D_Slot(mat_sp);
	return;
      }

    const MeshNumbering<Dimension2>& mesh_num = var_problem.GetMeshNumbering(0);
    Real_wp FirstMode1, FirstMode2; R2 ptA, ptB;
    Real_wp omega = var_problem.GetOmega();
    // loop over the slots
    for (int nslot = 0; nslot < this->list_slots.GetM(); nslot++)
      {
	// coefficients for the DtN (assuming that omega is not a resonant frequency)
	Real_wp L = this->list_slots(nslot).length;
	Real_wp epsilon = this->list_slots(nslot).epsilon;
	Real_wp coef_diag = omega/tan(omega*L)/epsilon;
	Real_wp coef_extra = -omega/sin(omega*L)/epsilon;
        //DISP(coef_diag); DISP(coef_extra);
        
        int order = mesh_num.GetOrder();    
        int nb_elt_int = this->list_slots(nslot).Slot_NumEdge.GetM()+2;
	IVect col_interac(nb_elt_int*order+2);
	Vector<Complex_wp> val_interac(nb_elt_int*order+2);
	int num_dof, c_interac, num_dof2;
	
	// loop on all the edges belonging to the slot
	for (int iedge = 0; iedge < this->list_slots(nslot).Slot_NumEdge.GetM(); iedge++)
	  {
	    // informations on the edge
	    int num_ref = this->list_slots(nslot).Slot_NumEdge(iedge);
	    ptA = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref).numVertex(0));
	    ptB = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref).numVertex(1));
	    Real_wp dsj = ptA.Distance(ptB);
	    Real_wp t0 = this->list_slots(nslot).Slot_LocalPos(iedge);
	    Real_wp t1 = this->list_slots(nslot).Slot_SecondPos(iedge);
	    //DISP(num_ref); DISP(ptA); DISP(ptB); DISP(coef_diag); DISP(coef_extra);
	    //DISP(t0); DISP(t1); DISP(dsj);
	    
	    // getting informations about the edge
            int num_edge = num_ref;
            int iquad = var_problem.mesh.Boundary(num_edge).numElement(0);
            int num_loc = var_problem.mesh.Element(iquad).GetPositionBoundary(num_edge);            
            order = mesh_num.GetOrderElement(iquad);
            const ElementReference<Dimension2, 1>& Fb = var_problem.GetReferenceElementH1(iquad);
            int nb_dof_edge = order+1;
            if (var_problem.mesh.Element(iquad).numVertex(num_loc)
                != var_problem.mesh.BoundaryRef(num_ref).numVertex(0))
              {
                Real_wp ti = 1.0-t1;
                t1 = 1.0-t0;
                t0 = ti;
              }

	    for(int k = 0; k < nb_dof_edge; k++)
	      {
                // value of \int \varphi_k
		FirstMode1 = FirstModeSlot2D(t0, t1, dsj, k, Fb);
		int num_dof_loc = Fb.GetLocalNumber(num_loc, k);
		num_dof = mesh_num.Element(iquad).GetNumberDof(num_dof_loc);
		
		// now loop on all the edges of the slot
		for (int iedge2 = 0; iedge2 < this->list_slots(nslot).Slot_NumEdge.GetM(); iedge2++)
		  {
		    // informations on the edge
		    int num_ref2 = this->list_slots(nslot).Slot_NumEdge(iedge2);
		    ptA = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref2).numVertex(0));
		    ptB = var_problem.mesh.Vertex(var_problem.mesh.BoundaryRef(num_ref2).numVertex(1));
		    Real_wp dsj2 = ptA.Distance(ptB);
		    Real_wp t0b = this->list_slots(nslot).Slot_LocalPos(iedge2);
		    Real_wp t1b = this->list_slots(nslot).Slot_SecondPos(iedge2);
		    
		    // getting informations about the edge
                    int num_edge2 = num_ref2;
                    int iquad2 = var_problem.mesh.Boundary(num_edge2).numElement(0);
                    int num_loc2 = var_problem.mesh.Element(iquad2).GetPositionBoundary(num_edge2);
                    int order2 = mesh_num.GetOrderElement(iquad2);
                    int nb_dof_edge2 = order2 + 1;
		    const ElementReference<Dimension2, 1>& Fb2 = var_problem.GetReferenceElementH1(iquad2);
                    if (var_problem.mesh.Element(iquad2).numVertex(num_loc2)
                        != var_problem.mesh.BoundaryRef(num_ref2).numVertex(0))
                      {
                        Real_wp ti = 1.0-t1b;
                        t1b = 1.0-t0b;
                        t0b = ti;
                      }

		    col_interac.Fill(0);
		    val_interac.Fill(0);
		    c_interac = 0;
		    for (int l = 0; l < nb_dof_edge2; l++)
		      {
			// value of \int \varphi_l
			FirstMode2 = FirstModeSlot2D(t0b, t1b, dsj2, l, Fb2);
			num_dof_loc = Fb2.GetLocalNumber(num_loc2, l);
			num_dof2 = mesh_num.Element(iquad2).GetNumberDof(num_dof_loc);
			col_interac(c_interac) = num_dof2;
                                                
			// left-left, right-right or left-right interaction ?
			if (this->list_slots(nslot).Slot_Extremity(iedge)
                            == this->list_slots(nslot).Slot_Extremity(iedge2))
			  {
			    // left-left or right-right
			    val_interac(c_interac) = FirstMode1*FirstMode2*coef_diag;
			  }
			else
			  {
			    // left-right
			    val_interac(c_interac) = FirstMode1*FirstMode2*coef_extra;
			  }
			
			c_interac++;
		      }
		    
		    // DISP(num_dof); DISP(c_interac); DISP(col_interac); DISP(val_interac);
		    // adding interactions to the matrix
		    mat_sp.AddInteractionRow(num_dof, c_interac, col_interac, val_interac);
		  } 
	      }
	  } 
      }
  }
    
  
  //! default constructor
  SlotModelParameters<Dimension2>::SlotModelParameters()
  {
    epsilon = 0; length = 0; order = 0; offset_nodl = 0; nb_points = 0;
  }
  
  //! we display informations about SlotModelParameters
  ostream& operator<<(ostream& out, const SlotModelParameters<Dimension2>& x)
  {
    out<<"Slot of extremities "<<x.ptA<<"   "<<x.ptB<<'\n';
    out<<"Extremities of each slot "<<x.ptC<<'\n'<<x.ptD<<'\n';
    out<<"Thickness "<<x.epsilon<<'\n';
    out<<"Length "<<x.length<<'\n';
    out<<"Order of approximation "<<x.order<<'\n';
    out<<"Number 1-D points "<<x.nb_points<<endl;
    return out;
  }
  
  //! computation of 1-D helmholtz problem related to the slot
  void SlotModelParameters<Dimension2>::ComputeVar1D()
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var1D.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if (rank_proc == 0)
      cout<<"on rentre dans ComputeVar1D "<<endl;
    
    var1D.InitIndices(PhysicalConstant::nb_max_indices);
    var1D.SetBoundaryCondition(BoundaryConditionEnum::LINE_NEUMANN, BoundaryConditionEnum::LINE_NEUMANN);
    var1D.mesh.CreateRegularMesh(0, ptA.Distance(ptB), nb_points, 1);
    var1D.mesh.SetOrder(order);
    var1D.ComputeMeshAndFiniteElement(string("EDGE_LOBATTO")); // DISP(var1D.mesh.Nodle);
    var1D.ComputeMassMatrix();
    var1D.TreatDirichletCondition();
    
    GlobalGenericMatrix<Complex_wp> nature_matrix;
    var1D.AddMatrixFEM(mat_direct1D, nature_matrix);
    var1D.AddBoundaryTerms(mat_direct1D, nature_matrix);
    // DISP(mat_direct1D);
  }
  
#endif

#ifdef MONTJOIE_WITH_TWO_DIM  
  template<>
  void VarHelmholtz_Base<Complex_wp, Dimension3>::TreatThinSlotCondition(VirtualMatrix<Complex_wp>& mat_sp)
  {
    // not implemented in 3-D
  }
#endif

#endif  
} // end namespace


#define FILE_MODEL_THIN_SLOT_HELMHOLTZ_CXX
#endif
