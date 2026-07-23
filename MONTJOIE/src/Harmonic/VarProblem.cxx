#ifndef MONTJOIE_FILE_VAR_PROBLEM_CXX

namespace Montjoie
{

  /**************
   * VarProblem *
   **************/
  
  
  //! adds memory used by the current object
  template<class Dimension>
  void VarProblem<Dimension>::GetMemoryUsed(map<string, size_t>& var) const
  {
    DistributedProblem<Dimension>::GetMemoryUsed(var);
    VarFiniteElement<Dimension>::GetMemoryUsed(var);
  }
  

  //! fills numbers of dofs on boundaries, inside and numbers of inside dofs
  template<class Dimension>
  void VarProblem<Dimension>::
  GetInternalNodesElement(int i, int nb_dof_loc,
			  int& nb_dof_edges, int& nb_dof_int,
			  Vector<int>& intern_node) const
  {
    if (this->FormulationDG() == ElementReference_Base::HDG)
      DistributedProblem<Dimension>::GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges,
							     nb_dof_int, intern_node);
    else
      VarComputationProblem::GetInternalNodesElement(i, nb_dof_loc, nb_dof_edges,
						     nb_dof_int, intern_node);
  }
  
  
  //! construction of finite elements
  template<class Dimension>
  void VarProblem<Dimension>::ConstructFiniteElement(const string& name_elt)
  {
    // previous finite elements are cleared if present
    this->ClearFiniteElement();
    
    // orders present in the mesh
    bool variable_order = this->mesh_num.IsOrderVariable();
    Vector<int> r_over(1 + this->other_mesh_num.GetM());
    r_over.Fill(this->order_over_integration);
    TinyVector<IVect, 4> order;    
    this->mesh_num.GetOrder(order);
    int rmax = 0;
    Vector<int> rn(this->other_mesh_num.GetM()+1);
    for (int k = 0; k < 4; k++)
      {
        if (order(k).GetM() > 0)
          {
            rn(0) = order(k)(0);
            rmax = max(rmax, order(k)(0));
          }
      }
    
    // orders for other numberings
    Vector<TinyVector<IVect, 4> > other_order(this->other_mesh_num.GetM());
    Vector<string> other_name(this->other_mesh_num.GetM());
    Vector<int> other_type(this->other_mesh_num.GetM());
    for (int i = 0; i < this->other_mesh_num.GetM(); i++)
      {
        this->other_mesh_num(i)->GetOrder(other_order(i));
        if (this->name_other_elements(i).size() > 0)
          other_name(i) = this->name_other_elements(i);
        else
          other_name(i) = name_elt;

        if (this->other_type_element.GetM() > i)
          other_type(i) = this->other_type_element(i);
        else
          other_type(i) = this->type_element;
        
        for (int k = 0; k < 4; k++)
          {
            if (other_order(i)(k).GetM() > 0)
              {
                rn(i+1) = other_order(i)(k)(0);
                rmax = max(rmax, other_order(i)(k)(0));
              }
          }
      }

    if (!variable_order && (this->other_mesh_num.GetM() > 0))
      {
        // using over-integration to have same quadrature points for all elements
        // (for uniform order)
        for (int k = 0; k < rn.GetM(); k++)
          r_over(k) = rmax - rn(k) + this->order_over_integration;
      }
    
    // adding orders of neighboring elements
    int nb = 0;
    for (int i = 0; i < this->NeighboringConnectivity.GetM(); i++)
      if (this->NeighboringConnectivity(i).GetM() > 0)
        nb++;

    this->mesh_num.ReallocateNeighborElements(nb);
    for (int i = 0; i < this->other_mesh_num.GetM(); i++)
      this->other_mesh_num(i)->ReallocateNeighborElements(nb);
    
    nb = 0; int nb_mesh_num = this->other_mesh_num.GetM() + 1;
    for (int i = 0; i < this->NeighboringConnectivity.GetM(); i++)
      if (this->NeighboringConnectivity(i).GetM() > 0)
        {
          int type_elt = this->GetTypeEltNeighboringFace(i);
          for (int nm = 0; nm < nb_mesh_num; nm++)
            {
              MeshNumbering<Dimension>& mesh_n = this->GetMeshNumbering(nm);
              int r = this->GetOrderEltNeighboringFace(i, nm);
              mesh_n.SetOrderNeighborElement(nb, r);
              mesh_n.SetTypeNeighborElement(nb, type_elt);
              mesh_n.SetEdgeNeighborElement(nb, i);

              TinyVector<IVect, 4>* order_ptr = NULL;
              if (nm == 0)
                order_ptr = &order;
              else
                order_ptr = &other_order(nm-1);
              
              TinyVector<IVect, 4>& order_n = *order_ptr;
              bool new_order = true;
              for (int k = 0; k < order_n(type_elt).GetM(); k++)
                if (order_n(type_elt)(k) == r)
                  new_order = false;
              
              if (new_order)
                order_n(type_elt).PushBack(r);
            }
          
          nb++;          
        }

    this->mesh_num.FinalizeNeighborElements(this->NeighboringConnectivity.GetM());
    for (int i = 0; i < this->other_mesh_num.GetM(); i++)
      this->other_mesh_num(i)->FinalizeNeighborElements(this->NeighboringConnectivity.GetM());
    
    // all the elements are affected the new finite element
    Vector<bool> change_elt(this->mesh.GetNbElt());
    change_elt.Fill(true);
    
    // new finite elements are added
    this->AddFiniteElement(name_elt, order, this->type_element,
                           other_name, other_order, other_type,
                           change_elt, this->mesh, this->mesh_num,
			   this->other_mesh_num, r_over, this->dg_formulation,
                           this->other_dg_formulation);
    
    // retrieves the type of quadrature for edges or triangles/quadrangles
    this->mesh_num.GetOrderQuadrature(order);
    RemoveDuplicate(order(0));

    // for hdg, we construct surface finite element
    if (this->FormulationDG() == ElementReference_Base::HDG)
      {
	change_elt.Resize(this->mesh.GetNbBoundary());
	change_elt.Fill(true);
	this->AddSurfaceFiniteElement(name_elt, order, change_elt, this->mesh, this->mesh_num,
				      this->order_over_integration, this->dg_formulation,
                                      this->type_element);
      }
    
    int type_integration_edge, type_integration_tri, type_integration_quad;
    this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, this->other_mesh_num,
                                     type_integration_edge,
				     type_integration_tri, type_integration_quad);
    
    if (Dimension::dim_N == 2)
      {
	this->mesh_num.number_map.ConstructQuadrature2D(order, type_integration_edge);
        for (int k = 0; k < this->other_mesh_num.GetM(); k++)
          this->other_mesh_num(k)->number_map.ConstructQuadrature2D(order, type_integration_edge);
      }
    else
      {
#ifdef MONTJOIE_WITH_THREE_DIM        
        this->mesh_num.number_map.ConstructQuadrature3D(order, type_integration_tri, type_integration_quad);
        for (int k = 0; k < this->other_mesh_num.GetM(); k++)
          this->other_mesh_num(k)->number_map.ConstructQuadrature3D(order, type_integration_tri, type_integration_quad);
#endif
      }
    
  }
  

  //! Updates finite elements
  template<class Dimension>
  void VarProblem<Dimension>::UpdateInterpolationElement()
  {
    VarFiniteElement<Dimension>::UpdateInterpolationElement(this->mesh, this->mesh_num);

    // retrieves the type of quadrature for edges or triangles/quadrangles
    TinyVector<IVect, 4> order;    
    this->mesh_num.GetOrderQuadrature(order);
    RemoveDuplicate(order(0));
    
    int type_integration_edge, type_integration_tri, type_integration_quad;
    this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, this->other_mesh_num,
                                     type_integration_edge,
				     type_integration_tri, type_integration_quad);
    
    if (Dimension::dim_N == 2)
      {
	this->mesh_num.number_map.ConstructQuadrature2D(order, type_integration_edge);
        for (int k = 0; k < this->other_mesh_num.GetM(); k++)
          this->other_mesh_num(k)->number_map.ConstructQuadrature2D(order, type_integration_edge);
      }
    else
      {
#ifdef MONTJOIE_WITH_THREE_DIM        
        this->mesh_num.number_map.ConstructQuadrature3D(order, type_integration_tri, type_integration_quad);
        for (int k = 0; k < this->other_mesh_num.GetM(); k++)
          this->other_mesh_num(k)->number_map.ConstructQuadrature3D(order, type_integration_tri, type_integration_quad);
#endif
      }
  
  }
  

  //! clears finite element objects
  template<class Dimension>
  void VarProblem<Dimension>::ClearFiniteElement()
  {
    VarFiniteElement<Dimension>::ClearFiniteElement();
  }
  

  
  //! checks if continuity is ensured between elements  
  template<class Dimension>
  void VarProblem<Dimension>::CheckContinuity()
  {
#ifdef SELDON_DEBUG_LEVEL_4
    // we check continuity of basis functions across boundaries    
    SetPoints<Dimension> PointsElem, PointsElem2;
    SetMatrices<Dimension> MatricesElem, MatricesElem2;
    for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
      {
	if (this->mesh.Boundary(i).GetNbElements() == 2)
	  {
	    int num_elem = this->mesh.Boundary(i).numElement(0);
	    int num_elem2 = this->mesh.Boundary(i).numElement(1);
            switch(this->type_element)
              {
              case 1:
                this->CheckLocalContinuity(i, num_elem, num_elem2,
                                           PointsElem, PointsElem2,
                                           MatricesElem, MatricesElem2,
                                           this->GetReferenceElementH1(num_elem),
                                           this->GetReferenceElementH1(num_elem2));
                break;
              case 2:
                this->CheckLocalContinuity(i, num_elem, num_elem2,
                                           PointsElem, PointsElem2,
                                           MatricesElem, MatricesElem2,
                                           this->GetReferenceElementHcurl(num_elem),
                                           this->GetReferenceElementHcurl(num_elem2));
                break;
              case 3:
                this->CheckLocalContinuity(i, num_elem, num_elem2,
                                           PointsElem, PointsElem2,
                                           MatricesElem, MatricesElem2,
                                           this->GetReferenceElementHdiv(num_elem),
                                           this->GetReferenceElementHdiv(num_elem2));
                break;
              }
            
          }     
      }
#endif
  }

  //! computes position of degrees of freedom (for interpolatory elements)
  template<class Dimension>
  void VarProblem<Dimension>::ComputeDofCoordinates(Vector<R_N>& PointsDof)
  {
    PointsDof.Reallocate(this->mesh_num.GetNbDof());
    SetPoints<Dimension> PointsElem;
    Vector<R_N> s;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
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
  }
  
  
  //! computes J_i, DF_i^-1, etc
  template<class Dimension>
  void VarProblem<Dimension>::ComputeGeometryQuantity(IVect& OrderFace)
  {
    // loop on the elements of the mesh
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
        int n = 0; int rmax = this->GetReferenceElement(i).GetOrder();
        for (int j = 0; j < this->other_mesh_num.GetM(); j++)
          {
            int r = this->GetReferenceElement(i, j+1).GetOrder();
            if (r > rmax)
              {
                rmax = r;
                n = j+1;
              }
          }
        
        const ElementReference_Dim<Dimension>& Fb = this->GetReferenceElement(i, n);
        this->ComputeLocalMassMatrix(i, Fb.GetNbPointsQuadratureInside(),
                                     Fb.LinearSparseMassMatrix(),
                                     PointsElem, MatricesElem, OrderFace, Fb.GetGeometricElement());
      }
  }
  

  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeReferenceGradientElement(int i, const VectReal_wp& JacobDof, const VectReal_wp& JacobQuadrature,
				    VectReal_wp& GradJacobQuad)
  {
    const ElementReference_Dim<Dimension>& Fb = this->GetReferenceElement(i);
    
    // computing gradient of Ji on reference element
    if (Fb.UseQuadraturePointsForSh())
      Fb.ApplyRhQuadratureTranspose(JacobQuadrature, GradJacobQuad);
    else
      {
	Fb.ApplyRhTranspose(JacobDof, GradJacobQuad);
	// computation of gradient of jacobian
	//typename Dimension::VectR_N GradJacobNodal(nb_points_nodal);
	//Fb.ComputeNodalGradientRef(JacobNodal, GradJacobNodal);
        
	// interpolation on quadrature points
	//GradJacobQuad.Reallocate(nb_points_quadrature);
	//VectReal_wp dJ_dx(nb_points_nodal), dJ_dx_quad(nb_points_quadrature);
	//for (int k = 0; k < Dimension::dim_N; k++)
	//{
	//  for (int j = 0; j < nb_points_nodal; j++)
	///    dJ_dx(j) = GradJacobNodal(j)(k);
	// 
	//  Fb.ApplyChNodalTranspose(dJ_dx, dJ_dx_quad);
	//  for (int j = 0; j < nb_points_quadrature; j++)
	//    GradJacobQuad(j*Dimension::dim_N + k) = dJ_dx_quad(j);
	//}
      }
  }
  
  
  //! computes geometrical values, physical indexes, etc
  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeMassMatrix(bool compute_rho, bool delete_points)
  {
    VarGeometryProblem<Dimension>::ComputeMassMatrix(compute_rho, delete_points);

#ifdef SELDON_WITH_MPI
    // taking the maximum of penalty coefficients among processors
    this->ReduceDistributedVectorFace(this->Glob_CoefPenalDG, MPI_MAX, 1);
#endif
    
  }

  
  //! returns the type of the mass matrix
  /*!
    the different types of matrices are in the class FemMassMatrix
    (diagonal, block-diagonal, etc)
   */
  template<class Dimension>
  int VarProblem<Dimension>::GetMassMatrixType(Vector<bool>& diag_elt) const
  {
    bool diagonal = true, semi_diagonal = false;
    bool block_diagonal = true, semi_block_diagonal = false;
    diag_elt.Reallocate(this->mesh.GetNbElt());
    diag_elt.Fill(false);
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb = this->GetReferenceElement(i);
	if (this->DiagonalScalarMassMatrix(Fb, i))
	  {
	    semi_diagonal = true;
	    diag_elt(i) = true;
	  }
	else
	  {
	    diagonal = false;
	    if (this->BlockDiagonalScalarMassMatrix(Fb, i))
	      {
		semi_block_diagonal = true;
		diag_elt(i) = true;
	      }
	    else
	      block_diagonal = false;
	  }
      }
    
#ifdef SELDON_WITH_MPI
    bool diag_tmp(diagonal), semi_diag_tmp(semi_diagonal);
    bool block_diag_tmp(block_diagonal), semi_block_diag_tmp(semi_block_diagonal);
    MPI_Allreduce(&diag_tmp, &diagonal, 1, MPI_CHAR, MPI_LAND, this->comm_group_mode);
    MPI_Allreduce(&semi_diag_tmp, &semi_diagonal, 1, MPI_CHAR, MPI_LOR, this->comm_group_mode);
    MPI_Allreduce(&block_diag_tmp, &block_diagonal, 1, MPI_CHAR, MPI_LAND, this->comm_group_mode);
    MPI_Allreduce(&semi_block_diag_tmp, &semi_block_diagonal, 1, MPI_CHAR, MPI_LOR, this->comm_group_mode);
#endif
    
    if (diagonal)
      return FemMassMatrix::DIAGONAL;

    if (block_diagonal)
      return FemMassMatrix::BLOCK_DIAGONAL_UNSYM;
    
    if (semi_block_diagonal)
      return FemMassMatrix::BLOCK_DIAG_MATRIX_FREE;

    if (semi_diagonal)
      return FemMassMatrix::DIAG_MATRIX_FREE;
    
    // last case, no element with mass lumping
    return FemMassMatrix::MATRIX_FREE;
  }


  //! returns the type of the elementary matrix (0 : dense, 1 : sparse)
  template<class Dimension>
  int VarProblem<Dimension>::GetElementaryMatrixType() const
  {
    if (Dimension::dim_N == 2)
      return 0;
    
    bool sparse_elem = true;
    for (int i = 0; i < this->mesh.GetNbElt(); i++)
      {
	const ElementReference_Dim<Dimension>& Fb = this->GetReferenceElement(i);
        if (!Fb.LumpedMassMatrix())
          sparse_elem = false;
      }

    if (sparse_elem)
      return 1;

    return 0;
  }


  //! returns true if the elementary matrix is symmetric
  template<>
  bool VarProblem<Dimension2>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const
  {    
    return symmetric_elementary_matrix;
  }

  //! returns true if the elementary matrix is symmetric
  template<>
  bool VarProblem<Dimension2>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const
  {    
    return symmetric_elementary_matrix;
  }


  //! returns true if the elementary matrix is symmetric
  template<>
  bool VarProblem<Dimension3>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const
  {
    for (int n = 0; n < this->GetNbMeshNumberings(); n++)
      {
        const MeshNumbering<Dimension3>& mesh_num = this->GetMeshNumbering(n);
        for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
          {
            int rf = mesh_num.GetOrderQuadrature(i);
            if (!mesh_num.number_map.DofInvariantByRotation(rf, this->mesh.Boundary(i)))
              return false;
          }
      }
    
    return symmetric_elementary_matrix;
  }


  //! returns true if the elementary matrix is symmetric
  template<>
  bool VarProblem<Dimension3>::IsSymmetricElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const
  {
    for (int n = 0; n < this->GetNbMeshNumberings(); n++)
      {
        const MeshNumbering<Dimension3>& mesh_num = this->GetMeshNumbering(n);
        for (int i = 0; i < this->mesh.GetNbBoundary(); i++)
          {
            int rf = mesh_num.GetOrderQuadrature(i);
            if (!mesh_num.number_map.DofInvariantByRotation(rf, this->mesh.Boundary(i)))
              return false;
          }
      }
    
    return symmetric_elementary_matrix;
  }

  
  //! returns true if the elementary mass matrix is diagonal
  template<class Dimension>
  bool VarProblem<Dimension>::IsSparseElementaryMatrix(const GlobalGenericMatrix<Real_wp>&) const
  {
    return false;
    if (this->GetElementaryMatrixType() == 1)
      return true;
    
    return this->sparse_elementary_matrix;
  }
  

  //! returns true if the elementary mass matrix is sparse
  template<class Dimension>
  bool VarProblem<Dimension>::IsSparseElementaryMatrix(const GlobalGenericMatrix<Complex_wp>&) const
  {
    return false;
    if (this->GetElementaryMatrixType() == 1)
      return true;
        
    return this->sparse_elementary_matrix;
  }

  
  //! Modifies the object with a line of the data file
  /*!
    \param[in] description_field keyword of the input file
    \param[in] parameters list of parameters associated to description_field
    \param[in] nb_param number of parameters
    \return 0 if the field has been found, -1 otherwise
   */
  template<class Dimension>
  void VarProblem<Dimension>::
  SetInputData(const string& description_field,
	       const VectString& parameters)
  {    
    VarComputationProblem::SetInputData(description_field, parameters);
    DistributedProblem<Dimension>::SetInputData(description_field, parameters);
    VarFiniteElement<Dimension>::SetInputData(description_field, parameters);
  }
    
  
  //! computes the number of degrees of freedom for the problem
  /*!
    This method should be overloaded if you want to add degrees
    of freedom for a new model or new boundary conditions for example
   */
  template<class Dimension>
  void VarProblem<Dimension>::ComputeNumberOfDofs()
  {
    if (this->dg_formulation == ElementReference_Base::HDG)
      {
	int nodl_scalar = this->mesh_num.GetNbDof();

	this->OffsetDofV.Reallocate(this->mesh.GetNbElt()+1);
        this->OffsetDofV(0) = 0;
	for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  this->OffsetDofV(i+1) = this->OffsetDofV(i) + this->GetReferenceElement(i).GetNbDof();
 	
	this->offset_dof_unknown.Reallocate(this->nb_unknowns + this->nb_unknowns_hdg + 1);
	this->offset_dof_condensed.Reallocate(this->nb_unknowns + this->nb_unknowns_hdg + 1);
	for (int m = 0; m <= this->nb_unknowns_hdg; m++)
	  {
            this->offset_dof_unknown(m) = m*nodl_scalar;
            this->offset_dof_condensed(m) = m*nodl_scalar;
          }
        
	int Nvol = this->OffsetDofV(this->mesh.GetNbElt());
	for (int m = 0; m < this->nb_unknowns; m++)
	  {
            this->offset_dof_unknown(m+this->nb_unknowns_hdg+1)
              = this->offset_dof_unknown(m+this->nb_unknowns_hdg) + Nvol;

            this->offset_dof_condensed(m+this->nb_unknowns_hdg+1) = this->offset_dof_condensed(this->nb_unknowns_hdg);
          }
        
	this->nodl = this->offset_dof_unknown(this->nb_unknowns + this->nb_unknowns_hdg);
      }
    else
      {
        this->offset_dof_unknown.Reallocate(this->nb_unknowns+1);
        this->offset_dof_condensed.Reallocate(this->nb_unknowns+1);
        this->offset_dof_unknown.Zero();
        this->offset_dof_condensed.Zero();
        
        IVect nb_ddl_cond(this->GetNbMeshNumberings());
        if (this->dg_formulation == ElementReference_Base::CONTINUOUS)
          {
            nb_ddl_cond(0) = this->mesh_num.OffsetDofVertexNumber.GetNormInf();
            nb_ddl_cond(0) = max(nb_ddl_cond(0), this->mesh_num.OffsetDofEdgeNumber.GetNormInf());
            nb_ddl_cond(0) = max(nb_ddl_cond(0), this->mesh_num.OffsetDofFaceNumber.GetNormInf());
            
            for (int k = 0; k < this->other_mesh_num.GetM(); k++)
              {
                nb_ddl_cond(k+1) = this->other_mesh_num(k)->OffsetDofVertexNumber.GetNormInf();
                nb_ddl_cond(k+1) = max(nb_ddl_cond(k+1), this->other_mesh_num(k)->OffsetDofEdgeNumber.GetNormInf());
                nb_ddl_cond(k+1) = max(nb_ddl_cond(k+1), this->other_mesh_num(k)->OffsetDofFaceNumber.GetNormInf());
              }
            
            this->offset_dof_condensed(0) = 0;
            for (int i = 0; i < this->nb_unknowns; i++)
              this->offset_dof_condensed(i+1) = this->offset_dof_condensed(i)
                + nb_ddl_cond(this->mesh_num_unknown(i));
          }
        
        this->nodl = 0;
        for (int i = 0; i < this->nb_unknowns; i++)
          {
            int p = this->mesh_num_unknown(i);
            if (p == 0)
              this->nodl += this->mesh_num.GetNbDof();
            else
              this->nodl += this->other_mesh_num(p-1)->GetNbDof();
            
            this->offset_dof_unknown(i+1) = this->nodl;
          }
      }
    
    if ((this->FirstOrderFormulation()) && (this->dg_formulation == ElementReference_Base::CONTINUOUS))
      {
        this->OffsetDofV.Reallocate(this->mesh.GetNbElt()+1);
        this->OffsetDofV(0) = 0;
	for (int i = 0; i < this->mesh.GetNbElt(); i++)
	  this->OffsetDofV(i+1) = this->OffsetDofV(i)
            + this->GetNbPointsQuadratureInside(i)*this->nb_unknowns_vec;
        
        this->nodl += this->OffsetDofV(this->mesh.GetNbElt());
      }        
  }


  //! Adding parallel dofs 
  /*!
    In the derived class, this method will update MatchingDofOrig_Subdomain
    and GlobDofNumber_Subdomain in order to have a consistent definition
    of dofs in parallel.
    This method is only called in parallel.
   */
  template<class Dimension>
  void VarProblem<Dimension>::PutOtherGlobalDofs()
  {
    
  }


  //! Computes E \times n and H \times on quadrature points
  /*!
    For Maxwell's equations, it computes E \times n and H \times, for Helmholtz u and du/dn
    and for other equations it should computed what is needed for integral representation
    \param[in] u_quadrature values of u (E for Maxwell) on quadrature points
    \param[in] grad_quadrature values of \grad u (curl E for Maxwell) on quadrature points
    \param[in] num_elem element number
    \param[in] pts quadrature points
    \param[in] normale normale for each quadrature point
    \param[in] compute_H if true H \times n must be computed as well
    \param[out] En_quad values of E \times n on quadrature points
    \param[out] Hn_quad values of H \times n on quadrature points
    
    u_quadrature(m)(j)(:) is the solution u for the unknown m and the quadrature point j
    it can be vectorial for H(curl)/  H(div) elements
    
    grad_quadrature(m)(j)(:) is the gradient (or curl or divergence) of u for the
    unknown m and the quadrature point j
    it is the gradient for H^1 elements (type=1), the curl for H(curl) elements (type=2)
    the divergence for H(div) elements (type=3)
    
    The values of E \times n and H \times n are stored in simple vectors En_quad, Hn_quad
    If E \times n contains three components for example, these values are stored
    in En_quad(3*j), En_quad(3*j+1), En_quad(3*j+2) where j is the quadrature point number
   */
  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeEnHnQuadrature(Vector<VectReal_wp>& u_quadrature,
			  Vector<VectReal_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, VectReal_wp& En_quad, VectReal_wp& Hn_quad) const
  {
    // this method should be overloaded for each considered equation
    cout << "ComputeEnHnQuadrature not available for this equation" << endl;
    abort();
  }

  

  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeEnHnQuadrature(Vector<VectComplex_wp>& u_quadrature,
			  Vector<VectComplex_wp>& grad_quadrature,
			  int num_elem, const Vector<R_N>& pts, const Vector<R_N>& normale,
			  bool compute_H, VectComplex_wp& En_quad, VectComplex_wp& Hn_quad) const
  {
    // this method should be overloaded for each considered equation
    cout << "ComputeEnHnQuadrature not available for this equation" << endl;
    abort();
  }


  //! computes local prolongation
  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeLocalProlongation(FiniteElementInterpolator& proj,
			     DistributedProblem<Dimension>& var_coarse_base, int rc, 
			     TinyVector<Matrix<Real_wp>, 4>& ProlongationElement) const
  {
    VarProblem<Dimension>& var_coarse
      = static_cast<VarProblem<Dimension>& >(var_coarse_base);
    
    if (this->FormulationDG() == ElementReference_Base::HDG)
      {
        typedef typename Dimension::DimensionBoundary DimensionB;
	Vector<const ElementReference_Dim<DimensionB>* > elt_fine;
	Vector<const ElementReference_Dim<DimensionB>* > elt_coarse;
	
	this->GetSurfaceFiniteElement(elt_fine);
	var_coarse.GetSurfaceFiniteElement(elt_coarse);
	
	proj.ComputeLocalProlongation(elt_fine, this->mesh_num.GetOrder(),
                                      elt_coarse, rc, ProlongationElement);
      }
    else
      {
        Vector<const ElementReference_Dim<Dimension>* > elt_fine;
        Vector<const ElementReference_Dim<Dimension>* > elt_coarse;
	
        this->GetReferenceElement(elt_fine);
	var_coarse.GetReferenceElement(elt_coarse);
	
	proj.ComputeLocalProlongation(elt_fine, this->mesh_num.GetOrder(),
				      elt_coarse, rc, ProlongationElement);
      }
  }  


  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Real_wp>&,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>&)
  {
    cout << "The method ComputeElementaryMatrix should be overloaded for each equation" << endl;
    cout << "You can call ComputeElementaryMatrixH1 for H^1 or DG elements" << endl;
    abort();
  }
  
  
  template<class Dimension>
  void VarProblem<Dimension>
  ::ComputeElementaryMatrix(int, IVect&, VirtualMatrix<Complex_wp>&,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>&)
  {
    cout << "The method ComputeElementaryMatrix should be overloaded for each equation" << endl;
    cout << "You can call ComputeElementaryMatrixH1 for H^1 or DG elements" << endl;
    abort(); 
  }
  
  
  //! copies input parameters of another similar problem
  template<class Dimension>
  void VarProblem<Dimension>::CopyInputData(const VarProblem_Base& var_base)
  {
    const VarProblem<Dimension>& var = static_cast<const VarProblem<Dimension>& >(var_base);
    
    DistributedProblem<Dimension>::CopyInputData(var);
    VarComputationProblem::CopyInputData(var);
  }
  
}
  
#define MONTJOIE_FILE_VAR_PROBLEM_CXX
#endif
