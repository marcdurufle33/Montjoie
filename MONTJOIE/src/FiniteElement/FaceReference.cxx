#ifndef MONTJOIE_FILE_FACE_REFERENCE_CXX

namespace Montjoie
{
  
  //! Destructor
  template<int type>
  FaceReference<type>::~FaceReference()
  {
    if (element_surface != NULL)
      {
	delete element_surface;
	element_surface = NULL;
      }
  }
  

  //! returns 1-D finite element associated with traces of the 2-D finite element space on any edge
  template<int type>
  const ElementReference<Dimension1, 1>&
  FaceReference<type>::GetSurfaceFiniteElement(int n) const
  {
    if (element_surface == NULL)
      {
	cout << "Null pointer for GetSurfaceFiniteElement" << endl;
	abort();
      }
    
    return *element_surface;
  }
  
  
  //! returns size of memory used by the object
  template<int type>
  size_t FaceReference<type>::GetMemorySize() const
  {
    size_t taille = ElementReference_Dim<Dimension2>::GetMemorySize();
    taille += ElementReferenceType<Dimension2, type>::GetMemorySize();
    taille += sizeof(int)*(EdgesDof.GetDataSize());
    taille += sizeof(int)*(EdgesNodal.GetDataSize());
    for (int i = 0; i < ProjOperatorOrder.GetM(); i++)
      taille += sizeof(Real_wp)*ProjOperatorOrder(i).GetDataSize();
    
    taille += sizeof(*this) - sizeof(ElementReference_Dim<Dimension2>);
    if (element_surface != NULL)
      taille += element_surface->GetMemorySize();
    
    return taille;
  }
  
  
  //! modification of the finite element with a line of the input file
  template<int type>
  void FaceReference<type>::SetInputData(const string& keyword, const Vector<string>& param)
  {
    // nothing (sigh)
  }
    

  //! computation of projection between finite element of different orders
  template<int type>
  void FaceReference<type>::
  SetVariableOrder(const Mesh<Dimension2>& mesh,
		   const MeshNumbering<Dimension2>& mesh_num)
  {
    if (!mesh_num.IsOrderVariable())
      return;
    
    int rmax = 0;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      rmax = max(rmax, mesh_num.GetOrderElement(i));
    
    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      rmax = max(rmax, mesh_num.GetOrderQuadrature(i));
    
    Vector<bool> order_present(rmax + 1);
    order_present.Fill(false);
    int type_elt = this->GetHybridType();
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	int r = mesh_num.GetOrderElement(i);
        if ((r == this->order) && (type_elt == mesh.Element(i).GetHybridType()))
	  {
	    // we are on an element with the order of this finite element
	    for (int j = 0; j < mesh.Element(i).GetNbEdges(); j++)
	      {
		int ne = mesh.Element(i).numEdge(j);
		int re = mesh_num.GetOrderQuadrature(ne);
                if (re != r)
		  order_present(re) = true;
	      }
	  }
      }
    
    // elements of other processors
    for (int i = 0; i < mesh_num.GetNbNeighborElt(); i++)
      {
        int r = mesh_num.GetOrderNeighborElement(i);
        if ( (r == this->order) && (type_elt = mesh_num.GetTypeNeighborElement(i)))
          {
            int ne = mesh_num.GetEdgeNeighborElement(i);
            int re = mesh_num.GetOrderQuadrature(ne);
            
            if (re != r)
              order_present(re) = true;
          }
      }
    
    this->ComputeVariableInterpolation(rmax, order_present);
  }
  
  
  //! forces computation of projection with orders present in array order_present
  template<int type>
  void FaceReference<type>::ComputeVariableInterpolation(int rmax,
							 Vector<bool>& order_present)
  {
    int nb_order = 0;
    for (int r = 1; r <= rmax; r++)
      if (order_present(r))
	nb_order++;
    
    if (nb_order > 0)
      {
	IVect order_elt(nb_order);
	Vector<VectReal_wp> points_quad(nb_order);
	int type_quad = this->GetTypeIntegrationEdge();
	nb_order = 0;
	for (int r = 1; r <= rmax; r++)
	  if (order_present(r))
	    {
	      Globatto<Real_wp> gauss;
	      gauss.ConstructQuadrature(r, type_quad);
	      points_quad(nb_order) = gauss.Points();
	      order_elt(nb_order) = r;
	      nb_order++;
	    }
	
	this->ComputeInterpolationProjectorOrder(order_elt, this->GetSurfaceFiniteElement(0),
						 points_quad);
      }
  }


  //! method to compute projector on other quadrature points
  template<int type>
  void FaceReference<type>
  ::ComputeInterpolationProjectorOrder(const IVect& order_elt,
				       const ElementReference<Dimension1, 1>& Fb,
				       const Vector<VectReal_wp>& Pts)
  {
    cout << "Method not implemented for type = " << type << endl;
    abort();
  }
  

  //! computation of u on nodal points of the element
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points
   */
  template<int type> template<class Vector1, class Vector2>
  void FaceReference<type>
  ::ComputeNodalValuesGen(const Vector1& Un, Vector2 & Unode) const
  {
    Unode.Reallocate(this->ValuePhi_Nodal.GetN());
    Mlt(SeldonTrans, this->ValuePhi_Nodal, Un, Unode);
  }
  
            
  //! displays informations about the class FaceReference<type>
  template<int type>
  ostream& operator <<(ostream& out, const FaceReference<type>& e)
  {
    out<<"Order of approximation "<<e.order<<endl;
    return out;
  }


  /****************
   * H^1 elements *
   ****************/
  

#ifdef MONTJOIE_WITH_NODAL_H1
  //! computation of values of u on nodal points of an edge
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points of edge num_loc
    \param[in] num_loc edge number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(this->EdgesNodal.GetM());
    for (int j = 0; j < this->EdgesNodal.GetM(); j++)
      {
        Unode(j) = 0;
        int node = EdgesNodal(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Unode(j) += Un(k)*this->ValuePhi_Nodal(k, node);
      }
  }
  

  //! computation of gradient of u on nodal points of an edge
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode gradient of u on nodal points of edge num_loc
    \param[in] num_loc edge number
  */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>
  ::ComputeGradientBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(2*this->EdgesNodal.GetM());
    for (int j = 0; j < this->EdgesNodal.GetM(); j++)
      {
        Unode(2*j) = 0; Unode(2*j+1) = 0;
        int node = this->EdgesNodal(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            Unode(2*j) += Un(k)*this->GradPhi_Nodal(k, node)(0);
            Unode(2*j+1) += Un(k)*this->GradPhi_Nodal(k, node)(1);
          }
      }
  }

  //! Integration against gradient of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
    where omega_k is the weight of integration
    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \nabla \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = \nabla \varphi_i(\xi_j)
  */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val = Real_wp(0);
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Gradient_Phi(i, j)(0)*Vh(2*j) + Gradient_Phi(i, j)(1)*Vh(2*j+1);
	
	Uh(i) = val;
      }
  }
  
  
  //! Integration against derivatives of basis functions
  /*!
    \param[in] Uh vector containing values \omega_k f(\xi_k)
                  where omega_k is the weight of integration
                  and \xi_k the point of integration
    \param[out] Vx Vx_i = \int_K f d/dx(\varphi_i)  dx
    \param[out] Vy Vy_i = \int_K f d/dy(\varphi_i)  dx
    This operation is equivalent to a matrix vector product
    Vx = Rh^1 Uh, Vy = Rh^2 Uh
    where (Rh)^1_{i,j} = d/dx(\varphi_i(\xi_j) )
    where (Rh)^2_{i,j} = d/dy(\varphi_i(\xi_j) )
  */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>
  ::ApplyRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy) const
  {
    typename Vector2::value_type val_x, val_y;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
        val_x = Real_wp(0);
	val_y = Real_wp(0);
        for (int j = 0; j < nb_points_quadrature_inside; j++)
          {
	    val_x += Gradient_Phi(i, j)(0)*Uh(2*j);
	    val_y += Gradient_Phi(i, j)(1)*Uh(2*j+1);
	  }
	
        Vx(i) = val_x;
        Vy(i) = val_y;
      }
  }


  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>
  ::ApplyRhQuadratureSplitGen(const Vector1& Uh, Vector2& Vh, Vector2&) const
  {
    cout << "not implemented for this finite element" << endl;
    abort();
  }    
  

  //! computation of gradient of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh gradient of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j grad phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = grad phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    TinyVector<typename Vector2::value_type, 2> grad;
    typename Vector2::value_type zero; SetComplexZero(zero);
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        grad.Fill(zero);
        for (int j = 0; j < this->nb_dof_loc; j++)
          Add(Uh(j), Gradient_Phi(j, i), grad);
        
        Vh(2*i) = grad(0); Vh(2*i+1) = grad(1);
      }
  }
  

  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class T1>
  void ElementReference<Dimension2, 1>::ApplyChGen(const Vector<T1>& Vh, Vector<T1>& Uh) const
  {
    T1 val;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	SetComplexZero(val);
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Value_Phi(i, j)*Vh(j);
        
        Uh(i) = val;
      }
  }


  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class T1>
  void ElementReference<Dimension2, 1>
  ::ApplyChTransposeGen(const Vector<T1>& Uh, Vector<T1>& Vh) const
  {
    T1 val;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	SetComplexZero(val);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  val += Value_Phi(j, i)*Uh(j);
	
	Vh(i) = val;
      }
  }
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh Uh
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>
  ::ApplyConstantRhGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(const_grad_matrix, Uh, Vh);
  }
  

  //! Multiplication by gradient matrix 
  /*!
    If we denote (Rh)_{i, j} = \int_K \varphi_j grad(\varphi_i) dx
    then Vh = Rh* Uh
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>
  ::ApplyConstantRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    Mlt(SeldonTrans, const_grad_matrix, Uh, Vh);
  }
  

  //! Multiplication by derivative matrices
  /*!
    If we denote (Rh^1)_{i, j} = \int_K \varphi_j d/dx(\varphi_i) dx
    (Rh^2)_{i, j} = \int_K \varphi_j d/dy(\varphi_i) dx
    then Vx = Rh^1 Uh, Vy = Rh^2 Uh
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 1>::
  ApplyConstantRhSplitGen(const Vector1& Uh, Vector2& Vx, Vector2& Vy) const
  {
    int m = const_grad_matrix.GetM();
    int n = const_grad_matrix.GetN();
    typename Vector2::value_type val;
    for (int i = 0; i < m; i++)
      {
	val = 0;
	for (int j = 0; j < n; j += 2)
	  val += const_grad_matrix(i, j)*Uh(j);
	
	Vx(i) = val;

	val = 0;
	for (int j = 1; j < n; j += 2)
	  val += const_grad_matrix(i, j)*Uh(j);
	
	Vy(i) = val;
      }
  }    
    

  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T>
  void ElementReference<Dimension2, 1>::
  AddConstantMassMatrixGen(int m, int n, const T& mass, VirtualMatrix<T>& A) const
  {
    T zero; SetComplexZero(zero);
    if (mass != zero)
      {
	Vector<int> num(this->nb_dof_loc);
	Vector<T> val(this->nb_dof_loc);
	for (int i = 0; i < this->nb_dof_loc; i++)
	  num(i) = n + i;
	
	for (int i = 0; i < this->nb_dof_loc; i++)
	  {
	    for (int j = 0; j < this->nb_dof_loc; j++)
	      val(j) = mass*mass_matrix(i, j);
	    
	    A.AddInteractionRow(m + i, this->nb_dof_loc, num, val);
	  }
      }
  }


  //! we add constant elementary matrix to the matrix A
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + mass M + \sum_p  E(p) (R^p)* + \sum_p  D(p) R^p + \sum_p  C(p, q) S^{p,q}
    where M is the mass matrix equal to
    M_{i, j} = \int_K \phi_j \phi_i dx
    R^p is a gradient matrix equal to :
    (R^p)_{i, j} = \int_K \phi_j d(\phi_i)/dx_p  dx
    S^{p,q} is a stiffness matrix equal to :
    (S^{p,q })_{i, j} = \int_K d(\phi_j)/dx_q d(\phi_i)/dx_p  dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension2, 1>::
  AddConstantElemMatrixGen(int m, int n, const T& mass,
			   const TinyMatrix<T, Prop, 2, 2>& C,
			   const TinyVector<T, 2>& D, 
			   const TinyVector<T, 2>& E,
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val.Zero();
	if (!null_term(0))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += mass*mass_matrix(i, j);
	
	if (!null_term(1))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += C(0, 0)*const_stiff_matrix(2*j, 2*i)
	      + C(0, 1)*const_stiff_matrix(2*j+1, 2*i)
	      + C(1, 0)*const_stiff_matrix(2*j, 2*i+1)
	      + C(1, 1)*const_stiff_matrix(2*j+1, 2*i+1);
	
	if (!null_term(2))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += D(0)*const_grad_matrix(i, j*2)
	      + D(1)*const_grad_matrix(i, j*2+1);
	
	if (!null_term(3))
	  for (int j = 0; j < this->nb_dof_loc; j++)
	    val(j) += E(0)*const_grad_matrix(j, i*2)
	      + E(1)*const_grad_matrix(j, i*2+1);
	
	A.AddInteractionRow(m + i, this->nb_dof_loc, num, val);
      }
  }
   
  
  //! adds \int C \nabla varphi_j \nabla varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] mass mass coefficient
    \param[in] C stiffness coefficient
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, mass, C, D or E is null (i=0, 1, 2, 3)
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + R + Rt + S
    M_ij = \int mass \varphi_j \varphi_i
    R_ij = \int D \cdot \nabla \varphi_i \varphi_j
    Rt_ij = \int E \cdot \nabla \varphi_j \varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
   */
  template<class T, class Prop>
  void ElementReference<Dimension2, 1>
  ::AddVariableElemMatrixGen(int off_row, int off_col, const Vector<T>& mass,
			     const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			     const Vector<TinyVector<T, 2> >& D,
			     const Vector<TinyVector<T, 2> >& E, 
			     const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc;
    TinyVector<T, 2> vec_u;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    SetComplexZero(vloc);
	    if (!null_term(0))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		vloc += this->Value_Phi(i, k)*this->Value_Phi(j, k)*mass(k);
	    
	    if (!null_term(1))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		{
		  Mlt(C(k), this->Gradient_Phi(j, k), vec_u);
		  vloc += DotProd(vec_u, this->Gradient_Phi(i, k));
		}
	    
	    if (!null_term(2))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		vloc += this->Value_Phi(j, k)*DotProd(D(k), this->Gradient_Phi(i, k));

	    if (!null_term(3))
	      for (int k = 0; k < this->nb_points_quadrature_inside; k++)
		vloc += this->Value_Phi(i, k)*DotProd(E(k), this->Gradient_Phi(j, k));
	    
	    val(j) = vloc;
	  }
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
   */
  template<class T>
  void ElementReference<Dimension2, 1>
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val.Zero();
	for (int j = 0; j < this->nb_dof_loc; j++)
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	    val(j) += this->Value_Phi(i, k)*this->Value_Phi(j, k)*A(k);
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }


  //! adds \int C \nabla varphi_j \nabla varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A mass coefficients
    \param[in] C tensor C on each quadrature point
    \param[in] D gradient coefficients
    \param[in] E gradient coefficients
    \param[in] null_term if null_term(i) is true, A, C, D or E is null (i=0, 1, 2, 3)
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + M + R + Rt + S
    M_ij = \int mass \varphi_j \varphi_i
    R_ij = \int D \cdot \nabla \varphi_i \varphi_j
    Rt_ij = \int E \cdot \nabla \varphi_j \varphi_i
    where S_ij = \int C \nabla varphi_j \nabla varphi_i
    This function uses function ApplyCh, ApplyRh, considering that
    they are optimized
   */
  template<class T, class Prop>
  void ElementReference<Dimension2, 1>::
  AddVariableElemMatrixOpt(int off_row, int off_col, const Vector<T>& A,
			   const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			   const Vector<TinyVector<T, 2> >& D,
			   const Vector<TinyVector<T, 2> >& E, 
			   const TinyVector<bool, 4>& null_term, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    int Nquad = this->nb_points_quadrature_inside;
    Vector<Real_wp> val_phi(Nquad); VectR2 grad_phi(Nquad);
    Vector<T> contrib(this->nb_dof_loc), feval_grad(2*Nquad), feval(Nquad), feval_bis(Nquad);
    R2 vec_u; TinyVector<T, 2> vec_v;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	this->GetGradientSinglePhiQuadrature(i, val_phi, grad_phi);
	
	feval.Zero();
        for (int k = 0; k < Nquad; k++)
          {
            vec_u = grad_phi(k);
            
	    if (!null_term(0))
	      feval(k) = A(k)*val_phi(k);

            if (!null_term(1))
	      MltTrans(C(k), vec_u, vec_v);
	    else
	      vec_v.Zero();
	    
	    if (!null_term(2))
	      feval(k) += DotProd(D(k), vec_u);

	    if (!null_term(3))
	      {
		vec_v(0) += E(k)(0)*val_phi(k);
		vec_v(1) += E(k)(1)*val_phi(k);
	      }
	    
            feval_grad(2*k) = vec_v(0);
            feval_grad(2*k+1) = vec_v(1);

          }
	
	if (this->UseQuadraturePointsForRh())
	  {
	    this->ApplyRhQuadrature(feval_grad, feval_bis);
	    feval += feval_bis;
	    this->ApplyCh(feval, val);	
	  }
	else
	  {
	    this->ApplyRh(feval_grad, contrib);	
	    this->ApplyCh(feval, val);	
	    val += contrib;
	  }	    
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }
        

  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A coefficient A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
    This function uses function ApplyCh, ApplyRh, considering that
    they are optimized
   */
  template<class T>
  void ElementReference<Dimension2, 1>::
  AddVariableMassMatrixOpt(int off_row, int off_col,
                           const Vector<T>& A, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    int Nquad = this->nb_points_quadrature_inside;
    Vector<Real_wp> val_phi(Nquad);
    Vector<T> feval(Nquad);
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	this->GetValueSinglePhiQuadrature(i, val_phi);
        for (int k = 0; k < Nquad; k++)
          feval(k) = val_phi(k)*A(k);
        
        this->ApplyCh(feval, val);
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }

  
  //! Integration against gradient of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension2, 1>
  ::ComputeIntegralGradientRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyRh(feval, res);
  }
  
  
  //! Surface integrals against gradient of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension2, 1>
  ::ComputeIntegralSurfaceGradientRef(const Vector1 & feval, Vector2& res, int n) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    ApplyNablaSh(one, n, feval, res);
  }
  
  
  //! computation of projection between finite element of different orders  
  void ElementReference<Dimension2, 1>
  ::ComputeInterpolationProjectorOrder(const IVect& order_elt,
				       const ElementReference<Dimension1, 1>& Fb,
				       const Vector<VectReal_wp>& Pts)
  {
    if (order_elt.GetM() <= 0)
      return;
    
    int rmax = 0;
    for (int i = 0; i < order_elt.GetM(); i++)
      rmax = max(order_elt(i), rmax);
    
    VectReal_wp phi(this->order+1);
    ProjOperatorOrder.Reallocate(rmax+1);
    for (int i = 0; i < order_elt.GetM(); i++)
      {
	int r = order_elt(i);
	ProjOperatorOrder(r).Reallocate(this->order+1, Pts(i).GetM());
	for (int j = 0; j < Pts(i).GetM(); j++)
	  {
	    Fb.ComputeValuesPhiRef(Pts(i)(j), phi);
	    for (int k = 0; k <= this->order; k++)
	      ProjOperatorOrder(r)(k, j) = phi(k);
	  }
      }
  }

  //! Computes the effect of rotation on a triangular face on degrees of freedom
  /*!
    \param[in] rf order of approximation for the face
    \param[in] Points2D quadrature points on the face
    \param[in] Weights2D quadrature weights on the face
    \param[in] ValuePhi values of L2 basis functions on quadrature points
    \param[out] nmap numbering scheme of the mesh
    We compute which operators should be applied to degrees of freedom
    when a triangular face is rotated
   */
  void ElementReference<Dimension2, 1>::FindH1RotationTri(int rf, const VectR2& Points2D,
                                                          const VectReal_wp& Weights2D,
                                                          const Matrix<Real_wp>& ValuePhi,
                                                          NumberMap& nmap)
  {
    int nb_dof = ValuePhi.GetM();
    
    // loop on degrees of freedom
    Matrix<int> FaceRot(6, nb_dof);
    Matrix<bool> OppositeSigneDof(6, nb_dof);
    for (int i = 0; i < nb_dof; i++)
      for (int j = 0; j < 6; j++)
	{
	  FaceRot(j, i) = i;
	  OppositeSigneDof(j, i) = false;
	}
    
    int nb_points = Points2D.GetM();
    
    Vector<Matrix<Real_wp> > CoefCombination(6);
    Matrix<Real_wp> Mh(nb_dof, nb_dof), Bh(nb_dof, nb_dof);
    VectR2 Points(nb_points); IVect permut(nb_points);
    Matrix<int> permut_rot(6, nb_points);
    Mh.Fill(0);
    for (int rot = 0; rot < 6; rot++)
      {
	CoefCombination(rot).Reallocate(nb_dof, nb_dof);
	CoefCombination(rot).Fill(0);
	for (int i = 0; i < nb_points; i++)
	  {
	    Real_wp l1 = Points2D(i)(0), l2 = Points2D(i)(1);
	    Real_wp l0 = 1.0 - l1 - l2, x, y;
	    switch(rot)
	      {
	      case 0:
		x = l1; y = l2; break;
	      case 1 :
		x = l0; y = l1; break;
	      case 2 :
		x = l2; y = l0; break;
	      case 3 :
		x = l2; y = l1; break;
	      case 4 :
		x = l0; y = l2; break;
	      case 5 :
		x = l1; y = l0; break;
	      }
	    
	    Points(i).Init(x, y);
	  }
	
	permut.Fill();
	Sort(nb_points, Points, permut);
	for (int i = 0; i < nb_points; i++)
	  permut_rot(rot, i) = permut(i);
      }
    
    Vector<int> rot_dof(6), rot_quad(6);
    Real_wp phi, phi2;
    rot_quad.Fill(-1);    
    for (int i = 0; i < nb_dof; i++)
      {
        rot_dof.Fill(-1);
        for (int j = 0; j < nb_dof; j++)
          for (int k = 0; k < nb_points; k++)
            {
              rot_quad(0) = permut_rot(0, k);
              phi = ValuePhi(i, rot_quad(0));
              Mh(i, j) += Weights2D(k)*phi*ValuePhi(j, rot_quad(0));
              
              for (int rot = 1; rot <= 5; rot++)
                {
                  phi2 = ValuePhi(j, permut_rot(rot, k));
                  CoefCombination(rot)(i, j) += Weights2D(k)*phi*phi2;
                }
            }
      }
    
    GetInverse(Mh);
    
    Matrix<Real_wp> Ah(nb_dof, nb_dof);
    Ah.Fill(0);
    bool linear_combination = false;
    CoefCombination(0).SetIdentity();
    Vector<Matrix<Real_wp> > CoefCombin(CoefCombination);
    for (int rot = 1; rot < 6; rot++)
      {
	Bh = CoefCombin(rot);
	Mlt(Mh, Bh, Ah);
        
	int rot2 = rot;
	if (rot == 2)
	  rot2 = 1;
	else if (rot == 1)
	  rot2 = 2;
	
	//DISP(rot); DISP(Ah);
	Real_wp threshold = 1000*epsilon_machine;
	// pour tous les ddls de la face
	Vector<bool> DofUsed(nb_dof); DofUsed.Fill(false);            
	for (int i = 0; i < nb_dof; i++)
	  {
	    int nb_coef = 0;
	    for (int j = 0; j < nb_dof; j++)
	      if (abs(Ah(j, i)) > threshold)
		nb_coef++;
	    
	    if (nb_coef == 0)
	      {
		cout << "Cas impossible " << endl;
		abort();
	      }
	    else if (nb_coef == 1)
	      {
		// on est tombe sur un autre ddl, avec eventuellement un signe
		for (int j = 0; j < nb_dof; j++)
		  if (abs(Ah(j, i)) > threshold)
		    {
		      FaceRot(rot2, i) = j;
		      if ( abs(Ah(j, i) -Real_wp(1)) <= threshold)
			OppositeSigneDof(rot2, i) = false;
		      else if ( abs(Ah(j, i) + Real_wp(1)) <= threshold)
			OppositeSigneDof(rot2, i) = true;
		      else
			{
			  // this dof can't be handled with a sign
			  linear_combination = true;
			}
                                      
		      DofUsed(j) = true;
		    }                       
	      }
	    else
	      {
		// no direct relation between dofs
		linear_combination = true;
	      }
	  }
	
	CoefCombination(rot2) = Ah;
      }

    // DISP(linear_combination); DISP(FaceRot); DISP(OppositeSigneDof);
    if (linear_combination)
      nmap.SetFacesDofRotationTri(rf, CoefCombination);
    else
      {
	nmap.SetFacesDofRotationTri(rf, FaceRot);
	nmap.SetSignDofRotationTri(rf, OppositeSigneDof);
      }
  }
#endif
  
  
  /********************
   * H(curl) elements *
   ********************/
  

#ifdef MONTJOIE_WITH_NODAL_HCURL
  void ElementReference<Dimension2, 2>
  ::ComputeValuesPhiHRef(const R2& point_loc, VectReal_wp& res) const
  {
    cout << "Not implemented for this finite element" << endl;
    abort();
  }
  

  //! construction of mass, stiffness matrices for edge elements
  /*!
    \param[in] Fb leaf finite element class
   */
  void ElementReference<Dimension2, 2>::ConstructHcurlElementaryMatrix()
  {
    this->Value_PhiVec.Reallocate(this->nb_dof_loc, this->GetNbPointsQuadrature());
    this->Curl_Phi2D.Reallocate(this->nb_dof_loc, this->GetNbPointsQuadrature());
    VectR2 phi; VectReal_wp curl_phi;
    for (int i = 0; i < this->GetNbPointsQuadrature(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	this->ComputeCurlPhiRef(this->PointsND(i), curl_phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    this->Value_PhiVec(j, i) = phi(j);
	    this->Curl_Phi2D(j, i) = curl_phi(j);
	  }
      }
    
    // mass matrix on reference element
    this->mass_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    Real_wp vloc;
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = i; j < this->nb_dof_loc; j++)
	{
	  vloc = 0.0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*DotProd(this->Value_PhiVec(i, k), this->Value_PhiVec(j, k));
          
	  this->mass_matrix(i,j) = vloc;
	}
    
    //this->mass_matrix.WriteText("Mh.dat");
    this->mass_matrix_chol = this->mass_matrix;
    GetCholesky(this->mass_matrix_chol);
    
    // projection on dofs
    this->MatProjectionDof.Reallocate(2*this->GetNbPointsDofInside(), this->nb_dof_loc);
    for (int i = 0; i < this->GetNbPointsDofInside(); i++)
      {
        this->ComputeValuesPhiRef(this->PointsDofND(i), phi);
        for (int j = 0; j < this->nb_dof_loc; j++)
          {
            this->MatProjectionDof(2*i, j) = phi(j)(0);
            this->MatProjectionDof(2*i+1, j) = phi(j)(1);
          }
      }
    
    GetQR(this->MatProjectionDof, this->tauProjectionDof);
    
    // values of basis functions on nodal points
    this->ValuePhi_Nodal.Reallocate(this->nb_dof_loc, 2*this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
        this->ComputeValuesPhiRef(this->PointsNodalND(i), phi);
        for (int j = 0; j < this->nb_dof_loc; j++)
          {
            this->ValuePhi_Nodal(j, 2*i) = phi(j)(0);
            this->ValuePhi_Nodal(j, 2*i+1) = phi(j)(1);
          }
      }
    
    // curl of basis functions on nodal points
    this->CurlPhi_Nodal.Reallocate(this->nb_dof_loc, this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
        this->ComputeCurlPhiRef(this->PointsNodalND(i), curl_phi);
        for (int j = 0; j < this->nb_dof_loc; j++)
          this->CurlPhi_Nodal(j, i) = curl_phi(j);
      }
    
    // mass and stiffness matrix for constant coefficients
    // \int_K C phi_i phi_j   and   \int_K curl(phi_i)  curl(phi_j)
    // where C is a 2x2 general matrix 
    this->const_mass_matrix(0, 0).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    this->const_mass_matrix(0, 1).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    this->const_mass_matrix(1, 0).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    this->const_mass_matrix(1, 1).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    
    this->const_curl_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
        
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < this->nb_dof_loc; j++)
        {
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i, k)(0)*this->Value_PhiVec(j, k)(0);
          
          this->const_mass_matrix(0, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(0)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(0, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(0);
          
          this->const_mass_matrix(1, 0)(i, j) = vloc;

          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(1, 1)(i, j) = vloc;
          
          // curl curl matrix
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Curl_Phi2D(i, k)*this->Curl_Phi2D(j, k);
          
          this->const_curl_matrix(i, j) = vloc;
        }

    if (this->nb_dof_H > 0)
      {
        Ch_Hloc.Reallocate(this->nb_dof_H, this->nb_points_quadrature_inside);
        for (int i = 0; i < this->nb_points_quadrature_inside; i++)
          {
            this->ComputeValuesPhiHRef(this->PointsND(i), curl_phi);
            for (int j = 0; j < this->nb_dof_H; j++)
              Ch_Hloc(j, i) = curl_phi(j);
          }
      }
  }


  //! fills array EdgesDof and check that dofs on edges are correctly numbered
  /*!
    \param[in] Fb leaf finite element class
    Assuming that basis functions have been constructed, we find which function
    is interior, which function is associated with edge 0, 1, 2, etc
    If the dofs are not correctly numbered (edge 0, edge 1, edge 2, then interior)
    the program is aborted
   */
  void ElementReference<Dimension2, 2>::FindDofsOnEdge()
  {
    int N = this->GetNbDof();
    if (N <= 0)
      {
        cout << "no degree of freedom "<<endl;
        abort();
      }
    
    Real_wp threshold = 1e4*epsilon_machine;
    // for each edge we are computing phi \times n 
    Matrix<bool> tangential_proj_null(N, this->GetNbBoundaries());
    tangential_proj_null.Fill(true);
    VectR2 phi; Real_wp phi_n;
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      for (int i = 0; i < this->GetNbQuadBoundary(num_loc); i++)
        {
          int k = this->GetQuadNumber(num_loc, i);
          this->GetValuePhiOnQuadraturePoint(k, phi);
          
          for (int j = 0; j < N; j++)
            {
              phi_n = phi(j)(0)*this->NormaleLoc(num_loc)(1) - phi(j)(1)*this->NormaleLoc(num_loc)(0);
              if (abs(phi_n) > threshold)
                {
                  tangential_proj_null(j, num_loc) = false;
                }
            }
        }
    
    // finding dofs inside element
    Vector<bool> DofInside(N); DofInside.Fill(true);
    int nb = 0;
    for (int i = 0; i < this->GetNbDof(); i++)
      {
        for (int j = 0; j < this->GetNbBoundaries(); j++)
          if (!tangential_proj_null(i, j))
            DofInside(i) = false;
            
        if (DofInside(i))
          nb++;
      }
    
    for (int i = N-nb; i < N; i++)
      if (!DofInside(i))
        {
          cout << "Internal dofs must be numbered at the end " << endl;
          abort();
        }
    
    if (this->nb_dof_boundaries <= 0)
      this->nb_dof_boundaries = N - nb;
    else
      {
        if (this->nb_dof_boundaries != N - nb)
          {
            cout << "Number of dofs on the boundary is noticed as "
                 << this->nb_dof_boundaries << endl;
            cout << "But we found " << N-nb << "dofs on the boundary " << endl;
            abort();
          }
      }
    
    IVect nb_dof_edge(this->GetNbBoundaries());
    nb_dof_edge.Fill(0);
    int nb_max = N - nb;
    for (int num_loc = this->GetNbBoundaries()-1; num_loc >= 0; num_loc--)
      {
        DofInside.Fill(true);
        nb = 0;
        for (int i = 0; i < nb_max; i++)
          {
            for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
              if (num_loc2 != num_loc)
                if (!tangential_proj_null(i, num_loc2))
                  DofInside(i) = false;
            
            if (DofInside(i))
              nb++;
          }
        
        for (int i = nb_max-nb; i < nb_max; i++)
          if (!DofInside(i))
            {
              cout << "dof number for the edge " << num_loc
                   <<" is incorrect " << endl;
              abort();
            }       
        
        nb_dof_edge(num_loc) = nb;
        nb_max -= nb;
      }
        
    // filling EdgesDof
    nb_max = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      nb_max = max(nb_max, nb_dof_edge(j));
    
    EdgesDof.Reallocate(nb_max, this->GetNbBoundaries());
    nb_max = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      {
        nb = nb_dof_edge(j);
        for (int k = 0; k < nb; k++)
          EdgesDof(k, j) = nb_max + k;
        
        nb_max += nb;
      }
    
    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);
  }


  //! computes the sign of dofs after symmetry of edges
  /*!
    This function is only a function to check that restriction of functions
    on edges are the same between the edges, and that a symmetry
    of an edge may change the sign only.
    This function does not modify the finite element class, neither nmap
   */
  void ElementReference<Dimension2, 2>::
  FindHcurlSignEdge(NumberMap& nmap) const
  {
    VectR2 s;
    switch (this->GetNbVertices())
      {
      case 3 :
        {
          s.Reallocate(3);
          s(0).Init(0, 0);
          s(1).Init(1, 0);
          s(2).Init(0, 1);
        }
        break;
      case 4 :
        {
          s.Reallocate(4);
          s(0).Init(0, 0);
          s(1).Init(1, 0);
          s(2).Init(1, 1);
          s(3).Init(0, 1);
        }
        break;
      }
    
    int nb_quad = this->elt_geom.Points1D().GetM();
    int nb_dof_edge = nmap.GetNbDofEdge(this->order);
    if (nb_dof_edge == 0)
      {
        cout << "You must have at least one dof per edge " << endl;
        abort();
      }
    
    Array3D<Real_wp> ValPhiEdge(this->GetNbEdges(), nb_dof_edge, nb_quad);
    ValPhiEdge.Fill(0);
    
    // computing tangent for each edge, and values of basis functions
    VectR2 tangente(this->GetNbEdges()), phi;
    int type_elt = this->GetHybridType();
    Real_wp threshold = 1e4*epsilon_machine;
    for (int num_loc = 0; num_loc < this->GetNbEdges(); num_loc++)
      {
        int n1 = MeshNumbering<Dimension2>::FirstExtremityEdge(type_elt, num_loc);
        int n2 = MeshNumbering<Dimension2>::SecondExtremityEdge(type_elt, num_loc);
        tangente(num_loc) = s(n2) - s(n1);
	for (int k = 0; k < nb_quad; k++)
          {
            R2 ptA = (1.0-this->Points1D(k))*s(n1) + this->Points1D(k)*s(n2);
            this->ComputeValuesPhiRef(ptA, phi);
            for (int j = 0; j < nb_dof_edge; j++)
              {
                int num_dof = num_loc*nb_dof_edge + j;
                ValPhiEdge(num_loc, j, k) = DotProd(phi(num_dof), tangente(num_loc));
                for (int num_loc2 = 0; num_loc2 < this->GetNbEdges(); num_loc2++)
                  if (num_loc2 != num_loc)
                    {
                      int num_dof = num_loc2*nb_dof_edge + j;
                      Real_wp val = DotProd(phi(num_dof), tangente(num_loc));
                      
                      if (abs(val) > threshold)
                        {
                          cout << "Basis function "<<num_dof << " is not null on edge "
                               << num_loc2 << endl;
                          abort();
                        }
                    }
              }
          }
      }
    
    // we check tangential traces on each edge
    for (int num_loc = 0; num_loc < this->GetNbEdges(); num_loc++)
      {
        Matrix<Real_wp> A(nb_dof_edge, nb_dof_edge), B = A, M = A;
        A.Fill(0); B.Fill(0); M.Fill(0);
        for (int j = 0; j < nb_dof_edge; j++)
          {
            // we check that it is the same as first edge
            for (int k = 0; k < nb_quad; k++)
              if (abs(ValPhiEdge(num_loc, j, k) - ValPhiEdge(0, j, k)) > threshold)
                {
                  cout << "values of tangential trace of function " << j 
                       << " on edge " <<  num_loc 
                       << " does not with function of first edge " << endl;
                  
                  DISP(ValPhiEdge(num_loc, j, k));
                  DISP(ValPhiEdge(0, j, k));
                  abort();
                }
            
            for (int k = 0; k < nb_quad; k++)
              {
                Real_wp val2 = ValPhiEdge(num_loc, j, nb_quad-1-k);
                Real_wp val = ValPhiEdge(num_loc, j, k);
                for (int j2 = 0; j2 < nb_dof_edge; j2++)
                  {
                    Real_wp val1 = ValPhiEdge(num_loc, j2, k);
                    B(j, j2) += val2*this->Weights1D(k)*val1;
                    M(j, j2) += val*this->Weights1D(k)*val1;
                  }
              }
          }
        
        GetInverse(M);
        Mlt(B, M, A);
        //DISP(num_loc); DISP(A);
        
        // we check that matrix A is a permutation matrix with sign
        for (int j = 0; j < nb_dof_edge; j++)
          {
            int nb = 0;
            for (int jb = 0; jb < nb_dof_edge; jb++)
              if (abs(A(j, jb)) > threshold)
                nb++;
            
            if (nb != 1)
              {
                abort();
              }
	    
            //int j2 = -1; int signe = 0;
            for (int jb = 0; jb < nb_dof_edge; jb++)
              {
                Real_wp val = A(j, jb);
                if (abs(val - 1.0) <= threshold)
                  {
                    //j2 = jb;
                    //signe = 1;
                  }
                else if (abs(val + 1.0) <= threshold)
                  {
                    //j2 = jb;
                    //signe = -1;
                  }
                else if (abs(A(j, jb)) > threshold)
                  {
                    abort();
                  }
              }
            // DISP(num_loc); DISP(j); DISP(j2); DISP(signe);
          }
      }
  }
  
  
  //! Computes the effect of rotation on a quadrilateral face on degrees of freedom
  /*!
    \param[out] nmap numbering scheme of the mesh
    We compute which operators should be applied to degrees of freedom
    when a quadrilateral face is rotated. For quadrilateral faces,
    this operator consists of permutations and change of signs.
    nmap is modified to collect this information
   */
  void ElementReference<Dimension2, 2>::FindHcurlRotationQuad(NumberMap& nmap, int offset) const
  {
    if (offset == -1)
      offset = this->GetNbDofBoundaries();
    
    const QuadrangleReference<2>& quad = dynamic_cast<const QuadrangleReference<2>& >(*this);
    const Matrix<int>& NumQuad2D = quad.NumQuad2D;
    int r = NumQuad2D.GetM()-1;
    // loop on interior degrees of freedom
    int nb_dof_inside = this->GetNbDof() - offset;
    Matrix<int> FaceRot(8, nb_dof_inside);
    Matrix<bool> OppositeSigneDof(8, nb_dof_inside);
    for (int i = 0; i < nb_dof_inside; i++)
      for (int j = 0; j < 8; j++)
        {
          FaceRot(j, i) = i;
          OppositeSigneDof(j, i) = false;
        }
    
    Vector<int> rot_dof(8), rot_quad(8);
    Vector<bool> sign_dof(8);
    R2 phi, phi2, phi_j; Real_wp threshold = 1e4*epsilon_machine;
    TinyVector<bool, 8> dof_found_p, dof_found_m;
    rot_quad.Fill(-1);
    for (int i = offset; i < this->GetNbDof(); i++)
      {
        rot_dof.Fill(-1);
        sign_dof.Fill(true);
        for (int j = offset; j < this->GetNbDof(); j++)
          {
            dof_found_p.Fill(true);
            dof_found_m.Fill(true);
            for (int k1 = 0; k1 <= r; k1++)
              for (int k2 = 0; k2 <= r; k2++)
                {
                  rot_quad(0) = NumQuad2D(k1, k2);
                  rot_quad(3) = NumQuad2D(r-k2, k1);
                  rot_quad(2) = NumQuad2D(r-k1, r-k2);
                  rot_quad(1) = NumQuad2D(k2, r-k1);
                  rot_quad(4) = NumQuad2D(k2, k1);
                  rot_quad(5) = NumQuad2D(r-k1, k2);
                  rot_quad(6) = NumQuad2D(r-k2, r-k1);
                  rot_quad(7) = NumQuad2D(k1, r-k2);
                  phi = this->Value_PhiVec(i, rot_quad(0));
                  
                  for (int rot = 1; rot <= 7; rot++)
                    {
                      phi_j = this->Value_PhiVec(j, rot_quad(rot));
                      switch (rot)
                        {
                        case 3 : phi2.Init(phi_j(1), -phi_j(0)); break;
                        case 2 : phi2.Init(-phi_j(0), -phi_j(1)); break;
                        case 1 : phi2.Init(-phi_j(1), phi_j(0)); break;
                        case 4 : phi2.Init(phi_j(1), phi_j(0)); break;
                        case 5 : phi2.Init(-phi_j(0), phi_j(1)); break;
                        case 6 : phi2.Init(-phi_j(1), -phi_j(0)); break;
                        case 7 : phi2.Init(phi_j(0), -phi_j(1)); break;
                        }
                      
                      if (phi.Distance(phi2) > threshold)
                        dof_found_p(rot) = false;
                      
                      Mlt(-Real_wp(1), phi2);
                      if (phi.Distance(phi2) > threshold)
                        dof_found_m(rot) = false;
                    }
                }
            
            for (int rot = 1; rot <= 7; rot++)
              {
                if (dof_found_p(rot))
                  {
                    rot_dof(rot) = j;
                    sign_dof(rot) = true;
                  }
                
                if (dof_found_m(rot))
                  {
                    rot_dof(rot) = j;
                    sign_dof(rot) = false;
                  }
              }
          }
        
        for (int rot = 1; rot <= 7; rot++)
          {
            if (rot_dof(rot) == -1)
              {
                cout << "Case not handled for quadrilateral faces" << endl;
                abort();
              }
            else
              {
                int i2 = i - offset;
                int j = rot_dof(rot) - offset;
                FaceRot(rot, i2) = j;
                OppositeSigneDof(rot, i2) = !sign_dof(rot);
              }
            
          }
      }
    
    nmap.SetFacesDofRotationQuad(this->GetOrder(), FaceRot);
    nmap.SetSignDofRotationQuad(this->GetOrder(), OppositeSigneDof);
  }


  //! Computes the effect of rotation on a triangular face on degrees of freedom
  /*!
    \param[in] Fb leaf finite element class
    \param[out] nmap numbering scheme of the mesh
    We compute which operators should be applied to degrees of freedom
    when a triangular face is rotated.
   */
  void ElementReference<Dimension2, 2>::FindHcurlRotationTri(NumberMap& nmap, int offset) const
  {
    if (offset == -1)
      offset = this->GetNbDofBoundaries();
    
    // loop on interior degrees of freedom
    int nb_dof = this->GetNbDof();
    int nb_dof_inside = this->GetNbDof() - offset;
    Matrix<int> FaceRot(6, nb_dof_inside);
    Matrix<bool> OppositeSigneDof(6, nb_dof_inside);
    for (int i = 0; i < nb_dof_inside; i++)
      for (int j = 0; j < 6; j++)
        {
          FaceRot(j, i) = i;
          OppositeSigneDof(j, i) = false;
        }
	
    // rotation of quadrature points
    Vector<Matrix<Real_wp> > CoefCombination(6);
    const VectR2& Points2D = this->PointsND();
    const VectReal_wp& Weights2D = this->WeightsND();
    int nb_points = this->GetNbPointsQuadratureInside();
    VectR2 Points(nb_points); IVect permut(nb_points);
    Matrix<int> permut_rot(6, nb_points);
    Vector<Matrix2_2> invDF(6);
    invDF(0).SetIdentity();
    invDF(1)(0, 0) = 0.0; invDF(1)(0, 1) = 1.0; invDF(1)(1, 0) = -1.0; invDF(1)(1, 1) = -1.0;
    invDF(2)(0, 0) = -1.0; invDF(2)(0, 1) = -1.0; invDF(2)(1, 0) = 1.0; invDF(2)(1, 1) = 0.0;
    invDF(3)(0, 0) = 0.0; invDF(3)(0, 1) = 1.0; invDF(3)(1, 0) = 1.0; invDF(3)(1, 1) = 0.0;
    invDF(4)(0, 0) = -1.0; invDF(4)(0, 1) = -1.0; invDF(4)(1, 0) = 0.0; invDF(4)(1, 1) = 1.0;
    invDF(5)(0, 0) = 1.0; invDF(5)(0, 1) = 0.0; invDF(5)(1, 0) = -1.0; invDF(5)(1, 1) = -1.0;
    for (int rot = 0; rot < 6; rot++)
      {
	CoefCombination(rot).Reallocate(nb_dof_inside, nb_dof_inside);
	CoefCombination(rot).Fill(0);
	for (int i = 0; i < nb_points; i++)
	  {
	    Real_wp l1 = Points2D(i)(0), l2 = Points2D(i)(1);
	    Real_wp l0 = 1.0 - l1 - l2, x, y;
	    switch(rot)
	      {
	      case 0:
		x = l1; y = l2; break;
	      case 1 :
		x = l0; y = l1; break;
	      case 2 :
		x = l2; y = l0; break;
	      case 3 :
		x = l2; y = l1; break;
	      case 4 :
		x = l0; y = l2; break;
	      case 5 :
		x = l1; y = l0; break;
	      }
	    
	    Points(i).Init(x, y);
	  }
	
	permut.Fill();
	Sort(nb_points, Points, permut);
	for (int i = 0; i < nb_points; i++)
	  permut_rot(rot, i) = permut(i);
      }

    // computation of operators for each rotation
    Vector<int> rot_dof(6);
    R2 phi, phi2, tmp;
    Matrix<Real_wp> Mh(nb_dof_inside, nb_dof_inside);
    Mh.Zero();
    for (int i = offset; i < nb_dof; i++)
      {
        rot_dof.Fill(-1);
        for (int j = offset; j < nb_dof; j++)
          for (int k = 0; k < nb_points; k++)
            {
              int krot = permut_rot(0, k);
              phi = this->Value_PhiVec(i, krot);
              Mh(i-offset, j-offset) += Weights2D(k)*DotProd(phi, this->Value_PhiVec(j, krot));
              
              for (int rot = 1; rot <= 5; rot++)
                {
                  tmp = this->Value_PhiVec(j, permut_rot(rot, k));
		  MltTrans(invDF(rot), tmp, phi2);
                  CoefCombination(rot)(i-offset, j-offset) += Weights2D(k)*DotProd(phi, phi2);
                }
            }
      }
    
    GetInverse(Mh);

    Matrix<Real_wp> Ah(nb_dof_inside, nb_dof_inside), Bh;
    Ah.Zero();
    bool linear_combination = false;
    CoefCombination(0).SetIdentity();
    Vector<Matrix<Real_wp> > CoefCombin(CoefCombination);
    for (int rot = 1; rot < 6; rot++)
      {
	Bh = CoefCombin(rot);
	Mlt(Mh, Bh, Ah);
        
	int rot2 = rot;
	if (rot == 2)
	  rot2 = 1;
	else if (rot == 1)
	  rot2 = 2;
	
	Real_wp threshold = 1000*epsilon_machine;
	// pour tous les ddls de la face
	Vector<bool> DofUsed(nb_dof_inside); DofUsed.Fill(false);            
	for (int i = offset; i < nb_dof; i++)
	  {
	    int nb_coef = 0;
	    for (int j = offset; j < nb_dof; j++)
	      if (abs(Ah(j-offset, i-offset)) > threshold)
		nb_coef++;
	    
	    if (nb_coef == 0)
	      {
		cout << "Cas impossible " << endl;
		abort();
	      }
	    else if (nb_coef == 1)
	      {
		// on est tombe sur un autre ddl, avec eventuellement un signe
		for (int j = offset; j < nb_dof; j++)
		  if (abs(Ah(j-offset, i-offset)) > threshold)
		    {
		      FaceRot(rot2, i-offset) = j-offset;
		      if ( abs(Ah(j-offset, i-offset) -Real_wp(1)) <= threshold)
			OppositeSigneDof(rot2, i-offset) = false;
		      else if ( abs(Ah(j-offset, i-offset) + Real_wp(1)) <= threshold)
			OppositeSigneDof(rot2, i-offset) = true;
		      else
			{
			  // this dof can't be handled with a sign
			  linear_combination = true;
			}
                                      
		      DofUsed(j-offset) = true;
		    }                       
	      }
	    else
	      {
		// no direct relation between dofs
		linear_combination = true;
	      }
	  }

	//DISP(rot2); DISP(Ah);
	CoefCombination(rot2) = Ah;
      }

    int rf = this->order;
    //DISP(linear_combination); DISP(FaceRot); DISP(OppositeSigneDof);
    if (linear_combination)
      nmap.SetFacesDofRotationTri(rf, CoefCombination);
    else
      {
	nmap.SetFacesDofRotationTri(rf, FaceRot);
	nmap.SetSignDofRotationTri(rf, OppositeSigneDof);
      }
  }
  
  
  //! computation of values of u on nodal points of an edge
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points of edge num_loc
    \param[in] num_loc edge number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 2>
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(2*this->EdgesNodal.GetM());
    typename Vector2::value_type vx, vy;
    for (int j = 0; j < this->EdgesNodal.GetM(); j++)
      {
        int node = this->EdgesNodal(j, num_loc);
        vx = 0; vy = 0;
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            vx += Un(k)*this->ValuePhi_Nodal(k, 2*node);
            vy += Un(k)*this->ValuePhi_Nodal(k, 2*node+1);
          }
        
        Unode(2*j) = vx;
        Unode(2*j+1) = vy;
      }
  }


  //! computation of curl of u on nodal points of an edge
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode curl of u on nodal points of edge num_loc
    \param[in] num_loc edge number
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 2>::
  ComputeCurlBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(this->EdgesNodal.GetM());
    for (int j = 0; j < this->EdgesNodal.GetM(); j++)
      {
        Unode(j) = 0;
        int node = this->EdgesNodal(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Unode(j) += Un(k)*this->CurlPhi_Nodal(k, node);
      }
  }

  
  //! Integration against curl of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f curl(\varphi_i) dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = curl(\varphi_i)(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 2>::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
        val = Real_wp(0);
        for (int j = 0; j < nb_points_quadrature_inside; j++)
          val += Curl_Phi2D(i, j)*Vh(j);
	
        Uh(i) = val;
      }
  }


  //! computation of curl of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh curl of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j curl phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = curl phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 2>::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type curl;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        curl = 0;
        for (int j = 0; j < this->nb_dof_loc; j++)
          curl += Uh(j)*Curl_Phi2D(j, i);
        
        Vh(i) = curl;
      }
  }


  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 2>::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val = 0;
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Value_PhiVec(i, j)(0)*Vh(2*j)
            + Value_PhiVec(i, j)(1)*Vh(2*j+1);        
	
	Uh(i) = val;
      }
  }


  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 2>
  ::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type vx, vy;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	vx = 0; vy = 0;
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
            vx += Value_PhiVec(j, i)(0)*Uh(j);
            vy += Value_PhiVec(j, i)(1)*Uh(j);
          }
	
	Vh(2*i) = vx; Vh(2*i+1) = vy;
      }
  }


  //! integration against basis functions on a boundary
  /*!
    \param[in] alpha coefficient
    \param[in] num_loc boundary number
    \param[in] Uh vector containing \omega_k f(\xi_k)
    \param[out] Vh Vh_i = Vh_i + alpha \int_{\partial K} f \psi_i ds    
    \param[in] r order of quadrature rules
    (\omega_k, \xi_k) is a quadrature rule of the boundary num_loc
   */
  template<class T>
  void ElementReference<Dimension2, 2>::
  ApplyShGen(const T& alpha, int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const
  {
    TinyVector<T, 2> f;
    for (int i = offset_faceSh(n); i < offset_faceSh(n+1); i++)
      {
        int j = i + nb_points_quadrature_inside;
        int p = i - offset_faceSh(n);
        CopyVector(Uh, p, f);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Vh(k) += alpha*DotProd(Value_PhiVec(k, j), f);
      }
  }


  //! computation of u on quadrature points of a boundary
  /*!
    \param[in] n boundary number
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points of the boundary
    \param[in] r order of quadrature rule of the boundary
    This operation is equivalent to a matrix vector product
    Vh = Sh* Uh
    where (Sh)_{i,j} = \varphi_i(\xi_j)  where \xi_j are quadrature points
   */
  template<class T>
  void ElementReference<Dimension2, 2>::
  ApplyShTransposeGen(int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const
  {
    abort();
  }

  
  //! Integration against curl of basis functions on boundaries
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f curl(\varphi_i) dx
    \param[in] n boundary number
   */  
  template<class T>
  void ElementReference<Dimension2, 2>
  ::ApplyNablaShGen(const T& alpha, int n, const Vector<T>& Uh, Vector<T>& Vh, int r) const
  {
    for (int i = offset_faceSh(n); i < offset_faceSh(n+1); i++)
      {
        int j = i + nb_points_quadrature_inside;
        int p = i - offset_faceSh(n);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Vh(k) += alpha*Curl_Phi2D(k, j)*Uh(p);
      }
  }

  
  //! Not implemented
  template<class T>
  void ElementReference<Dimension2, 2>
  ::ApplyNablaShTransposeGen(int num_loc, const Vector<T>& Uh, Vector<T>& Vh, int r) const
  {
    abort();
  }


  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension2, 2>::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 2, 2>& mass, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  val(j) = const_mass_matrix(0, 0)(i, j)*mass(0, 0)
	    + const_mass_matrix(0, 1)(i, j)*mass(0, 1)
	    + const_mass_matrix(1, 0)(i, j)*mass(1, 0)
	    + const_mass_matrix(1, 1)(i, j)*mass(1, 1);

	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }          
  }


  //! we add constant stiffness matrix C \int_K curl phi_j curl phi_i
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + C S
    where S is a stiffness matrix equal to :
    (S)_{i, j} = \int_K curl(\phi_j) curl(\phi_i)  dx
   */
  template<class T>
  void ElementReference<Dimension2, 2>::
  AddConstantStiffnessMatrixGen(int m, int n, const T& C, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  val(j) = const_curl_matrix(i, j)*C;
	
	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }
  }

    
  //! adds \int C curl(varphi_j) curl(varphi_i)
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] C coefficient C on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int C curl(varphi_j) curl(varphi_i)
   */
  template<class T>
  void ElementReference<Dimension2, 2>
  ::AddVariableStiffnessMatrixGen(int off_row, int off_col,
				  const Vector<T>& C, VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    SetComplexZero(vloc);
	    for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	      vloc += this->Curl_Phi2D(i, k)*this->Curl_Phi2D(j, k)*C(k);
          
	    val(j) = vloc;
	  }
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }
  
  
  //! adds \int A varphi_j \varphi_i
  /*!
    \param[in] off_row offset for row numbers when accessing to mat
    \param[in] off_col offset for column numbers when accessing to mat
    \param[in] A tensor A on each quadrature point
    \param[inout] mat modified matrix
    mat(off_row:, off_col:) = mat(off_row:, off_col:) + S
    where S_ij = \int A varphi_j  varphi_i
   */
  template<class T, class Prop>
  void ElementReference<Dimension2, 2>
  ::AddVariableMassMatrixGen(int off_row, int off_col,
			     const Vector<TinyMatrix<T, Prop, 2, 2> >& C,
			     VirtualMatrix<T>& mat) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = off_col + i;

    T vloc; TinyVector<T, 2> vec_u;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    SetComplexZero(vloc);
	    for (int k = 0; k < this->nb_points_quadrature_inside; k++)
	      {
		Mlt(C(k), this->Value_PhiVec(j, k), vec_u);
		vloc += DotProd(vec_u, this->Value_PhiVec(i, k));
	      }
	    
	    val(j) = vloc;
	  }
	
	mat.AddInteractionRow(off_row+i, this->nb_dof_loc, num, val);
      }
  }
  

  //! integrals against curl of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension2, 2>
  ::ComputeIntegralCurlRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyRh(feval, res);
  }
  
  
  //! surface integrals against curl of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension2, 2>
  ::ComputeIntegralSurfaceCurlRef(const Vector1 & feval, Vector2& res, int n) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    this->ApplyNablaSh(one, n, feval, res);
  }

  
  //! projection from values on dof points to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension2, 2>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib = feval;
    SolveQR(MatProjectionDof, tauProjectionDof, contrib);
  }

  
  //! projection from values on dof points to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension2, 2>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    VectReal_wp xreal(feval.GetM()), ximag(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      {
        xreal(i) = real(feval(i));
        ximag(i) = imag(feval(i));
      }
    
    SolveQR(MatProjectionDof, tauProjectionDof, xreal);
    SolveQR(MatProjectionDof, tauProjectionDof, ximag);
    
    contrib.Reallocate(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      contrib(i) = Complex_wp(xreal(i), ximag(i));
  }


  void ElementReference<Dimension2, 2>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    VectReal_wp feval_tangential(feval.GetM()/2);
    R2 vec_u;
    for (int i = 0; i < feval_tangential.GetM(); i++)
      {
	vec_u.Init(feval(2*i), feval(2*i+1));
	feval_tangential(i) = TimesProd(this->NormaleLoc(num_loc), vec_u);
      }
    
    this->element_surface->ComputeProjectionDofRef(feval_tangential, contrib);
  }


  void ElementReference<Dimension2, 2>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    VectComplex_wp feval_tangential(feval.GetM()/2);
    R2_Complex_wp vec_u;
    for (int i = 0; i < feval_tangential.GetM(); i++)
      {
	vec_u.Init(feval(2*i), feval(2*i+1));
	feval_tangential(i) = -TimesProd(vec_u, this->NormaleLoc(num_loc));
      }
    
    this->element_surface->ComputeProjectionDofRef(feval_tangential, contrib);
  }  

#endif


  /*******************
   * H(div) elements *
   *******************/
  
  
#ifdef MONTJOIE_WITH_NODAL_HDIV
  //! Computes the effect of rotation on a quadrilateral face on degrees of freedom
  /*!
    \param[in] rf order of approximation for the face
    \param[in] NumQuad2D numbering of quadrature points on the face
    \param[in] ValuePhi values of L2 basis functions on quadrature points
    \param[out] nmap numbering scheme of the mesh
    We compute which operators should be applied to degrees of freedom
    when a quadrilateral face is rotated. For quadrilateral faces,
    this operator consists of permutations and change of signs.
    nmap is modified to collect this information
   */
  void ElementReference<Dimension2, 3>::FindHdivRotationQuad(int rf, const Matrix<int>& NumQuad2D, 
							     const Matrix<Real_wp>& ValuePhi,
							     NumberMap& nmap)
  {
    int nb_dof = ValuePhi.GetM();
    int r = NumQuad2D.GetM()-1;
    
    // loop on interior degrees of freedom
    Matrix<int> FaceRot(8, nb_dof);
    Matrix<bool> OppositeSigneDof(8, nb_dof);
    for (int i = 0; i < nb_dof; i++)
      for (int j = 0; j < 8; j++)
	{
	  FaceRot(j, i) = i;
	  OppositeSigneDof(j, i) = false;
	}
    
    Vector<bool> sign_dof(8);
    Vector<int> rot_dof(8), rot_quad(8);
    Real_wp phi, phi2; Real_wp threshold = 1e4*epsilon_machine;
    TinyVector<bool, 8> dof_found_p, dof_found_m;
    rot_quad.Fill(-1);
    for (int i = 0; i < nb_dof; i++)
      {
        rot_dof.Fill(-1);
        sign_dof.Fill(true);
        for (int j = 0; j < nb_dof; j++)
          {
            dof_found_p.Fill(true);
            dof_found_m.Fill(true);
            for (int k1 = 0; k1 <= r; k1++)
              for (int k2 = 0; k2 <= r; k2++)
                {
                  rot_quad(0) = NumQuad2D(k1, k2);
                  rot_quad(3) = NumQuad2D(r-k2, k1);
                  rot_quad(2) = NumQuad2D(r-k1, r-k2);
                  rot_quad(1) = NumQuad2D(k2, r-k1);
                  rot_quad(4) = NumQuad2D(k2, k1);
                  rot_quad(5) = NumQuad2D(r-k1, k2);
                  rot_quad(6) = NumQuad2D(r-k2, r-k1);
                  rot_quad(7) = NumQuad2D(k1, r-k2);
                  phi = ValuePhi(i, rot_quad(0));
                  
                  for (int rot = 1; rot <= 7; rot++)
                    {
                      phi2 = ValuePhi(j, rot_quad(rot));
                      if (abs(phi - phi2) > threshold)
                        dof_found_p(rot) = false;
                      
		      phi2 *= -Real_wp(1);
                      if (abs(phi - phi2) > threshold)
                        dof_found_m(rot) = false;
                    }
                }
            
            for (int rot = 1; rot <= 7; rot++)
              {
                if (dof_found_p(rot))
                  {
                    rot_dof(rot) = j;
                    sign_dof(rot) = true;
                  }
                
                if (dof_found_m(rot))
                  {
                    rot_dof(rot) = j;
                    sign_dof(rot) = false;
                  }
              }
          }
        
        for (int rot = 1; rot <= 7; rot++)
          {
            if (rot_dof(rot) == -1)
              {
                cout << "Case not handled for quadrilateral faces" << endl;
                abort();
              }
            else
              {
                int i2 = i;
                int j = rot_dof(rot);
                FaceRot(rot, i2) = j;
                OppositeSigneDof(rot, i2) = !sign_dof(rot);
		if (rot >= 4)
		  OppositeSigneDof(rot, i2) = sign_dof(rot);
              }
            
          }
      }
    
    nmap.SetFacesDofRotationQuad(rf, FaceRot);
    nmap.SetSignDofRotationQuad(rf, OppositeSigneDof);
  }


  //! Computes the effect of rotation on a triangular face on degrees of freedom
  /*!
    \param[in] rf order of approximation for the face
    \param[in] Points2D quadrature points on the face
    \param[in] Weights2D quadrature weights on the face
    \param[in] ValuePhi values of L2 basis functions on quadrature points
    \param[out] nmap numbering scheme of the mesh
    We compute which operators should be applied to degrees of freedom
    when a quadrilateral face is rotated.
   */
  void ElementReference<Dimension2, 3>::FindHdivRotationTri(int rf, const VectR2& Points2D,
							    const VectReal_wp& Weights2D,
							    const Matrix<Real_wp>& ValuePhi,
							    NumberMap& nmap)
  {
    int nb_dof = ValuePhi.GetM();
    
    // loop on interior degrees of freedom
    Matrix<int> FaceRot(6, nb_dof);
    Matrix<bool> OppositeSigneDof(6, nb_dof);
    for (int i = 0; i < nb_dof; i++)
      for (int j = 0; j < 6; j++)
	{
	  FaceRot(j, i) = i;
	  OppositeSigneDof(j, i) = false;
	}
    
    int nb_points = Points2D.GetM();
    
    Vector<Matrix<Real_wp> > CoefCombination(6);
    Matrix<Real_wp> Mh(nb_dof, nb_dof), Bh(nb_dof, nb_dof);
    VectR2 Points(nb_points); IVect permut(nb_points);
    Matrix<int> permut_rot(6, nb_points);
    Mh.Fill(0);
    for (int rot = 0; rot < 6; rot++)
      {
	CoefCombination(rot).Reallocate(nb_dof, nb_dof);
	CoefCombination(rot).Fill(0);
	for (int i = 0; i < nb_points; i++)
	  {
	    Real_wp l1 = Points2D(i)(0), l2 = Points2D(i)(1);
	    Real_wp l0 = 1.0 - l1 - l2, x, y;
	    switch(rot)
	      {
	      case 0:
		x = l1; y = l2; break;
	      case 1 :
		x = l0; y = l1; break;
	      case 2 :
		x = l2; y = l0; break;
	      case 3 :
		x = l2; y = l1; break;
	      case 4 :
		x = l0; y = l2; break;
	      case 5 :
		x = l1; y = l0; break;
	      }
	    
	    Points(i).Init(x, y);
	  }
	
	permut.Fill();
	Sort(nb_points, Points, permut);
	for (int i = 0; i < nb_points; i++)
	  permut_rot(rot, i) = permut(i);
      }
    
    Vector<bool> sign_dof(6);
    Vector<int> rot_dof(6), rot_quad(6);
    Real_wp phi, phi2;
    TinyVector<bool, 6> dof_found_p, dof_found_m;
    rot_quad.Fill(-1);    
    for (int i = 0; i < nb_dof; i++)
      {
        rot_dof.Fill(-1);
        sign_dof.Fill(true);
        for (int j = 0; j < nb_dof; j++)
          {
            dof_found_p.Fill(true);
            dof_found_m.Fill(true);
            for (int k = 0; k < nb_points; k++)
	      {
		rot_quad(0) = permut_rot(0, k);
		phi = ValuePhi(i, rot_quad(0));
		Mh(i, j) += Weights2D(k)*phi*ValuePhi(j, rot_quad(0));
		
		for (int rot = 1; rot <= 5; rot++)
		  {
		    phi2 = ValuePhi(j, permut_rot(rot, k));
		    CoefCombination(rot)(i, j) += Weights2D(k)*phi*phi2;
		  }
	      }
          }
      }
    
    GetInverse(Mh);
    
    Matrix<Real_wp> Ah(nb_dof, nb_dof);
    Ah.Fill(0);
    bool linear_combination = false;
    CoefCombination(0).SetIdentity();
    Vector<Matrix<Real_wp> > CoefCombin(CoefCombination);
    for (int rot = 1; rot < 6; rot++)
      {
	Bh = CoefCombin(rot);
	Mlt(Mh, Bh, Ah);
	if (rot >= 3)
	  Mlt(-Real_wp(1), Ah);
	
	int rot2 = rot;
	if (rot == 2)
	  rot2 = 1;
	else if (rot == 1)
	  rot2 = 2;
	
	//DISP(rot); DISP(Ah);
	Real_wp threshold = 1000*epsilon_machine;
	// pour tous les ddls de la face
	Vector<bool> DofUsed(nb_dof); DofUsed.Fill(false);            
	for (int i = 0; i < nb_dof; i++)
	  {
	    int nb_coef = 0;
	    for (int j = 0; j < nb_dof; j++)
	      if (abs(Ah(j, i)) > threshold)
		nb_coef++;
	    
	    if (nb_coef == 0)
	      {
		cout << "Cas impossible " << endl;
		abort();
	      }
	    else if (nb_coef == 1)
	      {
		// on est tombe sur un autre ddl, avec eventuellement un signe
		for (int j = 0; j < nb_dof; j++)
		  if (abs(Ah(j, i)) > threshold)
		    {
		      FaceRot(rot2, i) = j;
		      if ( abs(Ah(j, i) -Real_wp(1)) <= threshold)
			OppositeSigneDof(rot2, i) = false;
		      else if ( abs(Ah(j, i) + Real_wp(1)) <= threshold)
			OppositeSigneDof(rot2, i) = true;
		      else
			{
			  // this dof can't be handled with a sign
			  linear_combination = true;
			}
                                      
		      DofUsed(j) = true;
		    }                       
	      }
	    else
	      {
		// no direct relation between dofs
		linear_combination = true;
	      }
	  }
	
	CoefCombination(rot2) = Ah;
      }
    
    if (linear_combination)
      nmap.SetFacesDofRotationTri(rf, CoefCombination);
    else
      {
	nmap.SetFacesDofRotationTri(rf, FaceRot);
	nmap.SetSignDofRotationTri(rf, OppositeSigneDof);
      }
  }

  
  //! computation of mass matrix, stiffness matrix, etc (for facet elements)  
  /*!
    \param[in] Fb leaf finite element class
   */
  void ElementReference<Dimension2, 3>::ConstructHdivElementaryMatrix()
  {
    int nb_points_quadrature = this->GetNbPointsQuadrature();
    Vector<R2> phi;
    Value_PhiVec.Reallocate(this->nb_dof_loc, nb_points_quadrature);        
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeValuesPhiRef(this->PointsND(i), phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  Value_PhiVec(j, i) = phi(j);
      }

    // mass matrix on reference element
    this->mass_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    Real_wp vloc;
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = i; j < this->nb_dof_loc; j++)
	{
	  vloc = 0.0;
	  for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*DotProd(this->Value_PhiVec(i, k), this->Value_PhiVec(j, k));
          
	  this->mass_matrix(i,j) = vloc;
	}
    
    //VectReal_wp lambda; 
    //GetEigenvalues(this->mass_matrix, lambda);
    //DISP(lambda);
	  
    this->mass_matrix_chol = this->mass_matrix;
    GetCholesky(this->mass_matrix_chol);
    
    // projection on dofs
    this->MatProjectionDof.Reallocate(2*this->GetNbPointsDofInside(), this->nb_dof_loc);
    for (int i = 0; i < this->GetNbPointsDofInside(); i++)
      {
        this->ComputeValuesPhiRef(this->PointsDofND(i), phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    this->MatProjectionDof(2*i, j) = phi(j)(0);
	    this->MatProjectionDof(2*i+1, j) = phi(j)(1);
	  }
      } 
    
    GetQR(this->MatProjectionDof, this->tauProjectionDof);
    
    ValuePhi_Nodal.Reallocate(this->nb_dof_loc, 2*this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
	this->ComputeValuesPhiRef(this->PointsNodalND(i), phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
	    ValuePhi_Nodal(j, 2*i) = phi(j)(0);
	    ValuePhi_Nodal(j, 2*i+1) = phi(j)(1);
	  }
      }
    
    Vector<Real_wp> div_phi;
    DivPhi_Nodal.Reallocate(this->nb_dof_loc, this->GetNbPointsNodalElt());
    for (int i = 0; i < this->GetNbPointsNodalElt(); i++)
      {
	this->ComputeDivPhiRef(this->PointsNodalND(i), div_phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  DivPhi_Nodal(j, i) = div_phi(j);
      }
    
    Div_Phi.Reallocate(this->nb_dof_loc, nb_points_quadrature);
    for (int i = 0; i < nb_points_quadrature; i++)
      {
	this->ComputeDivPhiRef(this->PointsND(i), div_phi);
	for (int j = 0; j < this->nb_dof_loc; j++)
	  Div_Phi(j, i) = div_phi(j);
      }
    
    // mass and stiffness matrix for constant coefficients
    // \int_K C phi_i phi_j   and   \int_K D div(phi_i)  div(phi_j)
    // where C is a 3x3 symmetric matrix
    for (int p = 0; p < 2; p++)
      for (int q = 0; q < 2; q++)
	this->const_mass_matrix(p, q).Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    
    this->const_div_matrix.Reallocate(this->nb_dof_loc, this->nb_dof_loc);
    
    for (int i = 0; i < this->nb_dof_loc; i++)
      for (int j = 0; j < this->nb_dof_loc; j++)
        {
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i, k)(0)*this->Value_PhiVec(j, k)(0);
          
          this->const_mass_matrix(0, 0)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(0)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(0, 1)(i, j) = vloc;
          
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(0);
          
          this->const_mass_matrix(1, 0)(i, j) = vloc;

	  vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Value_PhiVec(i,k)(1)*this->Value_PhiVec(j,k)(1);
          
          this->const_mass_matrix(1, 1)(i, j) = vloc;
	  
          // div div matrix
          vloc = 0;
          for (int k = 0; k < this->nb_points_quadrature_inside; k++)
            vloc += this->WeightsND(k)*this->Div_Phi(i, k)*this->Div_Phi(j, k);
          
          this->const_div_matrix(i, j) = vloc;
        }
  }


  //! fills array EdgesDof and check that dofs on edges are correctly numbered
  /*!
    \param[in] Fb leaf finite element class
    Assuming that basis functions have been constructed, we find which function
    is interior, which function is associated with edge 0, 1, 2, etc
    If the dofs are not correctly numbered (edge 0, edge 1, edge 2, then interior)
    the program is aborted
   */
  void ElementReference<Dimension2, 3>::FindDofsOnEdge()
  {
    int N = this->GetNbDof();
    if (N <= 0)
      {
        cout << "no degree of freedom "<<endl;
        abort();
      }
    
    Real_wp threshold = 1e4*epsilon_machine;
    // for each edge we are computing phi \cdot n 
    Matrix<bool> normal_proj_null(N, this->GetNbBoundaries());
    normal_proj_null.Fill(true);
    VectR2 phi; Real_wp phi_n;
    for (int num_loc = 0; num_loc < this->GetNbBoundaries(); num_loc++)
      for (int i = 0; i < this->GetNbQuadBoundary(num_loc); i++)
        {
          int k = this->GetQuadNumber(num_loc, i);
          this->GetValuePhiOnQuadraturePoint(k, phi);
          
          for (int j = 0; j < N; j++)
            {
              phi_n = phi(j)(0)*this->NormaleLoc(num_loc)(0) + phi(j)(1)*this->NormaleLoc(num_loc)(1);
              if (abs(phi_n) > threshold)
                {
                  normal_proj_null(j, num_loc) = false;
                }
            }
        }
    
    // finding dofs inside element
    Vector<bool> DofInside(N); DofInside.Fill(true);
    int nb = 0;
    for (int i = 0; i < this->GetNbDof(); i++)
      {
        for (int j = 0; j < this->GetNbBoundaries(); j++)
          if (!normal_proj_null(i, j))
            DofInside(i) = false;
            
        if (DofInside(i))
          nb++;
      }
    
    for (int i = N-nb; i < N; i++)
      if (!DofInside(i))
        {
          cout << "Internal dofs must be numbered at the end " << endl;
          abort();
        }
    
    if (this->nb_dof_boundaries <= 0)
      this->nb_dof_boundaries = N - nb;
    else
      {
        if (this->nb_dof_boundaries != N - nb)
          {
            cout << "Number of dofs on the boundary is noticed as "
                 << this->nb_dof_boundaries << endl;
            cout << "But we found " << N-nb << "dofs on the boundary " << endl;
            abort();
          }
      }
    
    IVect nb_dof_edge(this->GetNbBoundaries());
    nb_dof_edge.Fill(0);
    int nb_max = N - nb;
    for (int num_loc = this->GetNbBoundaries()-1; num_loc >= 0; num_loc--)
      {
        DofInside.Fill(true);
        nb = 0;
        for (int i = 0; i < nb_max; i++)
          {
            for (int num_loc2 = 0; num_loc2 < this->GetNbBoundaries(); num_loc2++)
              if (num_loc2 != num_loc)
                if (!normal_proj_null(i, num_loc2))
                  DofInside(i) = false;
            
            if (DofInside(i))
              nb++;
          }
        
        for (int i = nb_max-nb; i < nb_max; i++)
          if (!DofInside(i))
            {
              cout << "dof number for the edge " << num_loc
                   <<" is incorrect " << endl;
              abort();
            }       
        
        nb_dof_edge(num_loc) = nb;
        nb_max -= nb;
      }
        
    // filling EdgesDof
    nb_max = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      nb_max = max(nb_max, nb_dof_edge(j));
    
    EdgesDof.Reallocate(nb_max, this->GetNbBoundaries());
    nb_max = 0;
    for (int j = 0; j < this->GetNbBoundaries(); j++)
      {
        nb = nb_dof_edge(j);
        for (int k = 0; k < nb; k++)
          EdgesDof(k, j) = nb_max + k;
        
        nb_max += nb;
      }

    FillPositionDofBoundaries(EdgesDof, this->power_two_face, this->PosDofOnFace);
  }


  //! computation of values of u on nodal points of an edge
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode values of u on nodal points of edge num_loc
    \param[in] num_loc edge number
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 3>
  ::ComputeValueBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(2*this->EdgesNodal.GetM());
    typename Vector2::value_type vx, vy;
    for (int j = 0; j < this->EdgesNodal.GetM(); j++)
      {
        int node = this->EdgesNodal(j, num_loc);
        vx = 0; vy = 0;
        for (int k = 0; k < this->nb_dof_loc; k++)
          {
            vx += Un(k)*this->ValuePhi_Nodal(k, 2*node);
            vy += Un(k)*this->ValuePhi_Nodal(k, 2*node+1);
          }
        
        Unode(2*j) = vx;
        Unode(2*j+1) = vy;
      }
  }


  //! computation of divergence of u on nodal points of an edge
  /*!
    \param[in] Un components of u on degrees of freedom
    \param[out] Unode divergence of u on nodal points of edge num_loc
    \param[in] num_loc edge number
  */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 3>::
  ComputeDivBoundaryGen(const Vector1& Un, Vector2 & Unode, int num_loc) const
  {
    Unode.Reallocate(this->EdgesNodal.GetM());
    for (int j = 0; j < this->EdgesNodal.GetM(); j++)
      {
        Unode(j) = 0;
        int node = this->EdgesNodal(j, num_loc);
        for (int k = 0; k < this->nb_dof_loc; k++)
          Unode(j) += Un(k)*this->DivPhi_Nodal(k, node);
      }
  }


  //! Integration against divergence of basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f div(\varphi_i) dx
    This operation is equivalent to a matrix vector product
    Uh = Rh Vh
    where (Rh)_{i,j} = div(\varphi_i)(\xi_j)
   */  
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 3>::ApplyRhGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
        val = Real_wp(0);
        for (int j = 0; j < nb_points_quadrature_inside; j++)
          val += Div_Phi(i, j)*Vh(j);
	
        Uh(i) = val;
      }
  }


  //! computation of divergence of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh divergence of u on quadrature points
    This operation can be written as
    v_i \, = \sum_j div phi_j(xi_i) u_j
    that is to say Vh = Rh* Uh
    where (Rh)_{i, j} = div phi_i(xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 3>::ApplyRhTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type curl;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
        curl = 0;
        for (int j = 0; j < this->nb_dof_loc; j++)
          curl += Uh(j)*Div_Phi(j, i);
        
        Vh(i) = curl;
      }
  }


  //! Integration against basis functions
  /*!
    \param[in] Vh vector containing values \omega_k f(\xi_k)
               where omega_k is the weight of integration
               and \xi_k the point of integration
    \param[out] Uh Uh_i = \int_K f \varphi_i dx
    This operation is equivalent to a matrix vector product
    Uh = Ch Vh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 3>::ApplyChGen(const Vector1& Vh, Vector2& Uh) const
  {
    typename Vector2::value_type val;
    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	val = 0;
	for (int j = 0; j < nb_points_quadrature_inside; j++)
	  val += Value_PhiVec(i, j)(0)*Vh(2*j)
            + Value_PhiVec(i, j)(1)*Vh(2*j+1);        
	
	Uh(i) = val;
      }
  }


  //! computation of u on quadrature points
  /*!
    \param[in] Uh components of u on degrees of freedom
    \param[out] Vh values of u on quadrature points
    This operation is equivalent to a matrix vector product
    Vh = Ch* Uh
    where (Ch)_{i,j} = \varphi_i(\xi_j)
   */
  template<class Vector1, class Vector2>
  void ElementReference<Dimension2, 3>
  ::ApplyChTransposeGen(const Vector1& Uh, Vector2& Vh) const
  {
    typename Vector2::value_type vx, vy;
    for (int i = 0; i < nb_points_quadrature_inside; i++)
      {
	vx = 0; vy = 0;
	for (int j = 0; j < this->nb_dof_loc; j++)
	  {
            vx += Value_PhiVec(j, i)(0)*Uh(j);
            vy += Value_PhiVec(j, i)(1)*Uh(j);
          }
	
	Vh(2*i) = vx; Vh(2*i+1) = vy;
      }
  }


  //! Integration against basis functions on boundaries
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f \varphi_i dx
    \param[in] n boundary number
   */  
  template<class T>
  void ElementReference<Dimension2, 3>::
  ApplyShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    TinyVector<T, 2> f;
    for (int i = offset_faceSh(n); i < offset_faceSh(n+1); i++)
      {
        int p = i - offset_faceSh(n);
	int j = this->GetQuadNumber(n, p);
        CopyVector(feval, p, f);
        for (int k = 0; k < this->nb_dof_loc; k++)
          res(k) += alpha*DotProd(Value_PhiVec(k, j), f);
      }
  }

  
  //! Integration against divergence of basis functions on boundaries
  /*!
    \param[in] feval vector containing values \omega_k f(\xi_k)
                    where omega_k is the weight of integration
                    and \xi_k the point of integration
    \param[out] res res_i = \int_{\partial K} f div(\varphi_i) dx
    \param[in] n boundary number
   */  
  template<class T>
  void ElementReference<Dimension2, 3>::
  ApplyNablaShGen(const T& alpha, int n, const Vector<T>& feval, Vector<T>& res, int r) const
  {
    for (int i = offset_faceSh(n); i < offset_faceSh(n+1); i++)
      {
        int p = i - offset_faceSh(n);
	int j = this->GetQuadNumber(n, p);
        for (int k = 0; k < this->nb_dof_loc; k++)
          res(k) += alpha*Div_Phi(k, j)*feval(p);
      }
  }


  //! we add constant mass matrix
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] mass coefficient
    \param[inout] A matrix modified
    A(m:, n:) = A(m:, n:) + mass M
    where M is the mass matrix, M_ij = \int \varphi_j \varphi_i dx
   */
  template<class T, class Prop>
  void ElementReference<Dimension2, 3>::
  AddConstantMassMatrixGen(int m, int n, const TinyMatrix<T, Prop, 2, 2>& mass, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {	
	for (int j = 0; j < this->nb_dof_loc; j++)
	  val(j) = const_mass_matrix(0, 0)(i, j)*mass(0, 0)
	    + const_mass_matrix(0, 1)(i, j)*mass(0, 1)
	    + const_mass_matrix(1, 0)(i, j)*mass(1, 0)
	    + const_mass_matrix(1, 1)(i, j)*mass(1, 1);
	
	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }          
  }


  //! we add constant stiffness matrix C \int_K div phi_j div phi_i
  /*!
    \param[in] m offset for row numbers when accessing to A
    \param[in] n offset for column numbers when accessing to A
    \param[in] C coefficients
    \param[out] A modified matrix
    A(m:, n:) = A(m:, n:) + C S
    where S is a stiffness matrix equal to :
    (S)_{i, j} = \int_K div(\phi_j) div(\phi_i)  dx
   */
  template<class T>
  void ElementReference<Dimension2, 3>::
  AddConstantStiffnessMatrixGen(int m, int n, const T& C, VirtualMatrix<T>& A) const
  {
    Vector<int> num(this->nb_dof_loc);
    Vector<T> val(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      num(i) = n + i;

    for (int i = 0; i < this->nb_dof_loc; i++)
      {
	for (int j = 0; j < this->nb_dof_loc; j++)
	  val(j) = const_div_matrix(i, j)*C;
	
	A.AddInteractionRow(m+i, this->nb_dof_loc, num, val);
      }          
  }


  //! integrals against divergence of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension2, 3>
  ::ComputeIntegralDivRef(const Vector1 & feval, Vector2& res) const
  {
    this->ApplyRh(feval, res);
  }
  
  
  //! surface integrals against divergence of basis functions
  template<class Vector1,class Vector2>
  void ElementReference<Dimension2, 3>
  ::ComputeIntegralSurfaceDivRef(const Vector1 & feval, Vector2& res, int n) const
  {
    typename Vector2::value_type one; SetComplexOne(one);
    res.Fill(0);
    this->ApplyNablaSh(one, n, feval, res);
  }


  //! projection from values on dof points to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension2, 3>
  ::ComputeProjectionDofRef(const VectReal_wp& feval, VectReal_wp& contrib) const
  {
    contrib = feval;
    SolveQR(MatProjectionDof, tauProjectionDof, contrib);
  }

  
  //! projection from values on dof points to components on degrees of freedom
  /*!
    \param[in] feval values of u on dof points
    \param[out] contrib components of u on degrees of freedom
   */
  void ElementReference<Dimension2, 3>
  ::ComputeProjectionDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib) const
  {
    VectReal_wp xreal(feval.GetM()), ximag(feval.GetM());
    for (int i = 0; i < feval.GetM(); i++)
      {
        xreal(i) = real(feval(i));
        ximag(i) = imag(feval(i));
      }
    
    SolveQR(MatProjectionDof, tauProjectionDof, xreal);
    SolveQR(MatProjectionDof, tauProjectionDof, ximag);
    
    contrib.Reallocate(this->nb_dof_loc);
    for (int i = 0; i < this->nb_dof_loc; i++)
      contrib(i) = Complex_wp(xreal(i), ximag(i));
  }


  void ElementReference<Dimension2, 3>
  ::ComputeProjectionSurfaceDofRef(const VectReal_wp& feval, VectReal_wp& contrib, int num_loc) const
  {
    VectReal_wp feval_normal(feval.GetM()/2);
    R2 vec_u;
    for (int i = 0; i < feval_normal.GetM(); i++)
      {
	vec_u.Init(feval(2*i), feval(2*i+1));
	feval_normal(i) = DotProd(this->NormaleLoc(num_loc), vec_u);
      }
    
    this->element_surface->ComputeProjectionDofRef(feval_normal, contrib);
  }


  void ElementReference<Dimension2, 3>
  ::ComputeProjectionSurfaceDofRef(const VectComplex_wp& feval, VectComplex_wp& contrib, int num_loc) const
  {
    VectComplex_wp feval_normal(feval.GetM()/2);
    R2_Complex_wp vec_u;
    for (int i = 0; i < feval_normal.GetM(); i++)
      {
	vec_u.Init(feval(2*i), feval(2*i+1));
	feval_normal(i) = DotProd(this->NormaleLoc(num_loc), vec_u);
      }
    
    this->element_surface->ComputeProjectionDofRef(feval_normal, contrib);
  }  
#endif
  
}
  
#define MONTJOIE_FILE_FACE_REFERENCE_CXX
#endif
