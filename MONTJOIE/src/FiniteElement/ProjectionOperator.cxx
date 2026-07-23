#ifndef MONTJOIE_FILE_PROJECTION_OPERATOR_CXX

namespace Montjoie
{
  
  //! default constructor
  FiniteElementProjector::FiniteElementProjector()
  {
    this->type_projector = FiniteElementProjector::SPARSE;
    nb_row = 0;
    nb_nodes_u = 0;
    nb_nodes_v = 0;
  }
  

  //! setting operator to a permutation matrix
  /*!
    \param[in] icol column numbers
    \param[in] irow row numbers
    \param[in] n number of rows
    V = Ch U such that V(irow) = U(icol)
   */
  void FiniteElementProjector
  ::SetIdentity(const IVect& icol, const IVect& irow, int n)
  {
    type_projector = IDENTITY;
    num_row.Reallocate(icol.GetM());
    nb_row = n;
    for (int i = 0; i < icol.GetM(); i++)
      num_row(icol(i)) = irow(i);
  }
  

  //! projection of Un on Vn
  /*!
    Vn = Ch Un with Ch equal to a permutation matrix or equal
    to Ch_2 Ch_1 in 2-D, and equal to Ch_3 Ch_2 Ch_1 in 3-D    
   */
  template<class T1, class T2, int nb_unknowns>
  void FiniteElementProjector
  ::Project(const Vector<TinyVector<T1, nb_unknowns> >& Un,
            Vector<TinyVector<T2, nb_unknowns> >& Vn) const
  {
    if (type_projector == IDENTITY)
      {
        Vn.Reallocate(nb_row);
        FillZero(Vn);
        for (int i = 0; i < num_row.GetM(); i++)
          Vn(num_row(i)) = Un(i);
        
        return;
      }
    
    Vn.Reallocate(nb_nodes_v);
    
    Vector<T1> Utmp(nb_nodes_u), Vtmp(nb_nodes_v);
    // loop on all unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we get Un related to the considered unknown
	for (int i = 0; i < nb_nodes_u; i++)
	  Utmp(i) = Un(i)(m);
	
	ProjectScalar(Utmp, Vtmp);
	
	// we fill Vn
	for (int i = 0; i < nb_nodes_v; i++)
	  Vn(i)(m) = Vtmp(i);
      } 
      
  }
  

  //! projection of Un on Vn
  /*!
    Vn = Ch Un with Ch equal to a permutation matrix or equal
    to Ch_2 Ch_1 in 2-D, and equal to Ch_3 Ch_2 Ch_1 in 3-D    
   */
  template<class T1, class T2, int p, int q>
  void FiniteElementProjector
  ::Project(const Vector<TinyMatrix<T1, General, p, q> >& Un,
            Vector<TinyMatrix<T2, General, p, q> >& Vn) const
  {
    if (type_projector == IDENTITY)
      {
        Vn.Reallocate(nb_row);
        FillZero(Vn);
        for (int i = 0; i < num_row.GetM(); i++)
          Vn(num_row(i)) = Un(i);
        
        return;
      }
    
    Vn.Reallocate(nb_nodes_v);
    
    Vector<T1> Utmp(nb_nodes_u), Vtmp(nb_nodes_v);
    // loop on all unknowns
    for (int m = 0; m < p; m++)
      for (int n = 0; n < q; n++)
        {
          // we get Un related to the considered unknown
          for (int i = 0; i < nb_nodes_u; i++)
            Utmp(i) = Un(i)(m, n);
          
          // application of the projection operator
	  ProjectScalar(Utmp, Vtmp);
	  
          // we fill Vn
          for (int i = 0; i < nb_nodes_v; i++)
            Vn(i)(m, n) = Vtmp(i);
        }     
  }
  
  
  //! transpose projection of Un on Vn
  /*!
    Vn = Ch* Un with Ch equal to a permutation matrix or equal
    to Ch_2 Ch_1 in 2-D, and equal to Ch_3 Ch_2 Ch_1 in 3-D    
    This function is the transpose of Project
   */
  template<class T1, class T2, int nb_unknowns>
  void FiniteElementProjector::
  TransposeProject(const Vector<TinyVector<T1, nb_unknowns> >& Un,
                   Vector<TinyVector<T2, nb_unknowns> >& Vn) const
  {
    if (type_projector == IDENTITY)
      {
        Vn.Reallocate(num_row.GetM());
        for (int i = 0; i < num_row.GetM(); i++)
          Vn(i) = Un(num_row(i));
        
        return;
      }
    
    Vn.Reallocate(nb_nodes_u);
    
    Vector<T1> Utmp(nb_nodes_v), Vtmp(nb_nodes_u);
    // loop on all unknowns
    for (int m = 0; m < nb_unknowns; m++)
      {
	// we get Un related to the considered unknown
	for (int i = 0; i < nb_nodes_v; i++)
	  Utmp(i) = Un(i)(m);
	
	// application of the projection operator
	TransposeProjectScalar(Utmp, Vtmp);
	
	// we fill Vn
	for (int i = 0; i < nb_nodes_u; i++)
	  Vn(i)(m) = Vtmp(i);
      } 
    
  }
  
  
  //! projection of Un on Vn
  /*!
    Vn = Ch Un with Ch equal to a permutation matrix or equal
    to Ch_2 Ch_1 in 2-D, and equal to Ch_3 Ch_2 Ch_1 in 3-D    
   */
  template<class T1, class T2>
  void FiniteElementProjector
  ::Project(const Vector<T1>& Un, Vector<T2>& Vn) const
  {
    if (type_projector == IDENTITY)
      {
        Vn.Reallocate(nb_row);
        FillZero(Vn);
        for (int i = 0; i < num_row.GetM(); i++)
          Vn(num_row(i)) = Un(i);
        
        return;
      }
    
    Vn.Reallocate(nb_nodes_v);
    ProjectScalar(Un, Vn);
  }
  
  
  //! transpose projection of Un on Vn
  /*!
    Vn = Ch* Un with Ch equal to a permutation matrix or equal
    to Ch_2 Ch_1 in 2-D, and equal to Ch_3 Ch_2 Ch_1 in 3-D    
    This function is the transpose of Project
   */
  template<class T1, class T2>
  void FiniteElementProjector
  ::TransposeProject(const Vector<T1>& Un, Vector<T2>& Vn) const
  {
    if (type_projector == IDENTITY)
      {
        Vn.Reallocate(num_row.GetM());
        for (int i = 0; i < num_row.GetM(); i++)
          Vn(i) = Un(num_row(i));
        
        return;
      }

    Vn.Reallocate(nb_nodes_u);
    TransposeProjectScalar(Un, Vn);
  }


  /*****************************
   *  TensorizedProjector_Base *
   *****************************/
  
  
  //! default constructor
  template<class Dimension>
  TensorizedProjector_Base<Dimension>::TensorizedProjector_Base()
  {
  }
  
  
  template<class Dimension>
  size_t TensorizedProjector_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = 0;
    for (int i = 0; i < ch_loc.GetM(); i++)
      taille += ch_loc(i).GetMemorySize();

    return taille;
  }
  

  //! setting directly operator Ch
  /*!
    \param[in] m number of the operator to modify
    \param[in] A sparse matrix part of Ch
    Ch is equal to Ch_2 Ch_1 in 2-D, and equal to Ch_3 Ch_2 Ch_1 in 3-D
    With this function, you set Ch_m
   */
  template<class Dimension> template<class MatrixSparse>
  void TensorizedProjector_Base<Dimension>::SetChOperator(int m, const MatrixSparse& A)
  {
    type_projector = SPARSE;
    Copy(A, ch_loc(m));
  }

  
  /******************
   * DenseProjector *
   ******************/

  
  //! default constructor
  template<class Dimension>
  DenseProjector<Dimension>::DenseProjector()
  {
  }
  
  
  template<class Dimension>
  size_t DenseProjector<Dimension>::GetMemorySize() const
  {
    return ch_loc.GetMemorySize();
  }
  

  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] points_div the 2-D or 3-D points where we want to know the projection 
    (interpolation)
    The function ComputeValuesPhiNodalRef is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  template<class Dimension>
  void DenseProjector<Dimension>::
  InitProjector(const ElementGeomReference<Dimension>& FaceBasis,
                const typename Dimension::VectR_N& points_div)
  {
    this->type_projector = FiniteElementProjector::SPARSE;
    VectReal_wp phi;
    ch_loc.Reallocate(points_div.GetM(), FaceBasis.GetNbPointsNodalElt());
    this->nb_nodes_u = ch_loc.GetN();
    this->nb_nodes_v = ch_loc.GetM();
    for (int i = 0; i < points_div.GetM(); i++)
      {	
	FaceBasis.ComputeValuesPhiNodalRef(points_div(i), phi);
	for (int j = 0; j < ch_loc.GetN(); j++)
	  ch_loc(i, j) = phi(j);
      }
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (basis functions are defined in this object)
    \param[in] step_subdiv the 1-D points where we want to project (not used)
    \param[in] points_div the 2-D or 3-D points where we want to know the projection
    (interpolation)
    The function ComputeValuesPhiNodalRef is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  template<class Dimension>
  void DenseProjector<Dimension>::
  InitProjector(const ElementGeomReference<Dimension>& FaceBasis, const VectReal_wp& step_subdiv,
                const typename Dimension::VectR_N& points_div)
  {
    Init(FaceBasis, points_div);
  }
  
  
  //! projection of Un on Vn
  /*!
    Vn = Ch Un with Ch = dense matrix ch_loc
   */
  template<class Dimension>
  void DenseProjector<Dimension>::
  ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const
  {
    Mlt(ch_loc, Un, Vn);
  }


  //! projection of Un on Vn
  /*!
    Vn = Ch Un with Ch = dense matrix ch_loc
   */
  template<class Dimension>
  void DenseProjector<Dimension>::
  ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const
  {
    Mlt(ch_loc, Un, Vn);
  }


  //! projection of Un on Vn
  /*!
    Vn = Ch^T Un with Ch = dense matrix ch_loc
   */
  template<class Dimension>
  void DenseProjector<Dimension>::
  TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const
  {
    Mlt(SeldonTrans, ch_loc, Un, Vn);
  }


  //! projection of Un on Vn
  /*!
    Vn = Ch^T Un with Ch = dense matrix ch_loc
   */
  template<class Dimension>
  void DenseProjector<Dimension>::
  TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const
  {
    Mlt(SeldonTrans, ch_loc, Un, Vn);
  }


  //! default constructor
  TensorizedProjector<Dimension2>::TensorizedProjector()
  {
  }
  

  size_t TensorizedProjector<Dimension2>::GetMemorySize() const
  {
    size_t taille = TensorizedProjector_Base<Dimension2>::GetMemorySize();
    taille += NumNodes.GetMemorySize();
    taille += Ux_real.GetMemorySize();
    taille += Ux_cplx.GetMemorySize();
    return taille;
  }

  
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] points_div the 2-D points where we want to know the projection (interpolation)
    The function ComputeValuesNodalPhi1D is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  void TensorizedProjector<Dimension2>::
  Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectR2& points_div)
  {
    int N = points_div.GetM();
    if (N <= 0)
      {
	return;
      }
    
    VectReal_wp step_x(N), step_y(N);
    for (int i = 0; i < N; i++)
      {
        step_x(i) = points_div(i)(0);
        step_y(i) = points_div(i)(1);
      }
    
    Sort(step_x); Sort(step_y);
    Real_wp xprec = step_x(0);
    int nb = 1;
    for (int i = 1; i < N; i++)
      {
        if (abs(step_x(i) - xprec) > R2::threshold)
          {
            // new point
            step_x(nb) = step_x(i);
            xprec = step_x(nb);
            nb++;
          }
      }
    step_x.Resize(nb);
    
    xprec = step_y(0);
    nb = 1;
    for (int i = 1; i < N; i++)
      {
        if (abs(step_y(i) - xprec) > R2::threshold)
          {
            // new point
            step_y(nb) = step_y(i);
            xprec = step_y(nb);
            nb++;
          }
      }
    step_y.Resize(nb);
    
    InitProjector(FaceBasis, step_x, step_y, points_div);
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] step_x 1-D points along x and y direction
    \param[in] points_div the 2-D points where we want to know the projection (interpolation)
    The function ComputeValuesNodalPhi1D is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  void TensorizedProjector<Dimension2>::
  Init(const ElementGeomReference<Dimension2>& FaceBasis, const VectReal_wp& step_x,
       const VectR2& points_div)
  {
    InitProjector(FaceBasis, step_x, step_x, points_div);
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] step_x 1-D points along x direction
    \param[in] step_y 1-D points along y direction
    \param[in] points_div the 2-D points where we want to know the projection (interpolation)
    The function ComputeValuesNodalPhi1D is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  void TensorizedProjector<Dimension2>::
  InitProjector(const ElementGeomReference<Dimension2>& FaceBasis, const VectReal_wp& step_x,
           const VectReal_wp& step_y, const VectR2& points_div)
  {
    // construct of NumNodes(i,j) to have matching array between 
    // 2-indices numbering and 1-index numbering
    // for target points  points _div
    int nx = step_x.GetM()-1;
    int ny = step_y.GetM()-1;
    int nb_points = points_div.GetM();
    int r = FaceBasis.GetGeometryOrder();
    this->type_projector = FiniteElementProjector::SPARSE;
    this->nb_nodes_u = (r+1)*(r+1);
    this->nb_nodes_v = nb_points;
    
    NumNodes.Reallocate(nx+1, ny+1); NumNodes.Fill(-1);
    IVect DoublePoint(nb_points); DoublePoint.Fill(-1); 
    for (int i = 0; i < points_div.GetM(); i++)
      {
	int i1 = -1, i2 = -1;
	for (int j = 0; j <= nx; j++)
	  if (abs(step_x(j) - points_div(i)(0)) <= R2::threshold)
	    i1 = j;
	
	for (int j = 0; j <= ny; j++)
	  if (abs(step_y(j) - points_div(i)(1)) <= R2::threshold)
	    i2 = j;
	
	if ((i1 == -1)||(i2 == -1))
	  {
	    cout<<"Problem for the computation of projector "<<endl;
	    cout<<"Check that the 2-D points are obtained by tensorization " 
		<<"of 1-D points "<<endl;
            
	    DISP(step_x); DISP(step_y); DISP(points_div);
	    abort();
	  }
        
        if (NumNodes(i1, i2) >= 0)
          DoublePoint(i) = NumNodes(i1, i2);
        else
          NumNodes(i1, i2) = i;
      }    
    
    // computation of C1 and C2
    Matrix<Real_wp, General, ArrayRowSparse> Ch1((nx+1)*(r+1), (r+1)*(r+1));
    Matrix<Real_wp, General, ArrayRowSparse> Ch2(nb_points, (nx+1)*(r+1));
    Matrix<Real_wp> val_phix(r+1, nx+1);
    VectReal_wp phi;
    for (int i = 0; i <= nx; i++)
      {
	FaceBasis.ComputeValuesNodalPhi1D(step_x(i), phi);
	for (int j = 0; j <= r; j++)
	  val_phix(j, i) = phi(j);
      }

    Matrix<Real_wp> val_phiy(r+1, ny+1);
    for (int i = 0; i <= ny; i++)
      {
	FaceBasis.ComputeValuesNodalPhi1D(step_y(i), phi);
	for (int j = 0; j <= r; j++)
	  val_phiy(j, i) = phi(j);
      }
    
    // interpolation along x
    for (int i1 = 0; i1 <= nx; i1++)
      for (int j2 = 0; j2 <= r; j2++)
	for (int j1 = 0; j1 <= r; j1++)
	  Ch1.AddInteraction(i1*(r+1) + j2, FaceBasis.GetNumNodes2D(j1,j2), val_phix(j1, i1));
    
    // interpolation along y
    for (int i1 = 0; i1 <= nx; i1++)
      for (int i2 = 0; i2 <= ny; i2++)
        {
          int node = NumNodes(i1, i2);
          if (node >= 0)
            for (int j2 = 0; j2 <= r; j2++)
              Ch2.AddInteraction(NumNodes(i1,i2), i1*(r+1)+j2, val_phiy(j2,i2));
        }
    
    // duplicate points
    for (int i = 0; i < nb_points; i++)
      if (DoublePoint(i) >= 0)
        {
          int i2 = DoublePoint(i);
          for (int j2 = 0; j2 < Ch2.GetRowSize(i2); j2++)
            Ch2.AddInteraction(i, Ch2.Index(i2, j2), Ch2.Value(i2, j2));
        }
    
    // conversion to CSR format for efficiency
    Seldon::Copy(Ch1, this->ch_loc(0));
    Seldon::Copy(Ch2, this->ch_loc(1));

    Ux_real.Reallocate(this->ch_loc(0).GetM());
    Ux_cplx.Reallocate(this->ch_loc(0).GetM());    
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] lob definition of source points
    \param[in] NumNodes2D tensorial-scalar matching array for source points
    \param[in] step_subdiv the 1-D target points where we want to project
    (the \f$ \zeta_i$ \f$ points)
    \param[in] points_div the 2-D target points where we want to know the projection
    (interpolation)
  */
  template<class TypeFunction>
  void TensorizedProjector<Dimension2>::
  InitProjector(const TypeFunction& lob, const Matrix<int>& NumNodes2D,
           const VectReal_wp& step_subdiv, const VectR2& points_div)
  {
    InitProjector(lob, NumNodes2D, step_subdiv, step_subdiv, points_div);
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] lob definition of source points
    \param[in] NumNodes2D tensorial-scalar matching array for source points
    \param[in] step_x the 1-D target points along x
    \param[in] step_y the 1-D target points along y
    \param[in] points_div the 2-D target points where we want to know the projection
    (interpolation)
  */
  template<class TypeFunction>
  void TensorizedProjector<Dimension2>::
  InitProjector(const TypeFunction& lob, const Matrix<int>& NumNodes2D,
           const VectReal_wp& step_x, const VectReal_wp& step_y, const VectR2& points_div)
  {
    // construct of NumNodes(i,j) to have matching array 
    // between 2-indices numbering and 1-index numbering
    // tensorial-scalar matching array for target points
    int nx = step_x.GetM()-1;
    int ny = step_y.GetM()-1;
    
    NumNodes.Reallocate(nx+1, ny+1); NumNodes.Fill(-1);
    for (int i = 0; i < points_div.GetM(); i++)
      {
	int i1 = -1, i2 = -1;
	for (int j = 0; j <= nx; j++)
	  if (abs(step_x(j) - points_div(i)(0)) <= R2::threshold)
	    i1 = j;
	
	for (int j = 0; j <= ny; j++)
	  if (abs(step_y(j) - points_div(i)(1)) <= R2::threshold)
	    i2 = j;
	
	if ((i1 == -1)||(i2 == -1))
	  {
	    cout<<"Problem for the computation of projector "<<endl;
	    cout<<"Check that the 2-D points are obtained by tensorization of 1-D points "<<endl;
	    DISP(step_x); DISP(step_y); DISP(points_div);
	    abort();
	  }
	NumNodes(i1,i2) = i;
      }
    
    int r = lob.GetGeometryOrder();
    this->type_projector = FiniteElementProjector::SPARSE;
    this->nb_nodes_u = (r+1)*(r+1);
    this->nb_nodes_v = (nx+1)*(ny+1);

    Matrix<Real_wp> val_phix(r+1, nx+1), val_phiy(r+1, ny+1);
    for (int i = 0; i <= nx; i++)
      for (int j = 0; j <= r; j++)
	val_phix(j, i) = lob.EvaluatePhi(j, step_x(i));

    for (int i = 0; i <= ny; i++)
      for (int j = 0; j <= r; j++)
	val_phiy(j, i) = lob.EvaluatePhi(j, step_y(i));
    
    // computation of C1 and C2
    Matrix<Real_wp, General, ArrayRowSparse> Ch1((nx+1)*(r+1), (r+1)*(r+1));
    Matrix<Real_wp, General, ArrayRowSparse> Ch2((nx+1)*(ny+1), (nx+1)*(r+1));
    
    // DISP(val_phi);
    for (int i1 = 0; i1 <= nx; i1++)
      for (int j2 = 0; j2 <= r; j2++)
	for (int j1 = 0; j1 <= r; j1++)
	  Ch1.AddInteraction(i1*(r+1)+j2, NumNodes2D(j1,j2), val_phix(j1,i1));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int i2 = 0; i2 <= ny; i2++)
	{
          int node = NumNodes(i1, i2);
          if (node >= 0)
            for (int j2 = 0; j2 <= r; j2++)
              Ch2.AddInteraction(node, i1*(r+1)+j2, val_phiy(j2,i2));
        }
    
    // conversion to CSR format for efficiency
    Seldon::Copy(Ch1, this->ch_loc(0));
    Seldon::Copy(Ch2, this->ch_loc(1));
    
    Ux_real.Reallocate(this->ch_loc(0).GetM());
    Ux_cplx.Reallocate(this->ch_loc(0).GetM());
  }
  

  //! projection of Un on Vn
  void TensorizedProjector<Dimension2>
  ::ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const
  {
    VectReal_wp& Ux = const_cast<VectReal_wp&>(Ux_real);
    
    Mlt(this->ch_loc(0), Un, Ux);
    Mlt(this->ch_loc(1), Ux, Vn);
  }
  
  
  void TensorizedProjector<Dimension2>
  ::ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const
  {
    VectComplex_wp& Ux = const_cast<VectComplex_wp&>(Ux_cplx);
    
    Mlt(this->ch_loc(0), Un, Ux);
    Mlt(this->ch_loc(1), Ux, Vn);
  }
  
  
  void TensorizedProjector<Dimension2>
  ::TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const
  {
    VectReal_wp& Ux = const_cast<VectReal_wp&>(Ux_real);
    
    Mlt(SeldonTrans, this->ch_loc(1), Un, Ux);
    Mlt(SeldonTrans, this->ch_loc(0), Ux, Vn);
  }
  
  
  void TensorizedProjector<Dimension2>
  ::TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const
  {
    VectComplex_wp& Ux = const_cast<VectComplex_wp&>(Ux_cplx);
    
    Mlt(SeldonTrans, this->ch_loc(1), Un, Ux);
    Mlt(SeldonTrans, this->ch_loc(0), Ux, Vn);
  }

#ifdef MONTJOIE_WITH_THREE_DIM
  
  //! default constructor
  TensorizedProjector<Dimension3>::TensorizedProjector()
  {
  }
  

  size_t TensorizedProjector<Dimension3>::GetMemorySize() const
  {
    size_t taille = TensorizedProjector_Base<Dimension3>::GetMemorySize();
    taille += NumNodes.GetMemorySize();
    taille += Ux_real.GetMemorySize() + Uy_real.GetMemorySize();
    taille += Ux_cplx.GetMemorySize() + Uy_cplx.GetMemorySize();
    return taille;
  }

    
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] points_div the 3-D points where we want to know the projection (interpolation)
    The function ComputeValuesNodalPhi1D is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  void TensorizedProjector<Dimension3>::
  Init(const ElementGeomReference<Dimension3>& FaceBasis, const VectR3& points_div)
  {
    int N = points_div.GetM();
    if (N <= 0)
      return;
    
    VectReal_wp step_x(N), step_y(N), step_z(N);
    for (int i = 0; i < N; i++)
      {
        step_x(i) = points_div(i)(0);
        step_y(i) = points_div(i)(1);
        step_z(i) = points_div(i)(2);
      }
    
    Sort(step_x); Sort(step_y); Sort(step_z);
    Real_wp xprec = step_x(0);
    int nb = 1;
    for (int i = 1; i < N; i++)
      {
        if (abs(step_x(i) - xprec) > R3::threshold)
          {
            // new point
            step_x(nb) = step_x(i);
            xprec = step_x(nb);
            nb++;
          }
      }
    step_x.Resize(nb);
    
    xprec = step_y(0);
    nb = 1;
    for (int i = 1; i < N; i++)
      {
        if (abs(step_y(i) - xprec) > R3::threshold)
          {
            // new point
            step_y(nb) = step_y(i);
            xprec = step_y(nb);
            nb++;
          }
      }
    step_y.Resize(nb);

    xprec = step_z(0);
    nb = 1;
    for (int i = 1; i < N; i++)
      {
        if (abs(step_z(i) - xprec) > R3::threshold)
          {
            // new point
            step_z(nb) = step_z(i);
            xprec = step_z(nb);
            nb++;
          }
      }
    step_z.Resize(nb);
    
    InitProjector(FaceBasis, step_x, step_y, step_z, points_div);
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] step_subdiv 1-D points along x and y direction
    \param[in] points_div the 3-D points where we want to know the projection (interpolation)
    The function ComputeValuesNodalPhi1D is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  void TensorizedProjector<Dimension3>::
  Init(const ElementGeomReference<Dimension3>& FaceBasis,
       const VectReal_wp& step_subdiv, const VectR3& points_div)
  {
    InitProjector(FaceBasis, step_subdiv, step_subdiv, step_subdiv, points_div);
  }
  
    
  //! construction of projection operator
  /*!
    \param[in] FaceBasis definition of finite element (shape functions are defined in this object)
    \param[in] step_x 1-D points along x direction
    \param[in] step_y 1-D points along y direction
    \param[in] step_z 1-D points along z direction
    \param[in] points_div the 3-D points where we want to know the projection (interpolation)
    The function ComputeValuesNodalPhi1D is used to evaluate shape functions on points_div
    Ch operator is equal to (Ch)_{ij} = phi_i(xi_j)
    where phi_i are shape functions of FaceBasis, and xi_j points_div
  */
  void TensorizedProjector<Dimension3>::
  InitProjector(const ElementGeomReference<Dimension3>& FaceBasis, const VectReal_wp& step_x,
                const VectReal_wp& step_y, const VectReal_wp& step_z, const VectR3& points_div)
  {
    // construct of NumNodes(i,j) to have matching array 
    // between 3-indices numbering and 1-index numbering
    int nx = step_x.GetM()-1;
    int ny = step_y.GetM()-1;
    int nz = step_z.GetM()-1;
    int nb_points = points_div.GetM();
    
    NumNodes.Reallocate(nx+1, ny+1, nz+1); NumNodes.Fill(-1);
    IVect DoublePoint(nb_points); DoublePoint.Fill(-1); 
    // DISP(Points); DISP(points_div);
    for (int i = 0; i < points_div.GetM(); i++)
      {
	int i1 = -1, i2 = -1, i3 = -1;
	for (int j = 0; j <= nx; j++)
	  if (abs(step_x(j) - points_div(i)(0)) <= R3::threshold)
	    i1 = j;
	
	for (int j = 0; j <= ny; j++)
	  if (abs(step_y(j) - points_div(i)(1)) <= R3::threshold)
	    i2 = j;
	
	for (int j = 0; j <= nz; j++)
	  if (abs(step_z(j) - points_div(i)(2)) <= R3::threshold)
	    i3 = j;
	
	if ((i1 == -1)||(i2 == -1)||(i3 == -1))
	  {
	    cout<<"Problem for the computation of projector "<<endl;
	    cout<<"Check that the 3-D points are obtained by tensorization of 1-D points "<<endl;
	    DISP(step_x); DISP(step_y); DISP(step_z); DISP(points_div);
	    abort();
	  }
        
        if (NumNodes(i1, i2, i3) >= 0)
          DoublePoint(i) = NumNodes(i1, i2, i3);
        else
          NumNodes(i1, i2, i3) = i;
      }
    
    int r = FaceBasis.GetGeometryOrder();
    this->type_projector = FiniteElementProjector::SPARSE;
    this->nb_nodes_u = (r+1)*(r+1)*(r+1);
    this->nb_nodes_v = nb_points;

    Matrix<Real_wp> val_phix(r+1, nx+1), val_phiy(r+1, ny+1), val_phiz(r+1, nz+1);
    VectReal_wp phi;
    for (int i = 0; i <= nx; i++)
      {
	FaceBasis.ComputeValuesNodalPhi1D(step_x(i), phi);
	for (int j = 0; j <= r; j++)
	  val_phix(j, i) = phi(j);
      }
    for (int i = 0; i <= ny; i++)
      {
	FaceBasis.ComputeValuesNodalPhi1D(step_y(i), phi);
	for (int j = 0; j <= r; j++)
	  val_phiy(j, i) = phi(j);
      }

    for (int i = 0; i <= nz; i++)
      {
	FaceBasis.ComputeValuesNodalPhi1D(step_z(i), phi);
	for (int j = 0; j <= r; j++)
	  val_phiz(j, i) = phi(j);
      }

    // computation of C1, C2 and C3
    Matrix<Real_wp, General, ArrayRowSparse> Ch1((nx+1)*(r+1)*(r+1), (r+1)*(r+1)*(r+1));
    Matrix<Real_wp, General, ArrayRowSparse> Ch2((nx+1)*(ny+1)*(r+1), (nx+1)*(r+1)*(r+1));
    Matrix<Real_wp, General, ArrayRowSparse> Ch3(nb_points, (nx+1)*(ny+1)*(r+1));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int j3 = 0; j3 <= r; j3++)
	for (int j2 = 0; j2 <= r; j2++)
	  for (int j1 = 0; j1 <= r; j1++)
	    Ch1.AddInteraction(i1*(r+1)*(r+1)+j2*(r+1)+j3,
                               FaceBasis.GetNumNodes3D(j1,j2,j3), val_phix(j1,i1));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int i2 = 0; i2 <= ny; i2++)
	for (int j2 = 0; j2 <= r; j2++)
	  for (int j3 = 0; j3 <= r; j3++)
	    Ch2.AddInteraction(i1*(r+1)*(ny+1)+i2*(r+1)+j3, i1*(r+1)*(r+1)+j2*(r+1)+j3,
                               val_phiy(j2,i2));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int i2 = 0; i2 <= ny; i2++)
	for (int i3 = 0; i3 <= nz; i3++)
	  {
            int node = NumNodes(i1, i2, i3);
            if (node >= 0)
              for (int j3 = 0; j3 <= r; j3++)
                Ch3.AddInteraction(node, i1*(r+1)*(ny+1)+i2*(r+1)+j3, val_phiz(j3,i3));
          }

    // duplicate points
    for (int i = 0; i < nb_points; i++)
      if (DoublePoint(i) >= 0)
        {
          int i2 = DoublePoint(i);
          for (int j2 = 0; j2 < Ch3.GetRowSize(i2); j2++)
            Ch3.AddInteraction(i, Ch3.Index(i2, j2), Ch3.Value(i2, j2));
        }
    
    // conversion in CSR format for efficiency
    Seldon::Copy(Ch1, this->ch_loc(0));
    Seldon::Copy(Ch2, this->ch_loc(1));
    Seldon::Copy(Ch3, this->ch_loc(2));

    Ux_real.Reallocate((nx+1)*(r+1)*(r+1)); Ux_cplx.Reallocate((nx+1)*(r+1)*(r+1));
    Uy_real.Reallocate((nx+1)*(ny+1)*(r+1)); Uy_cplx.Reallocate((nx+1)*(ny+1)*(r+1));    
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] lob definition of source points
    \param[in] NumNodes3D tensorial-scalar matching array for source points
    \param[in] step_subdiv points for each coordinate
    \param[in] points_div the 3-D target points where we want to know the projection
    (interpolation)
  */
  template<class TypeFunction>
  void TensorizedProjector<Dimension3>::
  InitProjector(const TypeFunction& lob, const Array3D<int>& NumNodes3D,
                const VectReal_wp& step_subdiv, const VectR3& points_div)
  {
    Init(lob, NumNodes3D, step_subdiv,
         step_subdiv, step_subdiv, points_div);
  }
  
  
  //! construction of projection operator
  /*!
    \param[in] lob definition of source points
    \param[in] NumNodes3D tensorial-scalar matching array for source points
    \param[in] step_x points for x-coordinate
    \param[in] step_y points for y-coordinate
    \param[in] step_z points for z-coordinate
    \param[in] points_div the 3-D target points where we want to know the projection
    (interpolation)
  */
  template<class TypeFunction>
  void TensorizedProjector<Dimension3>::
  InitProjector(const TypeFunction& lob, const Array3D<int>& NumNodes3D,
                const VectReal_wp& step_x, const VectReal_wp& step_y, const VectReal_wp& step_z,
                const VectR3& points_div)
  {
    // construct of NumNodes(i,j,k) to have matching array 
    // between 3-subscript numbering and 1-subscript numbering
    // tensorial-scalar matching array for target points
    int nx = step_x.GetM()-1;
    int ny = step_y.GetM()-1;
    int nz = step_z.GetM()-1;
    
    NumNodes.Reallocate(nx+1, ny+1, nz+1); NumNodes.Fill(-1);
    for (int i = 0; i < points_div.GetM(); i++)
      {
	int i1 = -1, i2 = -1, i3 = -1;
	for (int j = 0; j <= nx; j++)
	  if (abs(step_x(j) - points_div(i)(0)) <= R3::threshold)
	    i1 = j;
	
	for (int j = 0; j <= ny; j++)
	  if (abs(step_y(j) - points_div(i)(1)) <= R3::threshold)
	    i2 = j;
	
	for (int j = 0; j <= nz; j++)
	  if (abs(step_z(j) - points_div(i)(2)) <= R3::threshold)
	    i3 = j;
	
	if ((i1 == -1)||(i2 == -1)||(i3 == -1))
	  {
	    cout<<"Problem for the computation of projector "<<endl;
	    cout<<"Check that the 3-D points are obtained by tensorization of 1-D points "<<endl;
	    DISP(step_x); DISP(step_y); DISP(step_z); DISP(points_div);
	    abort();
	  }
        
	NumNodes(i1,i2,i3) = i;
      }
    
    int r = lob.GetGeometryOrder();
    this->type_projector = FiniteElementProjector::SPARSE;
    this->nb_nodes_u = (r+1)*(r+1)*(r+1);
    this->nb_nodes_v = (nx+1)*(ny+1)*(nz+1);

    // computation of C1, C2 and C3
    Matrix<Real_wp, General, ArrayRowSparse> Ch1((nx+1)*(r+1)*(r+1), (r+1)*(r+1)*(r+1));
    Matrix<Real_wp, General, ArrayRowSparse> Ch2((nx+1)*(ny+1)*(r+1), (nx+1)*(r+1)*(r+1));
    Matrix<Real_wp, General, ArrayRowSparse> Ch3((nx+1)*(ny+1)*(nz+1), (nx+1)*(ny+1)*(r+1));
    Matrix<Real_wp> val_phix(r+1, nx+1), val_phiy(r+1, ny+1), val_phiz(r+1, nz+1);
    for (int i = 0; i <= nx; i++)
      for (int j = 0; j <= r; j++)
	val_phix(j,i) = lob.EvaluatePhi(j, step_x(i));

    for (int i = 0; i <= ny; i++)
      for (int j = 0; j <= r; j++)
	val_phiy(j,i) = lob.EvaluatePhi(j, step_y(i));

    for (int i = 0; i <= nz; i++)
      for (int j = 0; j <= r; j++)
	val_phiz(j,i) = lob.EvaluatePhi(j, step_z(i));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int j3 = 0; j3 <= r; j3++)
	for (int j2 = 0; j2 <= r; j2++)
	  for (int j1 = 0; j1 <= r; j1++)
	    Ch1.AddInteraction(i1*(r+1)*(r+1)+j2*(r+1)+j3, NumNodes3D(j1,j2,j3), val_phix(j1,i1));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int i2 = 0; i2 <= ny; i2++)
	for (int j2 = 0; j2 <= r; j2++)
	  for (int j3 = 0; j3 <= r; j3++)
	    Ch2.AddInteraction(i1*(r+1)*(ny+1)+i2*(r+1)+j3, i1*(r+1)*(r+1)+j2*(r+1)+j3,
                               val_phiy(j2,i2));
    
    for (int i1 = 0; i1 <= nx; i1++)
      for (int i2 = 0; i2 <= ny; i2++)
	for (int i3 = 0; i3 <= nz; i3++)
	  {
            int node = NumNodes(i1, i2, i3);
            if (node >= 0)
              for (int j3 = 0; j3 <= r; j3++)
                Ch3.AddInteraction(node, i1*(r+1)*(ny+1)+i2*(r+1)+j3, val_phiz(j3,i3));
          }
    
    // conversion in CSR format for efficiency
    Seldon::Copy(Ch1, this->ch_loc(0));
    Seldon::Copy(Ch2, this->ch_loc(1));
    Seldon::Copy(Ch3, this->ch_loc(2));
    
    Ux_real.Reallocate((nx+1)*(r+1)*(r+1)); Ux_cplx.Reallocate((nx+1)*(r+1)*(r+1));
    Uy_real.Reallocate((nx+1)*(ny+1)*(r+1)); Uy_cplx.Reallocate((nx+1)*(ny+1)*(r+1));
  }


  void TensorizedProjector<Dimension3>::AllocateTemporaryVectors()
  {
    this->type_projector = FiniteElementProjector::SPARSE;
    this->nb_nodes_u = this->ch_loc(0).GetN();
    this->nb_nodes_v = this->ch_loc(2).GetM();
    Ux_real.Reallocate(this->ch_loc(0).GetM());
    Ux_cplx.Reallocate(this->ch_loc(0).GetM());
    Uy_real.Reallocate(this->ch_loc(1).GetM());
    Uy_cplx.Reallocate(this->ch_loc(1).GetM());
  }
  

  //! projection of Un on Vn
  void TensorizedProjector<Dimension3>
  ::ProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const
  {
    VectReal_wp& Ux = const_cast<VectReal_wp&>(Ux_real);
    VectReal_wp& Uy = const_cast<VectReal_wp&>(Uy_real);
    
    Mlt(this->ch_loc(0), Un, Ux);
    Mlt(this->ch_loc(1), Ux, Uy);
    Mlt(this->ch_loc(2), Uy, Vn);
  }
  
  
  void TensorizedProjector<Dimension3>
  ::ProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const
  {
    VectComplex_wp& Ux = const_cast<VectComplex_wp&>(Ux_cplx);
    VectComplex_wp& Uy = const_cast<VectComplex_wp&>(Uy_cplx);
    
    Mlt(this->ch_loc(0), Un, Ux);
    Mlt(this->ch_loc(1), Ux, Uy);
    Mlt(this->ch_loc(2), Uy, Vn);
  }
  
  
  void TensorizedProjector<Dimension3>
  ::TransposeProjectScalar(const VectReal_wp& Un, VectReal_wp& Vn) const
  {
    VectReal_wp& Ux = const_cast<VectReal_wp&>(Ux_real);
    VectReal_wp& Uy = const_cast<VectReal_wp&>(Uy_real);
    
    Mlt(SeldonTrans, this->ch_loc(2), Un, Uy);
    Mlt(SeldonTrans, this->ch_loc(1), Uy, Ux);
    Mlt(SeldonTrans, this->ch_loc(0), Ux, Vn);
  }
  
  
  void TensorizedProjector<Dimension3>
  ::TransposeProjectScalar(const VectComplex_wp& Un, VectComplex_wp& Vn) const
  {
    VectComplex_wp& Ux = const_cast<VectComplex_wp&>(Ux_cplx);
    VectComplex_wp& Uy = const_cast<VectComplex_wp&>(Uy_cplx);
    
    Mlt(SeldonTrans, this->ch_loc(2), Un, Uy);
    Mlt(SeldonTrans, this->ch_loc(1), Uy, Ux);
    Mlt(SeldonTrans, this->ch_loc(0), Ux, Vn);
  }
  
#endif

  /*****************************
   * FiniteElementInterpolator *
   *****************************/

  
  //! default destructor
  FiniteElementInterpolator::~FiniteElementInterpolator()
  {
    for (int i = 0; i < stored_proj.GetM(); i++)
      {
	if (stored_proj(i) != NULL)
	  {
	    delete stored_proj(i);
	    stored_proj(i) = NULL;
	  }
      }
  }
  

  size_t FiniteElementInterpolator::GetMemorySize() const
  {
    size_t taille = stored_proj.GetMemorySize();
    for (int i = 0; i < stored_proj.GetM(); i++)
      if (stored_proj(i) != NULL)
        taille += stored_proj(i)->GetMemorySize();
    
    return taille;
  }
  
  
  //! initialisation of the projector
  /*!
    \param[in] Fb finite elements for all the elements of the mesh
    \param[in] points_div 2-D or 3-D points where the solution should be interpolated
    Once the initialisation process is done, the interpolation to points_div
    can be performed by calling Project/ TransposeProject/ ProjectScalar, etc
   */
  template<class FeElement, class Vector1>
  void FiniteElementInterpolator
  ::InitProjection(const Vector<FeElement*>& Fb, const Vector1& points_div)
  {
    Vector<int> ComputedProj(4);
    stored_proj.Reallocate(4);
    FiniteElementProjector* null_ptr = NULL;
    stored_proj.Fill(null_ptr);
    ComputedProj.Fill(-1);
    for (int i = 0; i < Fb.GetM(); i++)
      {
	int type = Fb(i)->GetHybridType();
	if (ComputedProj(type) == -1)
	  {
	    ComputedProj(type) = type;
	    stored_proj(type) = Fb(i)->GetNewNodalInterpolation();
	    stored_proj(type)->Init(*Fb(i), points_div(type));
	  }
	
	//proj(i) = stored_proj(type);
      }
  }
  

  //! initialisation of the projector
  /*!
    \param[in] Fb finite elements for all the elements of the mesh
    \param[in] step_div 1-D points where the solution should be interpolated
    \param[in] points_div 2-D or 3-D points where the solution should be interpolated
    Once the initialisation process is done, the interpolation to points_div
    can be performed by calling Project/ TransposeProject/ ProjectScalar, etc
   */
  template<class FeElement, class Vector1, class Vector2>
  void FiniteElementInterpolator
  ::InitProjection(const Vector<FeElement*>& Fb,
		   const Vector1& step_div, const Vector2& points_div)
  {
    Vector<int> ComputedProj(4);
    stored_proj.Reallocate(4);
    FiniteElementProjector* null_ptr = NULL;
    stored_proj.Fill(null_ptr);
    ComputedProj.Fill(-1);
    for (int i = 0; i < Fb.GetM(); i++)
      {
	int type = Fb(i)->GetHybridType();
	if (ComputedProj(type) == -1)
	  {
	    ComputedProj(type) = type;
	    stored_proj(type) = Fb(i)->GetNewNodalInterpolation();
	    stored_proj(type)->Init(*Fb(i), step_div, points_div(type));
	  }
	
	//proj(i) = stored_prod(type);
      }
  }


  //! Computes the prolongation operator from a coarse finite element to a fine finite element
  /*!
    \param[in] elt_fine list of fine finite elements
    \param[in] rf fine order of approximation 
    \param[in] elt_coarse list of coarse finite elements
    \param[in] rc coarse order of approximation
    \param[out] ProlongationMatrix operator for each reference element (if stored)
   */
  template<class Dimension>
  void FiniteElementInterpolator
  ::ComputeLocalProlongationLowOrder(const Vector<const ElementReference_Dim<Dimension>* >& elt_fine,
				     int rf,
				     const Vector<const ElementReference_Dim<Dimension>* >& elt_coarse,
				     int rc,
				     TinyVector<Matrix<Real_wp>, 4>& ProlongationMatrix)
  {
    Vector<int> ComputedProj(4);
    FiniteElementProjector* null_ptr = NULL;
    stored_proj.Reallocate(4);
    stored_proj.Fill(null_ptr);
    ComputedProj.Fill(-1);
    for (int i = 0; i < elt_fine.GetM(); i++)
      {
	int type_elt = elt_fine(i)->GetHybridType();
	if ((ComputedProj(type_elt) == -1) && (elt_fine(i)->GetOrder() == rf))
	  {
	    if ((elt_coarse(i)->GetHybridType() != type_elt) || 
		(elt_coarse(i)->GetOrder() != rc))
	      {
		cout << "Incompatible coarse and fine elements" << endl;
		abort();
	      }
	    
	    ComputedProj(type_elt) = type_elt;
	    stored_proj(type_elt) = elt_fine(i)->GetNewNodalInterpolation();
	    elt_fine(i)->ComputeLocalProlongationLowOrder(*stored_proj(type_elt), ProlongationMatrix(type_elt),
							  *elt_coarse(i)); 
	  }
      }
  }
  

  //! Computes the prolongation operator from a coarse finite element to a fine finite element
  /*!
    \param[in] elt_fine list of fine finite elements
    \param[in] rf fine order of approximation 
    \param[in] elt_coarse list of coarse finite elements
    \param[in] rc coarse order of approximation
    \param[out] ProlongationMatrix operator for each reference element (if stored)
   */
  template<class Dimension>
  void FiniteElementInterpolator
  ::ComputeLocalProlongation(const Vector<const ElementReference_Dim<Dimension>* >& elt_fine,
			     int rf,
			     const Vector<const ElementReference_Dim<Dimension>* >& elt_coarse,
			     int rc,
			     TinyVector<Matrix<Real_wp>, 4>& ProlongationMatrix)
  {
    Vector<int> ComputedProj(4);
    FiniteElementProjector* null_ptr = NULL;
    stored_proj.Reallocate(4);
    stored_proj.Fill(null_ptr);
    ComputedProj.Fill(-1);
    for (int i = 0; i < elt_fine.GetM(); i++)
      {
	int type_elt = elt_fine(i)->GetHybridType();
	if ((ComputedProj(type_elt) == -1) && (elt_fine(i)->GetOrder() == rf))
	  {
	    if ((elt_coarse(i)->GetHybridType() != type_elt) || 
		(elt_coarse(i)->GetOrder() != rc))
	      {
		cout << "Incompatible coarse and fine elements" << endl;
		abort();
	      }
	    
	    ComputedProj(type_elt) = type_elt;
	    stored_proj(type_elt) = elt_fine(i)->GetNewNodalInterpolation();
	    elt_fine(i)->ComputeLocalProlongation(*stored_proj(type_elt), ProlongationMatrix(type_elt),
						  *elt_coarse(i), *elt_fine(i)); 
	  }
      }
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_PROJECTION_OPERATOR_CXX
#endif

  
