#ifndef MONTJOIE_FILE_DISCONTINUOUS_MASS_MATRIX_CXX

namespace Montjoie
{
  
  //! destructor
  DiscontinuousBlockMassMatrix::~DiscontinuousBlockMassMatrix()
  {
  }
  
  
  /***********************************
   * DiscontinuousDiagonalMassMatrix *
   ***********************************/
  
  
  //! constructor with the diagonal
  DiscontinuousDiagonalMassMatrix::DiscontinuousDiagonalMassMatrix(const VectReal_wp& d)
  {
    SetDiagonal(d);
  }

  
  //! replaces y by M y
  void DiscontinuousDiagonalMassMatrix::MltMass(VectReal_wp& y) const
  {
    if (diagonal.GetM() <= 0)
      {
	cout << "Diagonal not present" << endl;
	abort();
      }
    
    for (int i = 0; i < diagonal.GetM(); i++)
      y(i) *= diagonal(i);    
  }
  
  
  //! replaces y by M^{-1} y
  void DiscontinuousDiagonalMassMatrix::SolveMass(VectReal_wp& y) const
  {
    if (diagonal.GetM() <= 0)
      {
	cout << "Diagonal not present" << endl;
	abort();
      }
    
    for (int i = 0; i < diagonal.GetM(); i++)
      y(i) /= diagonal(i);
  }
  

  //! replaces y by L^-1 y (or L^-T y)
  void DiscontinuousDiagonalMassMatrix
  ::SolveCholesky(const SeldonTranspose& transA, VectReal_wp& y) const
  {
    if (diagonal.GetM() <= 0)
      {
	cout << "Diagonal not present" << endl;
	abort();
      }
    
    for (int i = 0; i < diagonal.GetM(); i++)
      y(i) /= sqrt(diagonal(i));
  }


  //! inits the diagonal of the mass matrix
  void DiscontinuousDiagonalMassMatrix::SetDiagonal(const VectReal_wp& d)
  {
    diagonal  = d;
  }


  //! returns the size used by the object in bytes
  size_t DiscontinuousDiagonalMassMatrix::GetMemorySize() const
  {
    return diagonal.GetMemorySize();
  }
  

  /************************************
   * DiscontinuousBlockDiagMassMatrix *
   ************************************/
  
  
  //! constructor with the block-diagonal
  template<int d>
  DiscontinuousBlockDiagMassMatrix<d>
  ::DiscontinuousBlockDiagMassMatrix(const Vector<int>& num, const Vector<TinyMatrix<Real_wp, Symmetric, d, d> >& B)
  {
    SetBlockDiagonal(num, B);
  }

  
  //! replaces y by M y
  template<int d>
  void DiscontinuousBlockDiagMassMatrix<d>::MltMass(VectReal_wp& Y) const
  {
    if (row_num.GetM() <= 0)
      {
	cout << "Block-diagonal not present" << endl;
	abort();
      }

    TinyVector<Real_wp, d> vec_u;
    int offset = 0;
    for (int i = 0; i < blocks.GetM(); i++)
      {
	ExtractVector(Y, row_num, offset, vec_u);
	MltCholesky(SeldonTrans, blocks(i), vec_u);
	MltCholesky(SeldonNoTrans, blocks(i), vec_u);
	ExtractVector(vec_u, row_num, offset, Y);
	offset += d;
      }
  }
  
  
  //! replaces y by M^{-1} y
  template<int d>
  void DiscontinuousBlockDiagMassMatrix<d>::SolveMass(VectReal_wp& Y) const
  {
    if (row_num.GetM() <= 0)
      {
	cout << "Block-diagonal not present" << endl;
	abort();
      }

    TinyVector<Real_wp, d> vec_u;
    int offset = 0;
    for (int i = 0; i < blocks.GetM(); i++)
      {
	ExtractVector(Y, row_num, offset, vec_u);
	Seldon::SolveCholesky(SeldonNoTrans, blocks(i), vec_u);
	Seldon::SolveCholesky(SeldonTrans, blocks(i), vec_u);
	ExtractVector(vec_u, row_num, offset, Y);
	offset += d;
      }
  }
  

  //! replaces y by L^-1 y (or L^-T y)
  template<int d>
  void DiscontinuousBlockDiagMassMatrix<d>
  ::SolveCholesky(const SeldonTranspose& transA, VectReal_wp& Y) const
  {
    if (row_num.GetM() <= 0)
      {
	cout << "Block-diagonal not present" << endl;
	abort();
      }

    TinyVector<Real_wp, d> vec_u;
    int offset = 0;
    if (transA.Trans())
      for (int i = 0; i < blocks.GetM(); i++)
	{
	  ExtractVector(Y, row_num, offset, vec_u);
	  Seldon::SolveCholesky(SeldonTrans, blocks(i), vec_u);
	  ExtractVector(vec_u, row_num, offset, Y);
	  offset += d;
	}
    else      
      for (int i = 0; i < blocks.GetM(); i++)
	{
	  ExtractVector(Y, row_num, offset, vec_u);
	  Seldon::SolveCholesky(SeldonNoTrans, blocks(i), vec_u);
	  ExtractVector(vec_u, row_num, offset, Y);
	  offset += d;
	}
  }


  //! inits the diagonal of the mass matrix
  template<int d>
  void DiscontinuousBlockDiagMassMatrix<d>
  ::SetBlockDiagonal(const Vector<int>& num, const Vector<TinyMatrix<Real_wp, Symmetric, d, d> >& B)
  {
    row_num = num;
    blocks = B;
    for (int i = 0; i < blocks.GetM(); i++)
      GetCholesky(blocks(i));
  }


  //! returns the size used by the object in bytes
  template<int d>
  size_t DiscontinuousBlockDiagMassMatrix<d>::GetMemorySize() const
  {
    return blocks.GetMemorySize() + row_num.GetMemorySize();
  }

  
  /*********************************
   * DiscontinuousSparseMassMatrix *
   *********************************/


  //! constructor with a dense matrix
  DiscontinuousSparseMassMatrix
  ::DiscontinuousSparseMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked>& A, const Real_wp& eps)
  {
    ConvertDense(A, eps);
  }

  
  //! replaces y by M y
  void DiscontinuousSparseMassMatrix::MltMass(VectReal_wp& y) const
  {
    if (sparse_Mh.GetM() <= 0)
      {
	cout << "Sparse matrix not present" << endl;
	abort();
      }
    
    MltCholesky(SeldonTrans, sparse_Mh, y);
    MltCholesky(SeldonNoTrans, sparse_Mh, y);    
  }


  //! replaces y by M^{-1} y
  void DiscontinuousSparseMassMatrix::SolveMass(VectReal_wp& y) const
  {
    if (sparse_Mh.GetM() < 0)
      {
	cout << "Sparse matrix not present" << endl;
	abort();
      }
    
    Seldon::SolveCholesky(SeldonNoTrans, sparse_Mh, y);
    Seldon::SolveCholesky(SeldonTrans, sparse_Mh, y);
  }
  
  
  //! replaces y by L^-1 y (or L^-T y)
  void DiscontinuousSparseMassMatrix::SolveCholesky(const SeldonTranspose& transA, VectReal_wp& y) const
  {
    if (sparse_Mh.GetM() < 0)
      {
	cout << "Sparse matrix not present" << endl;
	abort();
      }
    
    Seldon::SolveCholesky(transA, sparse_Mh, y);
  }
  

  //! converts a dense matrix into a sparse one
  void DiscontinuousSparseMassMatrix
  ::ConvertDense(Matrix<Real_wp, Symmetric, RowSymPacked>& dense_Mh, const Real_wp& eps)
  {
    // first Cholesky factor is computed
    GetCholesky(dense_Mh);

    // then conversion into a sparse matrix
    ConvertToSparse(dense_Mh, sparse_Mh, eps);
  }
  

  //! returns the size used by the object in bytes
  size_t DiscontinuousSparseMassMatrix::GetMemorySize() const
  {
    return sparse_Mh.GetMemorySize();
  }
  

  /********************************
   * DiscontinuousDenseMassMatrix *
   ********************************/


  //! constructor with a dense matrix
  DiscontinuousDenseMassMatrix::DiscontinuousDenseMassMatrix(const Matrix<Real_wp, Symmetric, RowSymPacked>& A)
  {
    SetDense(A);
  }

  
  //! replaces y by M y
  void DiscontinuousDenseMassMatrix::MltMass(VectReal_wp& y) const
  {
    if (dense_Mh.GetM() <= 0)
      {
	cout << "Dense matrix not present" << endl;
	abort();
      }
    
    MltCholesky(SeldonTrans, dense_Mh, y);
    MltCholesky(SeldonNoTrans, dense_Mh, y);    
  }
  

  //! replaces y by M^{-1} y
  void DiscontinuousDenseMassMatrix::SolveMass(VectReal_wp& y) const
  {
    if (dense_Mh.GetM() < 0)
      {
	cout << "Dense matrix not present" << endl;
	abort();
      }
    
    Seldon::SolveCholesky(SeldonNoTrans, dense_Mh, y);
    Seldon::SolveCholesky(SeldonTrans, dense_Mh, y);
  }


  //! replaces y by L^-1 y (or L^-T y)
  void DiscontinuousDenseMassMatrix::SolveCholesky(const SeldonTranspose& transA, VectReal_wp& y) const
  {
    if (dense_Mh.GetM() < 0)
      {
	cout << "Dense matrix not present" << endl;
	abort();
      }
    
    Seldon::SolveCholesky(transA, dense_Mh, y);    
  }
  

  //! sets the dense matrix
  void DiscontinuousDenseMassMatrix
  ::SetDense(const Matrix<Real_wp, Symmetric, RowSymPacked>& M)
  {
    dense_Mh = M;
    
    // the mass matrix is factorized such that both Solve and Mlt are fast
    GetCholesky(dense_Mh);
  }
  

  //! returns the size used by the object in bytes
  size_t DiscontinuousDenseMassMatrix::GetMemorySize() const
  {
    return dense_Mh.GetMemorySize();
  }
  

  /*************************************
   * DiscontinuousMatrixFreeMassMatrix *
   *************************************/


  //! constructor with finite element class
  DiscontinuousMatrixFreeMassMatrix::DiscontinuousMatrixFreeMassMatrix(ElementReference_Base& elt)
    : Fb(elt)
  {
  }

  
  //! replaces y by M y
  void DiscontinuousMatrixFreeMassMatrix::MltMass(VectReal_wp& y) const
  {
    int Nquad = Fb.GetNbPointsQuadratureInside();
    xn_tmp.Reallocate(Nquad); xn_tmp.Zero();
    
    Fb.ApplyChTranspose(y, xn_tmp);
    
    for (int i = 0; i < Nquad; i++)
      xn_tmp(i) *= diagonal(i);
    
    Fb.ApplyCh(xn_tmp, y);
  }
  
  
  //! replaces y by M^{-1} y
  void DiscontinuousMatrixFreeMassMatrix::SolveMass(VectReal_wp& y) const
  {
    // conjugate gradient to solve the system
    Iteration<Real_wp> iter(1000, stopping_criterion);

    iter.HideMessages();

    VectReal_wp xn_cg(y.GetM());
    xn_cg.Zero();

    DiscontinuousMatrixFreeMassMatrix& prec = const_cast<DiscontinuousMatrixFreeMassMatrix& >(*this);
    Cg(*this, y, xn_cg, prec, iter);
  }
  

  //! replaces y by L^-1 y (or L^-T y)
  void DiscontinuousMatrixFreeMassMatrix::SolveCholesky(const SeldonTranspose& transA, VectReal_wp& y) const
  {
    cout << "Cholesky not possible with matrix-free storage" << endl;
    abort();
  }
          

  //! sets the stopping criterion used by the iterative solver
  void DiscontinuousMatrixFreeMassMatrix::SetStoppingCriterion(const Real_wp& eps)
  {
    stopping_criterion = eps;
  }


  //! returns the stopping criterion used by the iterative solver
  Real_wp DiscontinuousMatrixFreeMassMatrix::GetStoppingCriterion() const
  {
    return stopping_criterion;
  }


  //! Application of the preconditioning for matrix-free storage
  void DiscontinuousMatrixFreeMassMatrix::
  Solve(const VirtualMatrix<Real_wp>& A, const VectReal_wp& r, VectReal_wp& z)
  {
    if (diagonal.GetM() == 0)
      Copy(r, z);
    else
      for (int i = 0; i < r.GetM(); i++)
	z(i) = r(i)*diagonal(i);
  }
  

  //! returns the size used by the object in bytes
  size_t DiscontinuousMatrixFreeMassMatrix::GetMemorySize() const
  {
    size_t taille = xn_tmp.GetMemorySize() + diagonal.GetMemorySize();
    return taille;
  }
  

  /**********************************
   * DiscontinuousMassMatrixProblem *
   **********************************/


  //! retrieves physical coefficient 
  void DiscontinuousMassMatrixProblem::GetMassDampingCoefficient(int n, int k, int nu, int ref,
								 Real_wp& rho, Real_wp& sigma)
  {
    cout << "Method not overloaded" << endl;
    abort();
  }


  /******************************
   * DiscontinuousMassMatrixVol *
   ******************************/
  

  //! constructor
  template<class Dimension>
  DiscontinuousMassMatrixVol<Dimension>
  ::DiscontinuousMassMatrixVol(VarProblem<Dimension>& var_p,
			       VarInstationary_Base& var_t)
    : mesh(var_p.mesh), var_problem(var_p), var_time(var_t)
  {    
  }
  

  //! returns the number of blocks of mass matrix
  template<class Dimension>
  int DiscontinuousMassMatrixVol<Dimension>::GetNbElt() const
  {
    return mesh.GetNbElt();
  }


  //! gives informations about block i
  template<class Dimension>
  void DiscontinuousMassMatrixVol<Dimension>
  ::GetElementInfo(int i, int num, bool geom_mass,
		   int& ref, bool& affine, bool& varying) const
  {
    ref = mesh.Element(i).GetReference();
    varying = false;
    if (!geom_mass)
      varying = var_problem.IsVaryingMedia(num, ref);
    
    affine = mesh.IsElementAffine(i);      
  }
  

  //! retrieves jacobian (affine and curved)
  template<class Dimension>
  void DiscontinuousMassMatrixVol<Dimension>
  ::GetJacobianPanel(int n, bool& affine, bool& curved,
		     Real_wp& jacob_affine, VectReal_wp& decomp_jacob,
		     VectReal_wp& jacob_quad)
  {
    affine = mesh.IsElementAffine(n);
    curved = mesh.Element(n).IsCurved();
    if (affine)
      jacob_affine = var_problem.Glob_jacobian(n)(0);
    else
      {
	decomp_jacob = var_problem.Glob_decomp_jacobian(n);
	jacob_quad = var_problem.Glob_jacobian(n);
      }
  }
  

  //! retrieves physical coefficient 
  template<class Dimension>
  void DiscontinuousMassMatrixVol<Dimension>
  ::GetMassDampingCoefficient(int n, int k, int nu, int ref,
			      Real_wp& rho, Real_wp& sigma)
  {
    var_time.GetMassDampingCoefficient(n, k, nu, ref, rho, sigma);
  }
  

  //! returns jacobian for an affine element
  template<class Dimension>
  Real_wp DiscontinuousMassMatrixVol<Dimension>::GetAffineJacobian(int n) const
  {
    return var_problem.Glob_jacobian(n)(0);
  }


  //! returns reference element of element n
  template<class Dimension>
  const ElementReference_Base& DiscontinuousMassMatrixVol<Dimension>::GetReferenceElement(int n) const
  {
    return var_problem.GetReferenceElementBase(n);
  }
  
  
  /*******************************
   * DiscontinuousMassMatrixSurf *
   *******************************/
  

  //! constructor
  template<class Dimension>
  DiscontinuousMassMatrixSurf<Dimension>
  ::DiscontinuousMassMatrixSurf(VarProblem<Dimension>& var_p,
				VarInstationary_Base& var_t,
				VectReal_wp& val_ds, IVect& offset_ds, bool geom)
    : mesh(var_p.mesh), var_problem(var_p), var_time(var_t),
      ValueDsj(val_ds), OffsetDsj(offset_ds)
  {
    if (geom)
      ConstructDs();
    else
      ConstructCoef();
  }


  template<class Dimension>
  void DiscontinuousMassMatrixSurf<Dimension>::ConstructDs()
  {
    Vector<bool> affine_face(mesh.GetNbBoundary());
    affine_face.Fill(false);
    int nb = 0;
    OffsetDsj.Reallocate(affine_face.GetM()+1);
    OffsetDsj(0) = 0;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	if (var_problem.IsNewFace(i)(num_loc))
	  {
	    int num_face = var_problem.mesh.Element(i).numBoundary(num_loc);
	    Real_wp ds0 = var_problem.Glob_dsj(num_face)(0);
	    Real_wp eps = 1e-12*ds0;
	    bool affine = true;
	    for (int k = 1; k < var_problem.Glob_dsj(num_face).GetM(); k++)
	      if (abs(ds0 - var_problem.Glob_dsj(num_face)(k)) > eps)
		affine = false;
	    
            affine_face(num_face) = affine;
	    int nb_points = 1;
	    if (!affine)
	      nb_points = var_problem.Glob_dsj(num_face).GetM();
	    
	    OffsetDsj(num_face+1) = nb_points;
	    nb += nb_points;
	  }
    
    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      OffsetDsj(i+1) += OffsetDsj(i);
    
    ValueDsj.Reallocate(nb);
    const MeshNumbering<Dimension>& mesh_num = var_problem.GetMeshNumbering(0);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	if (var_problem.IsNewFace(i)(num_loc))
	  {
	    int num_face = var_problem.mesh.Element(i).numBoundary(num_loc);
	    int rf = mesh_num.GetOrderQuadrature(num_face);
	    const Matrix<int>& FacesQuadRotation = mesh_num.number_map.
	      GetRotationQuadraturePoints(rf, mesh.Boundary(num_face));
	    
            const VectReal_wp& weight = mesh_num.number_map.
              GetFluxWeight(rf, mesh.Boundary(num_face));
            
	    int rot = mesh.Element(i).GetOrientationBoundary(num_loc);
	    if (affine_face(num_face))
	      ValueDsj(OffsetDsj(num_face)) = var_problem.Glob_dsj(num_face)(0);
	    else
	      for (int k = 0; k < var_problem.Glob_dsj(num_face).GetM(); k++)
		{
		  int krot = FacesQuadRotation(rot, k);
		  ValueDsj(OffsetDsj(num_face)+k) = 2*var_problem.Glob_dsj(num_face)(krot)
                    * weight(k);                  
		}
	  }
  }


  template<class Dimension>
  void DiscontinuousMassMatrixSurf<Dimension>::ConstructCoef()
  {
    coef_surface.Reallocate(mesh.GetNbBoundary());
    coef_surface.Zero();
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	//int ref = mesh.Element(i).GetReference();
	for (int num_loc = 0; num_loc < mesh.Element(i).GetNbBoundary(); num_loc++)
	  {
            int nf = mesh.Element(i).numBoundary(num_loc);
            Real_wp tau = var_time.GetCoefficientTauHDG(nf);
            coef_surface(mesh.Element(i).numBoundary(num_loc)) += tau;
          }
      }

    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int ref = mesh.BoundaryRef(i).GetReference();
	if (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_ABSORBING)
	  {
	    Real_wp tau = var_time.GetCoefficientTauHDG(i);
	    coef_surface(i) += tau;
	  }
	else if ( (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_DIRICHLET)
		  || (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEUMANN)
		  || (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_INSIDE)
		  || (mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEIGHBOR) )
	  {
	    // nothing to add
	  }
	else
	  {
	    cout << "Other boundary conditions not implemented for HDG" << endl;
	    abort();
	  }
      }

#ifdef SELDON_WITH_MPI
    var_problem.ReduceDistributedVectorFace(coef_surface, MPI_SUM, 1);
#endif
    
  }

  
  //! returns the number of blocks of mass matrix
  template<class Dimension>
  int DiscontinuousMassMatrixSurf<Dimension>::GetNbElt() const
  {
    return OffsetDsj.GetM()-1;
  }


  //! gives informations about block i
  template<class Dimension>
  void DiscontinuousMassMatrixSurf<Dimension>
  ::GetElementInfo(int i, int num, bool geom_mass,
		   int& ref, bool& affine, bool& varying) const
  {
    ref = -1;
    varying = false;
    affine = false;
    if (OffsetDsj(i+1) == OffsetDsj(i)+1)
      affine = true;
  }
  

  //! retrieves jacobian (affine and curved)
  template<class Dimension>
  void DiscontinuousMassMatrixSurf<Dimension>
  ::GetJacobianPanel(int n, bool& affine, bool& curved,
		     Real_wp& jacob_affine, VectReal_wp& decomp_jacob,
		     VectReal_wp& jacob_quad)
  {
    affine = false;
    if (OffsetDsj(n+1) == OffsetDsj(n)+1)
      affine = true;

    curved = false;
    if (!affine)
      curved = mesh.IsBoundaryCurved(n);
    
    if (affine)
      jacob_affine = ValueDsj(OffsetDsj(n));
    else
      {
	decomp_jacob.Clear();
	int nb_val = OffsetDsj(n+1) - OffsetDsj(n);
	jacob_quad.Reallocate(nb_val);
	for (int j = 0; j < nb_val; j++)
	  jacob_quad(j) = ValueDsj(OffsetDsj(n) + j);
      }
  }
  

  //! retrieves physical coefficient 
  template<class Dimension>
  void DiscontinuousMassMatrixSurf<Dimension>
  ::GetMassDampingCoefficient(int n, int k, int nu, int ref,
			      Real_wp& rho, Real_wp& sigma)
  {
    sigma = 0.0;
    rho = coef_surface(n);
  }
  

  //! returns jacobian for an affine element
  template<class Dimension>
  Real_wp DiscontinuousMassMatrixSurf<Dimension>::GetAffineJacobian(int n) const
  {
    return ValueDsj(OffsetDsj(n));
  }


  //! returns reference element of element n
  template<class Dimension>
  const ElementReference_Base& DiscontinuousMassMatrixSurf<Dimension>::GetReferenceElement(int n) const
  {
    return var_problem.GetSurfaceElementBase(n);
  }


  //! returns the size used by the object in bytes
  template<class Dimension>
  size_t DiscontinuousMassMatrixSurf<Dimension>::GetMemorySize() const
  {
    size_t taille = coef_surface.GetMemorySize();
    return taille;
  }
  
  
  /***************************
   * DiscontinuousMassMatrix *
   ***************************/
  
  
  //! default constructor
  DiscontinuousMassMatrix_Base
  ::DiscontinuousMassMatrix_Base(VarInstationary_Base& var_t,
				 DiscontinuousMassMatrix_Base* Mgeom,
				 DiscontinuousMassMatrixProblem& var_p,
				 int type_mat)
    : GeomMass(Mgeom), level_time(), type_matrix(type_mat), var_problem(var_p)
  {
    threshold = 10*epsilon_machine;
    stopping_criterion = 1e-12;
    num_unknown = -1;
    coef_mass_ = 0;
    coef_damp_ = 0;
  }


  //! destructor
  DiscontinuousMassMatrix_Base::~DiscontinuousMassMatrix_Base()
  {
    ClearBlocks();
  }

  
  //! stored blocks are released
  void DiscontinuousMassMatrix_Base::ClearBlocks()
  {
    for (int i = 0; i < block_mass_store.GetM(); i++)
      if (block_mass_store(i) != NULL)
	{
	  delete block_mass_store(i);
	  block_mass_store(i) = NULL;
	}

    block_mass_store.Clear();
  }


  //! returns the size used by the object in bytes
  size_t DiscontinuousMassMatrix_Base::GetMemorySize() const
  {
    size_t taille = diagonal.GetMemorySize() + invDiagonal.GetMemorySize()
      + block_diagonal.GetMemorySize() + block_num_element.GetMemorySize()
      + block_size.GetMemorySize() + Seldon::GetMemorySize(coef_block_mass)
      + sizeof(void*)*block_mass_store.GetM();

    for (int i = 0; i < block_mass_store.GetM(); i++)
      if (block_mass_store(i) != NULL)
	taille += block_mass_store(i)->GetMemorySize();

    return taille;
  }
  

  //! multiplication by mass matrix, Vh(offset:) = M Uh(offset:)
  void DiscontinuousMassMatrix_Base
  ::MltMass(int offset0, const VectReal_wp& Uh, VectReal_wp& Vh)
  {  
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (diagonal.GetM() <= 0)
	    {
	      cout << "Diagonal matrix not present" << endl;
	      abort();
	    }
	  
	  for (int i = 0; i < diagonal.GetM(); i++)
	    Vh(offset0+i) = Uh(offset0+i)*diagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
	      abort();
	    }

	  cout << "not implemented" << endl;
          abort();
	}
	break;
      case FemMassMatrix::BLOCK_ELT :
        {
          if (block_mass.GetM() <= 0)
            {
	      cout << "Blocks not present" << endl;
              abort();
            }
          
          VectReal_wp Uloc;
	  int nb_dof_loc(0);
	  for (int i0 = 0; i0 < level_time.GetNbElt(); i0++)
            {
	      int i = level_time.GetElementNumber(i0);
	      int n = block_num_element(i);	      
	      nb_dof_loc = block_size(i+1) - block_size(i);
	      int offset = offset0 + block_size(i);
	      MltLocalMass(i, n, offset, Uh, nb_dof_loc, Uloc);
              for (int j = 0; j < Uloc.GetM(); j++)
                Vh(offset + j) = Uloc(j);
	    }
        }
        break;
      default :
	cout << "Storage not handled" << endl;
	abort();
	break;
      }
  }
  
  
  //! multiplication by mass matrix, Vh(offset:) = Vh(offset:) + alpha M Uh
  void DiscontinuousMassMatrix_Base
  ::MltMass(int offset0, const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (diagonal.GetM() <= 0)
	    {
	      cout << "Diagonal matrix not present" << endl;
	      abort();
	    }
	  
          for (int i = 0; i < diagonal.GetM(); i++)
            Vh(offset0+i) += alpha*Uh(offset0+i)*diagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
	      abort();
	    }

	  cout << "not implemented" << endl;
          abort();
	}
	break;
      case FemMassMatrix::BLOCK_ELT :
        {          
          VectReal_wp Uloc;
	  int nb_dof_loc(0);
	  for (int i0 = 0; i0 < level_time.GetNbElt(); i0++)
            {
	      int i = level_time.GetElementNumber(i0);
	      int n = block_num_element(i);
	      nb_dof_loc = block_size(i+1) - block_size(i);
	      int offset = offset0 + block_size(i);
	      MltLocalMass(i, n, offset, Uh, nb_dof_loc, Uloc);
              for (int j = 0; j < Uloc.GetM(); j++)
                Vh(offset + j) += alpha*Uloc(j); 
            }
        }
	break;
      default :
	cout << "Storage not handled" << endl;
	abort();
	break;
      }
  }
  
  
  //! solves mass matrix, Vh(offset:) = M^{-1} Vh(offset:)
  void DiscontinuousMassMatrix_Base
  ::SolveMass(int offset0, VectReal_wp& Vh)
  {
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invDiagonal.GetM() <= 0)
	    {
	      cout << "Diagonal matrix not present" << endl;
	      abort();
	    }
	  
	  for (int i = 0; i < invDiagonal.GetM(); i++)
	    Vh(offset0 + i) *= invDiagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
	      abort();
	    }

	  cout << "not implemented" << endl;
          abort();
	}
	break;
      case FemMassMatrix::BLOCK_ELT :
	{
          VectReal_wp Uloc;
	  int nb_dof_loc(0);
	  for (int i0 = 0; i0 < level_time.GetNbElt(); i0++)
            {
	      int i = level_time.GetElementNumber(i0);
	      int n = block_num_element(i);
	      nb_dof_loc = block_size(i+1) - block_size(i);
	      int offset = offset0 + block_size(i);
	      SolveLocalMass(i, n, offset, Vh, nb_dof_loc, Uloc);
	      for (int j = 0; j < Uloc.GetM(); j++)
		Vh(offset + j) = Uloc(j);                   
            }
	}
	break;
      default :
	{
	  cout << "Storage not implemented" << endl;
          abort();
        }
	break;
      }
  }

  
  //! solves mass matrix, Vh(offset:) = Vh(offset:) + alpha M^{-1} Uh
  void DiscontinuousMassMatrix_Base
  ::SolveMass(int offset0, const Real_wp& alpha, const VectReal_wp& Uh, VectReal_wp& Vh)
  {  
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invDiagonal.GetM() <= 0)
	    {
	      cout << "Diagonal matrix not present" << endl;
	      abort();
	    }
	  
	  for (int i = 0; i < invDiagonal.GetM(); i++)
	    Vh(offset0+i) += alpha*Uh(offset0+i)*invDiagonal(i);
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
	      abort();
	    }

	  cout << "Storage not implemented" << endl;
          abort();
	}
	break;
      case FemMassMatrix::BLOCK_ELT :
	{
          VectReal_wp Uloc; int nb_dof_loc(0);
	  for (int i0 = 0; i0 < level_time.GetNbElt(); i0++)
            {
	      int i = level_time.GetElementNumber(i0);
	      int n = block_num_element(i);
	      nb_dof_loc = block_size(i+1) - block_size(i);
	      int offset = offset0 + block_size(i);
	      SolveLocalMass(i, n, offset, Uh, nb_dof_loc, Uloc);
	      
	      for (int j = 0; j < Uloc.GetM(); j++)
		Vh(offset + j) += alpha*Uloc(j); 
            }
	}
	break;
      default :
	cout << "Storage not implemented" << endl;
	abort();
	break;
      }
  }

  
  //! solves by a cholesky factor of mass matrix M = L L^T, Vh = L^{-1} Vh or Vh = L^{-T} Vh  
  void DiscontinuousMassMatrix_Base
  ::SolveCholesky(int offset, const SeldonTranspose& transA, VectReal_wp& Vh)
  {  
    switch (type_matrix)
      {
      case FemMassMatrix::DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (invDiagonal.GetM() <= 0)
	    {
	      cout << "Diagonal matrix not present" << endl;
	      abort();
	    }
	  
	  for (int i = 0; i < invDiagonal.GetM(); i++)
	    Vh(offset + i) *= sqrt(invDiagonal(i));
	}
	break;
      case FemMassMatrix::BLOCK_DIAGONAL :
	{
	  // exiting if the matrix has not been computed
	  if (block_diagonal.GetM() <= 0)
	    {
	      cout << "Block-diagonal matrix not present" << endl;
	      abort();
	    }

	  cout << "Storage not implemented" << endl;
          abort();
	}
	break;
      case FemMassMatrix::BLOCK_ELT :
	{
          VectReal_wp Uloc;
          for (int i = 0; i < level_time.GetNbElt(); i++)
            {
	      SolveLocalCholesky(i, offset, transA, Vh, Uloc);
              
              for (int j = 0; j < Uloc.GetM(); j++)
                Vh(offset + j) = Uloc(j); 
              
              offset += Uloc.GetM();
            }
	}
	break;
      default :
	{
	  cout << "Storage not implemented" << endl;
          abort();
        }
	break;
      }
  }


  //! computation of mass matrix for an element of the mesh
  /*!
    \param[in] num_elem element number
    \param[in] ref reference of the element
    \param[in] varying if true, the block should be computed numerically
    \param[in] nb
   */
  void DiscontinuousMassMatrix_Base
  ::ComputeLocalMass(int num_elem, int ref, bool varying, int& nb,
		     bool geom_mass, const ElementReference_Base& Fb)
  {
    int nb_dof_loc = Fb.GetNbDof();
    bool affine, curved;
    Real_wp Ji; VectReal_wp decomp_Ji, val_Ji;
    var_problem.GetJacobianPanel(num_elem, affine, curved, Ji, decomp_Ji, val_Ji);
    
    int Nquad = Fb.GetNbPointsQuadratureInside();

    if ((!affine) && (Fb.use_warburton_trick))
      {
        // nothing to do here, we use reference mass matrix
      }
    else if (!varying)
      {
        if (geom_mass)
          {
            // we compute this mass matrix only for geometry mass 
            if (affine)
              {
                // we don't compute the mass matrix
                // we will use the reference mass matrix
		this->block_mass(this->num_unknown)(num_elem) = NULL;
              }
            else if (Fb.LinearSparseMassMatrix() && !curved)
              {
                if (nb >= this->block_mass_store.GetM())
                  {
                    abort();
                  }
		
		// the block is a sparse matrix
		DiscontinuousSparseMassMatrix* Mh;
		
		// Mh is computed as a dense matrix then converted to a sparse one
		Matrix<Real_wp, Symmetric, RowSymPacked> dense_Mh;
                Fb.ComputeMassMatrix(dense_Mh, decomp_Ji);

		Mh = new DiscontinuousSparseMassMatrix(dense_Mh, threshold);
                this->block_mass_store(nb) = Mh;		
		this->block_mass(this->num_unknown)(num_elem) = Mh;
                nb++;
              }
            else
              {    
                if (nb >= this->block_mass_store.GetM())
                  {
                    abort();
                  }

		// the block is either a sparse or dense matrix
                DiscontinuousBlockMassMatrix* Mh;
                
                if (Fb.LumpedMassMatrix())
                  {
		    // diagonal mass matrix
		    VectReal_wp diagonal(val_Ji);
		    Mh = new DiscontinuousDiagonalMassMatrix(diagonal);
                  }
                else
		  {
		    // Mh is computed as a dense matrix
		    Matrix<Real_wp, Symmetric, RowSymPacked> dense_Mh;
		    if (Fb.OptimizedComputationMassMatrix())
		      Fb.IntegrateMassMatrix(dense_Mh, val_Ji);
		    else
		      {
			dense_Mh.Reallocate(nb_dof_loc, nb_dof_loc);
			dense_Mh.Zero();
			VectReal_wp contrib(nb_dof_loc), feval(Nquad);
			for (int i = 0; i < nb_dof_loc; i++)
			  {
			    contrib.Zero(); contrib(i) = 1.0;
			    Fb.ApplyChTranspose(contrib, feval);
			    
			    for (int k = 0; k < Nquad; k++)
			      feval(k) *= val_Ji(k);
			    
			    Fb.ApplyCh(feval, contrib);
			    
			    for (int j = 0; j < nb_dof_loc; j++)
			      dense_Mh(i, j) = contrib(j);
			  }
		      }

		    // then it might be converted into a sparse matrix
		    if (Fb.SparseMassMatrix())
                      Mh = new DiscontinuousSparseMassMatrix(dense_Mh, threshold);
		    else
		      Mh = new DiscontinuousDenseMassMatrix(dense_Mh);
		  }
		
		this->block_mass_store(nb) = Mh;
		this->block_mass(this->num_unknown)(num_elem) = Mh;
                nb++;
              }
          }
        else
          {
            // nothing to do, we will use the geometric mass matrix
	    this->block_mass(this->num_unknown)(num_elem) = GeomMass->block_mass(0)(num_elem);
          }
      }
    else
      {
	// case where the mass matrix involves a variable physical coefficient (like rho)
        if (nb >= this->block_mass_store.GetM() )
          {
            abort();
          }
                
        DiscontinuousBlockMassMatrix* Mh = NULL;
        const VectReal_wp& weights = Fb.WeightsND();
	  
        if (Fb.LumpedMassMatrix())
          {
	    // diagonal mass matrix
            VectReal_wp diagonal(nb_dof_loc);
	    Real_wp rho, sigma;
            for (int k = 0; k < Nquad; k++)
              {
		var_problem.GetMassDampingCoefficient(num_elem, k, this->num_unknown, ref, rho, sigma);
		Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;
		
                if (affine)
                  coef *= Ji*weights(k);
                else
                  coef *= val_Ji(k);
                
                diagonal(k) = coef;
              }
	    
	    Mh = new DiscontinuousDiagonalMassMatrix(diagonal);
          }
        else
	  {
	    Matrix<Real_wp, Symmetric, RowSymPacked> dense_Mh;

	    VectReal_wp coef(Nquad);
	    Real_wp rho, sigma;
            for (int k = 0; k < Nquad; k++)
              {
		var_problem.GetMassDampingCoefficient(num_elem, k, this->num_unknown, ref, rho, sigma);		
                coef(k) = rho*this->coef_mass_ + sigma*this->coef_damp_;
		
                if (affine)
                  coef(k) *= Ji*weights(k);
                else
                  coef(k) *= val_Ji(k);
              }
	    
	    if (Fb.OptimizedComputationMassMatrix())
	      Fb.IntegrateMassMatrix(dense_Mh, coef);
	    else
	      {
		dense_Mh.Reallocate(nb_dof_loc, nb_dof_loc);
		VectReal_wp contrib(nb_dof_loc), feval(Nquad);
		for (int i = 0; i < nb_dof_loc; i++)
		  {
		    contrib.Zero(); contrib(i) = 1.0;
		    Fb.ApplyChTranspose(contrib, feval);
		    
		    for (int k = 0; k < Nquad; k++)
		      feval(k) *= coef(k);
		    
		    Fb.ApplyCh(feval, contrib);
		    
		    for (int j = 0; j < nb_dof_loc; j++)
		      dense_Mh(i, j) = contrib(j);
		  }
	      }

	    if (Fb.SparseMassMatrix())
	      Mh = new DiscontinuousSparseMassMatrix(dense_Mh, threshold);
	    else
	      Mh = new DiscontinuousDenseMassMatrix(dense_Mh);	    
	  }
	
	this->block_mass_store(nb) = Mh;
	this->block_mass(this->num_unknown)(num_elem) = Mh;
        nb++;
      }
  }


  //! allocates blocks for each element (for the stored mass)
  void DiscontinuousMassMatrix_Base::ReallocateStoredMass()
  {
    // previous blocks are cleared
    ClearBlocks();

    int nb_elt = var_problem.GetNbElt();
    DiscontinuousBlockMassMatrix* null_ptr = NULL;
    this->block_mass_store.Reallocate(nb_elt);
    this->block_mass_store.Fill(null_ptr);
    
    this->block_num_element.Reallocate(nb_elt);
    this->block_num_element.Fill();
    this->block_size.Reallocate(nb_elt+1);
    this->level_time.SetNbElt(nb_elt, 0);
    
    this->block_size(0) = 0;
    for (int i = 0; i < nb_elt; i++)
      {
	const ElementReference_Base& Fb = var_problem.GetReferenceElement(i);
	this->block_size(i+1) = block_size(i) + Fb.GetNbDof();
      }
  }
  

  //! modifies a single block of the stored mass
  void DiscontinuousMassMatrix_Base
  ::SetLocalStoredBlock(int i, DiscontinuousBlockMassMatrix* b)
  {
    if (this->block_mass_store(i) != NULL)
      delete this->block_mass_store(i);
    
    this->block_mass_store(i) = b;
  }


  //! shallow copy of blocks stored in GeomMass, coefficients are set
  void DiscontinuousMassMatrix_Base
  ::CopyBlockCoef(const VectReal_wp& coef, int num_u)
  {
    this->num_unknown = num_u;
    this->block_num_element = GeomMass->block_num_element;
    this->block_size = GeomMass->block_size;
    this->level_time = GeomMass->level_time;

    this->block_mass.Reallocate(num_u+1);
    this->block_mass(num_u) = GeomMass->block_mass_store;
    
    coef_block_mass.Reallocate(num_u+1);
    coef_block_mass(num_u) = coef;
  }

  
  //! computing mass matrix for unknowns m with m1 <= m < m2
  void DiscontinuousMassMatrix_Base
  ::ComputeMass(bool geom_mass, int m1, int m2)
  {
    this->type_matrix = FemMassMatrix::BLOCK_ELT;
    
    // first step, counting the number of matrices to compute
    int nb = 0;
    int ref; bool varying, affine;
    int nb_elt = var_problem.GetNbElt();
    for (int i = 0; i < nb_elt; i++)
      {
	const ElementReference_Base& Fb = var_problem.GetReferenceElement(i);
	
	var_problem.GetElementInfo(i, this->num_unknown, geom_mass,
				   ref, affine, varying);
	
        if ((!affine) && (Fb.use_warburton_trick))
          {
            // no matrix to store
          }
        else if (varying)
          {
            // we have to store here a mass matrix due to variation of physical media
            nb += m2-m1;
          }
        else
          {
            if (geom_mass)
              if (!affine)
                nb++;
          }               
      }
    
    // allocating matrices and computing them
    ClearBlocks();
    if (nb > 0)
      {
	DiscontinuousBlockMassMatrix* null_ptr = NULL;
	this->block_mass_store.Reallocate(nb);
	this->block_mass_store.Fill(null_ptr);
      }
    
    nb = 0;
    this->block_num_element.Reallocate(nb_elt);
    this->block_num_element.Fill();
    this->block_size.Reallocate(nb_elt+1);
    this->block_mass.Reallocate(m2);
    this->level_time.SetNbElt(nb_elt, 0);
    for (int m = m1; m < m2; m++)
      {
	this->block_mass(m).Reallocate(nb_elt);
	for (int k = 0; k < this->block_mass(m).GetM(); k++)
	  this->block_mass(m)(k) = NULL;
      }

    this->block_size(0) = 0;
    for (int m = m1; m < m2; m++)
      for (int i = 0; i < nb_elt; i++)
	{
	  const ElementReference_Base& Fb = var_problem.GetReferenceElement(i);
	  this->block_size(i+1) = block_size(i) + Fb.GetNbDof();
	  
	  this->num_unknown = m;
	  var_problem.GetElementInfo(i, this->num_unknown, geom_mass,
				     ref, affine, varying);

	  ComputeLocalMass(i, ref, varying, nb, geom_mass, Fb);
	}
  }


  //! multiplication by mass matrix for an element of the mesh
  void DiscontinuousMassMatrix_Base
  ::MltLocalMass(int i, int n, int offset, const VectReal_wp& Uh,
		 int& nb_dof_loc, VectReal_wp& Uloc)
  {
    const ElementReference_Base& Fb = var_problem.GetReferenceElement(n);
    nb_dof_loc = Fb.GetNbDof();
    Uloc.Reallocate(nb_dof_loc);
    for (int j = 0; j < Uloc.GetM(); j++)
      Uloc(j) = Uh(offset + j);

    if (this->coef_block_mass.GetM() > 0)
      {
	if (this->block_mass(this->num_unknown)(i) == NULL)
	  Fb.MltMassMatrix(Uloc);
	else
	  this->block_mass(this->num_unknown)(i)->MltMass(Uloc);
	
	Mlt(this->coef_block_mass(this->num_unknown)(i), Uloc);
	
	return;
      }

    int ref; bool affine, varying;
    var_problem.GetElementInfo(n, this->num_unknown, false,
			       ref, affine, varying);

    if ((!affine) && (Fb.use_warburton_trick))
      {
	Real_wp rho, sigma;
	var_problem.GetMassDampingCoefficient(n, 0, this->num_unknown, ref, rho, sigma);
        Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;
        
        if (coef != 0)
          {            
            Fb.MltMassMatrix(Uloc);
            Mlt(coef, Uloc);
          }
        else
          Uloc.Zero();
      }
    else if (varying)
      {
        this->block_mass(this->num_unknown)(i)->MltMass(Uloc);
      }
    else
      {
	Real_wp rho, sigma;
	var_problem.GetMassDampingCoefficient(n, 0, this->num_unknown, ref, rho, sigma);
        Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;

        if (affine)
          {
            Fb.MltMassMatrix(Uloc);
            coef *= var_problem.GetAffineJacobian(n);
          }
        else
	  this->block_mass(this->num_unknown)(i)->MltMass(Uloc);
	
        Mlt(coef, Uloc);
      }
  }
  
  
  //! solving by mass matrix for an element of the mesh
  void DiscontinuousMassMatrix_Base
  ::SolveLocalMass(int i, int n, int offset, const VectReal_wp& Uh,
		   int &nb_dof_loc, VectReal_wp& Uloc)
  {
    const ElementReference_Base& Fb = var_problem.GetReferenceElement(n);
    nb_dof_loc = Fb.GetNbDof();    
    
    Uloc.Reallocate(nb_dof_loc);
    for (int j = 0; j < Uloc.GetM(); j++)
      Uloc(j) = Uh(offset + j);

    if (this->coef_block_mass.GetM() > 0)
      {
	if (this->block_mass(this->num_unknown)(i) == NULL)
	  Fb.SolveMassMatrix(Uloc);
	else
	  this->block_mass(this->num_unknown)(i)->SolveMass(Uloc);

	Real_wp coef = this->coef_block_mass(this->num_unknown)(i);
	if (coef != Real_wp(1))
	  Mlt(Real_wp(1)/coef, Uloc);
	
	return;
      }

    int ref; bool affine, varying;
    var_problem.GetElementInfo(n, this->num_unknown, false,
			       ref, affine, varying);

    if ((!affine) && (Fb.use_warburton_trick))
      {
	Real_wp rho, sigma;
	var_problem.GetMassDampingCoefficient(n, 0, this->num_unknown, ref, rho, sigma);
	Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;

	Fb.SolveMassMatrix(Uloc);
	Mlt(1.0/coef, Uloc);
      }
    else if (varying)
      {
	this->block_mass(this->num_unknown)(i)->SolveMass(Uloc);
      }
    else
      {
	Real_wp rho, sigma;
	var_problem.GetMassDampingCoefficient(n, 0, this->num_unknown, ref, rho, sigma);
	Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;
	if (affine)
	  {
	    Fb.SolveMassMatrix(Uloc);
	    coef *= var_problem.GetAffineJacobian(n);
	  }
	else
	  this->block_mass(this->num_unknown)(i)->SolveMass(Uloc);

	Mlt(1.0/coef, Uloc);
      }
  }
  
  
  //! solves by a cholesky factor of mass matrix M = L L^T
  void DiscontinuousMassMatrix_Base::
  SolveLocalCholesky(int i, int offset, const SeldonTranspose& transA,
                     const VectReal_wp& Uh, VectReal_wp& Uloc)
  {
    const ElementReference_Base& Fb = var_problem.GetReferenceElement(i);
    int nb_dof_loc = Fb.GetNbDof();

    int ref; bool affine, varying;
    var_problem.GetElementInfo(i, this->num_unknown, false,
			       ref, affine, varying);
    
    Uloc.Reallocate(nb_dof_loc);
    for (int j = 0; j < Uloc.GetM(); j++)
      Uloc(j) = Uh(offset + j);
    
    if ((!affine) && (Fb.use_warburton_trick))
      {
	Real_wp rho, sigma;
	var_problem.GetMassDampingCoefficient(i, 0, this->num_unknown, ref, rho, sigma);
        Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;
        
        Fb.SolveCholesky(transA, Uloc);
        Mlt(1.0/sqrt(coef), Uloc);
        
      }
    else if (varying)
      {
        this->block_mass(this->num_unknown)(i)->SolveCholesky(transA, Uloc);
      }
    else
      {
	Real_wp rho, sigma;
	var_problem.GetMassDampingCoefficient(i, 0, this->num_unknown, ref, rho, sigma);
        Real_wp coef = rho*this->coef_mass_ + sigma*this->coef_damp_;
        if (affine)
          {
            Fb.SolveCholesky(transA, Uloc);
            coef *= var_problem.GetAffineJacobian(i);
          }
        else
	  this->block_mass(this->num_unknown)(i)->SolveCholesky(transA, Uloc);
	
        Mlt(1.0/sqrt(coef), Uloc);
      }
  }


  //! method that can be overloaded
  void DiscontinuousMassMatrix_Base
  ::SetDiagonalPML(int nb_new, int m1, int m2, int nodl, VectReal_wp& diag_masse)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  /*******************************
   * DiscontinuousMassMatrix_Dim *
   *******************************/
  

  //! constructor
  template<class Dimension>
  DiscontinuousMassMatrix_Dim<Dimension>
  ::DiscontinuousMassMatrix_Dim(DistributedProblem<Dimension>& var,
				VarInstationary_Base& var_t,
				DiscontinuousMassMatrixProblem& mass_object,
				DiscontinuousMassMatrix_Base* var_geom,
				int type_mat)
    : DiscontinuousMassMatrix_Base(var_t, var_geom, mass_object, type_mat),
      mesh(var.mesh), mesh_num(var.GetMeshNumbering(0)),
      var_boundary(var.GetBoundaryConditionProblem()), var_problem(var), var_time(var_t)
  {
  }
  
  
  //! for PML layers
  template<class Dimension>
  void DiscontinuousMassMatrix_Dim<Dimension>
  ::SetDiagonalPML(int nb_new, int m1, int m2, int nodl, VectReal_wp& diag_masse)
  {
    int nb_blocks = mesh.GetNbElt() + nb_new;
    int nb_old = block_mass_store.GetM();
    coef_block_mass.Reallocate(m2);
    for (int m = m1; m < m2; m++)
      {
	coef_block_mass(m).Reallocate(nb_blocks);
	coef_block_mass(m).Fill(1.0);
      }

    int npml = var_boundary.GetNbEltPML();
    int nb_elt_pml = npml;
    block_mass_store.Resize(nb_old+(nb_new+npml)*(m2-m1));
    int nb = nb_old, nb_elt = mesh.GetNbElt();
    
    if (nb_new > 0)
      {
	block_num_element.Resize(nb_blocks);
	block_size.Resize(nb_blocks+1);
	level_time.SetNbElt(nb_blocks, 0);
	npml = 0;
	for (int i = 0; i < nb_elt; i++)
	  if (var_problem.InsidePML(i))
	    {
	      block_num_element(nb_elt+npml) = i;
	      int nb_dof_loc = block_size(i+1) - block_size(i);
	      block_size(nb_elt+npml+1) = block_size(nb_elt+npml) + nb_dof_loc;
	      npml++;
	    }

	if (Dimension::dim_N == 3)
	  for (int i = 0; i < nb_elt; i++)
	    if (var_problem.InsidePML(i))
	      {
		block_num_element(nb_elt+npml) = i;
		int nb_dof_loc = block_size(i+1) - block_size(i);
		block_size(nb_elt+npml+1) = block_size(nb_elt+npml) + nb_dof_loc;
		npml++;
	      }
	
	for (int m = m1; m < m2; m++)
	  block_mass(m).Resize(nb_blocks);
      }
    
    for (int m = m1; m < m2; m++)
      {
	npml = 0;
	for (int i = 0; i < nb_elt; i++)
	  {
	    num_unknown = m;
	    bool affine = mesh.IsElementAffine(i);
	    int nb_dof_loc = mesh_num.GetNbLocalDof(i);
	    int ref = mesh.Element(i).GetReference();
	    
	    if (var_problem.InsidePML(i))
	      {
		DiscontinuousDiagonalMassMatrix* Mh;
		
		VectReal_wp diagonal(nb_dof_loc);
		int offset = mesh_num.Element(i).GetNumberDof(0) + (m-m1)*nodl;
		for (int j = 0; j < nb_dof_loc; j++)
		  diagonal(j) = diag_masse(offset+j);
		
		Mh = new DiscontinuousDiagonalMassMatrix(diagonal);
		block_mass_store(nb) = Mh;
		block_mass(num_unknown)(i) = block_mass_store(nb);
		nb++;
		
		if (nb_new > 0)
		  {
		    offset = mesh_num.GetNbDof() + 
		      mesh_num.GetDofPML(mesh_num.Element(i).GetNumberDof(0)) + (m-m1)*nodl;
		    
		    for (int j = 0; j < nb_dof_loc; j++)
		      diagonal(j) = diag_masse(offset+j);

		    Mh = new DiscontinuousDiagonalMassMatrix(diagonal);
		    block_mass_store(nb) = Mh;
		    block_mass(num_unknown)(npml+nb_elt) = block_mass_store(nb);
		    nb++;
		    
		    if (Dimension::dim_N == 3)
		      {
			offset += mesh_num.GetNbDofPML();
			
			for (int j = 0; j < nb_dof_loc; j++)
			  diagonal(j) = diag_masse(offset+j);

			Mh = new DiscontinuousDiagonalMassMatrix(diagonal);
			block_mass_store(nb) = Mh;
			block_mass(num_unknown)(nb_elt+nb_elt_pml+npml) = block_mass_store(nb);
			nb++;
		      }
		  }
		
		npml++;
	      }
	    else
	      {
		Real_wp rho, sigma;
		var_time.GetMassDampingCoefficient(i, 0, num_unknown, ref, rho, sigma);
		Real_wp coef = rho*coef_mass_ + sigma*coef_damp_;
		
		if ((!affine) && ElementReference_Base::use_warburton_trick)
		  {
		  }
		else if (var_problem.IsVaryingMedia(num_unknown, ref))
		  {
		    coef = 1.0;
		  }
		else
		  {
		    if (affine)
		      coef *= var_problem.Glob_jacobian(i)(0);
		  }
		
		coef_block_mass(m)(i) = coef;
	      }
	  }
      }
  }

  
  /***************************
   * DiscontinuousMassMatrix *
   ***************************/

  
  //! constructor with a given problem
  template<class Dimension>
  DiscontinuousMassMatrix<Dimension>
  ::DiscontinuousMassMatrix(VarProblem<Dimension>& var,
			    VarInstationary_Fem<Dimension>& var_t,
			    DiscontinuousMassMatrix_Base* var_geom,
			    int type_mat)
    : DiscontinuousMassMatrix_Dim<Dimension>(var, var_t, mass_object, var_geom, type_mat),
    mass_object(var, var_t)
  {
  }


  /**********************************
   * DiscontinuousSurfaceMassMatrix *
   **********************************/


  template<class Dimension>
  DiscontinuousSurfaceMassMatrix<Dimension>
  ::DiscontinuousSurfaceMassMatrix(VarProblem<Dimension>& var,
				   VarInstationary_Base& var_t,
				   DiscontinuousMassMatrix_Base* var_geom, int type_mat,
				   VectReal_wp& val_ds, IVect& offset_ds)
    : DiscontinuousMassMatrix_Base(var_t, var_geom, mass_surf, type_mat),
      mass_surf(var, var_t, val_ds, offset_ds, false)
  {
  }


  template<class Dimension>
  DiscontinuousSurfaceMassMatrix<Dimension>
  ::DiscontinuousSurfaceMassMatrix(VarProblem<Dimension>& var,
				   VarInstationary_Base& var_t,
				   DiscontinuousMassMatrix_Base* var_geom, int type_mat)
    : DiscontinuousMassMatrix_Base(var_t, var_geom, mass_surf, type_mat),
      mass_surf(var, var_t, ValueDsj, OffsetDsj, true)
  {
  }


  template<class Dimension>
  VectReal_wp& DiscontinuousSurfaceMassMatrix<Dimension>::GetDs()
  {
    return ValueDsj;
  }


  template<class Dimension>
  IVect& DiscontinuousSurfaceMassMatrix<Dimension>::GetOffsetDs()
  {
    return OffsetDsj;
  }

  
  template<class Dimension>
  size_t DiscontinuousSurfaceMassMatrix<Dimension>::GetMemorySize() const
  {
    size_t taille = DiscontinuousMassMatrix_Base::GetMemorySize();
    taille += ValueDsj.GetMemorySize() + OffsetDsj.GetMemorySize() + mass_surf.GetMemorySize();
    return taille;
  }
  
}

#define MONTJOIE_FILE_DISCONTINUOUS_MASS_MATRIX_CXX
#endif
