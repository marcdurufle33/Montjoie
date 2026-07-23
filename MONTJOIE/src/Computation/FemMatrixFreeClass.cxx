#ifndef MONTJOIE_FILE_FEM_MATRIX_FREE_CLASS_CXX

#include "FemMatrixFreeClassInline.cxx"

namespace Montjoie
{
  
  /****************************
   * ExtrapolVariablesProduct *
   ****************************/
  
  
  //! Retrieves informations about the current face
  /*!
    \param[in] var considered problem
    \param[in] mesh mesh associated with the considered problem
    \param[in] num_face number of the current face
    \param[in] num_elem1 first element close to the face
    \param[in] num_pos1_face local position of the face in the element num_elem1
    \param[in] offset_face1 offset for accessing values of u on quadrature points of the face
                            (in extrapolU)
    \param[out] ref_line reference of the current face
    \param[out] cond boundary condition associated with the current face
    \param[out] new_face if true, it is the first time this face is treated
    \param[out] num_elem2 second element close to the face
    \param[out] num_pos2_face local position of the face in the second element
    \param[out] offset_face2 offset for accessing values of u on quadrature points of the face
                            (in extrapolU)
    \param[out] rot difference of orientation between faces of the two elements
    \param[out] phase phase used for quasi-periodic condition
    \param[out] phase_conj conjugate of phase (if complex)
   */
  template<class Complexe> template<class Dimension, class T0>
  void ExtrapolVariablesProduct_Base<Complexe>::
  GetFaceOrientationAndNumber(const DistributedProblem<Dimension>& var,
                              int num_face, int num_elem1, int num_pos1_face,
                              int offset_face1, int & ref_line, int& cond,
                              bool& new_face, int& num_elem2, int& num_pos2_face,
                              int& offset_face2, int& rot,
                              T0& phase, T0& phase_conj)
  {
    const Mesh<Dimension>& mesh = var.mesh;
    const MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
    const VarBoundaryCondition_Base& var_boundary = var.GetBoundaryConditionProblem();
    
    phase = 1.0; phase_conj = 1.0;
    int rot1 = mesh.Element(num_elem1).GetOrientationBoundary(num_pos1_face);
    num_elem2 = mesh.Boundary(num_face).numElement(0);
    new_face = var.IsNewFace(num_elem1)(num_pos1_face);
    ref_line = mesh.Boundary(num_face).GetReference();
    cond = mesh.GetBoundaryCondition(ref_line);
    
    num_pos2_face = -1;
    int nv = mesh.Boundary(num_face).GetNbVertices();
    rot = 0;
    int rot2;
    if (mesh.Boundary(num_face).GetNbElements() == 2)
      if (num_elem2 == num_elem1)
        num_elem2 = mesh.Boundary(num_face).numElement(1);
    
    offset_face2 = mesh_num.OffsetQuadElementNumber(num_elem2);
    if (num_elem2 != num_elem1)
      {
        int nb = 0;
        for (int k = 0; k < mesh.Element(num_elem2).GetNbBoundary(); k++)
          {
            int nf = mesh.Element(num_elem2).numBoundary(k);
            if (nf == num_face)
              {
                num_pos2_face = k;
                break;
              }
            else
              nb += mesh_num.GetNbPointsQuadratureBoundary(nf);
          }
	    
        if (num_pos2_face < 0)
          {
            int nf2 = mesh_num.GetPeriodicBoundary(num_face);
            var_boundary.GetPeriodicPhase(num_face, phase);
            var_boundary.GetPeriodicPhase(nf2, phase_conj);            
            nb = 0;
            for (int k = 0; k < mesh.Element(num_elem2).GetNbBoundary(); k++)
              {
                int nf = mesh.Element(num_elem2).numBoundary(k);
                if (nf == nf2)
                  {
                    num_pos2_face = k;
                    break;
                  }
                else
                  nb += mesh_num.GetNbPointsQuadratureBoundary(nf);
              }
          }
        
        offset_face2 += nb;
	
        rot2 = mesh.Element(num_elem2).GetOrientationBoundary(num_pos2_face);
        rot = mesh_num.GetRotationFace(rot1, rot2, nv);
      }
    else
      offset_face2 = offset_face1;
  }

  
  /***************************
   * FemMatrixFreeClass_Base *
   ***************************/
  
  
  //! initializes parameters with default values
  template<class T>
  void FemMatrixFreeClass_Base<T>::InitDefaultValues()
  {
    coef_dirichlet = var_dir.GetCoefficientDirichlet();
    prod_level = NULL;
    row_scaling = false; row_scale = NULL;
    column_scaling = false; col_scale = NULL;
    sym_matrix = false;
    sym_dirichlet_condition = var_dir.UseSymmetricDirichlet();
    
    ignore_dirichlet_dof = false;

    source_real_inhg_dirichlet.Reallocate(1);
    source_cplx_inhg_dirichlet.Reallocate(1);

    condensed_solver = NULL;
    use_symmetrization = false;

    this->m_ = var_problem.GetNbDof();
    this->n_ = var_problem.GetNbDof();
  }
  
  
  //! sets the size of the matrix
  template<class T>
  void FemMatrixFreeClass_Base<T>::Reallocate(int m, int n)
  {
    this->m_ = m;
    this->n_ = n;
  }
  

  //! Destruction of the matrix
  template<class T>
  void FemMatrixFreeClass_Base<T>::Clear()
  {
    mat_iterative_sym.Clear();
    mat_boundary_sym.Clear();
    mat_iterative_unsym.Clear();
    mat_boundary_unsym.Clear();
    matCSR_iterative_sym.Clear();
    matCSR_boundary_sym.Clear();
    matCSR_iterative_unsym.Clear();
    matCSR_boundary_unsym.Clear();
    
    for (int k = 0; k < source_real_inhg_dirichlet.GetM(); k++)
      source_real_inhg_dirichlet(k).Clear();

    for (int k = 0; k < source_cplx_inhg_dirichlet.GetM(); k++)
	source_cplx_inhg_dirichlet(k).Clear();
    
    column_dirichlet.Clear();
  }
  

  //! returns size of memory used by the object in bytes
  template<class T>
  size_t FemMatrixFreeClass_Base<T>::GetMemorySize() const
  {
    size_t taille = mat_boundary_sym.GetMemorySize() + mat_boundary_unsym.GetMemorySize();
    taille += mat_iterative_sym.GetMemorySize() + mat_iterative_unsym.GetMemorySize();
    taille += matCSR_boundary_sym.GetMemorySize() + matCSR_boundary_unsym.GetMemorySize();
    taille += matCSR_iterative_sym.GetMemorySize() + matCSR_iterative_unsym.GetMemorySize();
    taille += Seldon::GetMemorySize(source_real_inhg_dirichlet);
    taille += Seldon::GetMemorySize(source_cplx_inhg_dirichlet);
    taille += column_dirichlet.GetMemorySize();
    taille += sizeof(*this);
    return taille;
  }
  
    
  //! initialisation of the matrix as a symmetric matrix
  template<class T>
  void FemMatrixFreeClass_Base<T>::InitSymmetricMatrix()
  {
    sym_matrix = true;
    mat_iterative_unsym.Clear(); mat_boundary_unsym.Clear();
    matCSR_iterative_unsym.Clear(); matCSR_boundary_unsym.Clear();
  }
  
  
  //! initialisation of the matrix as an unsymmetric matrix
  template<class T>
  void FemMatrixFreeClass_Base<T>::InitUnsymmetricMatrix()
  {
    sym_matrix = false;
    mat_iterative_sym.Clear(); mat_boundary_sym.Clear();
    matCSR_iterative_sym.Clear(); matCSR_boundary_sym.Clear();
  }
  

  //! returns true if the pointer successfully adress the finite element matrix
  template<class T>
  bool FemMatrixFreeClass_Base<T>
  ::SucceedInAffectingPointer(typename GenericStorage<T>::SparseSymMatrix*& ptrA,
			      typename GenericStorage<T>::SparseSymMatrixCSR*& ptrB)
  {
    ptrA = NULL; ptrB = NULL;
    if (mat_iterative_sym.GetM() > 0)
      {
        ptrA = &mat_iterative_sym;
        return true;
      }

    if (matCSR_iterative_sym.GetM() > 0)
      {
        ptrB = &matCSR_iterative_sym;
        return true;
      }      
    
    return false;
  }
  

  //! returns true if the pointer successfully adress the finite element matrix
  template<class T>
  bool FemMatrixFreeClass_Base<T>
  ::SucceedInAffectingPointer(typename GenericStorage<T>::SparseMatrix*& ptrA,
			      typename GenericStorage<T>::SparseMatrixCSR*& ptrB)
  {
    ptrA = NULL; ptrB = NULL;
    if (mat_iterative_unsym.GetM() > 0)
      {
        ptrA = &mat_iterative_unsym;
        return true;
      }

    if (matCSR_iterative_unsym.GetM() > 0)
      {
        ptrB = &matCSR_iterative_unsym;
        return true;
      }
    
    return false;
  }
  

  //! applies scaling on columns
  template<class T> template<class Vector1>
  void FemMatrixFreeClass_Base<T>
  ::ApplyRightScaling(const Vector1& B2, Vector1& C2, Vector1& B, Vector1& C) const
  {
    if (column_scaling)
      {
	B = B2;
	for (int i = 0; i < col_scale->GetM(); i++)
	  B(i) *= (*col_scale)(i);
	
	C.Reallocate(this->GetM()); C.Fill(0);
      }
    else
      {
	B.SetData(B2.GetM(), B2.GetData());
	C.SetData(C2.GetM(), C2.GetData());
      }
  }
  

  //! applies scaling on rows
  template<class T> template<class Vector1>
  void FemMatrixFreeClass_Base<T>
  ::ApplyLeftScaling(const Vector1& B2, Vector1& C2, Vector1& B, Vector1& C) const
  {
    if (row_scaling)
      {
	for (int i = 0; i < row_scale->GetM(); i++)
	  C2(i) += (*row_scale)(i)*C(i);
      }
    else
      {
	B.Nullify(); C.Nullify();
      }
  }
  
      
  //! matrices are converted to csr format to reduce memory usage
  template<class T>
  void FemMatrixFreeClass_Base<T>::CompressMatrix()
  {
    // compressing matrices for the iterative resolution
    if (mat_boundary_sym.GetM() > 0)
      {
	Seldon::Copy(mat_boundary_sym, matCSR_boundary_sym);
	mat_boundary_sym.Clear();
      }
    
    if (mat_iterative_sym.GetM() > 0)
      {
	Seldon::Copy(mat_iterative_sym, matCSR_iterative_sym);
	mat_iterative_sym.Clear();
      }

    if (mat_boundary_unsym.GetM() > 0)
      {
	Seldon::Copy(mat_boundary_unsym, matCSR_boundary_unsym);
	mat_boundary_unsym.Clear();
      }
    
    if (mat_iterative_unsym.GetM() > 0)
      {
	Seldon::Copy(mat_iterative_unsym, matCSR_iterative_unsym);
	mat_iterative_unsym.Clear();
      }
  }

  
  //! writes stored matrix in the disk
  template<class T>
  void FemMatrixFreeClass_Base<T>::WriteText(const string& file_name) const
  {
    if (mat_iterative_sym.GetM() > 0)
      return mat_iterative_sym.WriteText(file_name);
    
    if (matCSR_iterative_sym.GetM() > 0)
      return matCSR_iterative_sym.WriteText(file_name);

    if (mat_iterative_unsym.GetM() > 0)
      return mat_iterative_unsym.WriteText(file_name);
    
    if (matCSR_iterative_unsym.GetM() > 0)
      return matCSR_iterative_unsym.WriteText(file_name);
  }
  
  
  //! Adds terms due to boundary conditions in the matrix vector product
  template<>
  void FemMatrixFreeClass_Base<Real_wp>::
  AddExtraBoundaryTerms(const Real_wp& alpha,
			const Vector<Real_wp>& B, Vector<Real_wp>& C) const
  {
    Real_wp one(1);
    // additional terms are added without assembling the final result
    if (sym_matrix)
      {
	// additional terms from boundary
	if (mat_boundary_sym.GetM() > 0)
	  Seldon::MltAdd(alpha, mat_boundary_sym, B, one, C, false);
	else if (matCSR_boundary_sym.GetM() > 0)
	  Seldon::MltAdd(alpha, matCSR_boundary_sym, B, one, C, false);
      }
    else
      {
	// additional terms from boundary
	if (mat_boundary_unsym.GetM() > 0)
	  Seldon::MltAdd(alpha, mat_boundary_unsym, B, one, C, false);
	else if (matCSR_boundary_unsym.GetM() > 0)
	  Seldon::MltAdd(alpha, matCSR_boundary_unsym, B, one, C, false);
      }
  }


  //! Adds terms due to boundary conditions in the matrix vector product
  template<>
  void FemMatrixFreeClass_Base<Real_wp>::
  AddExtraBoundaryTerms(const Complex_wp& alpha,
			const Vector<Complex_wp>& B, Vector<Complex_wp>& C) const
  {
    Complex_wp one(1, 0);
    // additional terms are added without assembling the final result
    if (sym_matrix)
      {
	// additional terms from boundary
	if (mat_boundary_sym.GetM() > 0)
	  Seldon::MltAdd(alpha, mat_boundary_sym, B, one, C, false);
	else if (matCSR_boundary_sym.GetM() > 0)
	  Seldon::MltAdd(alpha, matCSR_boundary_sym, B, one, C, false);
      }
    else
      {
	// additional terms from boundary
	if (mat_boundary_unsym.GetM() > 0)
	  Seldon::MltAdd(alpha, mat_boundary_unsym, B, one, C, false);
	else if (matCSR_boundary_unsym.GetM() > 0)
	  Seldon::MltAdd(alpha, matCSR_boundary_unsym, B, one, C, false);
      }
  }


  //! Adds terms due to boundary conditions in the matrix vector product
  template<>
  void FemMatrixFreeClass_Base<Complex_wp>::
  AddExtraBoundaryTerms(const Real_wp& alpha, const Vector<Real_wp>& B,
                        Vector<Real_wp>& C) const
  {
    cout << "We can not perform this operation with a real vector" << endl;
    abort();
  }
  
  
  //! Adds terms due to boundary conditions in the matrix vector product
  template<>
  void FemMatrixFreeClass_Base<Complex_wp>::
  AddExtraBoundaryTerms(const Complex_wp& alpha, const Vector<Complex_wp>& B,
                        Vector<Complex_wp>& C) const
  {
    Complex_wp one(1, 0);
    // additional terms are added without assembling the final result
    if (sym_matrix)
      {
	// additional terms from boundary
	if (mat_boundary_sym.GetM() > 0)
	  Seldon::MltAdd(alpha, mat_boundary_sym, B, one, C, false);
	else if (matCSR_boundary_sym.GetM() > 0)
	  Seldon::MltAdd(alpha, matCSR_boundary_sym, B, one, C, false);
      }
    else
      {
	// additional terms from boundary
	if (mat_boundary_unsym.GetM() > 0)
	  Seldon::MltAdd(alpha, mat_boundary_unsym, B, one, C, false);
	else if (matCSR_boundary_unsym.GetM() > 0)
	  Seldon::MltAdd(alpha, matCSR_boundary_unsym, B, one, C, false);
      }
  }
  

  //! sets the number of Dirichlet conditions to store
  template<class T>
  void FemMatrixFreeClass_Base<T>::SetNbDirichletCondition(int n)
  {
    source_real_inhg_dirichlet.Reallocate(n);
    source_cplx_inhg_dirichlet.Reallocate(n);
  }
  
  
  //! applying inhomogenous Dirichlet condition to right hand side b_rhs
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::ApplyDirichletCondition(const SeldonTranspose& trans, Vector<Real_wp>& b_rhs, int k, bool assemble)
  {  
    if (var_dir.GetNbGlobalDirichletDof() <= 0)
      return;

    Vector<Real_wp> source_dirichlet(var_problem.GetNbDof());
    source_dirichlet.Zero();
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      {
        source_dirichlet(var_dir.GetDirichletDofNumber(i))
          = b_rhs(var_dir.GetDirichletDofNumber(i));        
      }

    // source_dirichlet is assembled to have the correct value on all processors
    var_comm.AssembleDirichlet(source_dirichlet);
    
    if (var_dir.GetNbDirichletDof() > 0)
      {
        // storing values of source on Dirichlet dofs
        Vector<Real_wp>& source_inhg = source_real_inhg_dirichlet(k);
        source_inhg.Reallocate(var_dir.GetNbDirichletDof());
        for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
          source_inhg(i) = source_dirichlet(var_dir.GetDirichletDofNumber(i));
      }
    
    // for unsymmetric matrices, only rows associated with Dirichlet dofs
    // are modified to set the equation u = u_dir
    // therefore the right hand side is not modified
    if (!sym_matrix && !sym_dirichlet_condition)
      {
        // right hand side is modified by multiplying with the correct coefficient
        for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
          b_rhs(var_dir.GetDirichletDofNumber(i)) *= coef_dirichlet;

        return;
      }
    
    // for symmetric matrices, rows and columns associated with Dirichlet dofs
    // are cleared except on the diagonal such that the condition u = 0 is set
    // on Dirichlet dofs
    // The right hand side b is modified as follows :
    // b <- b - A source_dirichlet
    // where source_dirichlet contains the inhomogeneous Dirichlet condition
    // Therefore a product matrix vector involving Dirichlet columns is performed
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      SetComplexZero(b_rhs(var_dir.GetDirichletDofNumber(i)));

    if (!trans.NoTrans())
      return;
    
    // performing the matrix-vector product b = b - A source_dirichlet
    VectReal_wp res;

    if (assemble)
      {
	res.Reallocate(b_rhs.GetM());
	res.Zero();
      }
    else
      res.SetData(b_rhs.GetM(), b_rhs.GetData());
    
    if (column_dirichlet.GetM() > 0)
      {
	// part of the matrix stored :
	MltAdd(-Real_wp(1), SeldonNoTrans,
	       column_dirichlet, source_dirichlet, Real_wp(1), res, false);
      }
    else
      {	
	this->MltAddVector(-Real_wp(1), source_dirichlet,
			   Real_wp(1), res, false);
      }
    
    // assembling the result if needed
    if (assemble)
      {
	var_comm.AddDomains(res);
	Add(Real_wp(1), res, b_rhs);
      }
    else
      res.Nullify();
    
    // then Dirichlet dofs are set to 0
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      b_rhs(var_dir.GetDirichletDofNumber(i)) = 0.0;
  }  
  

  //! applying inhomogenous Dirichlet condition to right hand side b_rhs
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::ApplyDirichletCondition(const SeldonTranspose& trans, Vector<Complex_wp>& b_rhs, int k, bool assemble)
  {  
    if (var_dir.GetNbGlobalDirichletDof() <= 0)
      return;

    Vector<Complex_wp> source_dirichlet(var_problem.GetNbDof());
    source_dirichlet.Zero();
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      {
        source_dirichlet(var_dir.GetDirichletDofNumber(i))
          = b_rhs(var_dir.GetDirichletDofNumber(i));        
      }

    // source_dirichlet is assembled to have the correct value on all processors
    var_comm.AssembleDirichlet(source_dirichlet);
    
    if (var_dir.GetNbDirichletDof() > 0)
      {
        // storing values of source on Dirichlet dofs
        Vector<Complex_wp>& source_inhg = source_cplx_inhg_dirichlet(k);
        source_inhg.Reallocate(var_dir.GetNbDirichletDof());
        for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
          source_inhg(i) = source_dirichlet(var_dir.GetDirichletDofNumber(i));
      }
    
    // for unsymmetric matrices, only rows associated with Dirichlet dofs
    // are modified to set the equation u = u_dir
    // therefore the right hand side is not modified
    if (!sym_matrix && !sym_dirichlet_condition)
      {
        // right hand side is modified by multiplying with the correct coefficient
        for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
          b_rhs(var_dir.GetDirichletDofNumber(i)) *= coef_dirichlet;

        return;
      }
    
    // for symmetric matrices, rows and columns associated with Dirichlet dofs
    // are cleared except on the diagonal such that the condition u = 0 is set
    // on Dirichlet dofs
    // The right hand side b is modified as follows :
    // b <- b - A source_dirichlet
    // where source_dirichlet contains the inhomogeneous Dirichlet condition
    // Therefore a product matrix vector involving Dirichlet columns is performed

    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      SetComplexZero(b_rhs(var_dir.GetDirichletDofNumber(i)));
    
    if (!trans.NoTrans())
      return;
    
    // performing the matrix-vector product b = b - A source_dirichlet
    VectComplex_wp res;
    if (assemble)
      {
	res.Reallocate(b_rhs.GetM());
	res.Zero();
      }
    else
      res.SetData(b_rhs.GetM(), b_rhs.GetData());
    
    if (column_dirichlet.GetM() > 0)
      {
	// if the matrix is stored :
	MltAdd(-Complex_wp(1, 0), SeldonNoTrans,
	       column_dirichlet, source_dirichlet,
	       Complex_wp(1, 0), res, false);
      }
    else
      {
	// otherwise
	this->MltAddVector(-Complex_wp(1, 0), source_dirichlet,
			   Complex_wp(1, 0), res, false);
      }
    
    // assembling the result
    if (assemble)
      {
	var_comm.AddDomains(res);
	Add(Complex_wp(1, 0), res, b_rhs);
      }
    else
      res.Nullify();
    
    // then Dirichlet dofs are set to 0
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      b_rhs(var_dir.GetDirichletDofNumber(i)) = 0.0;
  }  

  
  //! computes C = C + alpha A B
  //! where A contains columns erased because of Dirichlet dofs
  template<class T> template<class T0>
  void FemMatrixFreeClass_Base<T>
  ::MltAddHetereogeneousDirichlet(const T0& alpha, const SeldonTranspose& trans,
				  const Vector<T0>& B, Vector<T0>& C) const
  {
    T0 one; SetComplexOne(one);
    if (column_dirichlet.GetM() > 0)
      MltAdd(alpha, trans, column_dirichlet, B, one, C, false);
  }
    

  //! initialization method to take into account inhomogeneous Dirichlet condition
  template<class T> template<class MatrixSparse>
  void FemMatrixFreeClass_Base<T>
  ::SetDirichletCondition(MatrixSparse& mat_sp, int offset_row, int offset_col,
                          bool erase_col)
  {
    coef_dirichlet = var_dir.GetCoefficientDirichlet();
    source_real_inhg_dirichlet.Clear();
    source_cplx_inhg_dirichlet.Clear();
    column_dirichlet.Clear();
    if (var_dir.GetNbGlobalDirichletDof() <= 0)
      return;
    
    Vector<int> row_dir = var_dir.GetDirichletDofNumber();
    Vector<int> col_dir = var_dir.GetDirichletDofNumber();
    for (int i = 0; i < row_dir.GetM(); i++)
      {
	row_dir(i) += offset_row;
	col_dir(i) += offset_col;
      }
    
    if (sym_matrix || sym_dirichlet_condition)
      {
        IVect other_dof(this->m_ - var_dir.GetNbDirichletDof());
        int nb = 0;
        if (var_dir.GetNbDirichletDof() > 0)
          for (int i = 0; i < this->m_; i++)
            if (!var_dir.IsDofDirichlet(i))
              other_dof(nb++) = offset_row + i;
        
        // we are extracting columns of the matrix
        if (!erase_col)
          CopySubMatrix(mat_sp, other_dof,
                        col_dir, column_dirichlet);
        
        // we remove rows and columns related to Dirichlet dofs
        EraseRow(row_dir, mat_sp);
        if (!sym_matrix || !mat_sp.IsSymmetric())
          EraseCol(col_dir, mat_sp);
      }
    else
      {        
        // we remove rows related to Dirichlet dofs
        EraseRow(row_dir, mat_sp);
        if (erase_col)
          EraseCol(col_dir, mat_sp);
      }
    
    source_real_inhg_dirichlet.Reallocate(1);
    source_cplx_inhg_dirichlet.Reallocate(1);
    
    Vector<bool> OverlappedRow(var_problem.GetNbDof());
    OverlappedRow.Fill(false);
    for (int i = 0; i < var_comm.GetNbOverlappedDof(); i++)
      OverlappedRow(var_comm.GetOverlappedDofNumber(i)) = true;
    
    // on the diagonal, we put Dirichlet coefficient (e.g. 1)
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      if (!OverlappedRow(var_dir.GetDirichletDofNumber(i)))
        mat_sp.AddInteraction(row_dir(i), col_dir(i), coef_dirichlet);        
  }  
  

  //! initialization method to take into account inhomogeneous Dirichlet condition  
  template<class T>
  void FemMatrixFreeClass_Base<T>::SetDirichletCondition(int offset_row, int offset_col,
                                                         bool erase_col)
  {
    if (sym_matrix)
      {
	if (mat_iterative_sym.GetM() > 0)
	  SetDirichletCondition(mat_iterative_sym, offset_row, offset_col, erase_col);
	else if (matCSR_iterative_sym.GetM() > 0)
	  SetDirichletCondition(matCSR_iterative_sym, offset_row, offset_col, erase_col);
      }
    else
      {
	if (mat_iterative_unsym.GetM() > 0)
	  SetDirichletCondition(mat_iterative_unsym, offset_row, offset_col, erase_col);
	else if (matCSR_iterative_unsym.GetM() > 0)
	  SetDirichletCondition(matCSR_iterative_unsym, offset_row, offset_col, erase_col);
      }
  }


  template<class T> template<class Prop>
  void FemMatrixFreeClass_Base<T>::SetDirichletCondition(Matrix<T, Prop, DiagonalRow>& A,
							 int offset_row, int offset_col,
                                                         bool erase_col)
  {
    if (var_dir.GetNbGlobalDirichletDof() <= 0)
      return;

    coef_dirichlet = var_dir.GetCoefficientDirichlet();
    Vector<int> row_dir = var_dir.GetDirichletDofNumber();
    Vector<int> col_dir = var_dir.GetDirichletDofNumber();
    for (int i = 0; i < row_dir.GetM(); i++)
      {
	row_dir(i) += offset_row;
	col_dir(i) += offset_col;
      }

    Vector<bool> OverlappedRow(var_problem.GetNbDof());
    OverlappedRow.Fill(false);
    for (int i = 0; i < var_comm.GetNbOverlappedDof(); i++)
      OverlappedRow(var_comm.GetOverlappedDofNumber(i)) = true;

    T zero; SetComplexZero(zero);
    // on the diagonal, we put Dirichlet coefficient (e.g. 1)
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      A.Get(row_dir(i), col_dir(i)) = zero;

    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      if (!OverlappedRow(var_dir.GetDirichletDofNumber(i)))
        A.Get(row_dir(i), col_dir(i)) =  coef_dirichlet;        
  }
  
  
  template<class T> template<class Prop>
  void FemMatrixFreeClass_Base<T>::SetDirichletCondition(DistributedMatrix<T, Prop, BlockDiagRow>& mat_sp,
							 int offset_row, int offset_col, bool erase_col)
  {
    if (var_dir.GetNbGlobalDirichletDof() <= 0)
      return;

    coef_dirichlet = var_dir.GetCoefficientDirichlet();    
    Vector<int> row_dir = var_dir.GetDirichletDofNumber();
    Vector<int> col_dir = var_dir.GetDirichletDofNumber();
    for (int i = 0; i < row_dir.GetM(); i++)
      {
	row_dir(i) += offset_row;
	col_dir(i) += offset_col;
      }

    if (sym_matrix || sym_dirichlet_condition)
      {
        // we remove rows and columns related to Dirichlet dofs
        EraseRow(row_dir, mat_sp);
        EraseCol(col_dir, mat_sp);
      }
    else
      {        
        // we remove rows related to Dirichlet dofs
        EraseRow(row_dir, mat_sp);
        if (erase_col)
          EraseCol(col_dir, mat_sp);
      }
    
    Vector<bool> OverlappedRow(var_problem.GetNbDof());
    OverlappedRow.Fill(false);
    for (int i = 0; i < var_comm.GetNbOverlappedDof(); i++)
      OverlappedRow(var_comm.GetOverlappedDofNumber(i)) = true;
    
    // on the diagonal, we put Dirichlet coefficient (e.g. 1)
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      if (!OverlappedRow(var_dir.GetDirichletDofNumber(i)))
        mat_sp.AddInteraction(row_dir(i), col_dir(i), coef_dirichlet);        
    
  }
  
  
  template<class T> template<class Prop>
  void FemMatrixFreeClass_Base<T>::SetDirichletCondition(DistributedMatrix<T, Prop, BlockDiagRowSym>& mat_sp,
							 int offset_row, int offset_col, bool erase_col)
  {
    if (var_dir.GetNbGlobalDirichletDof() <= 0)
      return;

    coef_dirichlet = var_dir.GetCoefficientDirichlet();    
    Vector<int> row_dir = var_dir.GetDirichletDofNumber();
    Vector<int> col_dir = var_dir.GetDirichletDofNumber();
    for (int i = 0; i < row_dir.GetM(); i++)
      {
	row_dir(i) += offset_row;
	col_dir(i) += offset_col;
      }

    if (sym_matrix || sym_dirichlet_condition)
      {
        // we remove rows and columns related to Dirichlet dofs
        EraseRow(row_dir, mat_sp);
      }
    else
      {        
        // we remove rows related to Dirichlet dofs
        EraseRow(row_dir, mat_sp);
        if (erase_col)
          EraseCol(col_dir, mat_sp);
      }
    
    Vector<bool> OverlappedRow(var_problem.GetNbDof());
    OverlappedRow.Fill(false);
    for (int i = 0; i < var_comm.GetNbOverlappedDof(); i++)
      OverlappedRow(var_comm.GetOverlappedDofNumber(i)) = true;
    
    // on the diagonal, we put Dirichlet coefficient (e.g. 1)
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      if (!OverlappedRow(var_dir.GetDirichletDofNumber(i)))
        mat_sp.AddInteraction(row_dir(i), col_dir(i), coef_dirichlet);        

  }
  
  
  //! storing values of b on Dirichlet dofs
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::InitDirichletCondition(VectReal_wp& b_rhs, int k)
  {
    source_real_inhg_dirichlet(k).Reallocate(var_dir.GetNbDirichletDof());
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      source_real_inhg_dirichlet(k)(i) = b_rhs(var_dir.GetDirichletDofNumber(i));
  }
  

  //! storing values of b on Dirichlet dofs
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::InitDirichletCondition(VectComplex_wp& b_rhs, int k)
  {
    source_cplx_inhg_dirichlet(k).Reallocate(var_dir.GetNbDirichletDof());
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      source_cplx_inhg_dirichlet(k)(i) = b_rhs(var_dir.GetDirichletDofNumber(i));
  }
  
  
  //! setting x = f on Dirichlet dofs (inhomogeneous Dirichlet condition)
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::ImposeDirichletCondition(const SeldonTranspose& trans, VectReal_wp& x, int k)
  {
    for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
      x(var_dir.GetDirichletDofNumber(i)) = source_real_inhg_dirichlet(k)(i);      
    
    if (trans.Trans() && (column_dirichlet.GetM() > 0))
      {
        VectReal_wp res(x.GetM());
        column_dirichlet.MltVector(trans, x, res);                
        for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
          x(var_dir.GetDirichletDofNumber(i)) -= res(var_dir.GetDirichletDofNumber(i));
      }
  }
  

  //! setting x = f on Dirichlet dofs (inhomogeneous Dirichlet condition)  
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::ImposeDirichletCondition(const SeldonTranspose& trans, VectComplex_wp& x, int k)
  {
     for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
        x(var_dir.GetDirichletDofNumber(i)) = source_cplx_inhg_dirichlet(k)(i);
     
     if (trans.Trans() && (column_dirichlet.GetM() > 0))
       {
         VectComplex_wp res(x.GetM());
         column_dirichlet.MltVector(trans, x, res);                
         for (int i = 0; i < var_dir.GetNbDirichletDof(); i++)
           x(var_dir.GetDirichletDofNumber(i)) -= res(var_dir.GetDirichletDofNumber(i));
       }
  }
  
  
  //! TO DO
  template<class T>
  void FemMatrixFreeClass_Base<T>::AddRowSum(VectReal_wp&)
  {
    cout << "TO DO" << endl;
    abort();
    /*
    IVect num_ddl;
    Matrix<typename NatureMatrix::value_type,
	   typename TypeEquation::Property_Elem_Matrix,
	   typename GenericStorage<typename TypeEquation::Property_Elem_Matrix>
	   ::DenseStorage> mat_elt;
    
    for (int i = 0; i < mat.var.mesh.GetNbElt(); i++)
      {
	var.ComputeElementaryMatrix(i, num_ddl, mat_elt, nat_mat);
	int nb_dof = mat_elt.GetM();
	for(int j = 0; j < nb_dof; j++)
	  for(int k = 0; k < nb_dof; k++)
	    if (num_ddl(j) >= 0)
	      diagonal_scale_left(num_ddl(j)) += abs(mat_elt(j,k));
      }
    */
  }


  //! TO DO
  template<class T>
  void FemMatrixFreeClass_Base<T>::AddRowColSum(VectReal_wp&, VectReal_wp&)
  {
    cout << "TO DO" << endl;
    abort();
    /*
    IVect num_ddl;
    typedef typename TypeEquation::Nature_Matrix NatureMatrix;
    NatureMatrix nat_mat;
    Matrix<typename NatureMatrix::value_type,
	   typename TypeEquation::Property_Elem_Matrix,
	   typename GenericStorage<typename TypeEquation
				   ::Property_Elem_Matrix>::DenseStorage> mat_elt;
    
    for (int i = 0; i < mat.var.mesh.GetNbElt(); i++)
      {
	ComputeElementaryMatrix(i, num_ddl, mat_elt, nat_mat, mat.var);
	int nb_dof = mat_elt.GetM();
	for(int j = 0; j < nb_dof; j++)
	  for(int k = 0; k < nb_dof; k++)
	    if ((num_ddl(j) >= 0) && (num_ddl(k) >= 0))
	      {
		diagonal_scale_left(num_ddl(j)) += abs(mat_elt(j, k));
		diagonal_scale_right(num_ddl(k)) += abs(mat_elt(j, k));
	      }
      }
    */
  }

  
  //! Y = beta Y + alpha A X
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& X, const Treal& beta,
		 Vector<Treal>& Y, bool assemble) const
  {
    MltAddVector(alpha, SeldonNoTrans, X, beta, Y, assemble);
  }


  //! Y = beta Y + alpha A X
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& X, const Tcplx& beta,
		 Vector<Tcplx>& Y, bool assemble) const
  {
    MltAddVector(alpha, SeldonNoTrans, X, beta, Y, assemble);
  }
  
  
  //! Y = beta Y + alpha A X
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& X, const Treal& beta, Vector<Treal>& Y) const
  {
    MltAddVector(alpha, SeldonNoTrans, X, beta, Y, true);
  }


  //! Y = beta Y + alpha A X
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& X, const Tcplx& beta, Vector<Tcplx>& Y) const
  {
    MltAddVector(alpha, SeldonNoTrans, X, beta, Y, true);
  }


  //! Y = beta Y + alpha A X or with transpose of A or transpose conjugate
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans, const Vector<Treal>& X,
		 const Treal& beta, Vector<Treal>& Y) const
  {
    MltAddVector(alpha, trans, X, beta, Y, true);
  }


  //! Y = beta Y + alpha A X or with transpose of A or transpose conjugate
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans, const Vector<Tcplx>& X,
		 const Tcplx& beta, Vector<Tcplx>& Y) const
  {
    MltAddVector(alpha, trans, X, beta, Y, true);
  }
  
  
  //! Y = A X
  template<class T>
  void FemMatrixFreeClass_Base<T>::MltVector(const Vector<Treal>& X, Vector<Treal>& Y) const
  {
    MltVector(SeldonNoTrans, X, Y, true);
  }


  //! Y = A X
  template<class T>
  void FemMatrixFreeClass_Base<T>::MltVector(const Vector<Tcplx>& X, Vector<Tcplx>& Y) const
  {
    MltVector(SeldonNoTrans, X, Y, true);
  }

  
  //! Y = A X (or transpose of A)
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltVector(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y) const
  {
    MltVector(trans, X, Y, true);
  }


  //! Y = A X (or transpose of A)
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltVector(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y) const
  {
    MltVector(trans, X, Y, true);
  }
  

  //! computes Y = A X or A^T X or A^H X
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltVector(const SeldonTranspose& trans, const Vector<Treal>& X, Vector<Treal>& Y, bool assemble) const
  {
    MltVectorGen(trans, X, Y, assemble);
  }


  //! computes Y = A X or A^T X or A^H X
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltVector(const SeldonTranspose& trans, const Vector<Tcplx>& X, Vector<Tcplx>& Y, bool assemble) const
  {
    MltVectorGen(trans, X, Y, assemble);
  }
  
  
  //! computes Y = A X or A^T X or A^H X
  template<class T> template<class T0>
  void FemMatrixFreeClass_Base<T>
  ::MltVectorGen(const SeldonTranspose& trans, const Vector<T0>& X, Vector<T0>& Y, bool assemble) const
  {
    bool matrix_not_stored = false;
    if (this->sym_matrix)
      {
	if (this->mat_iterative_sym.GetM() > 0)
	  Mlt(trans, this->mat_iterative_sym, X, Y, false);
	else if (this->matCSR_iterative_sym.GetM() > 0)
	  Mlt(trans, this->matCSR_iterative_sym, X, Y, false);
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (this->mat_iterative_unsym.GetM() > 0)
	  Mlt(trans, this->mat_iterative_unsym, X, Y, false);
	else if (this->matCSR_iterative_unsym.GetM() > 0)
	  Mlt(trans, this->matCSR_iterative_unsym, X, Y, false);
	else
	  matrix_not_stored = true;
      }
    
    if (matrix_not_stored)
      {
	int ALL_LEVELS = MatrixVectorProductLevel::ALL_LEVELS;
	
	Y.Zero();
	SeldonTranspose transA(trans);
	if (trans.ConjTrans())
	  {
	    transA = SeldonTrans;
	    Conjugate(const_cast<Vector<T0>&>(X));
	  }
	    
	// use of a matrix-vector product on the fly
	// implementations in ProdMatVectH1.cxx, ProdMatVectHcurl2D.cxx, etc
	this->MltAddFree(this->nature_matrix, transA, ALL_LEVELS, X, Y);

	if (trans.ConjTrans())
	  {
	    Conjugate(const_cast<Vector<T0>&>(X));
	    Conjugate(Y);
	  }
      }
    else
      {
	T0 one; SetComplexOne(one);
        if (this->DirichletDofIgnored())
          MltAddHetereogeneousDirichlet(one, trans, X, Y);
      }
    
    // assembling vector C
    if (assemble)
      {
	if (condensed_solver != NULL)
	  condensed_solver->AddDomains(Y);
	else
	  var_comm.AddDomains(Y);
      }
  }
  

  //! computes Y = beta Y + alpha A X, transpose of A may be used depending on trans
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans, const Vector<Treal>& X,
		 const Treal& beta, Vector<Treal>& Y, bool assemble) const
  {
    MltAddVectorGen(alpha, trans, X, beta, Y, assemble);
  }


  //! computes Y = beta Y + alpha A X, transpose of A may be used depending on trans
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans, const Vector<Tcplx>& X,
		 const Tcplx& beta, Vector<Tcplx>& Y, bool assemble) const
  {
    MltAddVectorGen(alpha, trans, X, beta, Y, assemble);
  }
  
  
  //! computes Y = beta Y + alpha A X, transpose of A may be used depending on trans
  template<class T> template<class T0>
  void FemMatrixFreeClass_Base<T>
  ::MltAddVectorGen(const T0& alpha, const SeldonTranspose& trans, const Vector<T0>& X,
		    const T0& beta, Vector<T0>& Y, bool assemble) const
  {
    T0 zero; SetComplexZero(zero);
    if (beta == zero)
      {
	MltVectorGen(trans, X, Y, assemble);    
	Mlt(alpha, Y);
      }
    else
      {
	Vector<T0> AX(Y.GetM());	
	MltVectorGen(trans, X, AX, assemble);
	Mlt(beta, Y);
	Add(alpha, AX, Y);
      }
  }


  //! function that can be overloaded
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddFree(const GlobalGenericMatrix<T>& nat_mat, const SeldonTranspose& trans,
	       int level, const Vector<Treal>& X, Vector<Treal>& Y) const
  {
    cout << "No generic function for MltAddFree " << endl;
    abort();
  }


  //! function that can be overloaded
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::MltAddFree(const GlobalGenericMatrix<T>& nat_mat, const SeldonTranspose& trans,
	       int level, const Vector<Tcplx>& X, Vector<Tcplx>& Y) const
  {
    cout << "No generic function for MltAddFree " << endl;
    abort();
  }

  
  template<class T>
  void FemMatrixFreeClass_Base<T>
  ::ApplySor(const SeldonTranspose& trans, Vector<Treal>& x, const Vector<Treal>& r,
             const typename ClassComplexType<T>::Treal& omega,
             int nb_iter, int stage_ssor) const
  {
    cout << "SOR matrice free" << endl;
    bool matrix_not_stored = false;
    if (this->sym_matrix)
      {
	if (this->mat_iterative_sym.GetM() > 0)
	  this->mat_iterative_sym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else if (this->matCSR_iterative_sym.GetM() > 0)
          this->matCSR_iterative_sym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (this->mat_iterative_unsym.GetM() > 0)
          this->mat_iterative_unsym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else if (this->matCSR_iterative_unsym.GetM() > 0)
          this->matCSR_iterative_unsym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else
	  matrix_not_stored = true;
      }
    
    if (matrix_not_stored)
      cout << "Not implemented " << endl;    
  }
  

  template<class T> 
  void FemMatrixFreeClass_Base<T>
  ::ApplySor(const SeldonTranspose& trans, Vector<Tcplx>& x, const Vector<Tcplx>& r,
             const typename ClassComplexType<T>::Treal& omega,
             int nb_iter, int stage_ssor) const
  {
    bool matrix_not_stored = false;
    if (this->sym_matrix)
      {
	if (this->mat_iterative_sym.GetM() > 0)
	  this->mat_iterative_sym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else if (this->matCSR_iterative_sym.GetM() > 0)
          this->matCSR_iterative_sym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (this->mat_iterative_unsym.GetM() > 0)
          this->mat_iterative_unsym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else if (this->matCSR_iterative_unsym.GetM() > 0)
          this->matCSR_iterative_unsym.ApplySor(trans, x, r, omega, nb_iter, stage_ssor);
	else
	  matrix_not_stored = true;
      }
    
    if (matrix_not_stored)
      cout << "Not implemented " << endl;        
  }


  //! copies the real part of A in B
  template<class T>
  void CopyReal(const FemMatrixFreeClass_Base<complex<T> >& A, FemMatrixFreeClass_Base<T>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    
    CopyReal(A.mat_boundary_sym, B.mat_boundary_sym);
    CopyReal(A.mat_boundary_unsym, B.mat_boundary_unsym);
    CopyReal(A.mat_iterative_sym, B.mat_iterative_sym);
    CopyReal(A.mat_iterative_unsym, B.mat_iterative_unsym);
    CopyReal(A.matCSR_boundary_sym, B.matCSR_boundary_sym);
    CopyReal(A.matCSR_boundary_unsym, B.matCSR_boundary_unsym);
    CopyReal(A.matCSR_iterative_sym, B.matCSR_iterative_sym);
    CopyReal(A.matCSR_iterative_unsym, B.matCSR_iterative_unsym);

    B.row_scaling = A.row_scaling;
    B.column_scaling = A.column_scaling;
    B.row_scale = A.row_scale;
    B.col_scale = A.col_scale;

    B.ignore_dirichlet_dof = A.ignore_dirichlet_dof;
    B.sym_matrix = A.sym_matrix;
    B.sym_dirichlet_condition = A.sym_dirichlet_condition;

    B.coef_dirichlet = A.coef_dirichlet;
    B.prod_level = A.prod_level;

    CopyReal(A.column_dirichlet, B.column_dirichlet);
    
    B.nature_matrix.SetCoefMass(realpart(A.nature_matrix.GetCoefMass()));
    B.nature_matrix.SetCoefDamping(realpart(A.nature_matrix.GetCoefDamping()));
    B.nature_matrix.SetCoefStiffness(realpart(A.nature_matrix.GetCoefStiffness()));
    
  }
  

  template<class T>
  void CopyReal(const FemMatrixFreeClass_Base<T>& A, FemMatrixFreeClass_Base<T>& B)
  {
    cout << "not implemented" << endl;
    abort();
    //B = A;
  }


  //! returns a type of matrix adapted to a symmetric matrix
  int FemMassMatrix::GetSymmetryType(int type)
  {
    switch (type)
      {
      case BLOCK_DIAGONAL_UNSYM : return BLOCK_DIAGONAL;
      case MATRIX_SPARSE_UNSYM : return MATRIX_SPARSE;
      case DIAG_SPARSE_UNSYM : return DIAG_SPARSE;
      case BLOCK_DIAG_SPARSE_UNSYM : return BLOCK_DIAG_SPARSE;
      }

    return type;	
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_FEM_MATRIX_FREE_CLASS_CXX
#endif
