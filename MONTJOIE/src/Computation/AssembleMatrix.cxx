#ifndef MONTJOIE_FILE_ASSEMBLE_MATRIX_CXX

namespace Montjoie
{

  //! default constructor
  template<class T>
  CondensationBlockSolver_Base<T>::CondensationBlockSolver_Base()
  {
    treatment_stiff_inside = true;
    num_elem_condensed = -1;
    num_elem_global = -1;
    nb_elt_condensed = 0;
  }
  
  
  /************************
   * GlobalGenericMatrix *
   ************************/
  
  
  //! default constructor
  template<class T>
  GlobalGenericMatrix<T>::GlobalGenericMatrix()
  {
    SetComplexOne(coef_mass);
    SetComplexOne(coef_stiff);
    SetComplexOne(coef_sigma);
  }

  
  //! constructor with mass, damping and stiffness coefficient
  template<class T>
  GlobalGenericMatrix<T>::GlobalGenericMatrix(const T& a, const T& b, const T& c)
  {
    coef_mass = a; coef_sigma = b; coef_stiff = c;
  }
  
  
  //! display informations about GlobalGenericMatrix object
  template<class T>
  ostream& operator <<(ostream& out, const GlobalGenericMatrix<T>& var)
  {
    out << "Mass and stiffness coefficients : " << var.GetCoefMass()
        << ", " << var.GetCoefStiffness() << endl;
    
    out << "Damping coefficient : " << var.GetCoefDamping() << endl;
    
    return out;
  }
  

  /****************************
   * MatrixVectorProductLevel *
   ****************************/
  
  
  //! default constructor
  MatrixVectorProductLevel::MatrixVectorProductLevel()
  {
    current_elem = 0;
    nb_elt_outside_PML = 0;
    nb_elt = 0;
    nelem = 0;
    current_level = ALL_LEVELS;
    num_elem = NULL;
  }
  
  
  //! returns true if the element i is involved in the matrix vector product
  bool MatrixVectorProductLevel::TreatElement(int i)
  {
    switch (current_level)
      {
      case ALL_LEVELS : return true;
      case LVL_PML: return (i >= nb_elt_outside_PML);
      case LVL_NOPML : return (i < nb_elt_outside_PML);
      default:
        if (i == current_elem)
          {
	    // we select next element
            nelem++;
	    if (nelem >= (*num_elem)(current_level).GetM())
	      current_elem = nb_elt;
	    else
	      current_elem = (*num_elem)(current_level)(nelem);
	    
            return true;
          }
      }
    
    return false;
  }
  
  
  //! matrix vector product will be performed with elements of a given level
  void MatrixVectorProductLevel::SetLevel(int level)
  {
    nelem = 0;
    current_elem = -1;
    if (level >= 0)
      {
        if (level >= (*num_elem).GetM())
          {
            cout << "Unknown level " << endl;
            abort();
          }
        else
          {
            if ((*num_elem)(level).GetM() > 0)
              current_elem = (*num_elem)(level)(0);
          }
      }
    
    current_level = level;
  }
  
  
  //! sets the number of elements in elements in PML
  void MatrixVectorProductLevel::SetNbElt(int nelt, int nelt_pml)
  {
    nb_elt = nelt;
    nb_elt_outside_PML = nelt - nelt_pml;
  }
  

  //! returns the memory used by the objects in bytes
  size_t MatrixVectorProductLevel::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    return taille;
  }
 
  
  /*************
   * Functions *
   *************/


  //! Assembles finite element matrix 
  /*!
    This methods calls function ComputeElementaryMatrix to generate values
    In the case of static condensation, internal nodes are eliminated
    \param[in,out] mat_sp assembled sparse matrix, which can contain previous entries
    \param[in] mat_elem elementary matrix
    \param[in] nat_mat type of matrix to assemble
    \param[inout] vars given problem computing the elementary matrices
    \param[inout] solver how blocks are treated when a static condensation is used
   */
  template<class Complexe>
  void AssembleMatrix(VirtualMatrix<Complexe>& mat_sp, VirtualMatrix<Complexe>& mat_interac,
		      const GlobalGenericMatrix<Complexe>& nat_mat,
		      VarComputationProblem_Base& vars,
		      CondensationBlockSolver_Base<Complexe>& solver,
                      int offset_row, int offset_col)
  {
    // for null matrix, we allocate it
    int nb_rows = vars.GetNbRows();
    if (mat_sp.GetM() == 0)
      mat_sp.Reallocate(nb_rows, nb_rows);
    
    int num_dof, nb_dof_loc;
    IVect num_ddl;
    int c_interac; IVect col_interac;
    Vector<Complexe> val_interac, col_mat;
    Vector<Complexe, VectSparse> vec_mat;

    // loop over elements
    Vector<int> permut;
    Real_wp threshold_matrix = vars.GetThresholdMatrix();
    int nb_elt = vars.GetNbElt();
    solver.SetNbCondensedElt(nb_elt);
    for (int i = 0; i < nb_elt; i++) 
      {
	solver.SetElementNumber(i, i);
	
	// we compute elementary matrix for element i
	vars.ComputeElementaryMatrix(i, num_ddl, mat_interac, solver, nat_mat);
	solver.ModifyElementaryMatrix(i, num_ddl, mat_interac, nat_mat);

	nb_dof_loc = num_ddl.GetM();
        
        if (mat_interac.IsSparse())
          {
            // loop over local dofs
            for (int j = 0; j < nb_dof_loc; j++)
              if (num_ddl(j) >= 0)
                {
                  num_dof = num_ddl(j);
                  mat_interac.GetSparseRow(j, vec_mat);
                  
                  c_interac = 0;
                  for (int k0 = 0; k0 < vec_mat.GetM(); k0++)
                    if (num_ddl(vec_mat.Index(k0)) >= 0)
                      if (abs(vec_mat.Value(k0)) > threshold_matrix)
                        c_interac++;
                  
                  col_interac.Reallocate(c_interac);
                  val_interac.Reallocate(c_interac);
                  c_interac = 0;
                  for (int k0 = 0; k0 < vec_mat.GetM(); k0++)
                    if (num_ddl(vec_mat.Index(k0)) >= 0)
                      if (abs(vec_mat.Value(k0)) > threshold_matrix)
                        {
                          int k = vec_mat.Index(k0);
                          col_interac(c_interac) = offset_col + num_ddl(k);
                          val_interac(c_interac) = vec_mat.Value(k0);
                          c_interac++;
                        }
                  
                  // we add interactions
                  mat_sp.AddInteractionRow(offset_row + num_dof, c_interac,
                                           col_interac, val_interac);
                }
          }
        else
          {
            col_interac.Reallocate(nb_dof_loc);
            val_interac.Reallocate(nb_dof_loc);
            col_mat.Reallocate(nb_dof_loc);
            
            // dofs are sorted
            permut.Reallocate(nb_dof_loc);
            permut.Fill();
            Sort(num_ddl, permut);
            
            // loop over local dofs
            for (int j = 0; j < nb_dof_loc; j++)
              if (num_ddl(j) >= 0)
                {
                  num_dof = num_ddl(j); 
                  mat_interac.GetDenseRow(permut(j), col_mat);
                
                  // for each degree of freedom, we add interactions
                  // to the global matrix
                  c_interac = 0;
                  for (int k = 0; k < nb_dof_loc; k++)
                    if (num_ddl(k) >= 0)
                      if (abs(col_mat(permut(k))) > threshold_matrix)
                        {
                          col_interac(c_interac) = offset_col + num_ddl(k);
                          val_interac(c_interac) = col_mat(permut(k));
                          c_interac++;
                        }
                  
                  // we add interactions
                  mat_sp.AddInteractionRow(offset_row + num_dof, c_interac,
                                           col_interac, val_interac, true);
                }
	  }        
      }
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    \param[out]  diagonal_scale_left vector containing the sum of each row
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowSum(VectReal_wp& diagonal_scale_left,
		 FemMatrixFreeClass_Base<T>& mat)
  {
    bool matrix_not_stored = false;
    if (mat.IsSymmetric())
      {
	if (mat.mat_iterative_sym.GetM() > 0)
	  GetRowSum(diagonal_scale_left, mat.mat_iterative_sym);
	else if (mat.matCSR_iterative_sym.GetM() > 0)
	  {
	    abort();
	  }
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (mat.mat_iterative_unsym.GetM() > 0)
	  GetRowSum(diagonal_scale_left, mat.mat_iterative_unsym);
	else if (mat.matCSR_iterative_unsym.GetM() > 0)
	  {
	    abort();
	  }
	else
	  matrix_not_stored = true;
      }
    
    if (matrix_not_stored)
      {
        diagonal_scale_left.Reallocate(mat.GetM());
        diagonal_scale_left.Fill(0);
        
	if (mat.IsSymmetric())
	  {
	    if (mat.mat_boundary_sym.GetM() > 0)
	      GetRowSum(diagonal_scale_left, mat.mat_boundary_sym);
	    else if (mat.matCSR_boundary_sym.GetM() > 0)
	      {
		abort();
	      }
	  }
	else
	  {
	    if (mat.mat_boundary_unsym.GetM() > 0)
	      GetRowSum(diagonal_scale_left, mat.mat_boundary_unsym);
	    else if (mat.matCSR_boundary_unsym.GetM() > 0)
	      {
		abort();
	      }
	  }
	
        mat.AddRowSum(diagonal_scale_left);
      }
  }
  
  
  //! For each row and column of the matrix, computation of the sum of absolute values
  /*!
    \param[out]  diagonal_scale_left vector containing the sum of each row
    \param[out]  diagonal_scale_right vector containing the sum of each column
    \param[in] mat given matrix
   */
  template<class T>
  void GetRowColSum(VectReal_wp& diagonal_scale_left,
		    VectReal_wp& diagonal_scale_right,
		    FemMatrixFreeClass_Base<T>& mat)
  {
    bool matrix_not_stored = false;
    if (mat.IsSymmetric())
      {
	if (mat.mat_iterative_sym.GetM() > 0)
	  GetRowColSum(diagonal_scale_left, diagonal_scale_right, mat.mat_iterative_sym);
	else if (mat.matCSR_iterative_sym.GetM() > 0)
	  {
	    abort();
	  }
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (mat.mat_iterative_unsym.GetM() > 0)
	  GetRowColSum(diagonal_scale_left, diagonal_scale_right,
		       mat.mat_iterative_unsym);
	else if (mat.matCSR_iterative_unsym.GetM() > 0)
	  {
	    abort();
	  }
	else
	  matrix_not_stored = true;
      }
    
    if (matrix_not_stored)
      {
        diagonal_scale_right.Reallocate(mat.GetM());
        diagonal_scale_right.Fill(0);
        
        diagonal_scale_left.Reallocate(mat.GetM());
        diagonal_scale_left.Fill(0);
        
	if (mat.IsSymmetric())
	  {
	    if (mat.mat_boundary_sym.GetM() > 0)
	      GetRowColSum(diagonal_scale_left, diagonal_scale_right,
			   mat.mat_boundary_sym);
	    else if (mat.matCSR_boundary_sym.GetM() > 0)
	      {
		abort();
	      }
	  }
	else
	  {
	    if (mat.mat_boundary_unsym.GetM() > 0)
	      GetRowColSum(diagonal_scale_left, diagonal_scale_right,
			   mat.mat_boundary_unsym);
	    else if (mat.matCSR_boundary_unsym.GetM() > 0)
	      {
		abort();
	      }
	  }
        
        mat.AddRowColSum(diagonal_scale_left, diagonal_scale_right);        
      }
  }
  

  //! Matrix is scaled (rows and columns)
  /*!    
    \param[in] mat finite element matrix
    \param[in] diagonal_scale_left scaling vector for rows
    \param[in] diagonal_scale_right scaling vector for columns
   */
  template<class T>
  void ScaleMatrix(FemMatrixFreeClass_Base<T>& mat,
		   VectReal_wp& diagonal_scale_left,
		   VectReal_wp& diagonal_scale_right)
  {
    bool matrix_not_stored = false;
    if (mat.IsSymmetric())
      {
	if (mat.mat_iterative_sym.GetM() > 0)
	  ScaleMatrix(mat.mat_iterative_sym, diagonal_scale_left, diagonal_scale_right);
	else if (mat.matCSR_iterative_sym.GetM() > 0)
	  {
	    abort();
	  }
	else
	  matrix_not_stored = true;
      }
    else
      {
	if (mat.mat_iterative_unsym.GetM() > 0)
	  ScaleMatrix(mat.mat_iterative_unsym, diagonal_scale_left, diagonal_scale_right);
	else if (mat.matCSR_iterative_unsym.GetM() > 0)
	  {
	    abort();
	  }
	else
	  matrix_not_stored = true;
      }
    
    if (matrix_not_stored)      
      mat.SetScaling(diagonal_scale_left, diagonal_scale_right);
    
  }
  
  
  //! matrix-vector product with a matrix extracted from A
  /*!
    \param[in] A given matrix
    \param[in] nb_edge_ number of columns of the extracted matrix
    \param[in] num_edge column numbers of the extracted matrix
    \param[in] IndexEdge if IndexEdge(i) = -1, the row i is not modified
    \param[in] Ph_Ones vector to multiply
    \param[out] Ah_Ones result
   */
  template<class T, class Allocator, class Allocator2, class Allocator3>
  void SubMltAdd(const Matrix<complex<T>, Symmetric, ArrayRowSymComplexSparse, Allocator>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<complex<T>, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<complex<T>, Vect_Full, Allocator3>& Ah_Ones)
  {
    complex<T> value;
    for (int e = 0; e < nb_edge_; e++)
      {
	int edge = num_edge(e);
	value = 0.0;
	for (int n = 0; n < A.GetRealRowSize(edge); n++)
	  {
	    int col = A.IndexReal(edge, n);
	    if (IndexEdge(col) != -1)
	      {
		value = A.ValueReal(edge, n);
		Ah_Ones(edge) += value*Ph_Ones(col);
		if (col != edge)
		  Ah_Ones(col) += value*Ph_Ones(edge);
	      }
	  }
	for (int n = 0; n < A.GetImagRowSize(edge); n++)
	  {
	    int col = A.IndexImag(edge, n);
	    if (IndexEdge(col) != -1)
	      {
		value = A.ValueImag(edge, n);
		value *= Iwp;
		Ah_Ones(edge) += value*Ph_Ones(col);
		if (col != edge)
		  Ah_Ones(col) += value*Ph_Ones(edge);
	      }
	  }
      }
  }

  
  //! matrix-vector product with a matrix extracted from A
  /*!
    \param[in] A given matrix
    \param[in] nb_edge_ number of columns of the extracted matrix
    \param[in] num_edge column numbers of the extracted matrix
    \param[in] IndexEdge if IndexEdge(i) = -1, the row i is not modified
    \param[in] Ph_Ones vector to multiply
    \param[out] Ah_Ones result
   */
  template<class T, class Allocator, class Allocator2, class Allocator3>
  void SubMltAdd(const Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<T, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<T, Vect_Full, Allocator3>& Ah_Ones)
  {
    T value;
    for (int e = 0; e < nb_edge_; e++)
      {
	int edge = num_edge(e);
	value = 0.0;
	for (int n = 0; n < A.GetRowSize(edge); n++)
	  {
	    int col = A.Index(edge, n);
	    if (IndexEdge(col) != -1)
	      {
		value = A.Value(edge, n);
		Ah_Ones(edge) += value*Ph_Ones(col);
		if (col != edge)
		  Ah_Ones(col) += value*Ph_Ones(edge);
	      }
	  }
      }
  }

  
  //! matrix-vector product with a matrix extracted from A
  /*!
    \param[in] A given matrix
    \param[in] nb_edge_ number of columns of the extracted matrix
    \param[in] num_edge column numbers of the extracted matrix
    \param[in] IndexEdge if IndexEdge(i) = -1, the row i is not modified
    \param[in] Ph_Ones vector to multiply
    \param[out] Ah_Ones result
   */
  template<class T, class T2, class Allocator1, class Allocator2, class Allocator3>
  void SubMltAdd(const Matrix<T, General, ArrayRowSparse, Allocator1>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<T2, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<T2, Vect_Full, Allocator3>& Ah_Ones)
  {
    T2 value;
    for (int e = 0; e < nb_edge_; e++)
      {
	int edge = num_edge(e);
	value = 0.0;
	for (int n = 0; n < A.GetRowSize(edge); n++)
	  {
	    int col = A.Index(edge, n);
	    value = A.Value(edge, n);
	    Ah_Ones(edge) += value*Ph_Ones(col);
	  }
      }
  }

  
  //! matrix-vector product with a matrix extracted from A
  /*!
    \param[in] A given matrix
    \param[in] Trans product with A or transpose(A)
    \param[in] nb_edge_ number of columns of the extracted matrix
    \param[in] num_edge column numbers of the extracted matrix
    \param[in] IndexEdge if IndexEdge(i) = -1, the row i is not modified
    \param[in] Ph_Ones vector to multiply
    \param[out] Ah_Ones result
   */
  template<class T, class T2, class Allocator1, class Allocator2, class Allocator3>
  void SubMltAdd(const class_SeldonTrans& Trans,
		 const Matrix<T, General, ArrayRowSparse, Allocator1>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<T2, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<T2, Vect_Full, Allocator3>& Ah_Ones)
  {
    T2 value;
    for (int e = 0; e < nb_edge_; e++)
      {
	int edge = num_edge(e);
	value = 0.0;
	for (int n = 0; n < A.GetRowSize(edge); n++)
	  {
	    int col = A.Index(edge, n);
	    value = A.Value(edge, n);
	    Ah_Ones(col) += value*Ph_Ones(edge);
	  }
      }
  }
      
}

#define MONTJOIE_FILE_ASSEMBLE_MATRIX_CXX
#endif
