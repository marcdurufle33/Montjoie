#ifndef MONTJOIE_FILE_FEM_MATRIX_FREE_CLASS_INLINE_CXX

namespace Montjoie
{
  
  /****************************
   * ExtrapolVariablesProduct *
   ****************************/

  
  //! allocates vectors reserved for exchange (MPI)
  template<class Complexe>
  inline void ExtrapolVariablesProduct_Base<Complexe>
  ::ReallocateExchangeVector(int n)
  {
#ifdef SELDON_WITH_MPI
    request.Reallocate(n);
    Xneighbor.Reallocate(n);
    Xneighbor_tmp.Reallocate(n);
    Xsend.Reallocate(n);
    Xsend_tmp.Reallocate(n);
#endif
  }
  
  
  //! returns the vector storing u on quadrature points of element i
  template<class Complexe, class TypeEquation> template<class Dimension, int type>
  inline TinyVector<Vector<Complexe>, TypeEquation::nb_unknowns>&
  ExtrapolVariablesProductFEM<Complexe, TypeEquation>
  ::GetUnQuad(int i, const ElementReference<Dimension, type>& Fb)
  {
    if (Fb.UseQuadraturePointsForSh())
      return Un_quad(i);
    else
      return Uquad;
  }

  template<class Complexe, class TypeEquation>
  inline size_t ExtrapolVariablesProductFEM<Complexe, TypeEquation>::GetMemorySize() const
  {
    size_t taille = Seldon::GetMemorySize(extrapolU) + Seldon::GetMemorySize(extrapolDun)
      + Seldon::GetMemorySize(extrapolCgradU) + Seldon::GetMemorySize(Un_quad) 
      + Seldon::GetMemorySize(Uloc) + Seldon::GetMemorySize(Uquad) +
      Seldon::GetMemorySize(Uneighbor) + Seldon::GetMemorySize(DunNeighbor)
      + Seldon::GetMemorySize(CgradUNeighbor) + sizeof(unknown_to_derive) 
      + sizeof(fct_test_to_derive);
    
#ifdef SELDON_WITH_MPI
    taille += Seldon::GetMemorySize(this->Xneighbor) + Seldon::GetMemorySize(this->Xsend)
      + Seldon::GetMemorySize(this->Xneighbor_tmp) + Seldon::GetMemorySize(this->Xsend_tmp);
#endif

    return taille;
  }
  

  /***************************
   * FemMatrixFreeClass_Base *
   ***************************/


  //! constructor with class EllipticProblem
  template<class T> template<class TypeEquation>
  inline FemMatrixFreeClass_Base<T>::FemMatrixFreeClass_Base(const EllipticProblem<TypeEquation>& var)
    : var_dir(var), var_problem(var), var_comm(var)
  {
    InitDefaultValues();
  }
  
  
  //! sets coefficient to put on diagonal of Dirichlet rows
  template<class T>
  inline void FemMatrixFreeClass_Base<T>
  ::SetCoefficientDirichlet(const Real_wp& coef)
  {
    coef_dirichlet = coef;
  }
  
  
  //! sets mass, stiffness and damping coefficients (complex coefficients)
  template<class T>
  inline void FemMatrixFreeClass_Base<T>
  ::SetCoefficientMatrix(const GlobalGenericMatrix<T>& nat_mat)
  {
    nature_matrix = nat_mat;
  }


  //! returns mass coefficient
  template<class T>
  inline T FemMatrixFreeClass_Base<T>::GetCoefMass() const
  {
    return nature_matrix.GetCoefMass();
  }
  
  
  //! returns true if the matrix is symmetric
  template<class T>
  inline bool FemMatrixFreeClass_Base<T>::IsSymmetric() const
  {
    return sym_matrix;
  }


  //! returns true if discontinuous formulation is used
  template<class T>
  inline int FemMatrixFreeClass_Base<T>::FormulationDG() const
  {
    return var_problem.FormulationDG();
  }
  

  template<class T>
  inline void FemMatrixFreeClass_Base<T>::SetCondensedSolver(CondensationBlockSolver_Fem<T>* solver)
  {
    condensed_solver = solver;
  }
  
  
  //! returns A_{i, j} 
  template<class T>
  inline const T FemMatrixFreeClass_Base<T>::operator()(int i, int j) const
  {
    // direct element access is forbidden for this kind of matrix
    abort(); return Real_wp(0);
  }
    
  
  //! returns true if Dirichlet condition is skipped
  template<class T>
  inline bool FemMatrixFreeClass_Base<T>::DirichletDofIgnored() const
  {
    return ignore_dirichlet_dof;
  }
  
  
  //! asks to skip Dirichlet condition
  template<class T>
  inline void FemMatrixFreeClass_Base<T>::IgnoreDirichletDof()
  {
    ignore_dirichlet_dof = true; 
  }

  
  //! sets left and/or right scaling
  template<class T>
  inline void FemMatrixFreeClass_Base<T>
  ::SetScaling(VectReal_wp& diagonal_scale_left, VectReal_wp& diagonal_scale_right)
  {
    row_scaling = true;
    column_scaling = true;
    row_scale = &diagonal_scale_left;
    col_scale = &diagonal_scale_right;
  }
  
  
  //! returns true if the matrix A is symmetric
  template<class T>
  inline bool IsSymmetricMatrix(const FemMatrixFreeClass_Base<T>& A)
  {
    return A.IsSymmetric();
  }  

  
  //! returns true if the matrix A is complex
  template<class T>
  inline bool IsComplexMatrix(const FemMatrixFreeClass_Base<T>& A)
  {
    return false;
  }  


  //! returns true if the matrix A is complex
  template<class T>
  inline bool IsComplexMatrix(const FemMatrixFreeClass_Base<complex<T> >& A)
  {
    return true;
  }  
  
  
  /**********************
   * FemMatrixFreeClass *
   **********************/
  

  //! constructor with class EllipticProblem
  template<class T, class TypeEquation>
  inline FemMatrixFreeClass_Eq<T, TypeEquation>::
  FemMatrixFreeClass_Eq(const EllipticProblem<TypeEquation>& var_)
    : FemMatrixFreeClass_Base<T>(var_), var(var_)
  {
  }

  
  //! returns object useful for matrix-vector product
  template<class T, class TypeEquation>
  inline ExtrapolVariablesProductFEM<T, TypeEquation>&
  FemMatrixFreeClass_Eq<T, TypeEquation>::GetExtrapolVariables()
  {
    return var_extra;
  }


  template<class T, class TypeEquation>
  inline size_t FemMatrixFreeClass_Eq<T, TypeEquation>::GetMemorySize() const
  {
    size_t taille = FemMatrixFreeClass_Base<T>::GetMemorySize();
    taille += var_extra.GetMemorySize();
    return taille;
  }


  //! constructor with a given problem
  template<class T, class TypeEquation>
  inline FemMatrixFreeClass<T, TypeEquation>::
  FemMatrixFreeClass(const EllipticProblem<TypeEquation>& var_)
    : FemMatrixFreeClass_Eq<T, TypeEquation>(var_)
  {
  }

  
  
} // namespace Montjoie

#define MONTJOIE_FILE_FEM_MATRIX_FREE_CLASS_INLINE_CXX
#endif
