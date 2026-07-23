#ifndef MONTJOIE_FILE_ASSEMBLE_MATRIX_HXX

namespace Montjoie
{

  //! base class handling static condensation
  template<class T>
  class CondensationBlockSolver_Base
  {
    bool treatment_stiff_inside;
    int num_elem_condensed, nb_elt_condensed, num_elem_global;
    
  public:
    CondensationBlockSolver_Base();
    
    virtual ~CondensationBlockSolver_Base();
    virtual void ModifyElementaryMatrix(int i, IVect& num_ddl, VirtualMatrix<T>& mat_interac,
					const GlobalGenericMatrix<T>& nat_mat);

    void SetTreatmentStiffnessInside(bool t);
    bool TreatInsideStiffness() const;
    void SetElementNumber(int n, int ng);
    int GetCondensedElementNumber() const;
    int GetGlobalElementNumber() const;
    int GetNbCondensedElt() const;
    void SetNbCondensedElt(int n);

    virtual inline size_t GetMemorySize() const { return 0; }
    
  };

  
  //! generic "Nature_Matrix" class containing mass and stiffness coefficients
  /*!
    We are often considering the following evolution system :
    M d^2 U/dt^2 + S dU/dt + K U = 0
    Therefore, for some time-schemes, it can be useful to construct matrix :
    alpha M + sigma S + beta K
    where alpha is the "mass coefficient", beta the "stiffness coefficient"
    and sigma the "damping coefficient"
   */
  template<class T>
  class GlobalGenericMatrix
  {
  protected :
    T coef_mass; //!< mass coefficient
    T coef_stiff; //!< stiffness coefficient
    T coef_sigma; //!< damping coefficient
    
  public :
    typedef T value_type;
    typedef GlobalGenericMatrix<typename ClassComplexType<T>::Tcplx> Nature_Matrix_Complex;
    
    GlobalGenericMatrix();
    GlobalGenericMatrix(const T& a, const T& b, const T& c);
    
    const T& GetCoefMass() const;
    const T& GetCoefStiffness() const;
    const T& GetCoefDamping() const;

    void SetCoefMass(const T& alpha);
    void SetCoefStiffness(const T& alpha);
    void SetCoefDamping(const T& alpha);
    
  };

  
  //! Class used to select elements for which matrix vector product is performed
  class MatrixVectorProductLevel
  {
  public :
    //! matrix-vector product done on all elements or specific elements
    /*!
      ALL_LEVELS : matrix-vector product for all elements
      LVL_PML : matrix-vector product for elements in PML elements
      LVL_NOPML : matrix-vector product for all elements except PML elements
    */
    enum {ALL_LEVELS = -1, LVL_PML = -2, LVL_NOPML = -3};
    
  protected :
    //! element numbers for each level
    Vector<IVect>* num_elem;
    //! number of elements not in PML
    int nb_elt_outside_PML;
    //! number of elements in the mesh
    int nb_elt;
    //! current element to treat for the matrix vector product
    int current_elem;
    //! the number of elements already treated
    int nelem;
    //! current level to use for the matrix vector product
    int current_level;
    
  public :
    MatrixVectorProductLevel();
    
    // Inline methods
    void SetLevel();        
    Vector<IVect>& GetLevelArray();
    void SetLevelArray(Vector<IVect>& lvl);
    
    int GetNbElt() const;
    int GetNbElt(int lvl) const;
    int GetElementNumber(int i) const;
    int GetLocalElementNumber() const;
    
    // Other methods
    bool TreatElement(int i);
    void SetLevel(int level);
    void SetNbElt(int nelt, int nelt_pml);
    
    size_t GetMemorySize() const;
    
  };
  
  
  //! display informations about GlobalGenericMatrix object
  template<class T>
  ostream& operator <<(ostream& out, const GlobalGenericMatrix<T>& var);

  template<class Complexe>
  void AssembleMatrix(VirtualMatrix<Complexe>& mat_sp, VirtualMatrix<Complexe>& mat_elem,
		      const GlobalGenericMatrix<Complexe>& nat_mat,
		      VarComputationProblem_Base& vars,
		      CondensationBlockSolver_Base<Complexe>& solver,
		      int offset_row, int offset_col);
  
  template<class T>
  void GetRowSum(VectReal_wp& diagonal_scale_left,
		 FemMatrixFreeClass_Base<T>& mat);

  template<class T>
  void GetRowColSum(VectReal_wp& diagonal_scale_left,
		    VectReal_wp& diagonal_scale_right,
		    FemMatrixFreeClass_Base<T>& mat);

  template<class T>
  void ScaleMatrix(FemMatrixFreeClass_Base<T>& mat,
		   VectReal_wp& diagonal_scale_left,
		   VectReal_wp& diagonal_scale_right);
  
  template<class T, class Allocator, class Allocator2, class Allocator3>
  void SubMltAdd(const Matrix<complex<T>, Symmetric, ArrayRowSymComplexSparse, Allocator>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<complex<T>, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<complex<T>, Vect_Full, Allocator3>& Ah_Ones);
  
  template<class T, class Allocator, class Allocator2, class Allocator3>
  void SubMltAdd(const Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<T, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<T, Vect_Full, Allocator3>& Ah_Ones);
  
  template<class T, class T2, class Allocator1, class Allocator2, class Allocator3>
  void SubMltAdd(const Matrix<T, General, ArrayRowSparse, Allocator1>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<T2, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<T2, Vect_Full, Allocator3>& Ah_Ones);
  
  template<class T, class T2, class Allocator1, class Allocator2, class Allocator3>
  void SubMltAdd(const class_SeldonTrans& Trans,
		 const Matrix<T, General, ArrayRowSparse, Allocator1>& A,
		 int nb_edge_, const IVect& num_edge, const IVect& IndexEdge,
		 const Vector<T2, Vect_Full, Allocator2>& Ph_Ones,
		 Vector<T2, Vect_Full, Allocator3>& Ah_Ones);
  
} // namespace Montjoie

#define MONTJOIE_FILE_ASSEMBLE_MATRIX_HXX
#endif
