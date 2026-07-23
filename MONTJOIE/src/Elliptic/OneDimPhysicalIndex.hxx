#ifndef MONTJOIE_FILE_ONE_DIM_PHYSICAL_INDEX_HXX

namespace Montjoie
{

#ifdef MONTJOIE_WITH_ONE_DIM
  class VarProblem_1D;
#endif
  
  //! base class for problems involving physical indices (rho, mu, epsilon, sigma, etc)
  class VarPhysicalProblem
  {
  public :    
    //! numbering of each referenced media
    IVect ElementRho;
    
  };
  
  
  //! Base class for a 1-D variable parameter
  template<class T>
  class VariableParameter1D_Base
  {
  protected:
    // number of degrees of freedom for the 1-D parameter
    int nb_points_grid_c;
    // position of sub-intervals
    VectReal_wp points;
    // position of points used to project a function
    VectReal_wp all_points;
    // order of approximation
    int order;
    
  public:
    // Inline methods
    VariableParameter1D_Base();
    virtual ~VariableParameter1D_Base();
    
    int GetNbDof() const;
    
    // Other methods
    static VariableParameter1D_Base<T>* GetNewApproximation(const Real_wp& a, const Real_wp& b, const Vector<string>& param);
    
    int FindInterval(const Real_wp& xinterp) const;
    
    // virtual methods
    virtual size_t GetMemorySize() const = 0;
    virtual void ComputeValuesPhiRef(const Real_wp& x, VectReal_wp& phi) const = 0;
    virtual void Init(const Vector<T>& cn) = 0;
    virtual void Mlt(const T& coef) = 0;
    
    virtual T Evaluate(const Real_wp& x) const = 0;
    virtual void EvaluateDerivative(const Real_wp& x, T&, T&) const = 0;
    virtual void EvaluateSecondDerivative(const Real_wp& x, T&, T&, T&) const;
    virtual VariableParameter1D_Base<T>* GetDuplicate() const = 0;
    
    virtual void ProjectFunction(T (&f)(Real_wp), Vector<T>& projF);
    
  };
  
  
  //! variable index defined with a cubic spline
  template<class T>
  class SplineInterpolationIndex : public VariableParameter1D_Base<T>
  {
  protected:
    SplineInterpolation<T> spline;
    
  public:
    // Inline methods
    size_t GetMemorySize() const;
    void ComputeValuesPhiRef(const Real_wp& x, VectReal_wp& phi) const;
    void Init(const Vector<T>& cn);
    void Mlt(const T& coef);
    
    T Evaluate(const Real_wp& x) const;
    void EvaluateDerivative(const Real_wp& x, T&, T&) const;
    void EvaluateSecondDerivative(const Real_wp& x, T&, T&, T&) const;
    VariableParameter1D_Base<T>* GetDuplicate() const;
    
    // Other methods
    void Init(const string& data_file, T& cte_infty);
    void Init(const Vector<Real_wp>& x, const Vector<T>& y, T& cte,
	      bool first_order = false);

    void SetApproximation(const VectReal_wp& points_reg, int r);
    
  };


  //! variable index defined with a B-spline
  template<class T>
  class BSplineInterpolationIndex : public VariableParameter1D_Base<T>
  {
  protected:
    BSplineInterpolation<T> bspline;
    VectReal_wp knots;
    
  public:
    // Inline methods
    size_t GetMemorySize() const;
    void ComputeValuesPhiRef(const Real_wp& x, VectReal_wp& phi) const;
    void Init(const Vector<T>& cn);
    void Mlt(const T& coef);
    
    T Evaluate(const Real_wp& x) const;
    void EvaluateDerivative(const Real_wp& x, T&, T&) const;
    VariableParameter1D_Base<T>* GetDuplicate() const;
    
    // Other methods
    void Init(const string& data_file, int r, T& cte_infty);
    void ProjectFunction(T (&f)(Real_wp), Vector<T>& projF);
    
    void SetApproximation(const Vector<string>& param, const VectReal_wp& points_reg);
    
  };


  //! variable index defined with piecewise polynomials
  template<class T>
  class PolynomialInterpolationIndex : public VariableParameter1D_Base<T>
  {
  protected:
    //! if true, the coefficients are discontinuous
    bool discontinuous;
    //! number of dofs per element
    int offset_per_elt;
    // values of the indices on the interpolation points
    Vector<T> yn;
    // basis functions on the unit interval
    Globatto<Real_wp> lob_basis;

  public:
    PolynomialInterpolationIndex();

    size_t GetMemorySize() const;
    void ComputeValuesPhiRef(const Real_wp& x, VectReal_wp& phi) const;
    void Init(const Vector<T>& cn);
    void Mlt(const T& coef);
    
    T Evaluate(const Real_wp& x) const;
    void EvaluateDerivative(const Real_wp& x, T&, T&) const;
    
    VariableParameter1D_Base<T>* GetDuplicate() const;

    void SetApproximation(const Vector<string>& param, const VectReal_wp& points_reg);
    
  };

  
  //! varying index in 1-D 
  template<class T>
  class PhysicalVaryingMedia<Dimension1, T>
  {
  protected:
    //! type of field (CONSTANT, VARIABLE or USER)
    int type;
    //! coefficient on quadrature points
    Vector<T> eval_coef, grad_coef;
    //! constant coefficient for a constant field
    T cte_coef;
    //! offset to use for accessing to the array eval_coef
    Vector<int> offset_eval;
    //! if true the gradient is also computed
    bool compute_grad;
    //! pointer to the object representing the variable index
    VariableParameter1D_Base<T>* index; 
    
  public:    
    enum {CONSTANT, VARIABLE, USER};
    
    PhysicalVaryingMedia();
    ~PhysicalVaryingMedia();

    PhysicalVaryingMedia(const PhysicalVaryingMedia<Dimension1, T>&);
    PhysicalVaryingMedia<Dimension1, T>& operator=(const PhysicalVaryingMedia<Dimension1, T>&);
    
    // inline methods
    void SetConstant(const T& coef);     
    const T& GetConstant() const;    
    void Mlt(const T& coef);
    
    void EnableComputationGradient(bool enable_grad = true);
    const T& GetCoefficient(const VarPhysicalProblem& var, int i, int j) const;
    void GetCoefGradient(const VarPhysicalProblem& var, int i, int j, T&, T&) const;
    
    bool IsVarying() const;
    int GetNbPoints() const;
    int GetOffsetPoint(const VarPhysicalProblem& var, int i) const;
    
    // other methods        
    size_t GetMemorySize() const;
    void SetInputData(const VectString& param, int& nb);

#ifdef MONTJOIE_WITH_ONE_DIM
    void InitArray(int ref, const VarProblem_1D& var);
#endif
    
    void ComputeCoefficient(const VarPhysicalProblem& var, int i, int j, const Real_wp& x);
    const T EvaluateCoefficient(const Real_wp& x);
    
    void SetUserIndex(const Vector<T>& user_coef, const T& cte);
    
  };


}

#define MONTJOIE_FILE_ONE_DIM_PHYSICAL_INDEX_HXX
#endif

