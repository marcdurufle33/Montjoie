#ifndef MONTJOIE_FILE_ONE_DIM_PHYSICAL_INDEX_INLINE_CXX

namespace Montjoie
{

  /****************************
   * VariableParameter1D_Base *
   ****************************/


  //! Default constructor
  template<class T>
  inline VariableParameter1D_Base<T>::VariableParameter1D_Base()
  {
    nb_points_grid_c = 0;
    order = 1;
  }
  

  //! Destructor
  template<class T>
  inline VariableParameter1D_Base<T>::~VariableParameter1D_Base()
  {
  }
  
  
  //! returns the number of degrees of freedom representing the parameter
  template<class T>
  inline int VariableParameter1D_Base<T>::GetNbDof() const
  {
    return nb_points_grid_c;
  }

  
  /****************************
   * SplineInterpolationIndex *
   ****************************/
  

  //! returns the size of the object in bytes
  template<class T>
  inline size_t SplineInterpolationIndex<T>::GetMemorySize() const
  {
    return spline.GetMemorySize();
  }

  
  //! Evaluates basis functions at a given point
  template<class T>
  inline void SplineInterpolationIndex<T>
  ::ComputeValuesPhiRef(const Real_wp& point, VectReal_wp & vec_phi) const
  {
    spline.EvaluateFctBasis(point, vec_phi);
  }
  

  //! Sets values of the index
  template<class T>
  inline void SplineInterpolationIndex<T>::Init(const Vector<T>& cn)
  {
    bool first_order = (this->order == 1);
    spline.Init(this->points, cn, first_order);
  }


  //! Multiplies values of the index by a coefficient
  template<class T>
  inline void SplineInterpolationIndex<T>::Mlt(const T& coef)
  {
    spline.Mlt(coef);
  }
  
  
  //! Evaluates the index at a given point
  template<class T>
  inline T SplineInterpolationIndex<T>::Evaluate(const Real_wp& x) const
  {
    return spline.Evaluate(x);
  }
  
  
  //! Evaluates the index and its derivative at a given point
  template<class T>
  inline void SplineInterpolationIndex<T>::EvaluateDerivative(const Real_wp& x, T& f, T& df) const
  {
    T d2f;
    spline.EvaluateDerivative(x, f, df, d2f);
  }


  //! evaluates f, df and d^2 f
  template<class T>
  inline void SplineInterpolationIndex<T>::EvaluateSecondDerivative(const Real_wp& x, T& f, T& df, T& d2f) const
  {
    spline.EvaluateDerivative(x, f, df, d2f);
  }


  //! returns a duplicate of the variable parameter
  template<class T>
  VariableParameter1D_Base<T>* SplineInterpolationIndex<T>::GetDuplicate() const
  {
    return new SplineInterpolationIndex<T>(*this);
  }
  
  
  /*****************************
   * BSplineInterpolationIndex *
   *****************************/
  

  //! returns the size of the object in bytes
  template<class T>
  inline size_t BSplineInterpolationIndex<T>::GetMemorySize() const
  {
    return bspline.GetMemorySize() + knots.GetMemorySize();
  }

    
  //! Evaluates basis functions at a given point
  template<class T>
  inline void BSplineInterpolationIndex<T>
  ::ComputeValuesPhiRef(const Real_wp& point, VectReal_wp & vec_phi) const
  {
    bspline.EvaluateFctBasis(point, vec_phi);
  }
  

  //! Sets values of the index
  template<class T>
  inline void BSplineInterpolationIndex<T>::Init(const Vector<T>& cn)
  {
    bspline.SetData(knots, cn, this->order);
  }
  

  //! Multiplies values of the index by a coefficient
  template<class T>
  inline void BSplineInterpolationIndex<T>::Mlt(const T& coef)
  {
    bspline.Mlt(coef);
  }
  
  
  //! Evaluates the index at a given point
  template<class T>
  inline T BSplineInterpolationIndex<T>::Evaluate(const Real_wp& x) const
  {
    return bspline.Evaluate(x);
  }
  
  
  //! Evaluates the index and its derivative at a given point
  template<class T>
  inline void BSplineInterpolationIndex<T>::EvaluateDerivative(const Real_wp& x, T& f, T& df) const
  {
    bspline.EvaluateDerivative(x, f, df);
  }


  //! returns a duplicate of the variable parameter
  template<class T>
  VariableParameter1D_Base<T>* BSplineInterpolationIndex<T>::GetDuplicate() const
  {
    return new BSplineInterpolationIndex<T>(*this);
  }
  
  
  /************************
   * PhysicalVaryingMedia *
   ************************/


  //! sets the index to be constant
  template<class T>
  inline void PhysicalVaryingMedia<Dimension1, T>::SetConstant(const T& coef)
  {
    type = CONSTANT;
    cte_coef = coef;
  }
  

  //! returns the constant associated with the current index
  template<class T>
  inline const T& PhysicalVaryingMedia<Dimension1, T>::GetConstant() const
  {
    return cte_coef;
  }


  //! multiplies the index by a constant
  template<class T>
  inline void PhysicalVaryingMedia<Dimension1, T>::Mlt(const T& coef)
  {
    index->Mlt(coef);
    cte_coef *= coef;
    eval_coef *= coef;
    grad_coef *= coef;
  }
  

  //! enables the computation of gradients
  template<class T>
  inline void PhysicalVaryingMedia<Dimension1, T>::EnableComputationGradient(bool enable_grad)
  {
    compute_grad = enable_grad;
  }


  //! returns the coefficient stored at element i and quadrature point j
  template<class T>
  inline const T& PhysicalVaryingMedia<Dimension1, T>
  ::GetCoefficient(const VarPhysicalProblem& var, int i, int j) const
  {
    if (type == CONSTANT)
      return cte_coef;
    
    return eval_coef(offset_eval(var.ElementRho(i)) + j);
  }


  //! retrieves the physical index (in coef) at element i and quadrature point j
  //! and its gradient (in dcoef)
  template<class T>
  inline void PhysicalVaryingMedia<Dimension1, T>
  ::GetCoefGradient(const VarPhysicalProblem& var, int i, int j, T& coef, T& dcoef) const
  {
    if (type == CONSTANT)
      {
	coef = cte_coef;
	SetComplexZero(dcoef);
      }
    else
      {
	coef = eval_coef(offset_eval(var.ElementRho(i)) + j);
	dcoef = grad_coef(offset_eval(var.ElementRho(i)) + j);
      }
  }
  

  //! returns true if the physical index is variable
  template<class T>
  inline bool PhysicalVaryingMedia<Dimension1, T>::IsVarying() const
  {
    if (type == CONSTANT)
      return false;
    
    return true;
  }


  //! returns the number of coefficients stored
  template<class T>
  inline int PhysicalVaryingMedia<Dimension1, T>::GetNbPoints() const
  {
    return offset_eval(offset_eval.GetM()-1);
  }


  //! returns the offset associated with element i
  template<class T>
  inline int PhysicalVaryingMedia<Dimension1, T>
  ::GetOffsetPoint(const VarPhysicalProblem& var, int i) const
  {
    return offset_eval(var.ElementRho(i));
  }
  
}

#define MONTJOIE_FILE_ONE_DIM_PHYSICAL_INDEX_INLINE_CXX
#endif
