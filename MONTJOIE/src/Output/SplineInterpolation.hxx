#ifndef MONTJOIE_FILE_SPLINE_INTERPOLATION_HXX

namespace Montjoie
{

  template<class T>
  class SplineInterpolation
  {
  protected :
    //! x and y values
    Vector<Real_wp> x; Vector<T> y;
    //! y'' as evaluted by cubic spline algorithm
    Vector<T> ypp; VectReal_wp coef_h;
    //! threshold
    Real_wp threshold;
    //! using first-order functions ?
    bool first;

    void ComputeTridiagonalMatrix(const VectReal_wp& h,
				  TinyBandMatrix<Real_wp, 1>& A) const;
    
  public :
    SplineInterpolation();
    
    void SetThreshold(const Real_wp& eps);
    const Vector<T>& GetSecondDerivatives() const;

    size_t GetMemorySize() const;    
    void Init(const Vector<Real_wp>& xdiv, const Vector<T>& ydiv, bool first_ = false);

    void FindInterval(const Real_wp& xinterp, int& a, int& b) const;
    
    T Evaluate(const Real_wp& xinterp) const;
    void EvaluateDerivative(const Real_wp& xinterp, T& yinterp,
                            T& dy_interp, T& d2y_interp) const;
    
    void EvaluateFctBasis(const Real_wp& point, VectReal_wp & vec_phi) const;
    //void FillPolynomialBasis(Vector<UnivariatePolynomial<T> >& P);
    
    void Mlt(const T& coef);
    void Clear();
    
  };


  class BSplineInterpolation_Base
  {
  protected:
    //! order of approximation
    int order;
    //! number of basis functions
    int nb_basis_fct;
    //! knots
    Vector<Real_wp> t;
    //! temporary vector that will store the value of basis functions
    mutable Vector<Real_wp> value_phi, value_dphi;
    //! threshold
    Real_wp threshold;
    
  public:
    BSplineInterpolation_Base();
    
    int GetOrder() const;
    int GetNbDof() const;
    const VectReal_wp& GetPoints() const;
    
    void SetThreshold(const Real_wp& eps);
    void FindInterval(const Real_wp& xinterp, int& a, int& b) const;
    
    int EvaluateFctBasis(const Real_wp& point, VectReal_wp& vec_phi,
			  bool init_vec = true) const;

    int EvaluateDeriveFctBasis(const Real_wp& point, VectReal_wp& vec_phi,
			       VectReal_wp& der_phi, bool init_vec = true) const;
    
    void Init(const Vector<Real_wp>& xdiv, int r,
              const Vector<Real_wp>& knots, Matrix<Real_wp>& VDM, VectReal_wp& tau);
    
  };


  template<class T>
  class BSplineInterpolation : public BSplineInterpolation_Base
  {
  protected:
    Vector<T> weights;
    
  public:
    const Vector<T>& GetWeights() const;
    
    size_t GetMemorySize() const;

    void Init(const Vector<Real_wp>& xdiv, const Vector<T>& ydiv, int r,
              const Vector<Real_wp>& knots);
    
    void SetData(const VectReal_wp& knots, const Vector<T>& weights, int r);

    T Evaluate(const Real_wp& point) const;
    void EvaluateDerivative(const Real_wp& point, T& y, T& dy) const;

    void Mlt(const T& coef);
    void Clear();
    
  };
  
}

#define MONTJOIE_FILE_SPLINE_INTERPOLATION_HXX
#endif
