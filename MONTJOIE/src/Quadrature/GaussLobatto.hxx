#ifndef MONTJOIE_FILE_GAUSS_LOBATTO_HXX

namespace Montjoie
{
  // declaration of classes to compute integration 1-D points
  
  //! base class for 1-D integration
  template<class T = Real_wp>
  class Globatto
  {
  public :
    // quadrature rules proposed in the object
    enum {QUADRATURE_GAUSS, QUADRATURE_LOBATTO, QUADRATURE_RADAU, QUADRATURE_RADAU_RIGHT,
	  QUADRATURE_JACOBI, QUADRATURE_RADAU_JACOBI, QUADRATURE_LOBATTO_JACOBI,
	  QUADRATURE_GAUSS_SQUARED, QUADRATURE_GAUSS_BLENDED, TCHEBYSHEV};
    
  protected :
    int order; // order of approximation
    int ndloc; // number of quadrature points on a edge
    int nb_points_quadrature; // number of quadrature points over element
    T threshold;
    
    Vector<T> points_; // Points of integration on the edge
    Vector<T> weights_; // Weights of integration
    //! Constants of basis functions
    /*
      cte_phi(i)= \f$ 1/ \prod_{k=1..(r+1) ,k \ne i} ( \xi_k- \xi_i) \f$
    */
    Vector<T> cte_phi;
    
    //!< matrix so that \f$ \frac{\partial \varphi_j}{\partial x}
    //! (\xi_k) = \mbox{valGrad}(j,k) \f$
    Matrix<T> val_grad;  
    
    // method to compute cte_phi from points
    void ComputeFactorPhi();
  
  public :
    static T blending_default;
    
    Globatto();

    size_t GetMemorySize() const;
    int GetOrder() const;
    int GetGeometryOrder() const;
    int GetNbPointsQuad() const;
    
    // Gauss or Gauss-Lobatto quadrature according the value of type_quadrature
    void ConstructQuadrature(int, int type_quadrature = QUADRATURE_GAUSS,
                             T alpha = 0, T beta = 0);
    
    // set array Points equal to Points1D, and recompute cte_phi
    template<class Vector1>
    void AffectPoints(const Vector1&);

    template<class Vector1>
    void AffectWeights(const Vector1&);
    
    // direct access to quadrature formulas
    const T& Points(int i) const;
    const T& Weights(int i) const;

    const Vector<T>& Weights() const;
    const Vector<T>& Points() const;
    
    // returns \phi_i(x) 
    T EvaluatePhi(int i, const T& x) const;  
    
    // returns d\phi_i / dx (x) 
    T EvaluatePhiGrad(int i, const T& x) const ;
    
    const T& GradPhi(int, int) const;
    const Matrix<T>& GradPhi() const;
    
    void ComputeValuesPhiRef(const T&, Vector<T>& phi) const;
    
    // computes val_grad
    void ComputeGradPhi(T eps = 0); 
    
    // prints object on stream out
    template<class T2>
    friend ostream& operator <<(ostream& out, const Globatto<T2>& e);
    
  };
  
  
  //! Interpolation with n subdivisions and r+1 Gauss-Lobatto points on each interval
  /*!
    The number of shape functions is equal to n r + 1
   */
  class SubdivGlobatto
  {
  protected : 
    Vector<VectReal_wp> points, CteG;
    VectReal_wp xpos, points1d;
    
  public :
    
    const VectReal_wp& Points() const;
    const Real_wp& Points(int i) const;
    
    void Init(bool regular, int nb_subdiv, int r,
              int type_quad = Globatto<Real_wp>::QUADRATURE_LOBATTO);
    
    void InitPoints(int nb_subdiv, const VectReal_wp& pts);
    
    int GetOrder() const;
    int GetGeometryOrder() const;
    
    Real_wp EvaluatePhi(int i, const Real_wp& x) const;
    
  };
  
} // end namespace

#define MONTJOIE_FILE_GAUSS_LOBATTO_HXX
#endif


