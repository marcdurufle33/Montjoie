#ifndef MONTJOIE_FILE_GAUSS_LOBATTO_INLINE_CXX

namespace Montjoie
{
  
  /************
   * Globatto *
   ************/
  
  
  //! default constructor
  template<class T>
  inline Globatto<T>::Globatto()
  {
    order = 0; ndloc = 0; nb_points_quadrature = 0;
    threshold = GetPrecisionMachine(T(1));
  }
  
  
  //! returns the order of quadrature formula
  template<class T>
  inline int Globatto<T>::GetOrder() const
  {
    return this->order;
  }
  

  //! returns the order of quadrature formula
  template<class T>
  inline int Globatto<T>::GetGeometryOrder() const
  {
    return this->order;
  }

  
  //! returns the number of quadrature points
  template<class T>
  inline int Globatto<T>::GetNbPointsQuad() const
  {
    return this->nb_points_quadrature;
  }
  
    
  //! returns 1-D weights
  template<class T>
  inline const Vector<T>& Globatto<T>::Weights() const
  {
    return weights_;
  }
  

  //! returns 1-D points
  template<class T>
  inline const Vector<T>& Globatto<T>::Points() const
  {
    return points_;
  }


  //! returns 1-D weights
  template<class T>
  inline const T& Globatto<T>::Weights(int i) const
  {
    return weights_(i);
  }
  

  //! returns 1-D points
  template<class T>
  inline const T& Globatto<T>::Points(int i) const
  {
    return points_(i);
  }
  
  
  //! returns Phi_i(x)
  /*!
    basis functions are constructed on integration points
    \param[in] i basis function number
    \param[in] x argument
   */
  template<class T>
  inline T Globatto<T>::EvaluatePhi(int i, const T& x) const
  {  
    T p(1);
    for (int k = 0; k <= order; k++)
      if (k != i)
	p *= x - points_(k);
    
    p *= cte_phi(i);
    return p;
  }	
  
  //! returns dPhi_i / dx (x)
  /*!
    \param[in] i basis function number
    \param[in] x argument
    basis functions are constructed on integration points
   */
  template<class T>
  inline T Globatto<T>::EvaluatePhiGrad(int i, const T& x) const
  {
    // we compute derivative
    T p(1);
    T s(0);
    for (int l = 0; l <= order; l++)
      if (l != i)
	{
	  p = T(1);
	  for (int k = 0; k <= order; k++)
	    if ((k != i) && (k != l))
	      p *= x - points_(k);
	  
	  s += p;
	}
    
    s *= cte_phi(i);
    return s;
  }
  
  
  //! returns dphi_i(xi_j)
  template<class T>
  inline const T& Globatto<T>::GradPhi(int i, int j) const
  {
    return val_grad(i, j);
  }
  
  
  //! returns matrix contanining dphi_i(xi_j)
  template<class T>
  inline const Matrix<T>& Globatto<T>::GradPhi() const
  {
    return val_grad;
  }
    
} // end namespace

#define MONTJOIE_FILE_GAUSS_LOBATTO_INLINE_CXX
#endif
