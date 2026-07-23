#ifndef SELDON_FILE_TINY_SYMMETRIC_TENSOR_CXX

#include "TinySymmetricTensor.hxx"

namespace Seldon
{
  
  //! initialisation to an isotropic tensor
  /*!
    C_ijkl = lambda (delta_{i, j} + delta_{k, l} - delta_{i, j} delta_{k, l})
             + mu (delta_{i, k} + delta_{j, l})
  */
  template<class T, int m>
  void TinySymmetricTensor<T, m>::FillIsotrope(const T& lambda, const T& mu)
  {
    C.Fill(0);
    for (int i = 0; i < m; i++)
      {
        C(id(i, i), id(i, i)) = lambda + 2.0*mu;
        for (int j = (i+1); j < m; j++)
          {
            C(id(i, j), id(i, j)) = mu;
	    C(id(i, i), id(j, j)) = lambda;
          }
      }
  }


  //! Adds an isotropic tensor
  /*!
    C_ijkl += lambda (delta_{i, j} + delta_{k, l} - delta_{i, j} delta_{k, l})
             + mu (delta_{i, k} + delta_{j, l})
  */
  template<class T, int m>
  void TinySymmetricTensor<T, m>::AddIsotrope(const T& lambda, const T& mu)
  {
    for (int i = 0; i < m; i++)
      {
        C(id(i, i), id(i, i)) += lambda + 2.0*mu;
        for (int j = (i+1); j < m; j++)
          {
            C(id(i, j), id(i, j)) += mu;
	    C(id(i, i), id(j, j)) += lambda;
          }
      }
  }
  
  
  //! Applies a rotation to the tensor
  /*!
    The matrix matQ is assumed to be an orthogonal matrix (representing a rotation)
    For example in 2-D, it should be equal to [cos theta, -sin theta; sin theta, cos theta]
   */
  template<class T, int m> template<class T0>
  void TinySymmetricTensor<T, m>
  ::ApplyRotation(const TinyMatrix<T0, General, m, m>& matQ)
  {
    TinyMatrix<T, Symmetric, N, N> cte_D;
    
    for (int p = 0; p < m; p++)
      for (int q = p; q < m; q++)
        for (int r = 0; r < m; r++)
          for (int s = r; s < m; s++)
            {
              T tmp; SetComplexZero(tmp);
              for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                  for (int k = 0; k < m; k++)
                    for (int l = 0; l < m; l++)
                      tmp += matQ(p, i) * matQ(q, j) * matQ(r, k)
                        * matQ(s, l) * C(id(i, j), id(k, l));
              
              cte_D(id(p, q), id(r, s)) = tmp;
            }
    
    C = cte_D;            
  }
  
  
  //! replaces the tensor by its inverse
  template<class T, int m>
  void TinySymmetricTensor<T, m>::GetInverse()
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal sqrt_two = sqrt(Treal(2));
    
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        for (int k = 0; k < m; k++)
          for (int l = k+1; l < m; l++)
            C(id(i, j), id(k, l)) *= sqrt_two;
    
    for (int k = 0; k < m; k++)
      for (int l = k+1; l < m; l++)
        C(id(k, l), id(k, l)) *= sqrt_two;
    
    Seldon::GetInverse(C);
    
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        for (int k = 0; k < m; k++)
          for (int l = k+1; l < m; l++)
            C(id(i, j), id(k, l)) *= 0.5*sqrt_two;    
    
    for (int k = 0; k < m; k++)
      for (int l = k+1; l < m; l++)
        C(id(k, l), id(k, l)) *= 0.5*sqrt_two;
  }
  
  
  //! Computation of eigenvalues of the tensor
  template<class T, class T0, int p>
  void GetEigenvalues(const TinySymmetricTensor<T, p>& C, TinyVector<T0, p>& k,
		      TinyVector<T, p>& lambda)
  {
    // computing Chrystoffel's tensor
    TinyMatrix<T, Symmetric, p, p> Gamma;
    T val;
    for (int i = 0; i < p; i++)
      for (int j = i; j < p; j++)
        {
          SetComplexZero(val);
          for (int l = 0; l < p; l++)
            for (int m = 0; m < p; m++)
              val += C(i, l, j, m)*k(l)*k(m);
          
          Gamma(i, j) = val;
        }

    GetEigenvalues(Gamma, lambda);
  }
  
} // end namespace

#define SELDON_FILE_TINY_SYMMETRIC_TENSOR_CXX
#endif
