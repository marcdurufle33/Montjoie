#ifndef MONTJOIE_FILE_GENERIC_EQUATION_INLINE_CXX

namespace Montjoie
{
  
  /************************
   * GenericEquation_Base *
   ************************/
  
  
  //! default constructor
  template<class T>
  inline GenericEquation_Base<T>::GenericEquation_Base()
  {
  }

  
  //! returns true if the global finite element matrix is symmetric
  template<class T>
  inline bool GenericEquation_Base<T>::SymmetricGlobalMatrix()
  {
    return false; 
  }

  
  //! returns true if the elementary finite element matrix is symmetric
  template<class T>
  inline bool GenericEquation_Base<T>::SymmetricElementaryMatrix()
  {
    return false;
  }
  
  
  //! modification of equation parameters with data file
  /*template<class T>
  inline void GenericEquation_Base<T>::
  SetInputData(const string& keyword, const Vector<string>& param)
  {
  }*/
  
  
  //! if true, DF_i^*-1 will be stored 
  template<class T>
  inline bool GenericEquation_Base<T>::ComputeDFjm1()
  {
    return true;
  }
  
  
  /////////////////
  // MASS MATRIX //
  
  
  //! computation of geometrical transformations that need to be computed 
  //! for mass and stiffness matrix
  /*!
    \param[inout] var problem to be considered
    \param[in] i element number
    \param[in] Fb finite element class
    This method is used for any type of element
  */
  template<class T> template<class TypeEquation, class Dimension>
  inline void GenericEquation_Base<T>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
		      int i, const ElementReference_Dim<Dimension>& Fb) 
  {
  }
  
  
  //! returns the matrix M_K, in the integral \f$ \int M_K u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] num_elem number of the element where M_K needs to be evaluated
    \param[in] jloc number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] mass matrix to be provided
    this method is used for nodal finite element
    The integral can be decomposed by the \f$ \sum_{m,n} \int mass(m,n) u^n v^m \f$
    where u^n is the n-th unknown of the vector u
    By default, M_K is set to zero matrix.
    You will have to overload method if you want to define a new equation
   */
  template<class T>
  template<class TypeEquation, class T0, class MatStiff>
  inline void GenericEquation_Base<T>
  ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
		  int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
		  int ref, MatStiff& mass)
  {
    mass.Fill(0);
  }

  
  //! Applying the tensor M_K to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] U unknown vector U
    \param[out] V result to M_K u
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    This method is used only for nodal H1 finite element method
    The user has to return V = M_K U, with \f$ v^m =  \sum M_K^{m,n} u^n \f$
    You will have to overload method if you want to define
    a new equation solved by H^1 nodal finite element method.
  */
  template<class T>
  template<class TypeEquation, class T0, class Vector1>
  inline void GenericEquation_Base<T>
  ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var,
		    int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
		    int ref, Vector1& U, Vector1& V)
  {
    V.Fill(0);
  }
  
  
  // MASS MATRIX //
  /////////////////
  
  
  ///////////////////////////////////////////////
  // METHODS TO OVERLOAD FOR H1 TRANSFORMATION //
  
  
  //! what are the unknowns/fct tests whose derivatives have to be computed ?
  /*!
    \param[in] vars considered problem to be solved
    \param[in] nat_mat mass and stiffness coefficients
    \param[out] unknown_to_derive if unknown_to_derive(i) is true, the unknown i is derivated
    \param[out] fct_test_to_derive if fct_to_derive(i) is true, the test function i is derivated
  */
  template<class T>
  template<class TypeEquation, class NatureMat, class Vector1>
  inline void GenericEquation_Base<T>
  ::GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
			const NatureMat& nat_mat,
			Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
  {
    // default  : all unknowns and test functions are derived
    unknown_to_derive.Fill(true);
    fct_test_to_derive.Fill(true);
  }

  
  //! returns the tensor C, in the integral \f$ \int C \nabla u \nabla v \f$
  /*!
    \param[in] vars given problem
    \param[in] num_elem number of the element where C needs to be evaluated
    \param[in] jloc number of the local quadrature point in the element
    \param[out] Cgrad_grad tensor to be provided
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    The integral can be decomposed by the \f$ \sum_{i,j,m,n} \int Cgrad_grad(m,n)(i,j)
    \frac{\partial u^n}{\partial x_j} \frac{\partial v^m}{\partial x_i} \f$
    where u^n is the n-th unknown of the vector u
    By default, Cgrad_grad is set to zero
    You will have to overload method if you want to define a new equation
   */
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff>
  inline void GenericEquation_Base<T>
  ::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars, int num_elem,
		      int jloc, const NatureMat& nat_mat, int ref, MatStiff& Cgrad_grad)
  {
    FillZero(Cgrad_grad);
  }
  
    
  //! returns the tensors D and E
  //! in the integrals \f$ \int D grad(u) v + \int E u grad(v) \f$
  /*!
    \param[in] vars given problem
    \param[in] num_elem number of the element where C needs to be evaluated
    \param[in] jloc number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Dgrad_phi tensor D to be provided
    \param[out] Ephi_grad tensor E to be provided
    The integral can be rewrittent as \f$ \sum_{m,n,i} \int E(m,n)(i) 
    \frac{\partial u^n}{\partial x_i} v^m 
    + \int D(m,n)(i) u^n \frac{\partial v^m}{\partial x_i} \f$
    where u^n is the n-th unknown of the vector u
    By default, the two tensors are equal to zero
    These tensors are usually associated with the stiffness matrix
    You will have to overload method if you want to define a new equation
   */
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff>
  inline void GenericEquation_Base<T>
  ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
		     int num_elem, int jloc, const NatureMat& nat_mat, int ref,
		     MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
  }
  
  
  //! Applying the tensor C to grad(u)
    /*!
      \param[in] var problem to be solved
      \param[in] i element number
      \param[in] j quadrature point number inside the element
      \param[in] nat_mat mass and stiffness coefficients
      \param[in] ref reference of the physical domain
      \param[in] dU gradient of the unknown vector U
      \param[out] dV result to C grad(u)
      The user has to return dV = C dU, with \f$ \frac{\partial v^m}{\partial x_i}
      =  \sum C^{m,n}_{i,j} \frac{\partial u^n}{\partial x_j} \f$
      You will have to overload method if you want to define a new equation
     */
  template<class T>
  template<class TypeEquation, class NatureMat, class Vector1>
  inline void GenericEquation_Base<T>
  ::ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var, int i, int j,
			 const NatureMat& nat_mat, int ref, Vector1& dU, Vector1& dV)
  {
    FillZero(dV);
  }
    
  
  //! Applying the tensor E to grad(u)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dU gradient of the unknown vector U
    \param[out] V result to E grad(u)
    The user has to return V = E dU,
    with \f$ v^m =  \sum E^{m,n}_i \frac{\partial u^n}{\partial x_i} \f$
    You will have to overload method if you want to define a new equation
  */
  template<class T>
  template<class TypeEquation, class NatureMat, class Vector1, class Vector2>
  inline void GenericEquation_Base<T>
  ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var, int i, int j,
			 const NatureMat& nat_mat, int ref, Vector1& dU, Vector2& V)
  {
    V.Fill(0);
  }
    
  
  //! Applying the tensor D to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result D u
    The user has to return dV = D U,
    with \f$ \frac{\partial v^m}{\partial x_i} =  \sum D^{m,n}_i u^n \f$
    You will have to overload method if you want to define a new equation
  */
  template<class T>
  template<class TypeEquation, class NatureMat, class Vector1, class Vector2>
  inline void GenericEquation_Base<T>
  ::ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var, int i, int j,
		       const NatureMat& nat_mat, int ref, Vector1& U, Vector2& dV)
  {
    FillZero(dV);
  }
  
  
  // METHODS TO OVERLOAD FOR H1 TRANSFORMATION //
  ///////////////////////////////////////////////
    
  
    
  /////////////////////////////////////////////////
  // METHODS TO OVERLOAD FOR BOUNDARY CONDITIONS //
  
  
  //! computation of the impedance matrix N
  /*!
    \param[out] Nabc matrix N
    \param[in] ref reference of the boundary
    \param[in] normale normale 
    \param[in] iquad element number
    \param[in] k quadrature point number on the element
    \param[in] ref_d reference of the physical domain
    \param[in] vars considered problem to be solved
    \param[in] Fb finite element class
    N is associated to the term \f$ \int_\Gamma N u \varphi \f$
  */
  template<class T>
  template<class Matrix1, class TypeEquation, class NatureMat, class Dimension>
  inline void GenericEquation_Base<T>::
  GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale, int ref, int iquad, int k,
          const NatureMat& nat_mat, int ref_d,
	  const EllipticProblem<TypeEquation>& vars, const ElementReference_Dim<Dimension>& Fb)
  {
    Nabc.Fill(T(0));
  }
  
  
  //! multiplication the boundary condition matrix N
  /*!
    \param[in] normale normale n_i to the boundary 
    \param[in] ref reference of the boundary
    \param[in,out] Un vector you have to apply N
    \param[in] num_elem1 element number
    \param[in] k quadrature point number on the element
    \param[in] ref_d reference of the physical domain
    \param[in] vars considered problem to be solved
    \param[in] Fb finite element class
    N is associated to the term \f$ \int_\Gamma N u \varphi \f$
  */
  template<class T>
  template<class Vector1, class TypeEquation, class NatureMat, class Dimension>
  inline void GenericEquation_Base<T>::
  MltNabc(const typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
	  int num_elem1, int k, const NatureMat& nat_mat, int ref_d,
	  const EllipticProblem<TypeEquation>& vars, const ElementReference_Dim<Dimension>& Fb)
  {
    Un.Fill(T(0));
  }
  
  
  // METHODS TO OVERLOAD FOR BOUNDARY CONDITIONS //
  /////////////////////////////////////////////////
  
  
  ////////////////////////////////////////////////////////
  // METHODS TO OVERLOAD FOR PENALIZATION/STABILIZATION //
  
  
  //! Computation of the "penalization" matrix 
  /*!
    \param[in] Nabc the penalization matrix
    \param[in] normale normale n_i to the boundary
    \param[in] iquad element number
    \param[in] k quadrature point number on the element
    \param[in] vars considered problem to be solved
  */
  template<class T>
  template<class Matrix1, class TypeEquation, class NatureMat, class Dimension>
  inline void GenericEquation_Base<T>::
  GetPenalDG(Matrix1& Nabc, const typename Dimension::R_N& normale,
	     int iquad, int k, int nf, const NatureMat& nat_mat, int ref, int ref2,
	     const EllipticProblem<TypeEquation>& vars,
	     const ElementReference_Dim<Dimension>& Fb)
  {
    Nabc.Fill(T(0));
  }
  
  
  //! Multiplication by the "penalization" matrix 
  /*!
    \param[in] normale normale n_i to the boundary
    \param[in] Vn vector to be multiplied
    \param[out] Un output vector
    \param[in] i element number
    \param[in] n quadrature point number on the element
    \param[in] vars considered problem to be solved
  */
  template<class T>
  template<class Vector1, class TypeEquation, class NatureMat, class Dimension>
  inline void GenericEquation_Base<T>::
  MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector1& Un, 
	     int i, int k, int nf, const NatureMat& nat_mat, int ref, int ref2,
	     const EllipticProblem<TypeEquation>& vars,
	     const ElementReference_Dim<Dimension>& Fb)
  {
    Un.Fill(T(0));
  }
  
  
  // METHODS TO OVERLOAD FOR PENALIZATION/STABILIZATION //
  ////////////////////////////////////////////////////////  
  

} // namespace Montjoie

#define MONTJOIE_FILE_GENERIC_EQUATION_INLINE_CXX
#endif


