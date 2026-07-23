#ifndef MONTJOIE_FILE_GENERIC_EQUATION_CXX

namespace Montjoie
{  
  
  /************************
   * GenericEquation_Base *
   ************************/
      
    
  //////////////////////////////////
  // PML METHODS, DO NOT OVERLOAD //
    
  
  //! Modification of the tensor C in order to use PML layers
  /*!
    \param[in] vars problem to be solved
    \param[in] num_elem number of the element where C is expressed
    \param[in] jloc local number of the quadrature point
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Cgrad_grad the new C tensor
    This method is used for nodal H^1 finite element
    The pml layers consist in replacing d/dx_k -> 1/(1+i sigma/omega) d/dx_k
    where k is the direction of the layer
    This modification induces a modification of the C tensor.
    DO NOT OVERLOAD THIS FUNCTION (NO NEED)
   */
  template<class T>
  template<class TypeEquation, class NatureMat,
           class Prop3, class Prop4, int p, int q>
  void GenericEquation_Base<T>::
  ApplyPmlGradGrad(const EllipticProblem<TypeEquation>& vars, int num_elem,
		   int jloc, const NatureMat& nat_mat, int ref,
		   TinyMatrix<TinyMatrix<Complex_wp, Prop3, p, p>, Prop4, q, q>& Cgrad_grad)
  {
    int i1 = num_elem - vars.mesh.GetNbElt() + vars.GetNbEltPML();
    typedef typename TypeEquation::Dimension Dimension;
    TinyVector<Complex_wp, Dimension::dim_N> tau;
    Complex_wp prod(1, 0);
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        tau(i) = vars.GetTauPML(i1, jloc, i);
        prod *= tau(i);
      }
    
    prod = Real_wp(1)/prod;
    for (int m = 0; m < Cgrad_grad.GetM(); m++)
      for (int n = 0; n < Cgrad_grad.GetN(); n++)
	for (int i = 0; i < Dimension::dim_N; i++)
	  for (int j = 0; j < Dimension::dim_N; j++)
	    Cgrad_grad(m, n)(i, j) *= tau(i)*tau(j)*prod;
    
  }
  

  //! Modification of the tensors D and E in order to use PML layers
  /*!
    \param[in] vars problem to be solved
    \param[in] num_elem number of the element where C is expressed
    \param[in] jloc local number of the quadrature point
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Dgrad_phi the new C tensor
    \param[out] Ephi_grad the new C tensor
    This method is used for nodal H^1 finite element
    The pml layers consist in replacing d/dx_k -> 1/(1+i sigma/omega) d/dx_k
    where k is the direction of the layer
    This modification induces a modification of the D and E tensors.
    DO NOT OVERLOAD THIS FUNCTION (NO NEED)
   */
  template<class T>
  template<class TypeEquation, class NatureMat, class Prop3, int p, int q>
  void GenericEquation_Base<T>::
  ApplyPmlGradPhi(const EllipticProblem<TypeEquation>& vars,
		  int num_elem, int jloc, const NatureMat& nat_mat, int ref,
		  TinyMatrix<TinyVector<Complex_wp, p>, Prop3, q, q>& Dgrad_phi,
		  TinyMatrix<TinyVector<Complex_wp, p>, Prop3, q, q>& Ephi_grad)
  {
    int i1 = num_elem - vars.mesh.GetNbElt() + vars.GetNbEltPML();
    typedef typename TypeEquation::Dimension Dimension;
    TinyVector<Complex_wp, Dimension::dim_N> tau;
    Complex_wp prod(1, 0);    
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        tau(i) = vars.GetTauPML(i1, jloc, i);
        prod *= tau(i);
      }

    prod = Real_wp(1)/prod;
    for (int m = 0; m < Dgrad_phi.GetM(); m++)
      for (int n = 0; n < Dgrad_phi.GetN(); n++)
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Dgrad_phi(m,n)(i) *= tau(i)*prod;
	    Ephi_grad(m,n)(i) *= tau(i)*prod;
	  }
    
  }
  

  //! Modification of tensor A because of PMLs
  template<class T>
  template<class TypeEquation, class NatureMat, class Prop3, int p>
  void GenericEquation_Base<T>::
  ApplyPmlPhi(const EllipticProblem<TypeEquation>& vars,
              int num_elem, int jloc, const NatureMat& nat_mat, int ref,
              TinyMatrix<Complex_wp, Prop3, p, p>& Aphi)
  {
    int i1 = num_elem - vars.mesh.GetNbElt() + vars.GetNbEltPML();
    typedef typename TypeEquation::Dimension Dimension;
    Complex_wp prod(1, 0);    
    for (int i = 0; i < Dimension::dim_N; i++)
      prod *= vars.GetTauPML(i1, jloc, i);
          
    Aphi *= Real_wp(1)/prod;
  }
  
  
  //! no PML in real numbers
  template<class T>
  template<class TypeEquation,
           class NatureMat, class Prop3, class Prop4, int p, int q>
  void GenericEquation_Base<T>::
  ApplyPmlGradGrad(const EllipticProblem<TypeEquation>& vars, int num_elem,
		    int jloc, const NatureMat& nat_mat, int ref, 
		    TinyMatrix<TinyMatrix<Real_wp, Prop3, p, p>, Prop4, q, q >& Cgrad_grad)
  {
  }
    
  
  //! no PML in real numbers
  template<class T>
  template<class TypeEquation, class NatureMat, class Prop3, int p, int q>
  void GenericEquation_Base<T>::
  ApplyPmlGradPhi(const EllipticProblem<TypeEquation>& vars,
		   int num_elem, int jloc, const NatureMat& nat_mat, int ref,
		   TinyMatrix<TinyVector<Real_wp, p>, Prop3, q, q>& Dgrad,
		   TinyMatrix<TinyVector<Real_wp, p>, Prop3, q, q>& Egrad)
  {
  }

  
  //! no PML in real numbers
  template<class T>
  template<class TypeEquation, class NatureMat, class Prop3, int p>
  void GenericEquation_Base<T>::
  ApplyPmlPhi(const EllipticProblem<TypeEquation>& vars,
              int num_elem, int jloc, const NatureMat& nat_mat, int ref,
              TinyMatrix<Real_wp, Prop3, p, p>& Aphi)
  {
  }
  
  
  // PML METHODS, DO NOT OVERLOAD //
  //////////////////////////////////
  
  
  ///////////////////////////////////////////////
  // MODIFICATION OF TENSOR ACCORDING GEOMETRY //
  // DO NOT OVERLOAD                           //
  
  
  //! Returns tensor C, after change of variables to unit reference    
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff, class TypeElt>
  void GenericEquation_Base<T>::
  GetMassPhiDFiTensor(const EllipticProblem<TypeEquation>& vars,
                      int num_elem, int jloc,
                      const NatureMat& nat_mat, int ref, MatStiff& Amass,
                      bool variable, bool affine, const TypeElt& Fb)
  {
    // we get the tensor independent from the geometry
    MatStiff Aphys;
    TypeEquation::GetTensorMass(vars, num_elem, jloc, nat_mat, ref, Aphys);
    if (vars.InsidePML(num_elem))
      ApplyPmlPhi(vars, num_elem, jloc, nat_mat, ref, Aphys);

    Amass = Aphys;    
    Real_wp jacob_weight(0);
    if (variable)
      {
        if (affine)
          jacob_weight = vars.Glob_jacobian(num_elem)(0)*Fb.WeightsND(jloc);
        else
          jacob_weight = vars.Glob_jacobian(num_elem)(jloc);
      }
    else
      jacob_weight = vars.Glob_jacobian(num_elem)(0);
    
    Amass *= jacob_weight;
  }

  
  //! Returns tensor C, after change of variables to unit reference    
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff, class TypeElt>
  void GenericEquation_Base<T>::
  GetGradGradDFiTensor(const EllipticProblem<TypeEquation>& vars,
		       int num_elem, int jloc,
		       const NatureMat& nat_mat, int ref, MatStiff& Cgrad_grad,
                       bool variable, bool affine, const TypeElt& Fb)
  {
    // we get the tensor non-dependant of the geometry
    MatStiff Cphys;
    TypeEquation::GetGradGradTensor(vars, num_elem, jloc, nat_mat, ref, Cphys);
    
    if (vars.InsidePML(num_elem))
      ApplyPmlGradGrad(vars, num_elem, jloc, nat_mat, ref, Cphys);
    
    Cgrad_grad = Cphys;
    typedef typename TypeEquation::Dimension Dimension;
    typename Dimension::MatrixN_N dfjm1;
    if (vars.Glob_DFjm1(num_elem).GetM() > 1)
      dfjm1 = vars.Glob_DFjm1(num_elem)(jloc); 
    else
      dfjm1 = vars.Glob_DFjm1(num_elem)(0);
    
    Real_wp jacob = Det(dfjm1);
    if (Dimension::dim_N == 3)
      jacob = sqrt(jacob);
    
    Real_wp jacob_weight(0);
    if (variable)
      {
        if (affine)
          jacob_weight = vars.Glob_jacobian(num_elem)(0)*Fb.WeightsND(jloc);
        else
          jacob_weight = vars.Glob_jacobian(num_elem)(jloc);
      }
    else
      jacob_weight = vars.Glob_jacobian(num_elem)(jloc);
    
    Mlt(1.0/jacob, dfjm1);
    // and we apply the matrix DF-1
    for (int i = 0; i < TypeEquation::nb_unknowns; i++)
      for (int m = 0; m < Dimension::dim_N; m++)
	for (int k = 0; k < TypeEquation::nb_unknowns; k++)
	  for (int n = 0; n < Dimension::dim_N; n++)
	    {
	      Cgrad_grad(i,k)(m,n) = T(0);
	      for (int j = 0; j < Dimension::dim_N; j++)
		for (int l = 0; l < Dimension::dim_N; l++)
		  Cgrad_grad(i,k)(m,n) += dfjm1(n,l)*dfjm1(m,j)*Cphys(i,k)(j,l);
	      
	      Cgrad_grad(i,k)(m,n) *= jacob_weight;
	    }
  }
  
  
  //! Returns tensors D and E, after change of variables to unit reference  
  template<class T>
  template<class TypeEquation, class NatureMat, class MatStiff, class TypeElt>
  void GenericEquation_Base<T>::
  GetGradPhiDFiTensor(const EllipticProblem<TypeEquation>& vars, int num_elem,
		      int jloc, const NatureMat& nat_mat, int ref,
		      MatStiff& Dgrad_phi, MatStiff& Ephi_grad,
                      bool variable, bool affine, const TypeElt& Fb)
  {
    // we get the tensor non-dependant of the geometry
    MatStiff Dphys, Ephys;
    TypeEquation::GetGradPhiTensor(vars, num_elem, jloc, nat_mat, ref, Dphys, Ephys);
    
    if (vars.InsidePML(num_elem))
      ApplyPmlGradPhi(vars, num_elem, jloc, nat_mat, ref, Dphys, Ephys);
    
    Dgrad_phi = Dphys; Ephi_grad = Ephys;
    typedef typename TypeEquation::Dimension Dimension;
    typename Dimension::MatrixN_N dfjm1;
    if (vars.Glob_DFjm1(num_elem).GetM() > 1)
      dfjm1 = vars.Glob_DFjm1(num_elem)(jloc);
    else
      dfjm1 = vars.Glob_DFjm1(num_elem)(0);
    
    Real_wp jacob = Det(dfjm1);
    if (Dimension::dim_N == 3)
      jacob = sqrt(jacob);
    
    Real_wp jacob_weight(0);
    if (variable)
      {
        if (affine)
          jacob_weight = vars.Glob_jacobian(num_elem)(0)*Fb.WeightsND(jloc);
        else
          jacob_weight = vars.Glob_jacobian(num_elem)(jloc);
      }
    else
      jacob_weight = vars.Glob_jacobian(num_elem)(jloc);
    
    Real_wp poids = jacob_weight/jacob;
    // and we apply the matrix DF-1
    for (int i = 0; i < TypeEquation::nb_unknowns; i++)
      for (int m = 0; m < Dimension::dim_N; m++)
	for (int k = 0; k < TypeEquation::nb_unknowns; k++)
	  {
	    Dgrad_phi(i,k)(m) = T(0); Ephi_grad(i,k)(m) = T(0);
	    for (int j = 0; j < Dimension::dim_N; j++)
	      {
		Dgrad_phi(i,k)(m) += dfjm1(m,j)*Dphys(i,k)(j);
		Ephi_grad(i,k)(m) += dfjm1(m,j)*Ephys(i,k)(j);
	      }
	    
	    Dgrad_phi(i,k)(m) *= poids;
	    Ephi_grad(i,k)(m) *= poids;
	  }
    
  }

  
  // MODIFICATION OF TENSOR ACCORDING GEOMETRY //
  // DO NOT OVERLOAD                           //
  ///////////////////////////////////////////////
  
} // namespace Montjoie

#define MONTJOIE_FILE_GENERIC_EQUATION_CXX
#endif


