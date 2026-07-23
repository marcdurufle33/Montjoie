#ifndef MONTJOIE_FILE_REISSNER_MINDLIN_CXX

namespace Montjoie
{

  /********************************
   * ReissnerMindlinEquation_Base *
   ********************************/
  
  
  //! Applying the stiffness tensor C to grad(u)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dU gradient of the unknown vector U
    \param[out] dV result to C grad(u)
    The function returns dV = C dU
    The stiffness tensor C is involved in the term \int C grad(u) grad(v) dx
  */
  template<class T>
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void ReissnerMindlinEquation_Base<T>::
  ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var, int i, int j,
                       const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& dU, Vector2& dV)
  {
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    Real_wp delta3_div12 = pow(delta, 3)/12;
    const TinySymmetricTensor<Real_wp, 2>& C = var.ref_tensorC(ref);    
    const TinyMatrix<Real_wp, Symmetric, 2, 2>& G = var.ref_Grot(ref);
    
    T0 Eps_xy = dU(0)(1) + dU(1)(0);
    dV(0)(0) = C(0, 0) * dU(0)(0) + C(0, 1)*Eps_xy + C(0, 2)*dU(1)(1);
    dV(0)(1) = C(1, 0) * dU(0)(0) + C(1, 1)*Eps_xy + C(1, 2)*dU(1)(1);
    dV(1)(0) = dV(0)(1);
    dV(1)(1) = C(2, 0) * dU(0)(0) + C(2, 1)*Eps_xy + C(2, 2)*dU(1)(1);
    
    dV *= delta3_div12;
    
    dV(2)(0) = delta*(G(0, 0)*dU(2)(0) + G(0, 1)*dU(2)(1));
    dV(2)(1) = delta*(G(1, 0)*dU(2)(0) + G(1, 1)*dU(2)(1));
    
    Mlt(nat_mat.GetCoefStiffness(), dV);
  }
    
  
  //! Applies tensor mass M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U vector to be multiplied by M
    \param[out] V result vector V = M U
    the tensor M is involved in term \int M u v dx
  */
  template<class T>
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void ReissnerMindlinEquation_Base<T>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& V)
  {
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    const TinyMatrix<Real_wp, Symmetric, 2, 2>& G = var.ref_Grot(ref);
    
    V(0) = delta*(G(0, 0)*U(0) + G(0, 1)*U(1));
    V(1) = delta*(G(1, 0)*U(0) + G(1, 1)*U(1));
    V(2) = 0;
    Mlt(nat_mat.GetCoefStiffness(), V);
  }
  
  
  //! Applying the tensor D to grad(u)
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] dU gradient of the unknown vector V
    \param[out] V result to D grad(u)
    The tensor D is involved in term \int D grad(u) v dx
  */
  template<class T> 
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void ReissnerMindlinEquation_Base<T>::
  ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var, int i, int j,
                       const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& dU, Vector2& V)
  {
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    const TinyMatrix<Real_wp, Symmetric, 2, 2>& G = var.ref_Grot(ref);

    V(0) = delta*(G(0, 0)*dU(2)(0) + G(0, 1)*dU(2)(1));
    V(1) = delta*(G(1, 0)*dU(2)(0) + G(1, 1)*dU(2)(1));
    V(2) = 0;
    Mlt(nat_mat.GetCoefStiffness(), V);
  }
    

  //! Applying the tensor E to u
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U unknown vector U
    \param[out] dV result E u
    The tensor E is involved in term \int E u grad(v) dx
  */
  template<class T>
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void ReissnerMindlinEquation_Base<T>::
  ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var, int i, int j,
                       const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& dV)
  {
    FillZero(dV);
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    const TinyMatrix<Real_wp, Symmetric, 2, 2>& G = var.ref_Grot(ref);
    
    dV(2)(0) = delta*(G(0, 0)*U(0) + G(0, 1)*U(1));
    dV(2)(1) = delta*(G(1, 0)*U(0) + G(1, 1)*U(1));
    Mlt(nat_mat.GetCoefStiffness(),dV);
  }
  
  
  //! Retrieving the stiffness tensor C
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Cgrad_grad stiffness tensor C
    The stiffness tensor C is involved in the term \int C grad(u) grad(v) dx
  */
  template<class T>
  template<class TypeEquation, class T0, class MatStiff>
  void ReissnerMindlinEquation_Base<T>::
  GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                    int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, MatStiff& Cgrad_grad)
  {
    Real_wp delta = vars.ref_delta(ref).GetCoefficient(vars, i, j);
    Real_wp delta3_div12 = pow(delta, 3)/12;
    
    const TinySymmetricTensor<Real_wp, 2>& C = vars.ref_tensorC(ref);    
    TinyMatrix<Real_wp, Symmetric, 2, 2> G = vars.ref_Grot(ref);
    FillZero(Cgrad_grad);
    
    Cgrad_grad(0, 0)(0, 0) = C(0, 0, 0, 0);
    Cgrad_grad(0, 0)(0, 1) = C(0, 0, 0, 1);
    Cgrad_grad(0, 0)(1, 0) = C(0, 1, 0, 0);
    Cgrad_grad(0, 0)(1, 1) = C(0, 1, 0, 1);

    Cgrad_grad(0, 1)(0, 0) = C(0, 0, 1, 0);
    Cgrad_grad(0, 1)(0, 1) = C(0, 0, 1, 1);
    Cgrad_grad(0, 1)(1, 0) = C(0, 1, 1, 0);
    Cgrad_grad(0, 1)(1, 1) = C(0, 1, 1, 1);
    
    Cgrad_grad(1, 0)(0, 0) = C(1, 0, 0, 0);
    Cgrad_grad(1, 0)(0, 1) = C(1, 0, 0, 1);
    Cgrad_grad(1, 0)(1, 0) = C(1, 1, 0, 0);
    Cgrad_grad(1, 0)(1, 1) = C(1, 1, 0, 1);

    Cgrad_grad(1, 1)(0, 0) = C(1, 0, 1, 0);    
    Cgrad_grad(1, 1)(0, 1) = C(1, 0, 1, 1);    
    Cgrad_grad(1, 1)(1, 0) = C(1, 1, 1, 0);
    Cgrad_grad(1, 1)(1, 1) = C(1, 1, 1, 1);
    
    Cgrad_grad *= delta3_div12;
    
    Cgrad_grad(2, 2)(0, 0) = delta*G(0, 0);
    Cgrad_grad(2, 2)(0, 1) = delta*G(0, 1);
    Cgrad_grad(2, 2)(1, 0) = delta*G(1, 0);
    Cgrad_grad(2, 2)(1, 1) = delta*G(1, 1);
    
    Mlt(nat_mat.GetCoefStiffness(), Cgrad_grad);
    //DISP(Cgrad_grad);
  }  
  
  
  //! Retrieves tensor mass M
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] mass tensor M
    the tensor M is involved in term \int M u v dx
  */
  template<class T>
  template<class TypeEquation, class T0, class Matrix1>
  void ReissnerMindlinEquation_Base<T>::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Matrix1& mass)
  {
    Real_wp delta = vars.ref_delta(ref).GetCoefficient(vars, i, j);
    TinyMatrix<Real_wp, Symmetric, 2, 2> G = vars.ref_Grot(ref);
    
    FillZero(mass);
    mass(0, 0) = delta*G(0, 0);
    mass(0, 1) = delta*G(0, 1);
    mass(1, 0) = delta*G(1, 0);
    mass(1, 1) = delta*G(1, 1);

    //DISP(j); DISP(jacob);  DISP(mass);
    
    Mlt(nat_mat.GetCoefStiffness(), mass);
    //DISP(mass);
  }
  
  
  //! Retrieves tensors D and E
  /*!
    \param[in] var problem to be solved
    \param[in] i element number
    \param[in] j quadrature point number inside the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] Dgrad_phi tensor D
    \param[out] Ephi_grad tensor E
    The tensor D is involved in term \int D grad(u) v dx
    The tensor E is involved in term \int E u grad(v) dx
  */
  template<class T>
  template<class TypeEquation, class T0, class MatStiff>
  void ReissnerMindlinEquation_Base<T>::
  GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars, int i, int j,
                   const GlobalGenericMatrix<T0>& nat_mat,
		   int ref, MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
  {
    Real_wp delta = vars.ref_delta(ref).GetCoefficient(vars, i, j);
    TinyMatrix<Real_wp, Symmetric, 2, 2> G = vars.ref_Grot(ref);
    
    FillZero(Dgrad_phi);
    FillZero(Ephi_grad);
    
    Dgrad_phi(2, 0)(0) = G(0, 0);
    Dgrad_phi(2, 0)(1) = G(0, 1);
    Dgrad_phi(2, 1)(0) = G(1, 0);
    Dgrad_phi(2, 1)(1) = G(1, 1);
    
    Ephi_grad(0, 2)(0) = G(0, 0);
    Ephi_grad(0, 2)(1) = G(0, 1);
    Ephi_grad(1, 2)(0) = G(1, 0);
    Ephi_grad(1, 2)(1) = G(1, 1);
    
    Mlt(nat_mat.GetCoefStiffness()*delta, Dgrad_phi);
    Mlt(nat_mat.GetCoefStiffness()*delta, Ephi_grad);
    //DISP(Ephi_grad); DISP(Dgrad_phi);
  }
  
  
  //! Retrieves tensor mass M
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] mass tensor M
    the tensor M is involved in term \int M u v dx
  */
  template<class TypeEquation, class T0, class Matrix1>
  void ReissnerMindlinEquation::
  GetTensorMass(const EllipticProblem<TypeEquation>& vars,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Matrix1& mass)
  {
    ReissnerMindlinEquation_Base<Real_wp>::GetTensorMass(vars, i, j, nat_mat, ref, mass);
    
    T0 rho = vars.ref_rho(ref).GetCoefficient(vars, i, j);
    Real_wp delta = vars.ref_delta(ref).GetCoefficient(vars, i, j);
    Real_wp delta3_div12 = pow(delta, 3)/12;
    rho *= nat_mat.GetCoefMass();
    mass(2, 2) += rho*delta;
    
    rho *= delta3_div12;
    mass(0, 0) += rho;
    mass(1, 1) += rho;
    
  }
  

  //! Applies tensor mass M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U vector to be multiplied by M
    \param[out] V result vector V = M U
    the tensor M is involved in term \int M u v dx
  */
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void ReissnerMindlinEquation::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& V)
  {
    ReissnerMindlinEquation_Base<Real_wp>::ApplyTensorMass(var, i, j, nat_mat, ref, U, V);
    
    T0 rho = var.ref_rho(ref).GetCoefficient(var, i, j);
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    Real_wp delta3_div12 = pow(delta, 3)/12;
    rho *= nat_mat.GetCoefMass();
    
    V(2) += rho*delta*U(2);
    
    rho *= delta3_div12;
    V(0) += rho*U(0);
    V(1) += rho*U(1);
  }
  
  
  //! Retrieves tensor mass M
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] mass tensor M
    the tensor M is involved in term \int M u v dx
  */
  template<class TypeEquation, class T0, class Matrix1>
  void HarmonicReissnerMindlinEquation::
  GetTensorMass(const EllipticProblem<TypeEquation>& var,
                int i, int j, const GlobalGenericMatrix<T0>& nat_mat, int ref, Matrix1& mass)
  {
    ReissnerMindlinEquation_Base<Complex_wp>::GetTensorMass(var, i, j, nat_mat, ref, mass);
    
    Complex_wp rho = var.ref_rho(ref).GetCoefficient(var, i, j);
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    Real_wp delta3_div12 = pow(delta, 3)/12;
    rho *= -nat_mat.GetCoefMass()*var.GetSquareOmega();
    mass(2, 2) += rho*delta;
    
    rho *= delta3_div12;
    mass(0, 0) += rho;
    mass(1, 1) += rho;
  }
  
  
  //! Applies tensor mass M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] U vector to be multiplied by M
    \param[out] V result vector V = M U
    the tensor M is involved in term \int M u v dx
  */
  template<class TypeEquation, class T0, class Vector1, class Vector2>
  void HarmonicReissnerMindlinEquation::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
                  const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector2& V)
  {
    ReissnerMindlinEquation_Base<Complex_wp>::ApplyTensorMass(var, i, j, nat_mat, ref, U, V);
    
    Complex_wp rho = var.ref_rho(ref).GetCoefficient(var, i, j);
    Real_wp delta = var.ref_delta(ref).GetCoefficient(var, i, j);
    Real_wp delta3_div12 = pow(delta, 3)/12;
    rho *= -var.GetSquareOmega()*nat_mat.GetCoefMass();
    
    V(2) += rho*delta*U(2);
    
    rho *= delta3_div12;
    V(0) += rho*U(0);
    V(1) += rho*U(1);
  }


  /************************
   * ReissnerMindlin_Base *
   ************************/
  
  
  //! default constructor
  template<class TypeEquation>
  ReissnerMindlin_Base::ReissnerMindlin_Base(EllipticProblem<TypeEquation>& var)
    : var_problem(var)
  {
  }
  
  
  //! allocation of arrays containing physical properties
  void ReissnerMindlin_Base::InitIndices(int n)
  {
    ref_rho.Reallocate(n);
    ref_delta.Reallocate(n);
    ref_Eyoung.Reallocate(n);
    ref_nu_poisson.Reallocate(n);
    ref_shear_modulus.Reallocate(n);
    ref_timoshenko.Reallocate(n);
    ref_theta_orthotrope.Reallocate(n);
    for (int i = 0; i < n; i++)
      {
        ref_rho(i).SetIdentity();
        ref_delta(i).SetConstant(0.1);
        ref_Eyoung(i).SetIdentity();
	R2 nu_cte(0.25, 0.25);
        ref_nu_poisson(i).SetConstant(nu_cte);
        ref_shear_modulus(i).SetIdentity();
        ref_timoshenko(i).SetIdentity();
        ref_theta_orthotrope(i).SetConstant(0);
      }
  }
  
  
  //! return number of physical domains
  int ReissnerMindlin_Base::GetNbPhysicalIndices() const
  {
    return ref_rho.GetM();
  }
  
  
  //! modification of physical indexes according to data file
  void ReissnerMindlin_Base::SetIndices(int i, const VectString& parameters)
  {
    int nb = 1;
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }

#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif
    
    if (!parameters(0).compare("ISOTROPE"))
      {
        // MateriauDielec = 1 ISOTROPE rho lambda mu
        ref_rho(i).SetInputData(nb, parameters, parameters(0));
        
        Real_wp lambda, mu;
        PhysicalVaryingMedia<Dimension2, Real_wp> fct_lambda, fct_mu;
        GenericPhysicalIndice<Dimension2, Real_wp>::
          SetInputVaryingMedia(nb, fct_lambda, lambda, parameters);
        
        GenericPhysicalIndice<Dimension2, Real_wp>::
          SetInputVaryingMedia(nb, fct_mu, mu, parameters);
        
        if ((fct_lambda.IsVarying()) || (fct_mu.IsVarying()))
          {
            cout << " Not treated " << endl;
            abort();
          }
        
        // we get E and nu from Lame coefficients lambda and mu
        Real_wp nu = 0.5*lambda/(lambda+mu);
        Real_wp E = (3.0*lambda + 2.0*mu)*mu/(lambda+mu);
        //DISP(lambda); DISP(mu); DISP(E); DISP(nu);
        
        ref_Eyoung(i).SetConstant(R2(E, E));
        ref_nu_poisson(i).SetConstant(R2(nu, nu));
        ref_shear_modulus(i).SetConstant(R3(mu, mu, mu));
        ref_timoshenko(i).SetConstant(R2(1, 1));
      }
    else if (!parameters(0).compare("YOUNG_POISSON"))
      {
        // MateriauDielec = 1 YOUNG_POISSON rho E nu
        ref_rho(i).SetInputData(nb, parameters, parameters(0));
        
        PhysicalVaryingMedia<Dimension2, Real_wp> fct_E, fct_nu;
        Real_wp E, nu;
        GenericPhysicalIndice<Dimension2, Real_wp>::
          SetInputVaryingMedia(nb, fct_E, E, parameters);
        
        GenericPhysicalIndice<Dimension2, Real_wp>::
          SetInputVaryingMedia(nb, fct_nu, nu, parameters);
        
        if ((fct_E.IsVarying()) || (fct_nu.IsVarying()))
          {
            cout << " Not treated " << endl;
            abort();
          }
        
        Real_wp G;
        G = 0.5*E/(1.0+nu);
        ref_Eyoung(i).SetConstant(R2(E, E));
        ref_nu_poisson(i).SetConstant(R2(nu, nu));
        ref_shear_modulus(i).SetConstant(R3(G, G, G));
        ref_timoshenko(i).SetConstant(R2(1, 1));
      }
    else if (!parameters(0).compare("ORTHOTROPE"))
      {
        // MateriauDielec = 1 ORTHOTROPE rho Exx Eyy nu_xy nu_yx Gxy Gxz Gyz Kx Ky theta
        ref_rho(i).SetInputData(nb, parameters, parameters(0));
        ref_Eyoung(i).SetInputData(nb, parameters, parameters(0));
        ref_nu_poisson(i).SetInputData(nb, parameters, parameters(0));
        ref_shear_modulus(i).SetInputData(nb, parameters, parameters(0));
        ref_timoshenko(i).SetInputData(nb, parameters, parameters(0));
        ref_theta_orthotrope(i).SetInputData(nb, parameters, parameters(0));
        
        if (ref_nu_poisson(i).IsVarying())
          {
            cout << "Not handled" << endl;
            abort();
          }
        
        R2 E = ref_Eyoung(i).GetConstant();
        R2 nu = ref_nu_poisson(i).GetConstant();
        
        if (abs(nu(0)*E(1)/E(0) - nu(1)) > 1e-12)
          {
            if (rank_proc == 0)
              {
                cout << "WARNING : " << endl;
                cout << "Poisson's coefficient should satisfy the relation "<< endl;
                cout << "nu_xy E_y / Ex = nu_yx" << endl;
                
                nu(1) = nu(0)*E(1)/E(0);
                cout << "True value of nu_yx is now " << nu(1) << endl;
              }
          }            
      }
  }
  

  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
   */
  void ReissnerMindlin_Base
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    if (i >= ref_rho.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho.GetM() << endl;
        abort();
      }
    
    cout << "Not implemented" << endl;
    abort();
  }


  //! returns the name associated with the physical index num
  string ReissnerMindlin_Base::GetPhysicalIndexName(int m) const
  {
    switch(m)
      {
      case 0: return string("delta");
      }

    return string();
  }  

  
  //! computation of physical properties of the problem
  void ReissnerMindlin_Base::ComputePhysicalCoefficients()
  {    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    ref_tensorC.Reallocate(ref_rho.GetM());
    ref_Grot.Reallocate(ref_rho.GetM());
    TinySymmetricTensor<Real_wp, 2> C;
    TinyMatrix<Real_wp, General, 2, 2> Q, A;
    TinyMatrix<Real_wp, Symmetric, 2, 2> Gmat, Gdiag;
    for (int ref = 1; ref < ref_rho.GetM(); ref++)
      {
        R2 E = ref_Eyoung(ref).GetConstant();
        R2 nu = ref_nu_poisson(ref).GetConstant();
        R3 G = ref_shear_modulus(ref).GetConstant();
        Real_wp theta = ref_theta_orthotrope(ref).GetConstant();
        
        // Timoshenko coefficients
        R2 Ks = ref_timoshenko(ref).GetConstant();
        Real_wp D1 = E(0)/(1.0 - nu(0)*nu(1));
        Real_wp D2 = D1*nu(1);
        Real_wp D3 = E(1)/(1.0 - nu(0)*nu(1));
        
        C.Fill(0);
        C(0, 0, 0, 0) = D1;
        C(0, 0, 1, 1) = D2;
        
        C(0, 1, 0, 1) = G(0);
        C(0, 1, 1, 0) = G(0);
        C(1, 0, 0, 1) = G(0);
        C(1, 0, 1, 0) = G(0);
        
        C(1, 1, 1, 1) = D3;
        C(1, 1, 0, 0) = D2;
        
        Q(0, 0) = cos(theta); Q(0, 1) = -sin(theta);
        Q(1, 0) = sin(theta); Q(1, 1) = Q(0, 0);
        
	//DISP(ref); DISP(C.C);
        C.ApplyRotation(Q);
	//DISP(C.C);
        ref_tensorC(ref) = C;
        
        Gdiag.Fill(0);
        Gdiag(0, 0) = G(1)*Ks(0)*Ks(0);
        Gdiag(1, 1) = G(2)*Ks(1)*Ks(1);
        
        // calcul de Q Gdiag Q^T
        MltTrans(Gdiag, Q, A);
        Mlt(Q, A, Gmat);
        
        ref_Grot(ref) = Gmat;
      }

    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        Vector<bool> IsRefUsed(ref_rho.GetM());
        IsRefUsed.Fill(false);
        for (int i = 0; i < var_problem.mesh.GetNbElt(); i++)
          {
            int ref = var_problem.mesh.Element(i).GetReference();
            IsRefUsed(ref) = true;
          }
        
        int nb_ref_used = 0;
        for (int i = 0; i < ref_rho.GetM(); i++)
          if (IsRefUsed(i))
            nb_ref_used++;
        
        IVect RefUsed(nb_ref_used);
        nb_ref_used = 0;
        for (int i = 0; i < ref_rho.GetM(); i++)
          if (IsRefUsed(i))
            RefUsed(nb_ref_used++) = i;
                
        // retrieving the maximal values for rho and C
        Real_wp rho_max(0), Gmax(0), delta_max(0);
        TinyMatrix<Real_wp, Symmetric, 2, 2> cteG;
        for (int n = 0; n < RefUsed.GetM(); n++)
          {
            int ref = RefUsed(n);
            delta_max = max(delta_max, this->ref_delta(ref).GetConstant());
            rho_max = max(rho_max, abs(this->ref_rho(ref).GetConstant()));
            cteG = this->ref_Grot(ref);
            Gmax = max(Gmax, abs(cteG(0, 0)));
            Gmax = max(Gmax, abs(cteG(0, 1)));
            Gmax = max(Gmax, abs(cteG(1, 1)));
          }

#ifdef SELDON_WITH_MPI
        R3 vec_u(delta_max, rho_max, Gmax), vec_v;
        Vector<int64_t> vec_tmp;
        MpiAllreduce(var_problem.comm_group_mode, &vec_u(0), vec_tmp, &vec_v(0), 3, MPI_MAX);
        delta_max = vec_v(0); rho_max = vec_v(1); Gmax = vec_v(2);
#endif

        
        u_bar = 1.0/(delta_max*Gmax);
        omega_bar = sqrt(Gmax/rho_max);
        Real_wp omegap = var_problem.GetOmega()/omega_bar;
        if ((rank_proc == 0) && (var_problem.print_level >= 2))
          {
            cout << "Value of displacement u used in adimensionalization = " << u_bar << endl;
            cout << "Value of pulsation omega used in adimensionalization = " << omega_bar << endl;
            cout << "New value of pulsation = " << omegap << endl;
          }
        
        Real_wp coef_rho = square(omega_bar)*u_bar;
        Real_wp coef_C = u_bar;
        for (int n = 0; n < RefUsed.GetM(); n++)
          {
            int ref = RefUsed(n);
            this->ref_rho(ref).Mlt(coef_rho);
            this->ref_tensorC(ref) *= coef_C;
            this->ref_Grot(ref) *= coef_C;
          }
        
        var_problem.SetOmega(omegap);
      }
    
  }
  
  
  //! modification of shift due to adimensionalization
  template<class Complexe>
  void ReissnerMindlin_Base
  ::UpdateShiftAdimensionalization(Complexe& shift, Complexe& shift_imag)
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        if (var_problem.FirstOrderFormulation())
          {
            shift /= omega_bar;
            shift_imag /= omega_bar;
          }
        else
          {
            shift /= omega_bar*omega_bar;
            shift_imag /= omega_bar*omega_bar;
          }
      }
  }
    

  //! modification of eigenvalues and eigenvectors due to adimensionalization  
  template<class Complexe>
  void ReissnerMindlin_Base::
  UpdateEigenvaluesAdimensionalization(Vector<Complexe>& lambda,
                                       Vector<Complexe>& lambda_imag,
				       Matrix<Complexe, General, ColMajor>& eigen_vec)
  {
    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_YES)
      {
        Complexe coef; SetComplexOne(coef);
        if (var_problem.FirstOrderFormulation())
          coef = omega_bar;
        else
          coef = omega_bar*omega_bar;
        
        Mlt(coef, lambda);
        Mlt(coef, lambda_imag);
      }
  }

  
  //! additional inputs to read in data file
  void ReissnerMindlin_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    
    if (!description_field.compare("ThicknessPlate"))
      {
        if (!parameters(0).compare("VARIABLE"))
          {
            int ref = to_num<int>(parameters(1));
            int nb = 2;
            ref_delta(ref).SetInputData(nb, parameters, string("ISOTROPE"));
          }
        else
          {
            Real_wp delta = to_num<Real_wp>(parameters(0));
            for (int ref = 0; ref < ref_delta.GetM(); ref++)
              ref_delta(ref).SetConstant(delta);
          }
      }
  }
  
  
  //! returns true if physical indexes are varying inside domain i
  bool ReissnerMindlin_Base::IsVaryingMedia(int i) const
  {
    return ref_delta(i).IsVarying();
  }

  
  //! returns the maximal velocity of waves in a media
  Real_wp ReissnerMindlin_Base::GetVelocityOfMedia(int ref) const
  {
    return 1.0;
  }
  
  
  //! returns the velocity of waves at infinity
  Real_wp ReissnerMindlin_Base::GetVelocityOfInfinity() const
  {
    return 1.0;
  }
  

  //! filling varying indices of the problem
  void ReissnerMindlin_Base::
  GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                    IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
		    Vector<bool>& compute_hess)
  {
    // to do
    int nb = 0;
    for (int i = 0; i < ref_delta.GetM(); i++)
      nb += ref_delta(i).GetNbVaryingMedia();
    
    rho_real.Reallocate(nb);
    num_ref.Reallocate(nb);
    num_index.Reallocate(nb);
    num_component.Reallocate(nb);
    compute_grad.Reallocate(nb);
    compute_hess.Reallocate(nb);
    compute_grad.Fill(false);
    compute_hess.Fill(false);
    nb = 0;
    for (int i = 0; i < ref_delta.GetM(); i++)
      {
        int nb0 = nb;
        ref_delta(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
          {
            num_index(j) = 0;
            num_ref(j) = i;
          }
      }    
  }

  
  /**********************
   * ReissnerMindlin_Eq *
   **********************/

  
  template<class TypeEquation>
  ReissnerMindlin_Eq<TypeEquation>::ReissnerMindlin_Eq()
    : ReissnerMindlin_Base(this->GetLeafClass())
  {
  }
  
  
  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>::InitIndices(int n)
  {
    ReissnerMindlin_Base::InitIndices(n);
  }
  
  
  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    ReissnerMindlin_Base::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    ReissnerMindlin_Base::SetPhysicalIndex(name_media, i, parameters);
  }
   

  template<class TypeEquation>
  string ReissnerMindlin_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return ReissnerMindlin_Base::GetPhysicalIndexName(m);
  }
 
  
  template<class TypeEquation>
  int ReissnerMindlin_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return ReissnerMindlin_Base::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  bool ReissnerMindlin_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return ReissnerMindlin_Base::IsVaryingMedia(i);
  }


  template<class TypeEquation>
  bool ReissnerMindlin_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return ReissnerMindlin_Base::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  Real_wp ReissnerMindlin_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return ReissnerMindlin_Base::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  Real_wp ReissnerMindlin_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return ReissnerMindlin_Base::GetVelocityOfInfinity();
  }
    
  
  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_cplx,
		      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real,
		      IVect& num_ref, IVect& num_index, IVect& num_component,
		      Vector<bool>& compute_grad, Vector<bool>& compute_hess)
  {
    ReissnerMindlin_Base::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
					    compute_grad, compute_hess);
  }


  template<class TypeEquation>    
  void ReissnerMindlin_Eq<TypeEquation>::UpdateShiftAdimensionalization(Real_wp& sr, Real_wp& si)
  {
    ReissnerMindlin_Base::UpdateShiftAdimensionalization(sr, si);
  }
  

  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>
  ::UpdateEigenvaluesAdimensionalization(Vector<Real_wp>& lambda, Vector<Real_wp>& lambda_imag,
					 Matrix<Real_wp, General, ColMajor>& eigen_vec)
  {
    ReissnerMindlin_Base::UpdateEigenvaluesAdimensionalization(lambda, lambda_imag, eigen_vec);
  }


  template<class TypeEquation>    
  void ReissnerMindlin_Eq<TypeEquation>::UpdateShiftAdimensionalization(Complex_wp& sr, Complex_wp& si)
  {
    ReissnerMindlin_Base::UpdateShiftAdimensionalization(sr, si);
  }
  

  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>
  ::UpdateEigenvaluesAdimensionalization(Vector<Complex_wp>& lambda, Vector<Complex_wp>& lambda_imag,
					 Matrix<Complex_wp, General, ColMajor>& eigen_vec)
  {
    ReissnerMindlin_Base::UpdateEigenvaluesAdimensionalization(lambda, lambda_imag, eigen_vec);
  }
  
  
  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();
    ReissnerMindlin_Base::ComputePhysicalCoefficients();
  }
  

  template<class TypeEquation>
  void ReissnerMindlin_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    ReissnerMindlin_Base::SetInputData(description_field, parameters);
  }


  /*******************
   * EllipticProblem *
   *******************/
  
  
  void EllipticProblem<ReissnerMindlinEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));
    
  }


  void EllipticProblem<ReissnerMindlinEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  
  FemMatrixFreeClass<Real_wp, ReissnerMindlinEquation>
  ::FemMatrixFreeClass(const EllipticProblem<ReissnerMindlinEquation>& var_)
    : FemMatrixFreeClass_Eq<Real_wp, ReissnerMindlinEquation>(var_)
  {
  }
  
  
  void FemMatrixFreeClass<Real_wp, ReissnerMindlinEquation>
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Real_wp>& X, Vector<Real_wp>& Y) const
  {
    MltAddVectorH1(Real_wp(1), nat_mat, trans, lvl, *this, 
		   X, Real_wp(1), Y, false);
  }


  void FemMatrixFreeClass<Real_wp, ReissnerMindlinEquation>
  ::MltAddFree(const GlobalGenericMatrix<Real_wp>& nat_mat,
	       const SeldonTranspose& trans, int lvl, 
	       const Vector<Complex_wp>& X, Vector<Complex_wp>& Y) const
  {
    cout << "Not implemented" << endl;
    abort();
  }
  

  void EllipticProblem<HarmonicReissnerMindlinEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Real_wp>& mat_elem,
			    CondensationBlockSolver_Base<Real_wp>&,
			    const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }


  void EllipticProblem<HarmonicReissnerMindlinEquation>
  ::ComputeElementaryMatrix(int i, IVect& num_dof, VirtualMatrix<Complex_wp>& mat_elem,
			    CondensationBlockSolver_Base<Complex_wp>&,
			    const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(i, num_dof, mat_elem, nat_mat, *this,
				      this->GetReferenceElementH1(i));

  }
  
}

#define MONTJOIE_FILE_REISSNER_MINDLIN_CXX
#endif
