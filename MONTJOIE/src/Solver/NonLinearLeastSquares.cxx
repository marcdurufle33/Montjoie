#ifndef MONTJOIE_FILE_NON_LINEAR_LEAST_SQUARES_CXX

namespace Montjoie
{

  //! Implementation of Levenberg-Marquardt algorithm by storing the jacobian
  template<class T, class VectorSol, class VectorRhs, class MatrixJac, class MatrixAAt>
  T SolveLeastSquaresLvm(VirtualLeastSquaresFunction<T, VectorSol, VectorRhs, MatrixJac, MatrixAAt>& fct,
			 VectorSol& xsol, Real_wp epsilon1, Real_wp epsilon2,
                         unsigned nb_max_iter, int print_level, Real_wp tau)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = fct.GetM();
    int n = fct.GetN();
    if (m < n)
      {
        //cout << "The number of equations must be larger than the number of unknowns" << endl;
        //abort();
      }

    if (epsilon2 < Real_wp(0))
      epsilon2 = epsilon1;
    
    T one, zero;
    SetComplexOne(one); SetComplexZero(zero);
    
    VectorRhs fvec(m);
    VectorSol g(n), x_new(n), h(n);
    MatrixJac fjac(m, n);
    fct.InitJacobian(fjac);
    
    // initial guess
    fct.FindInitGuess(xsol);
    
    // f and J are evaluated
    fct.EvaluateJacobian(xsol, fvec, fjac);

    Treal nu = 2.0;
    // storing A = J^T J
    MatrixAAt A(n, n), B(n, n);;
    A.Zero();
    
    MltAdd(one, SeldonTrans, fjac, SeldonNoTrans, fjac, zero, A);
    
    // we compute g = J^T f(x)
    Mlt(SeldonTrans, fjac, fvec, g);
    
    Treal residual = fct.GetNorm2_Rhs(fvec);
    Treal norm_g = Norm2(g);
    if (norm_g <= epsilon1)
      {
	// initial guess already good
	return residual;
      }
    
    Treal mu = tau*MaxAbs(A);
    unsigned num_iter = 0;
    while (num_iter < nb_max_iter)
      {	
	num_iter++;

	// we solve (A + mu I) h = -g
	B = A;
	for (int i = 0; i < A.GetM(); i++)
	  B(i, i) += mu;
	
	h = g;
	GetCholesky(B);	
	SolveCholesky(SeldonNoTrans, B, h);
	SolveCholesky(SeldonTrans, B, h);
	Mlt(-one, h);
	
	Treal norm_h = Norm2(h);
	if (norm_h <= epsilon2*(Norm2(xsol) + epsilon2))
	  {
	    // we have reached a local minimum
	    return residual;
	  }
	
	if (norm_h >= (residual+epsilon2) / epsilon_machine)
	  {
	    // singular matrix ?
	    //return residual;
	  }
	
	// next iterate x_new = x + h
	x_new = xsol + h;
	
	fct.EvaluateF(x_new, fvec);
	Treal new_residual = fct.GetNorm2_Rhs(fvec);
	Treal dF = (square(residual) - square(new_residual));
	// dL = <h, mu h - g>
	Treal dL = 0;
	for (int i = 0; i < n; i++)
	  dL += h(i)*(mu*h(i) - g(i));
	
	Treal rho = dF / dL;
	if (rho > 0)
	  {
	    // acceptable step
	    xsol = x_new;
	    fct.GiveIterate(num_iter, xsol, new_residual);
	    
	    // A is updated and g as well
	    fct.EvaluateJacobian(xsol, fvec, fjac);
	    MltAdd(one, SeldonTrans, fjac, SeldonNoTrans, fjac, zero, A);
	    
	    Mlt(SeldonTrans, fjac, fvec, g);
	    
	    residual = new_residual;
	    norm_g = Norm2(g);
	    if (print_level >= 2)
	      cout << "Residu at iteration " << num_iter << " : " << residual << endl;
	    else if (print_level >= 1)
	      if (num_iter%10 == 0)
		cout << "Residu at iteration " << num_iter << " : " << residual << endl;
	    
	    if (norm_g <= epsilon1)
	      return new_residual;
	    
	    mu *= max(Treal(1)/3, Treal(1) - pow(2*rho-1, 3));
	    nu = 2.0;
	  }
	else
	  {
	    // mu is increased
	    mu *= nu; nu *= Treal(2);
	  }

      }

    return residual;
  }

  
#ifdef SELDON_WITH_MKL
  template<class T>
  T SolveLeastSquaresMkl(VirtualLeastSquaresFunction<T>& fct, Vector<T>& xsol,
                         Real_wp epsilon, unsigned nb_max_iterations)
  {
    int m = fct.GetM();
    int n = fct.GetN();
    if (m < n)
      {
        cout << "The number of equations must be larger than the number of unknowns" << endl;
        abort();
      }
    
    Vector<T> fvec(m);
    Matrix<T, General, ColMajor> fjac(m, n);
    fvec.Zero();
    fjac.Zero();
    
    // initial guess
    fct.FindInitGuess(xsol);
    
    // initialize solver (allocate mamory, set initial values)
    // iter1   in:     maximum number of iterations
    // iter2   in:     maximum number of iterations of calculation of trial-step
    // rs      in:     initial step bound
    _TRNSP_HANDLE_t handle;
    int iter1, iter2;
    double rs = 1.0;
    double eps[6];
    eps[0] = epsilon;
    eps[1] = epsilon;
    eps[2] = epsilon;
    eps[3] = epsilon_machine;
    eps[4] = epsilon_machine;
    eps[5] = epsilon;
    iter1 = nb_max_iterations;
    iter2 = nb_max_iterations/10;
    dtrnlsp_init(&handle, &n, &m, xsol.GetData(), eps, &iter1, &iter2, &rs);
    
    int info[6];
    dtrnlsp_check(&handle, &n, &m, fjac.GetData(), fvec.GetData(), eps, info);

    // main loop
    bool test_loop = true;
    int rci_request;
    while (test_loop)
      {
        dtrnlsp_solve(&handle, fvec.GetData(), fjac.GetData(), &rci_request);
        
        if ((rci_request <= -1) && (rci_request >= -6))
          test_loop = false;
        
        if (rci_request == 1)
          fct.EvaluateF(xsol, fvec);
        else if (rci_request == 2)
          fct.EvaluateJacobian(xsol, fvec, fjac);
      }
    
    int st_cr, nb_iter;
    double r1, r2;
    dtrnlsp_get(&handle, &nb_iter, &st_cr, &r1, &r2);    
    dtrnlsp_delete(&handle);
    
    fct.EvaluateF(xsol, fvec);
    return Norm2(fvec);
  }
#endif

}

#define MONTJOIE_FILE_NON_LINEAR_LEAST_SQUARES_CXX
#endif
