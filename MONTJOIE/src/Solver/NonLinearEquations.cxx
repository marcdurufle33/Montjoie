#ifndef MONTJOIE_FILE_NONLINEAR_EQUATIONS_CXX

namespace Montjoie

{
  //! solves a non-linear set of equations, knowing gradient of system
  /*!
    \param[in,out] eqns definition of the non-linear system
    \param[in,out] x_sol on input, initial guess, on output computed solution
    \param[out] fvec residual vector
    \param[out] fjac jacobian of system
    \param[out] scale_eqn scaling factor
    \param[in] Control integer parameters
    \param[in] RControl real parameters
   */
  template<class T, class Vector1, class Matrix1>
  int SolveMinpack(NonLinearEquations_Base<T, Vector1, Matrix1>& eqns,
		   Vector1& x_sol, Vector1& fvec, Matrix1& fjac, Vector1& scale_eqn,
		   const IVect& Control, const VectReal_wp& RControl)
  {
    T epsilon_machine, xtol, factor;
    T zero(0), one(1), p1, p5, p001, p0001;
    p1 = T(1)/T(10); p5 = 5*p1; p001 = one/1000; p0001 = p001/10; 
    int max_eval_function, mode_scale;
    int nb_eval_function = 0, nb_eval_jacobian = 0, n, nb_iter;
    int ncsuc, ncfail, nslow1, nslow2;
    IVect iwa; Vector1 wa1, wa2, wa3, wa4, qtf; Vector<T> r;
    T xnorm(0), delta, temp, sum, prered, actred, fnorm, norm_rhs, ratio, pnorm;
  
    n = x_sol.GetM();
    
    // initialization
    wa1 = x_sol; wa2 = x_sol; wa3 = x_sol;
    qtf = x_sol; wa4 = x_sol;
    iwa.Reallocate(n); r.Reallocate(n*(n+1)/2); scale_eqn = x_sol;
    // fjac.Reallocate(n,n); fvec.Reallocate(n);
    wa1.Fill(zero); wa2.Fill(zero); wa3.Fill(zero); wa4.Fill(zero); qtf.Fill(zero);
    iwa.Fill(0); r.Fill(zero); fvec.Fill(zero); fjac.Fill(zero);
    
    epsilon_machine = GetPrecisionMachine(one);
    int info=0;
  
    // we get control parameters
    xtol = RControl(0);
    factor = RControl(1);
    
    max_eval_function = Control(0);
    mode_scale = Control(3);
    //print_level = Control(4);
    
    // check the input parameters for errors.
    if ((n <= 0)||(xtol < 0)||(max_eval_function <= 0)||(factor <= zero))
      {
	return info;
      }
    if (mode_scale==2)
      {
	for (int j=0;j<n;j++)
	  if (scale_eqn(j)<=zero)
	    return info;
      }
  
    // evaluate the function at the starting point
    //  and calculate its norm.
    eqns.EvaluateFunction(x_sol,fvec); nb_eval_function++;
    // DISP(fvec);
    norm_rhs = Norm2(fvec); // DISP(norm_rhs);
  
    // initialize iteration counter and monitors.
  
    nb_iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;
  
    // beginning of the outer loop.
    bool test_loop=true,jeval;
    while (test_loop)
      {
	jeval = true;
	// calculate the jacobian matrix
	eqns.EvaluateJacobian(x_sol,fjac);
	nb_eval_jacobian++; // DISP(fjac);
      
	// compute the qr factorization of the jacobian.
	Compute_QR_Factorisation(fjac,false,iwa,wa1,wa2,wa3);
	// DISP(iwa); DISP(wa1); DISP(wa2); DISP(wa3);
	// on the first iteration and if mode is 1, scale according
	// to the norms of the columns of the initial jacobian.
	if (nb_iter==1)
	  {
	    if (mode_scale==1)
	      {
		for (int j = 0; j < n; j++)
		  {
		    if (wa2(j)==zero)
		      scale_eqn(j)=one;
		    else
		      scale_eqn(j)=wa2(j);
		  }
	      }
	    // on the first iteration, calculate the norm of the scaled x
	    // and initialize the step bound delta.
	    for (int j = 0; j < n; j++)
	      wa3(j)=scale_eqn(j)*x_sol(j);
	  
	    xnorm = Norm2(wa3);
	    delta = factor*xnorm;
	    if (delta == zero)
	      delta = factor;
	  }
      
	// form (q transpose)*fvec and store in qtf.
	qtf = fvec;
	for (int j = 0; j < n; j++)
	  {
	    if (fjac(j,j) != zero)
	      {
		sum = zero;
		for (int i = j; i < n; i++)
		  sum += fjac(i,j)*qtf(i);
	      
		temp = -sum/fjac(j,j);
		for (int i = j; i < n; i++)
		  qtf(i)+=fjac(i,j)*temp;
	      }
	  }
	// DISP(qtf); DISP(delta); DISP(xnorm);
	// copy the triangular factor of the qr factorization into Qr.
	bool sing = false;
	for (int j = 0; j < n; j++)
	  {
	    int l = j, jm1 = j-1;
	    if (jm1 >= 0)
	      for (int i = 0; i <= jm1; i++)
		{
		  r(l) = fjac(i,j);
		  l += n-1-i;
		}
	    r(l) = wa1(j);
	    if (wa1(j) == zero) 
	      sing = true;
	  }
      
	// accumulate the orthogonal factor in fjac
      
	qform(fjac, wa1); // DISP(r); DISP(wa1); DISP(fjac);
      
	// rescale if necessary.
	if (mode_scale == 1)
	  {
	    for (int j=0;j<n;j++)
	      scale_eqn(j) = max(scale_eqn(j),wa2(j));
	  }
	bool test_inner_loop = true;
	// beginning of the inner loop.
	while (test_inner_loop)
	  {
	    // determine the direction p.
	    dogleg(r, scale_eqn, qtf, delta, wa1, wa2, wa3);
	    // store the direction p and x + p. calculate the norm of p.
	    for (int j = 0; j < n; j++)
	      {
		wa1(j) = -wa1(j);
		wa2(j) = x_sol(j) + wa1(j);
		wa3(j) = scale_eqn(j)*wa1(j);
	      }
	    pnorm = Norm2(wa3);
	  
	    // on the first iteration, adjust the initial step bound.
	    if (nb_iter == 1)
	      delta = min(delta, pnorm);
	  
	    // evaluate the function at x + p and calculate its norm.
	    eqns.EvaluateFunction(wa2, wa4); // DISP(wa2); DISP(wa4);
	    nb_eval_function++;
	    fnorm = Norm2(wa4);
	  
	    // compute the scaled actual reduction.
	    actred = -one;
	    if (fnorm < norm_rhs)
	      actred = one-square(fnorm/norm_rhs);
	  
	    // compute the scaled predicted reduction.
	    //DISP(norm_rhs); DISP(fnorm); DISP(actred);
	    int l = 0;
	    for (int i = 0 ; i < n; i++)
	      {
		sum = zero;
		for (int j = i; j < n; j++)
		  {
		    sum += r(l)*wa1(j);
		    l++;
		  }
		wa3(i) = qtf(i) + sum;
	      }
	    temp = Norm2(wa3);
	    prered = zero;
	    if (temp < norm_rhs)
	      prered = one - square(temp/norm_rhs);
	  
	    // compute the ratio of the actual to the predicted
	    // reduction.
	  
	    ratio = zero;
	    if (prered > zero)
	      ratio = actred/prered;
	  
	    // update the step bound.
	  
	    if (ratio < p1)
	      {
		ncsuc = 0;
		ncfail++;
		delta *= p5;
	      }
	    else
	      {
		ncfail = 0;
		ncsuc++;
		if ((ratio >= p5)||(ncsuc > 1))
		  delta = max(delta, pnorm/p5);
	      
		if (abs(ratio-one)<=p1)
		  delta = pnorm/p5;
	      }
	  
	    // DISP(ratio); DISP(pnorm); DISP(prered); DISP(wa1);
	    // DISP(wa2); DISP(wa3);
	    // test for successful iteration
	    if (ratio >= p0001)
	      {
		// successful iteration. update x, fvec, and their norms.
		for (int j = 0; j < n; j++)
		  {
		    x_sol(j) = wa2(j);
		    wa2(j) = scale_eqn(j)*x_sol(j);
		    fvec(j) = wa4(j);
		  }
		xnorm = Norm2(wa2);
		norm_rhs = fnorm;
		nb_iter++;
	      }
	  
	    // DISP(xnorm); DISP(norm_rhs);
	    // determine the progress of the iteration.
	  
	    nslow1++;
	    if (actred >= p001) 
	      nslow1 = 0;
	    if (jeval) 
	      nslow2++;
	    if (actred >= p1) 
	      nslow2 = 0;
	  
	    // test for convergence.
            //DISP(delta); DISP(norm_rhs);
	    if ((delta <= (xtol*xnorm))|| norm_rhs==zero) 
	      {
		info = 0;
		return info;
	      }
	    // tests for termination and stringent tolerances.
	  
	    if (nb_eval_function >= max_eval_function) 
	      return 2;
	    if (p1*max(p1*delta, pnorm)<= epsilon_machine*xnorm) 
	      return 0;
	    if (nslow2==5) 
	      return 4;
	    if (nslow1==10)
	      return 5;
	    //
	    // criterion for recalculating jacobian.
	    // 
	    if (ncfail==2)
	      test_inner_loop=false;
	    else
	      {
		// calculate the rank one modification to the jacobian
		// and update qtf if necessary.
	      
		for (int j = 0; j < n; j++)
		  {
		    sum = zero;
		    for (int i = 0; i < n; i++)
		      sum += fjac(i,j)*wa4(i);
		  
		    wa2(j) = (sum - wa3(j))/pnorm;
		    wa1(j) = scale_eqn(j)*((scale_eqn(j)*wa1(j))/pnorm);
		    if (ratio >= p0001) 
		      qtf(j) = sum;
		  }
	      
		//
		// compute the qr factorization of the updated jacobian.
	      
		// DISP(wa2);DISP(wa1);
		r1updt(r, wa1, wa2, wa3, sing);
		// DISP(r); DISP(wa1); DISP(wa2); DISP(wa3);
		r1mpyq(fjac, wa2, wa3);
		// DISP(fjac); DISP(wa2); DISP(wa3);
		r1mpyq(qtf, wa2, wa3);
		// DISP(qtf); DISP(wa2); DISP(wa3);

		//  end of the inner loop.
	      
		jeval = false;
	      }
	  }
	// end of the outer loop.
      }
  
    return 0;
  } 
  

  //! intermediary function used to solve non-linear system
  template<class Vector, class Matrix>
  void Compute_QR_Factorisation(Matrix& A, bool pivot, IVect& ipvt, 
				Vector& rdiag, Vector& acnorm, Vector& wa)
  {
    //     **********
    //
    //    method qrfac
    //
    //     this subroutine uses householder transformations with column
    //     pivoting (optional) to compute a qr factorization of the
    //     m by n matrix a. that is, qrfac determines an orthogonal
    //     matrix q, a permutation matrix p, and an upper trapezoidal
    //     matrix r with diagonal elements of nonincreasing magnitude,
    //     such that a*p = q*r. the householder transformation for
    //     column k, k = 1,2,...,min(m,n), is of the form
    //
    //                           t
    //           i - (1/u(k))*u*u
    //
    //     where u has zeros in the first k-1 positions. the form of
    //     this transformation and the method of pivoting first
    //     appeared in the corresponding linpack subroutine.
    //
    //     the subroutine statement is
    //
    //     void Compute_QR_Factorisation(A,pivot,ipvt,rdiag,acnorm,wa)
    //
    //     where
    //
    //       A is an m by n array. on input a contains the matrix for
    //         which the qr factorization is to be computed. on output
    //        the strict upper trapezoidal part of a contains the strict
    //         upper trapezoidal part of r, and the lower trapezoidal
    //        part of a contains a factored form of q (the non-trivial
    //         elements of the u vectors described above).
    //
    //       pivot is a logical input variable. if pivot is set true,
    //         then column pivoting is enforced. if pivot is set false,
    //        then no column pivoting is done.
    //
    //       ipvt is an integer output array of length lipvt. ipvt
    //         defines the permutation matrix p such that a*p = q*r.
    //         column j of p is column ipvt(j) of the identity matrix.
    //         if pivot is false, ipvt is not referenced.
    //
    //       rdiag is an output array of length n which contains the
    //         diagonal elements of r.
    //
    //       acnorm is an output array of length n which contains the
    //         norms of the corresponding columns of the input matrix a.
    //         if this information is not needed, then acnorm can coincide
    //         with rdiag.
    //
    //       wa is a work array of length n. if pivot is false, then wa
    //         can coincide with rdiag.
    //
    typedef typename Vector::value_type real;
    real epsilon_machine,zero(0),one(1),temp,sum,p05,ajnorm;
    p05 = real(5)/real(100);
    epsilon_machine = GetPrecisionMachine(one);
    
    int m = A.GetM(), n = A.GetN();
    if (pivot)
      ipvt.Reallocate(n);
  
    // allocate if it is not already done
    // rdiag.Reallocate(n); acnorm.Reallocate(n); wa.Reallocate(n);
    // compute the initial column norms and initialize several arrays.
  
    for (int j = 0; j < n; j++)
      {
	acnorm(j) = Norm2_Column(A, 0, j);
	rdiag(j) = acnorm(j);
	wa(j) = rdiag(j);
	if (pivot) 
	  ipvt(j) = j;
      }
  
    // reduce a to r with householder transformations.
  
    int minmn = min(m,n);
    for (int j = 0; j < minmn; j++)
      {
	if (pivot)
	  {
	    // bring the column of largest norm into the pivot position.
	  
	    int kmax = j;
	    for (int k = j; k < n; k++)
	      if (rdiag(k) > rdiag(kmax)) 
		kmax = k;

	    if (kmax != j) 
	      {
		// we permute columns j and kmax
		for (int i = 0; i < m;i++)
		  {
		    temp = A(i,j);
		    A(i,j) = A(i,kmax);
		    A(i,kmax) = temp;
		  }
		rdiag(kmax) = rdiag(j);
		wa(kmax) = wa(j);
		int k = ipvt(j);
		ipvt(j) = ipvt(kmax);
		ipvt(kmax) = k;
	      }
	  } // end of search pivot

	//        compute the householder transformation to reduce the
	//        j-th column of a to a multiple of the j-th unit vector.
      
	ajnorm = Norm2_Column(A,j,j);
	if (ajnorm!=zero)
	  {
	    if (A(j,j) < zero) 
	      ajnorm = -ajnorm;
	  
	    for (int i = j; i < m;i++)
	      A(i,j) = A(i,j)/ajnorm;
	  
	    A(j,j) += one;

	    //        apply the transformation to the remaining columns
	    //        and update the norms.
	    //
	  
	    int jp1 = j + 1;
	    for (int k = jp1; k < n;k++)
	      {
		sum = zero;
		for (int i = j; i < m;i++)
		  sum += A(i,j)*A(i,k);
	      
		temp = sum/A(j,j);
		for (int i = j; i < m;i++)
		  A(i,k) -= temp*A(i,j);
	      
		if (pivot&& (rdiag(k)!=zero))
		  {
		    temp = A(j,k)/rdiag(k);
		    rdiag(k) = rdiag(k)*sqrt(max(zero,real(one-temp*temp)));
		    if (p05*square(real(rdiag(k)/wa(k)))<= epsilon_machine)
		      {
			rdiag(k) = Norm2_Column(A,jp1,k);
			wa(k) = rdiag(k);
		      }
		  }
	      }
	  }
	rdiag(j) = -ajnorm;
      }
      
    return ;
      
    //     last card of subroutine qrfac.
  }
  
  //! intermediary function used to solve non-linear system
  template<class Vector, class Matrix>
  void QrSolve(Matrix& r, IVect& ipvt, Vector& diag,
	       Vector& qtb, Vector& x, Vector& sdiag, Vector& wa)
  {

    //
    // subroutine qrsolv
    //
    // given an m by n matrix a, an n by n diagonal matrix d,
    // and an m-vector b, the problem is to determine an x which
    // solves the system
  
    // a*x = b ,     d*x = 0 ,
  
    // in the least squares sense.
  
    // this subroutine completes the solution of the problem
    // if it is provided with the necessary information from the
    // qr factorization, with column pivoting, of a. that is, if
    // a*p = q*r, where p is a permutation matrix, q has orthogonal
    // columns, and r is an upper triangular matrix with diagonal
    // elements of nonincreasing magnitude, then qrsolv expects
    // the full upper triangle of r, the permutation matrix p,
    // and the first n components of (q transpose)*b. the system
    // a*x = b, d*x = 0, is then equivalent to
  
    // t                t
    // r*z = q *b ,  p *d*p*z = 0 ,
  
    // where x = p*z. if this system does not have full rank,
    // then a least squares solution is obtained. on output qrsolv
    // also provides an upper triangular matrix s such that
    //
    // t   t               t
    // p *(a *a + d*d)*p = s *s .
  
    // s is computed within qrsolv and may be of separate interest.
  
    // the subroutine statement is
  
    // void QrSolve(r,ipvt,diag,qtb,x,sdiag,wa)
  
    // where
  
    // r is an n by n array. on input the full upper triangle
    // must contain the full upper triangle of the matrix r.
    // on output the full upper triangle is unaltered, and the
    // strict lower triangle contains the strict upper triangle
    // (transposed) of the upper triangular matrix s.
  
    // ipvt is an integer input array of length n which defines the
    // permutation matrix p such that a*p = q*r. column j of p
    // is column ipvt(j) of the identity matrix.
  
    // diag is an input array of length n which must contain the
    // diagonal elements of the matrix d.
  
    // qtb is an input array of length n which must contain the first
    // n elements of the vector (q transpose)*b.
  
    // x is an output array of length n which contains the least
    // squares solution of the system a*x = b, d*x = 0.
  
    // sdiag is an output array of length n which contains the
    // diagonal elements of the upper triangular matrix s.
  
    // wa is a work array of length n.
  
    //     copy r and (q transpose)*b to preserve input and initialize s.
    // in particular, save the diagonal elements of r in x.
    typedef typename Vector::value_type real;
    real cos_teta, cotan_teta, p5, p25, qtbpj,
      sin_teta, sum, tan_teta, temp, zero(0);
    p5 = real(5)/real(10); p25 = p5/2;
    
    int n = x.GetM();
    for (int j = 0; j < n; j++)
      {
	for (int i = j; i < n; i++)
	  r(i,j) = r(j,i);

	x(j) = r(j,j);
	wa(j) = qtb(j);
      }
  
    //     eliminate the diagonal matrix d using a givens rotation.
  
    for (int j = 0; j < n; j++)
      {
      
	//        prepare the row of d to be eliminated, locating the
	//        diagonal element using p from the qr factorization.
      
	int l = ipvt(j);
	if (diag(l) != zero)
	  {
	    for (int k = j; k < n; k++)
	      sdiag(k) = zero;

	    sdiag(j) = diag(l);
	  
	    // the transformations to eliminate the row of d
	    // modify only a single element of (q transpose)*b
	    // beyond the first n, which is initially zero.
	  
	    qtbpj = zero;
	    for (int k = j; k < n; k++)
	      {
	      
		// determine a givens rotation which eliminates the
		// appropriate element in the current row of d.
	      
		if (sdiag(k) != zero)
		  {
		    if (abs(r(k,k)) < abs(sdiag(k)))
		      {
			cotan_teta = r(k,k)/sdiag(k);
			sin_teta = p5/sqrt(p25+p25*cotan_teta*cotan_teta);
			cos_teta = sin_teta*cotan_teta;
		      }
		    else
		      {		     
			tan_teta = sdiag(k)/r(k,k);
			cos_teta = p5/sqrt(p25+p25*tan_teta*tan_teta);
			sin_teta = cos_teta*tan_teta;
		      }
		  
		  
		    // compute the modified diagonal element of r and
		    // the modified element of ((q transpose)*b,0).
		  
		    r(k,k) = cos_teta*r(k,k) + sin_teta*sdiag(k);
		    temp = cos_teta*wa(k) + sin_teta*qtbpj;
		    qtbpj = -sin_teta*wa(k) + cos_teta*qtbpj;
		    wa(k) = temp;
		  
		    // accumulate the tranformation in the row of s.
		  
		    int kp1 = k + 1;
		    for (int i = kp1; i < n; i++)
		      {
			temp = cos_teta*r(i,k) + sin_teta*sdiag(i);
			sdiag(i) = -sin_teta*r(i,k) + cos_teta*sdiag(i);
			r(i,k) = temp;
		      }
		  }
	      }

	    // store the diagonal element of s and restore
	    // the corresponding diagonal element of r.
	  }
	sdiag(j) = r(j,j);
	r(j,j) = x(j);
      }
  
    // solve the triangular system for z. if the system is
    // singular, then obtain a least squares solution.
    
    int nsing = n-1;
    for (int j = 0; j < n; j++)
      {
	if ((sdiag(j) == zero)&& (nsing == (n-1)))
	  nsing = j - 1;
	if (nsing < (n-1)) 
	  wa(j) = zero;
      }
    if (nsing>=0)
      {
	for (int k = 0; k <= nsing; k++)
	  {
	    int j = nsing - k ;
	    sum = zero;
	    int jp1 = j + 1;
	    for (int i = jp1; i <= nsing; i++)
	      sum = sum + r(i,j)*wa(i);
	  
	    wa(j) = (wa(j) - sum)/sdiag(j);
	  }
      }

    // permute the components of z back to components of x.
  
    for (int j = 0; j < n; j++)
      {
	int l = ipvt(j);
	x(l) = wa(j);
      }
  
    // last card of subroutine qrsolv.
  }
  
  //! intermediary function used to solve non-linear system
  template<class Vector, class Matrix>
  void qform(Matrix& q, Vector& wa)
  {
    // subroutine qform
    //
    //     this subroutine proceeds from the computed qr factorization of
    //     an m by n matrix a to accumulate the m by m orthogonal matrix
    //     q from its factored form.
    //
    //     the subroutine statement is
    //
    //       void qform(q,wa)
    //
    //     where
    //
    //
    //       q is an m by m array. on input the full lower trapezoid in
    //         the first min(m,n) columns of q contains the factored form.
    //         on output q has been accumulated into a square matrix.
  
    //       wa is a work array of length m.
  

    // zero out upper triangle of q in the first min(m,n) columns.
    //
    typedef typename Vector::value_type real;
    real one(1), sum, temp, zero(0);
    int m = q.GetM(); int n = q.GetN();
  
    int minmn = min(m,n);
    for (int j = 1; j < minmn; j++)
      {
	int jm1 = j - 1;
	for (int  i = 0; i <= jm1; i++)
	  q(i,j) = zero;
      }
  
    // initialize remaining columns to those of the identity matrix.
    //
  
    for (int j = n; j < m; j++)
      {
	for (int i = 0; i < m; i++)
	  q(i,j) = zero;

	q(j,j) = one;
      }

  
    //     accumulate q from its factored form.
    //
    for (int l = 0; l < minmn; l++)
      {
	int k = minmn - l-1;
	for (int i = k; i < m; i++)
	  {
	    wa(i) = q(i,k);
	    q(i,k) = zero;
	  }
	q(k,k) = one;
	if (wa(k) != zero)
	  {
	    for (int j = k; j < m; j++)
	      {
		sum = zero;
		for (int i = k; i < m; i++)
		  sum = sum + q(i,j)*wa(i);

		temp = sum/wa(k);
		for (int i = k; i < m; i++)
		  q(i,j) = q(i,j) - temp*wa(i);
	      }
	  }
      }
  }
  
  //! intermediary function used to solve non-linear system
  template<class real, class Vector, class Vector2>
  void dogleg(Vector2& r, Vector& diag, Vector& qtb, real& delta,
	      Vector& x, Vector& wa1, Vector& wa2)
  {

    //     **********
    //
    //     subroutine dogleg
    //
    //    given an m by n matrix a, an n by n nonsingular diagonal
    //     matrix d, an m-vector b, and a positive number delta, the
    //     problem is to determine the convex combination x of the
    //     gauss-newton and scaled gradient directions that minimizes
    //     (a*x - b) in the least squares sense, subject to the
    //     restriction that the euclidean norm of d*x be at most delta.
    //
    //     this subroutine completes the solution of the problem
    //     if it is provided with the necessary information from the
    //     qr factorization of a. that is, if a = q*r, where q has
    //     orthogonal columns and r is an upper triangular matrix,
    //     then dogleg expects the full upper triangle of r and
    //     the first n components of (q transpose)*b.
    //
    //     the subroutine statement is
    //
    //       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
    //
    //     where
  
    //       r is an input array of length lr which must contain the upper
    //         triangular matrix r stored by rows.
    //
    //
    //       diag is an input array of length n which must contain the
    //         diagonal elements of the matrix d.
    //
    //       qtb is an input array of length n which must contain the first
    //         n elements of the vector (q transpose)*b.
    //
    //       delta is a positive input variable which specifies an upper
    //         bound on the euclidean norm of d*x.
    //
    //       x is an output array of length n which contains the desired
    //         convex combination of the gauss-newton direction and the
    //         scaled gradient direction.
    //
    //       wa1 and wa2 are work arrays of length n.
    //
    //
    //     first, calculate the gauss-newton direction.
    //
    int n=qtb.GetM();
    real alpha, bnorm, epsilon_machine, gnorm, one(1),
      qnorm, sgnorm, sum, temp, zero(0);
    
    epsilon_machine = GetPrecisionMachine(one);
  
    int  jj = (n*(n + 1))/2;
    for (int k = 0; k < n; k++)
      {
	int j = n - k - 1;
	int jp1 = j + 1;
	jj = jj - k-1;
	int l = jj + 1;
	sum = zero;
	for (int i = jp1; i < n; i++)
	  {
	    sum = sum + r(l)*x(i);
	    l = l + 1;
	  }
	temp = r(jj);
	if (temp == zero)
	  {
	    l = j;
	    for (int i = 0; i <= j; i++)
	      {
		temp = max(temp,real(abs(r(l))));
		l = l + n-1 - i;
	      }
	    temp = epsilon_machine*temp;
	    if (temp == zero) 
	      temp = epsilon_machine;
	  }
	x(j) = (qtb(j) - sum)/temp;
      }
    // DISP(qtb); DISP(x);
    //
    //     test whether the gauss-newton direction is acceptable.
    //
    for (int j = 0; j < n; j++)
      {
	wa1(j) = zero;
	wa2(j) = diag(j)*x(j);
      }
    // DISP(wa1); DISP(wa2); DISP(delta);
    qnorm = Norm2(wa2);// DISP(qnorm);
    if (qnorm > delta)
      {
	//
	//     the gauss-newton direction is not acceptable.
	//     next, calculate the scaled gradient direction.
	//
	int l = 0;
	for (int j = 0; j < n; j++)
	  {
	    temp = qtb(j);
	    for (int i = j; i < n; i++)
	      {
		wa1(i) = wa1(i) + r(l)*temp;
		l = l + 1;
	      }
	    wa1(j) = wa1(j)/diag(j);
	  }
	// DISP(wa1);
	//
	//     calculate the norm of the scaled gradient and test for
	//     the special case in which the scaled gradient is zero.
	//
      
	gnorm = Norm2(wa1);
	sgnorm = zero;
	alpha = delta/qnorm; // DISP(gnorm); DISP(sgnorm); DISP(alpha);
	if (gnorm != zero)
	  {
	  
	    //     calculate the point along the scaled gradient
	    //     at which the quadratic is minimized.
	  
	    for (int j = 0; j < n; j++)
	      wa1(j) = (wa1(j)/gnorm)/diag(j);

	    l = 0;
	    for (int j = 0; j < n; j++)
	      {
		sum = zero;
		for (int i = j; i < n; i++)
		  {
		    sum = sum + r(l)*wa1(i);
		    l = l + 1;
		  }
		wa2(j) = sum;
	      }
	  
	    temp = Norm2(wa2);
	    sgnorm = (gnorm/temp)/temp;
	    // DISP(wa1); DISP(wa2); DISP(sgnorm);
	    //
	    //     test whether the scaled gradient direction is acceptable.
	    //
	    alpha = zero;
	    if (sgnorm < delta)
	      {
		//
		//     the scaled gradient direction is not acceptable.
		//     finally, calculate the point along the dogleg
		//     at which the quadratic is minimized.
	      
		bnorm = Norm2(qtb);
		temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta);
		temp = temp - (delta/qnorm)*square(real(sgnorm/delta))
		  + sqrt(square(real(temp-(delta/qnorm))) + 
			 (one-square(real(delta/qnorm)))*(one-square(real(sgnorm/delta))));
		alpha = ((delta/qnorm)*(one - square(real(sgnorm/delta))))/temp;
	      }
	    // DISP(alpha); DISP(bnorm); 
	  }
      
	//     form appropriate convex combination of the gauss-newton
	//     direction and the scaled gradient direction.
      
	temp = (one - alpha)*min(sgnorm, delta);
	for (int j = 0; j < n; j++)
	  x(j) = temp*wa1(j) + alpha*x(j);
      
	// DISP(x);
      }
  }
  
  //! intermediary function used to solve non-linear system
  template<class Vector, class Matrix>
  void r1mpyq(Matrix& a, Vector& v, Vector& w)
  {
    // **********
    //
    //     subroutine r1mpyq
    //
    //     given an m by n matrix a, this subroutine computes a*q where
    //     q is the product of 2*(n - 1) transformations
    //
    //           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    //
    //     and gv(i), gw(i) are givens rotations in the (i,n) plane which
    //     eliminate elements in the i-th and n-th planes, respectively.
    //     q itself is not given, rather the information to recover the
    //     gv, gw rotations is supplied.
    //
    //     the subroutine statement is
    //
    //       subroutine r1mpyq(m,n,a,lda,v,w)
    //
    //     where
    //
    //       a is an m by n array. on input a must contain the matrix
    //         to be postmultiplied by the orthogonal matrix q
    //         described above. on output a*q has replaced a.

    //       v is an input array of length n. v(i) must contain the
    //         information necessary to recover the givens rotation gv(i)
    //         described above.
    //
    //       w is an input array of length n. w(i) must contain the
    //         information necessary to recover the givens rotation gw(i)
    //         described above.

    //
    //     apply the first set of givens rotations to a.
    //
    typedef typename Vector::value_type real;
    real cos_teta(0), one(1), sin_teta(0),temp;
    int m = a.GetM(); int n = a.GetN();
  
    int nm1 = n - 1;
    if (nm1 < 1) 
      return;
  
    for (int nmj =0; nmj < nm1; nmj++)
      {
	int j = n - nmj-2;
	if (abs(v(j)) > one) 
	  cos_teta = one/v(j);
	if (abs(v(j)) > one) 
	  sin_teta = sqrt(one-cos_teta*cos_teta);
	if (abs(v(j)) <= one) 
	  sin_teta = v(j);
	if (abs(v(j)) <= one) 
	  cos_teta = sqrt(one-sin_teta*sin_teta);
      
	for (int i = 0; i < m; i++)
	  {
	    temp = cos_teta*a(i,j) - sin_teta*a(i,nm1);
	    a(i,nm1) = sin_teta*a(i,j) + cos_teta*a(i,nm1);
	    a(i,j) = temp;
	  }
      }
    //
    //     apply the second set of givens rotations to a.
    //
    for (int j = 0; j < nm1; j++)
      {
	if (abs(w(j)) > one) 
	  {
	    cos_teta = one/w(j);
	    sin_teta = sqrt(one-cos_teta*cos_teta);
	  }
	else
	  {
	    sin_teta = w(j);
	    cos_teta = sqrt(one-sin_teta*sin_teta);
	  }
	for (int i = 0; i < m; i++)
	  {
	    temp = cos_teta*a(i,j) + sin_teta*a(i,nm1);
	    a(i,nm1) = -sin_teta*a(i,j) + cos_teta*a(i,nm1);
	    a(i,j) = temp;
	  }
      }
  }
  
  //! intermediary function used to solve non-linear system
  template<class Vector>
  void r1mpyq(Vector& a, Vector& v, Vector& w)
  {
    // **********
    //
    //     subroutine r1mpyq
    //
    //     given an m by n matrix a, this subroutine computes a*q where
    //     q is the product of 2*(n - 1) transformations
    //
    //           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    //
    //     and gv(i), gw(i) are givens rotations in the (i,n) plane which
    //     eliminate elements in the i-th and n-th planes, respectively.
    //     q itself is not given, rather the information to recover the
    //     gv, gw rotations is supplied.
    //
    //     the subroutine statement is
    //
    //       subroutine r1mpyq(m,n,a,lda,v,w)
    //
    //     where
    //
    //       a is an m by n array. on input a must contain the matrix
    //         to be postmultiplied by the orthogonal matrix q
    //         described above. on output a*q has replaced a.

    //       v is an input array of length n. v(i) must contain the
    //         information necessary to recover the givens rotation gv(i)
    //         described above.
    //
    //       w is an input array of length n. w(i) must contain the
    //         information necessary to recover the givens rotation gw(i)
    //         described above.

    //
    //     apply the first set of givens rotations to a.
    //
    typedef typename Vector::value_type real;
    real cos_teta(0),one(1),sin_teta(0),temp;
    int n = a.GetM();
  
    int nm1 = n - 1;
    if (nm1 < 1) 
      return;
  
    for (int nmj =0; nmj < nm1; nmj++)
      {
	int j = n - nmj-2;
	if (abs(v(j)) > one) 
	  cos_teta = one/v(j);
	if (abs(v(j)) > one) 
	  sin_teta = sqrt(one-cos_teta*cos_teta);
	if (abs(v(j)) <= one) 
	  sin_teta = v(j);
	if (abs(v(j)) <= one) 
	  cos_teta = sqrt(one-sin_teta*sin_teta);
      
	temp = cos_teta*a(j) - sin_teta*a(nm1);
	a(nm1) = sin_teta*a(j) + cos_teta*a(nm1);
	a(j) = temp;
      }
    //
    //     apply the second set of givens rotations to a.
    //
    for (int j = 0; j <nm1; j++)
      {
	if (abs(w(j)) > one) 
	  {
	    cos_teta = one/w(j);
	    sin_teta = sqrt(one-cos_teta*cos_teta);
	  }
	else
	  {
	    sin_teta = w(j);
	    cos_teta = sqrt(one-sin_teta*sin_teta);
	  }
	temp = cos_teta*a(j) + sin_teta*a(nm1);
	a(nm1) = -sin_teta*a(j) + cos_teta*a(nm1);
	a(j) = temp;
      }
  }
  
  //! intermediary function used to solve non-linear system
  template<class Vector2, class Vector>
  void r1updt(Vector2& s, Vector& u, Vector& v, Vector& w, bool& sing)
  {
    //     **********
    //
    //     subroutine r1updt
    //
    //     given an m by n lower trapezoidal matrix s, an m-vector u,
    //     and an n-vector v, the problem is to determine an
    //     orthogonal matrix q such that
    //
    //                   t
    //           (s + u*v )*q
    //
    //     is again lower trapezoidal.
    //
    //     this subroutine determines q as the product of 2*(n - 1)
    //     transformations
    //
    //           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    //
    //     where gv(i), gw(i) are givens rotations in the (i,n) plane
    //     which eliminate elements in the i-th and n-th planes,
    //     respectively. q itself is not accumulated, rather the
    //     information to recover the gv, gw rotations is returned.
    //
    //     the subroutine statement is
    //
    //       subroutine r1updt(m,n,s,ls,u,v,w,sing)
    //
    //     where
    //
    //
    //       s is an array of length ls. on input s must contain the lower
    //         trapezoidal matrix s stored by columns. on output s contains
    //         the lower trapezoidal matrix produced as described above.

    //       u is an input array of length m which must contain the
    //         vector u.
    //
    //       v is an array of length n. on input v must contain the vector
    //         v. on output v(i) contains the information necessary to
    //         recover the givens rotation gv(i) described above.
    //
    //       w is an output array of length m. w(i) contains information
    //         necessary to recover the givens rotation gw(i) described
    //         above.
    //
    //       sing is a logical output variable. sing is set true if any
    //         of the diagonal elements of the output s are zero. otherwise
    //         sing is set false.
    //

    //
    //     giant is the largest magnitude.
    //
    typedef typename Vector::value_type real;
    real cos_teta, cotan_teta, giant, one(1), p5, p25, sin_teta, tan_teta,
      tau, temp, zero(0);
    
    p5 = real(5)/10; p25 = p5/2;
    giant = 1e300;
    int n = v.GetM(); int m = u.GetM();
    //
    //     initialize the diagonal element pointer.
    //
    int jj = (n*(2*m - n + 1))/2 - (m - n)-1;
    //
    //     move the nontrivial part of the last column of s into w.
    //
    int l = jj;
    for (int i = (n-1); i < m; i++)
      {
	w(i) = s(l);
	l = l + 1;
      }

    //
    //     rotate the vector v into a multiple of the n-th unit vector
    //     in such a way that a spike is introduced into w.
    //
    int nm1 = n - 1;
    for (int nmj = 0; nmj < nm1; nmj++)
      {
	int j = nm1 - nmj-1;
	jj = jj - (m - j );
	w(j) = zero;
	if (v(j) != zero)
	  {
	    //	 
	    //        determine a givens rotation which eliminates the
	    //        j-th element of v.
	    //
	    if (abs(v(n-1)) < abs(v(j)))
	      {
		cotan_teta = v(nm1)/v(j);
		sin_teta = p5/sqrt(p25+p25*cotan_teta*cotan_teta);
		cos_teta = sin_teta*cotan_teta;
		tau = one;
		if (abs(cos_teta)*giant > one) 
		  tau = one/cos_teta;
	      }
	    else
	      {
		tan_teta = v(j)/v(nm1);
		cos_teta = p5/sqrt(p25+p25*tan_teta*tan_teta);
		sin_teta = cos_teta*tan_teta;
		tau = sin_teta;
	      }
	  
	    //        apply the transformation to v and store the information
	    //        necessary to recover the givens rotation.
	  
	    v(nm1) = sin_teta*v(j) + cos_teta*v(nm1);
	    v(j) = tau;
	  
	    //        apply the transformation to s and extend the spike in w.
	  
	    int l = jj;
	    for (int i = j; i < m; i++)
	      {
		temp = cos_teta*s(l) - sin_teta*w(i);
		w(i) = sin_teta*s(l) + cos_teta*w(i);
		s(l) = temp;
		l = l + 1;
	      }
	  }
      }
  
    // add the spike from the rank 1 update to w.
  
    for (int i = 0; i < m; i++)
      w(i) = w(i) + v(nm1)*u(i);

    //     eliminate the spike.
  
    sing = false;
    for (int j = 0; j < nm1; j++)
      {
	if (w(j) != zero)
	  {
	    //
	    //        determine a givens rotation which eliminates the
	    //        j-th element of the spike.
	    //
	    if (abs(s(jj)) < abs(w(j))) 
	      {
		cotan_teta = s(jj)/w(j);
		sin_teta = p5/sqrt(p25+p25*cotan_teta*cotan_teta);
		cos_teta = sin_teta*cotan_teta;
		tau = one;
		if (abs(cos_teta)*giant > one) 
		  tau = one/cos_teta;
	      }
	    else
	      {
		tan_teta = w(j)/s(jj);
		cos_teta = p5/sqrt(p25+p25*tan_teta*tan_teta);
		sin_teta = cos_teta*tan_teta;
		tau = sin_teta;
	      }
	    //
	    //        apply the transformation to s and reduce the spike in w.
	    //
	  
	    int l = jj;
	    for (int i = j; i < m; i++)
	      {
		temp = cos_teta*s(l) + sin_teta*w(i);
		w(i) = -sin_teta*s(l) + cos_teta*w(i);
		s(l) = temp;
		l = l + 1;
	      }
	    //
	    //        store the information necessary to recover the
	    //        givens rotation.
	    //
	    w(j) = tau;
	  }
	//
	//        test for zero diagonal elements in the output s.
      
	if (s(jj) == zero) 
	  sing = true;
      
	jj = jj + (m - j );
      }
  
    //     move w back into the last column of the output s.
      
    l = jj;
    for (int i = (n-1); i < m; i++)
      {
	s(l) = w(i);
	l = l + 1;
      }
  
    if (s(jj) == zero) 
      sing = true;
  }


#ifdef MONTJOIE_WITH_GSL
  //! evaluation of the function to minimize
  template<class T>
  double VirtualMinimizedFunction<T>::my_f(const gsl_vector *v, void *params)
  {
    VirtualMinimizedFunction<T>* fct = static_cast<VirtualMinimizedFunction<T>* >(params);
    T feval;
    Vector<T> x(fct->GetM());
    for (int i = 0; i < fct->GetM(); i++)
      x(i) = gsl_vector_get(v, i);
    
    fct->EvaluateFunction(x, feval);

    return toDouble(feval);
  }
  
  //! evaluation of the derivative of the function to minimize
  template<class T>
  void VirtualMinimizedFunction<T>::my_df(const gsl_vector *v, void *params, gsl_vector *df)
  {
    VirtualMinimizedFunction<T>* fct = static_cast<VirtualMinimizedFunction<T>* >(params);
    T feval;
    Vector<T> x(fct->GetM()), param, df_vec(fct->GetM());
    for (int i = 0; i < fct->GetM(); i++)
      x(i) = gsl_vector_get(v, i);
    
    fct->EvaluateFunctionGradient(x, feval, df_vec);
    
    for (int i = 0; i < fct->GetM(); i++)
      gsl_vector_set(df, i, toDouble(df_vec(i)));
    
  }
  
  //! evaluation of the function to minimize, and its derivatives
  template<class T>
  void VirtualMinimizedFunction<T>::my_fdf(const gsl_vector *v, void *params, double* f, gsl_vector *df)
  {
    VirtualMinimizedFunction<T>* fct = static_cast<VirtualMinimizedFunction<T>* >(params);
    T feval;
    Vector<T> x(fct->GetM()), param, df_vec(fct->GetM());
    for (int i = 0; i < fct->GetM(); i++)
      x(i) = gsl_vector_get(v, i);
    
    fct->EvaluateFunctionGradient(x, feval, df_vec);
    
    for (int i = 0; i < fct->GetM(); i++)
      gsl_vector_set(df, i, toDouble(df_vec(i)));
    
    *f = toDouble(feval);
  }
#endif

#ifdef MONTJOIE_WITH_ALGLIB
  template<class T>
  void VirtualMinimizedFunction<T>::alglib_fdf(const alglib::real_1d_array &x, double &func,
					       alglib::real_1d_array &grad, void *ptr)
  {
    VirtualMinimizedFunction<T>* fct = static_cast<VirtualMinimizedFunction<T>* >(ptr);

    T feval;
    Vector<T> df_vec, xsol;

    xsol.SetData(x.length(), const_cast<T*>(x.getcontent()));
    df_vec.SetData(grad.length(), grad.getcontent());
    
    fct->EvaluateFunctionGradient(xsol, feval, df_vec);

    xsol.Nullify();
    df_vec.Nullify();
    
    func = toDouble(feval);
  }

  template<class T>
  void VirtualMinimizedFunction<T>::alglib_rep(const alglib::real_1d_array &x, double func, void* ptr)
  {
    VirtualMinimizedFunction<T>* fct = static_cast<VirtualMinimizedFunction<T>* >(ptr);
    Vector<T> xsol;
    xsol.SetData(x.length(), const_cast<T*>(x.getcontent()));

    fct->step_number++;
    fct->GiveIntermediateSolution(fct->step_number, xsol);

    xsol.Nullify();
  }
#endif
  

#ifdef MONTJOIE_WITH_GSL
  //! Minimizing a general function with Gsl
  /*!
    \param[in] fct given problem
    \param[out] xsol solution parameters
   */
  template<class T>
  T MinimizeParametersGsl(VirtualMinimizedFunction<T>& fct, Vector<T>& xsol,
			  double epsilon, unsigned nb_max_iterations)
  {
    // GSL class for definition of function f and derivative df
    gsl_multimin_function_fdf my_func;
    
    my_func.f = &fct.my_f;
    my_func.df = &fct.my_df;
    my_func.fdf = &fct.my_fdf;
    my_func.n = fct.GetM();
    my_func.params = &fct;

    // case where only f is needed
    gsl_multimin_function my_func_f;
    my_func_f.f = &fct.my_f;
    my_func_f.n = fct.GetM();
    my_func_f.params = &fct;
    
    // minimized function
    const gsl_multimin_fdfminimizer_type *Tmin = NULL;
    gsl_multimin_fdfminimizer *s = NULL;

    const gsl_multimin_fminimizer_type *Tmin_f = NULL;
    gsl_multimin_fminimizer *s_f = NULL;
    
    // solution stored in a gsl vector
    gsl_vector *x, *step;
    //DISP(fct.GetM()); DISP(xsol.GetM());
    x = gsl_vector_alloc(fct.GetM());
    
    // initial guess
    fct.FindInitGuess(xsol);
    for (int i = 0; i < xsol.GetM(); i++)
      gsl_vector_set(x, i, toDouble(xsol(i)));
    
    // selecting the required algorithm
    double step_size_simplex = fct.GetStepSize();
    double step_size = fct.GetStepSize();
    double line_minim_param = fct.GetLineMinimizationParameter();
    bool algo_with_only_f = false;
    switch (fct.GetGslAlgorithm())
      {
      case VirtualMinimizedFunction<T>::BFGS :
	Tmin = gsl_multimin_fdfminimizer_vector_bfgs;
	break;
      case VirtualMinimizedFunction<T>::BFGS2 :
	Tmin = gsl_multimin_fdfminimizer_vector_bfgs2;
	break;
      case VirtualMinimizedFunction<T>::CG :
	Tmin = gsl_multimin_fdfminimizer_conjugate_pr;
	break;
      case VirtualMinimizedFunction<T>::CG_FR :
	Tmin = gsl_multimin_fdfminimizer_conjugate_fr;
	break;
      case VirtualMinimizedFunction<T>::STEEPEST_DESCENT :
	Tmin = gsl_multimin_fdfminimizer_steepest_descent;
	break;
      case VirtualMinimizedFunction<T>::SIMPLEX :
	algo_with_only_f = true;
	Tmin_f = gsl_multimin_fminimizer_nmsimplex;
	break;
      case VirtualMinimizedFunction<T>::SIMPLEX2 :
	algo_with_only_f = true;
	Tmin_f = gsl_multimin_fminimizer_nmsimplex2;
	break;
      case VirtualMinimizedFunction<T>::SIMPLEX2_RAND :
	algo_with_only_f = true;
	Tmin_f = gsl_multimin_fminimizer_nmsimplex2rand;
	break;
      default:
        {
          cout << "Unknown algorithm" << endl;
          abort();
        }
      }
    
    // allocating the minimizer
    if (algo_with_only_f)
      {
	step = gsl_vector_alloc(my_func.n);
	for (int i = 0; i < xsol.GetM(); i++)
	  gsl_vector_set(step, i, step_size_simplex);
	
	s_f = gsl_multimin_fminimizer_alloc (Tmin_f, my_func_f.n);      
	gsl_multimin_fminimizer_set(s_f, &my_func_f, x, step);
      }
    else
      {
	s = gsl_multimin_fdfminimizer_alloc (Tmin, my_func.n);	
	gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, line_minim_param);
      }
    
    // main loop of minimization
    size_t iter = 0;
    int status;

    do
      {
	iter++;
	if (algo_with_only_f)
	  {
	    status = gsl_multimin_fminimizer_iterate(s_f);
	    for (int i = 0; i < xsol.GetM(); i++)
	      xsol(i) = gsl_vector_get(s_f->x, i);	    
	  }
	else
	  {
	    status = gsl_multimin_fdfminimizer_iterate(s);
	    for (int i = 0; i < xsol.GetM(); i++)
	      xsol(i) = gsl_vector_get(s->x, i);
	  }
	
	fct.GiveIntermediateSolution(iter, xsol);
	if (status)
	  break;
	
	if (algo_with_only_f)
	  {
	    double size = gsl_multimin_fminimizer_size(s_f);
	    status = gsl_multimin_test_size(size, epsilon);
	  }
	else
	  status = gsl_multimin_test_gradient (s->gradient, epsilon);
      
      }
    while (status == GSL_CONTINUE && iter < nb_max_iterations);
    
    // the solution is stored in xsol
    if (algo_with_only_f)
      for (int i = 0; i < xsol.GetM(); i++)
	xsol(i) = gsl_vector_get(s_f->x, i);
    else
      for (int i = 0; i < xsol.GetM(); i++)
	xsol(i) = gsl_vector_get(s->x, i);

    // releasing memory needed by Gsl
    if (algo_with_only_f)
      {
	gsl_vector_free(step);
	gsl_multimin_fminimizer_free(s_f);
      }
    else
      gsl_multimin_fdfminimizer_free(s);
    
    gsl_vector_free(x);
    
    // we return the minimum reached
    T feval; fct.EvaluateFunction(xsol, feval);
    return feval;
  }
#endif


#ifdef MONTJOIE_WITH_ALGLIB
  //! Minimizing a general function with Alglib
  /*!
    \param[in] fct given problem
    \param[out] xsol solution parameters
   */
  template<class T>
  T MinimizeParametersAlglib(VirtualMinimizedFunction<T>& fct, Vector<T>& xsol,
			     double epsilon, unsigned nb_max_iterations)
  {
    alglib::real_1d_array x;

    // initial guess
    fct.FindInitGuess(xsol);
    
    x.setcontent(xsol.GetM(), xsol.GetData());

    double epsg = epsilon;
    double epsf = 0;
    double epsx = 0;
    alglib::ae_int_t maxits = nb_max_iterations;

    if ((fct.GetGslAlgorithm() == fct.BFGS) || (fct.GetGslAlgorithm() == fct.BFGS2))
      {
	alglib::minlbfgsstate state;
	alglib::minlbfgsreport rep;
	int m = min(int(xsol.GetM()), 3);
	
	alglib::minlbfgscreate(m, x, state);
	alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
	alglib::minlbfgssetxrep(state, true);
	
	alglib::minlbfgsoptimize(state, &fct.alglib_fdf, &fct.alglib_rep,
				 reinterpret_cast<void*>(&fct));
	
	alglib::minlbfgsresults(state, x, rep);
	for (int i = 0; i < xsol.GetM(); i++)
	  xsol(i) = x[i];
      }
    else if ((fct.GetGslAlgorithm() == fct.CG) || (fct.GetGslAlgorithm() == fct.CG_FR))
      {
	alglib::mincgstate state;
	alglib::mincgreport rep;

	alglib::mincgcreate(x, state);
	alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
	alglib::mincgsetxrep(state, true);
	
	alglib::mincgoptimize(state, &fct.alglib_fdf, &fct.alglib_rep,
			      reinterpret_cast<void*>(&fct));
	
	alglib::mincgresults(state, x, rep);
	for (int i = 0; i < xsol.GetM(); i++)
	  xsol(i) = x[i];
      }
    else
      {
	cout << "Unknown algorithm" << endl;
	abort();
      }
    
    // we return the minimum reached
    T feval; fct.EvaluateFunction(xsol, feval);
    return feval;
  }
#endif

} // end namespace

#define MONTJOIE_FILE_NONLINEAR_EQUATIONS_CXX
#endif
