#ifndef MONTJOIE_FILE_GAUSS_LOBATTO_POINTS_CXX

namespace Montjoie
{
  //! numerical computation of Gauss or Gauss-Lobatto points and weights
  /*!
    \param[in] nb_points_gauss number of Gauss or GL points to compute
    \param[out] points obtained Gauss points
    \param[out] weights obtained Gauss weights
    \param[in] gauss_lobatto if true, Gauss-Lobatto points are computed
   */
  template<class Vector>
  void ComputeGaussFormulas(int nb_points_gauss, Vector& points,
			    Vector& weights, bool gauss_lobatto)
  {
    // For details about that function, take a look at
    //
    //     P. N. swarztrauber, Computing the points and weights for 
    //     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
    //     24(2002) pp. 945-954.
    //     

    typedef typename Vector::value_type real;
    real unity(1), zero(0); 
    int nb_points_quadrature = nb_points_gauss;

    // Points(i) will be equal to cos(theta(i))
    Vector theta;
    
    // for Gauss-Lobatto, we compute internal nodes
    if (gauss_lobatto)
      nb_points_quadrature--;
    
    theta.Reallocate(nb_points_gauss);
    points.Reallocate(nb_points_gauss);
    weights.Reallocate(nb_points_gauss);
    
    // coefficient of Fourier transform of Legendre polynom and derivative
    // we need the second derivative for Gauss-Lobatto points   
    Vector coef(nb_points_quadrature), deriv_coef(nb_points_quadrature), 
      second_deriv_coef(nb_points_quadrature);
    
    coef.Fill(zero); deriv_coef.Fill(zero); second_deriv_coef.Fill(zero);
    
    real pi_;
    pi_ = pi_wp;
    if (pi_ < real(3))
      {
	cout<<"Is Montjoie initialized ? "<<endl;
	abort();
      }
    
    real pis2 = 0.5*pi_;
    if (nb_points_quadrature <= 0) 
      {
	cout << "Specify at least one point to compute " << endl;
	abort();
      }
    
    // exact expressions are used for order 1, 2, 3
    if (gauss_lobatto)
      {
	if (nb_points_quadrature == 1)
	  {
	    // [0,1] with weights 1/2, 1/2
	    points(0) = zero;
	    points(1) = unity;
	    
	    weights(0) = 0.5;
	    weights(1) = 0.5;
	    return;
	  }
	else if (nb_points_quadrature == 2)
	  {
	    // [0, 1/2, 1] with weights 1/6, 2/3, 1/6
	    points(0) = zero;
	    points(1) = unity/real(2);
	    points(2) = unity;
	    
	    weights(0) = unity/real(6);
	    weights(1) = unity*real(2)/real(3);
	    weights(2) = weights(0);
	    	    
	    return;
	  }
	else if (nb_points_quadrature == 3)
	  {
	    points(0) = zero;
	    points(1) = (real(5)-Sqrt(real(5)))/real(10);
	    points(2) = unity - points(1);
	    points(3) = unity;
	    
	    weights(0) = real(1)/real(12);
	    weights(1) = real(5)/real(12);
	    weights(2) = weights(1);
	    weights(3) = weights(0);
	  }
      }
    else
      {
	if (nb_points_quadrature == 1)
	  {
	    // 1/2 with weight 1
	    points(0) = real(0.5);
	    weights(0) = unity;
	    return;
	  }
	else if (nb_points_quadrature == 2)
	  {
	    // [(3 - sqrt(3))/6, (3 + sqrt(3))/6] with weights 1/2, 1/2
	    points(0) = (real(3) - Sqrt(real(3)))/real(6);
	    points(1) = unity - points(0);
	    
	    weights(0) = 0.5;
	    weights(1) = 0.5;
	    
	    return;
	  }
	else if (nb_points_quadrature == 3)
	  {
	    points(0) = real(0.5) - Sqrt(real(15))/real(10);
	    points(1) = real(0.5);
	    points(2) = unity - points(0);
	    
	    weights(0) = real(5)/real(18);
	    weights(1) = real(4)/real(9);
	    weights(2) = weights(0);
	    
	    return;
	  }
      }
    
    // for order greater than 3, numerical computation
    real eps = GetPrecisionMachine(unity);
    int parity_nb_points = nb_points_quadrature%2;
    int half_nb_points = nb_points_gauss/2;
    int nhalf = (nb_points_quadrature+1)/2;
    
    real cz = zero;
    ComputeFourierCoef(nb_points_quadrature, cz, coef,
		       deriv_coef, second_deriv_coef);
    
    real dtheta = zero;
    if (gauss_lobatto&&(parity_nb_points==1))
      dtheta = pis2/(nhalf-1);
    else
      dtheta = pis2/nhalf;
    
    real dthalf = dtheta/real(2);
    real cmax = dtheta/real(5);
    // the zeros of Legendre polynom
    real zprev = zero, zlast = zero, zhold = zero;
    real pb = zero, dpb = zero, ddpb = zero, dcor = zero, sgnd = zero;
    int nix;
    //
    //     estimate first point next to theta = pi/2
    //
    if (parity_nb_points!=0)
      {
	// nb_points_quadrature = 2 nhalf-1
	// if odd the first zero is at the middle pi/2 
	// and the following pi/2-pi/(2*nhalf)
	if (gauss_lobatto)
	  {
	    zero = pis2 - dthalf;
	    zprev = pis2 + dthalf;
	    nix = nhalf - 1;
	  }
	else
	  {
	    zero = pis2 - dtheta;
	    zprev = pis2;
	    nix = nhalf - 1; // index of the point
	  }
      }
    else
      {
	// if even, no zero on the middle, the first zero is on pi/2-pi/(4*nhalf)
	if (gauss_lobatto)
	  {
	    zero = pis2 - dtheta;
	    zprev = pis2;
	    nix = nhalf-1; // index of the point
	  }
	else
	  {
	    zero = pis2-dthalf; 
	    zprev = pis2+dthalf;
	    nix = nhalf;
	  }
      }
  
    bool each_point_not_computed = true;
    while (each_point_not_computed)
      {
	int nb_iter = 0;
	bool test_convergence_newton = true;
	real residu(1), residu_prec(1);
	while ((test_convergence_newton)&&(nb_iter < 100)) 
	  {
	    nb_iter++;
	    zlast = zero;
	    //
	    //     newton iterations
	    //
	    ComputeLegendrePol_and_Derivative(nb_points_quadrature, zero, cz, coef,
					      deriv_coef, second_deriv_coef,
					      pb, dpb, ddpb);
	    if (gauss_lobatto)
	      dcor = dpb/ddpb;
	    else
	      dcor = pb/dpb;
	    
	    sgnd = unity;
	    if (dcor != real(0)) 
	      sgnd = dcor/abs(dcor);
	    
	    // we don't move the point further than 0.2*delta_theta
	    real tmp = abs(dcor); 
	    dcor = sgnd*min(tmp, cmax);
	    zero = zero-dcor;
	    residu_prec = residu; 
	    residu = abs(zero-zlast);
	    // we check if the stopping criteria are reached
	    if ((abs(zero-zlast) <= eps*abs(zero))
		||((nb_iter>5)&&(residu_prec < residu)))
	      test_convergence_newton = false;
	  }
	
	theta(nix-1) = zero;
	zhold = zero;
	//      weights(nix) = (nb_points_quadrature+nb_points_quadrature+1)/(dpb*dpb)
	//    
	//     yakimiw's formula permits using old pb and dpb
	//
	if (gauss_lobatto)
	  {
	    real tmp = pb - dcor*dpb + 0.5*dcor*dcor*ddpb;
	    tmp*=tmp;
	    weights(nix-1) = unity/tmp;
	  }
	else
	  {
	    real tmp = dpb + pb*cos(zlast)/sin(zlast); tmp*=tmp;
	    weights(nix-1) = real(nb_points_quadrature+nb_points_quadrature+1)/tmp;
	  }
	nix--;
      
	if (nix == 0)
	  each_point_not_computed = false;
	else if (nix <= (nhalf-1))
	  zero += zero-zprev;
	
	zprev = zhold;
      }
  
    //
    //     extend points and weights via symmetries
    //
    if ((!gauss_lobatto)&&(parity_nb_points != 0))
      {
	theta(nhalf-1) = pis2;
	ComputeLegendrePol_and_Derivative(nb_points_quadrature, pis2, cz, coef,
					  deriv_coef, second_deriv_coef, 
					  pb, dpb, ddpb);
	
	weights(nhalf-1) = (nb_points_quadrature 
			    + nb_points_quadrature+1)/(dpb*dpb);
      }
    if ((gauss_lobatto)&&(parity_nb_points==0))
      {
	theta(nhalf-1) = pis2;
	ComputeLegendrePol_and_Derivative(nb_points_quadrature, pis2, cz, coef,
					  deriv_coef, second_deriv_coef,
					  pb, dpb, ddpb);
	weights(nhalf-1) = unity/(pb*pb);
      }
    
    if (gauss_lobatto)
      {
	for (int i = (nhalf-1);i >= 0; i--)
	  {
	    theta(i+1) = theta(i);
	    weights(i+1) = weights(i);
	  }
	theta(0) = 0.0;
	ComputeLegendrePol_and_Derivative(nb_points_quadrature, theta(0), cz, coef,
					  deriv_coef, second_deriv_coef, pb, dpb, ddpb);
	
	weights(0) = unity/(pb*pb);
      }

    for (int i = 0; i < half_nb_points; i++)
      {
	weights(nb_points_gauss-i-1) = weights(i);
	theta(nb_points_gauss-i-1) = pi_-theta(i);
      }

    real sum(0);
    for (int i = 0; i < nb_points_gauss; i++)
      {
	sum += weights(i);
      }
    for (int i = 0; i < nb_points_gauss; i++)
      {
	weights(i) = weights(i)/sum;
	points(i) = real(0.5)+0.5*cos(theta(i));
      }
    
    // if points are not sorted
    Sort(points, weights);
  }

  
  //! intermediary function used for the computation of Gauss points
  template<class real,class Vector>
  void ComputeFourierCoef(const int n, real& cz, Vector& cp,
			  Vector& dcp, Vector& ddcp)
  {
    //
    //     computes the fourier coefficients of the legendre
    //     polynomial p_n0 and its derivative. 
    //     n is the degree and n/2 or (n+1)/2
    //     coefficients are returned in cp depending on whether
    //     n is even or odd. The same number of coefficients
    //     are returned in dcp. For n even the constant 
    //     coefficient is returned in cz. 
    //
    real t1(0), t2(0), t3(0), t4(0), unity(1), number_two(2);
    int ncp = (n+1)/2;
    t1 = -unity;
    t2 = n+unity;
    t3 = 0.0;
    t4 = n+n+unity;
    if(n%2==0)
      {
	cp(ncp) = unity;
	for (int j = ncp; j >= 2; j--)
	  {
	    t1 += number_two;
	    t2 += -unity;
	    t3 += unity;
	    t4 += -number_two;
	    cp(j-1) = (t1*t2)/(t3*t4)*cp(j);
	  }
	
	t1 = t1+number_two;
	t2 = t2-unity;
	t3 = t3+unity;
	t4 = t4-number_two;
	cz = (t1*t2)/(t3*t4)*cp(1);
	for (int j = 1; j <= ncp; j++)
	  {
	    dcp(j) = (j+j)*cp(j);
	    ddcp(j) = (j+j)*dcp(j);
	  }
      }
    else
      {
	cp(ncp) = unity;
	for (int j = ncp-1; j >= 1; j--)
	  {
	    t1 = t1+number_two;
	    t2 = t2-unity;
	    t3 = t3+unity;
	    t4 = t4-number_two;
	    cp(j) = (t1*t2)/(t3*t4)*cp(j+1);
	  }
	for (int j = 1; j <= ncp; j++)
	  {
	    dcp(j) = (j+j-1)*cp(j);
	    ddcp(j) = (j+j-1)*dcp(j);
	  }
      }
    return;
  }

  
  //! intermediary function used for the computation of Gauss points
  template<class real, class Vector>
  void ComputeLegendrePol_and_Derivative
  (int n, real theta, const real& cz, const Vector& cp, const Vector& dcp,
   const Vector& ddcp, real& pb, real& dpb, real& ddpb)
  {
    //     computes pn(theta) and its derivative dpb(theta) with 
    //     respect to theta
    //

    real cdt = cos(theta+theta);
    real sdt = sin(theta+theta);
    real cth(0), sth(0), chh(0);
    if (n%2==0)
      {
      
	//     n even
      
	int kdo = n/2;
	pb = cz/2;
	dpb = real(0);
	ddpb = real(0);
	if (n > 0) 
	  {
	    cth = cdt;
	    sth = sdt;
	    for (int k = 1; k <= kdo; k++)
	      {
		//      pb = pb+cp(k)*cos(2*k*theta)
		pb = pb+cp(k)*cth;
		//      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
		dpb = dpb-dcp(k)*sth;
		ddpb = ddpb-ddcp(k)*cth;
		chh = cdt*cth-sdt*sth;
		sth = sdt*cth+cdt*sth;
		cth = chh;
	      }
	  }
      }
    else
      {
      
	//     n odd
      
	int kdo = (n+1)/2;
	pb = real(0);
	dpb = real(0);
	ddpb = real(0);
	cth = cos(theta);
	sth = sin(theta);
	for (int k = 1; k <= kdo; k++)
	  {
	    //      pb = pb+cp(k)*cos((2*k-1)*theta)
	    pb = pb+cp(k)*cth;
	    //      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta)
	    dpb = dpb-dcp(k)*sth;
	    ddpb = ddpb-ddcp(k)*cth;
	    chh = cdt*cth-sdt*sth;
	    sth = sdt*cth+cdt*sth;
	    cth = chh;
	  }
      }

    return;
  }
  
  
  //! computes blending of Gauss and Gauss-Lobatto rules as defined by Ainsworth
  template<class T>
  void ComputeGaussBlendedFormulas(int nb_points, Vector<T>& points,
				   Vector<T>& weights, const T& tau)
  {
    if (nb_points <= 1)
      {
	points.Clear();
	weights.Clear();
	return;
      }

    T x_init, x;
    Vector<T> Pn, dPn;
    LegendrePolynomial<T> Leg(nb_points+1);
    
    // using a basic Newton algorithm
    points.Reallocate(nb_points);
    weights.Reallocate(nb_points);
    // boucle sur les points (seulement la moitie, les autres etant obtenus par symetrie)
    for (int i = 0; i < (nb_points+1)/2; i++)
      {
        // initialisation avec le point de Tchebychev
	x_init = -cos(0.5*pi_wp*(2*i+1)/nb_points);
	// tolerance de l'algo de Newton
	T epsilon = 20*epsilon_machine;
	T test, test_prec = 1e100, feval, der_feval;
	if ((nb_points%2 == 1) && (i == (nb_points+1)/2))
	  {
	    // la solution est ici connue, c'est x = 0
	    x_init = T(0);
	    test_prec = T(0);
	  }
	
	x = x_init;
	Leg.EvaluateDerive(nb_points+1, x, Pn, dPn);
	feval = Pn(nb_points) - tau*Pn(nb_points-2);
	test = abs(feval);
        // boucle de l'algo de Newton
	while ((test < test_prec) && (test > epsilon))
	  {
	    // newton iteration
	    der_feval = dPn(nb_points) - tau*dPn(nb_points-2);
	    x -= feval/der_feval;

	    // on met a jour Pn et dPn
	    Leg.EvaluateDerive(nb_points+1, x, Pn, dPn);
	    feval = Pn(nb_points) - tau*Pn(nb_points-2);
	    test_prec = test; test = abs(feval);
	  }
	
        // point d'integration sur [0, 1]
	points(i) = 0.5 + 0.5*x;
		
        // poids d'integration
	int p = nb_points-1;
	weights(i) = (p*(tau+1.0) + tau) / (p*(p+1)*Pn(p)*(dPn(p+1) - tau*dPn(p-1)));

	// symetrie pour obtenir les autres points
	points(nb_points-1-i) = T(1) - points(i);
	weights(nb_points-1-i) = weights(i);
      }
  }  
  
} // end namespace

#define MONTJOIE_FILE_GAUSS_LOBATTO_POINTS_CXX
#endif
