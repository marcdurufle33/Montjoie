#ifndef SELDON_FILE_EIGENVALUE_CXX

namespace Seldon
{

#ifndef SELDON_WITH_LAPACK
  
  /***************************************
   * C++ translations of Lapack routines *
   ***************************************/
  
  ////////////////////////////
  // General complex matrix //
  ////////////////////////////
  
  
  template<class Matrix1, class Vector1, class reel>
  void TransformMatrix_to_HessenbergForm(int low, int high, Matrix1& A,
                                         Vector1& Q, const reel& ftmp)
  {
    typedef typename Matrix1::value_type complexe;
    complexe fr, zero(reel(0), reel(0)); fr = zero;
    int n = A.GetM();
    reel f(0), g(0), scale(0), h(0), zero_reel(0), one(1);
    
    for (int m = low+1; m <= high-1; m++)
      {
	h = zero_reel;
	Q(m) = zero;
	scale = zero_reel;
	for (int i = m; i <= high; i++)
	  scale += abs(real(A(i,m-1))) + abs(imag(A(i,m-1)));
	
	if (scale != zero_reel)
	  {
	    for (int i = high; i >= m; i--)
	      {
		Q(i) = A(i,m-1)/scale;
		h += real(Q(i)*conj(Q(i)));
	      }
	    	    
	    if (h > zero_reel)
	      g = sqrt(h);
	    else
	      g = zero_reel;
	    
	    if (Q(m) == zero)
	      f = zero_reel;
	    else
	      f = abs(Q(m));
	    
	    if (f==zero_reel)
	      {
		Q(m) = complexe(g, zero_reel);
		A(m,m-1) = complexe(scale, zero_reel);
	      }
	    else
	      {
		h += f*g;
		g /= f;
		Q(m) *= one + g;
	      }
	    
	    for (int j = m; j < n; j++)
	      {
		fr = zero;
		for (int i = high; i >= m; i--)
		  fr += conj(Q(i))*A(i,j);
		
		fr /= h;
		for (int i = m; i <= high; i++)
		  A(i,j) -= fr*Q(i);
	      }
	    	    
	    // .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
	    for (int i = 0; i <= high; i++)
	      {
		fr = zero;
		// .......... for j=igh step -1 until m do -- ..........
		for (int j = high; j >= m; j--)
		  fr += Q(j)*A(i,j);
		
		fr /= h;
		for (int j = m; j <= high; j++)
		  A(i,j) -= conj(Q(j))*fr;
	      }
	    
	    Q(m) *= scale;
	    A(m,m-1) *= -g;
	  }
      }
  }
  
  
  template<class Matrix1, class Vector1, class reel>
  int ComputeEigenValuesQR_Complex(int low, int high, Matrix1& H, Vector1& W,
                                   const reel& norm_tmp)
  {
    typedef typename Matrix1::value_type complexe;
    complexe zero(reel(0), reel(0)), y, s, t, x, zz, z;
    reel norm(0), tst1(0), tst2(0), zero_reel(0);
    
    y = zero; s = zero; t = zero; x = zero; zz = zero; z = zero;
    W.Fill(zero);
    int n = H.GetM();
    int ierr = 0, en, itn, its(0), enm1(0), lp1;
    
    if (low != high) 
      {
	
	//  .......... create real subdiagonal elements ..........
	int l = low + 1;
	
	for (int i = l; i <= high; i++)
	  {
	    int ll = min(i+1, high);
	    if (imag(H(i,i-1)) != zero_reel)
	      {	      
		if (H(i, i-1) == zero)
		  norm = zero_reel;
		else
		  norm = abs(H(i, i-1));
		y = H(i,i-1) / norm;
		H(i,i-1) = complexe(norm, zero_reel);
				
		for (int j = i; j <= high;j++)
		  {
		    H(i,j) = conj(y)*H(i,j);
		  }
		
		for (int j = low; j <= ll; j++)
		  {
		    H(j,i) = y*H(j,i);
		  }
	      }
	  }
	
	//     .......... store roots isolated by cbal ..........
	for (int i = 0; i < n; i++)
	  if ((i<low)|| (i>high))
	    {
	      W(i) = H(i,i);
	    }
	
	en = high;
	t = zero;
	itn = 30*n;
	bool init_iteration = true;
	//     .......... search for next eigenvalue ..........
	while (en >= low)
	  {
	    if (init_iteration)
	      {
		its = 0;
		enm1 = en - 1;
	      }
	    //     .......... look for single small sub-diagonal element
	    //                for l=en step -1 until low d0 -- ..........
	    for (int ll = low; ll <= en; ll++)
	      {
		l = en + low - ll;
		if (l != low)
		  {
		    if (H(l-1,l-1)==zero)
		      tst1 = zero_reel;
		    else
		      tst1 = abs(H(l-1,l-1));
		    
		    if (H(l,l)!=zero)
		      tst1 += abs(H(l,l));
		    
		    tst2 = tst1 + abs(real(H(l,l-1)));
		    if (tst2==tst1)
		      break;
		  }
		else
		  break;
	      }
	    
	    
	    //     .......... form shift ..........
	    if (l==en)
	      {
		W(en) = H(en,en) + t;
		en = enm1;
		init_iteration=true;
	      }
	    else
	      {
		if (itn==0) 
		  return en;
		
		if ( ((its==10)||(its==20)) && (en >= 2))
		  {
		    //     .......... form exceptional shift ..........
		    s = complexe(abs(real(H(en,enm1))) + abs(real(H(enm1,en-2))), zero_reel);
		  }
		else
		  {
		    s = H(en,en);
		    x = H(enm1,en) * real(H(en,enm1));
		    if (x != zero)
		      {
			y = (H(enm1,enm1) - s)/reel(2);
			z = y*y + x;
			zz = sqrt(z);
			tst1 = real(y*conj(zz));
			if (tst1 < zero_reel)
			  zz = -zz;
			
			x = x/(y+zz);
			s -= x;
		      }
		  }
		
		for (int i = low; i <= en; i++)
		  H(i,i) -= s;

		t += s;
		its ++;
		itn --;
		//     .......... reduce to triangle (rows) ..........
		lp1 = l + 1;
		
		for (int i = lp1; i <= en; i++)
		  {
		    s = complexe(real(H(i,i-1)), zero_reel);
		    norm = abs(complexe(abs(H(i-1, i-1)), real(s)));
		    x = H(i-1,i-1) / norm;
		    W(i-1) = x;
		    H(i-1, i-1) = complexe(norm, zero_reel);
		    H(i, i-1) = complexe(zero_reel, real(s)/norm);
		    
		    for (int j = i; j <= en;j++)
		      {
			y = H(i-1,j);
			zz = H(i,j);
			H(i-1, j) = conj(x)*y + zz*imag(H(i, i-1));
			H(i, j) = x*zz - y*imag(H(i,i-1));
		      }
		  }
		
		s = imag(H(en, en));
		if (s != zero)
		  {
		    if (H(en, en) == zero)
		      norm = zero_reel;
		    else
		      norm = abs(H(en, en));
		    s = H(en, en) / norm;
		    H(en,en) = complexe(norm, zero_reel);
		  }
		
		//     .......... inverse operation (columns) ..........
		for (int j = lp1; j <= en; j++)
		  {
		    x = W(j-1);
		    
		    for (int i = l; i <= j; i++)
		      {
			y = complexe(real(H(i,j-1)), zero_reel);
			zz = H(i,j);
			if (i != j)
			  {
			    y = H(i,j-1);
			    H(i,j-1) = complexe(real(H(i, j-1)),
                                                imag(x*y) + imag(H(j, j-1)) * imag(zz));
			  }
			
			H(i,j-1) = complexe(real(x*y) + imag(H(j, j-1))*real(zz), imag(H(i, j-1)));
			H(i,j) = conj(x)*zz - imag(H(j,j-1)) * y;
		      }
		  }
		
		if (imag(s) != zero_reel)
		  {
		    for (int i = l; i <= en; i++)
		      {
			y = H(i, en);
			H(i, en) = s*y;
		      }
		  }
		
		init_iteration = false;
	      }
	  }
      }
    
    return ierr;
  }
  

  template<class Matrix1, class Vector1, class Matrix2, class Matrix3, class reel>
  int ComputeEigenValuesQR_Complex(int low, int high, Matrix1& H, Matrix2& Q,
                                   Vector1& W, Matrix3& Z, const reel& norm_tmp)
  {
    typedef typename Matrix1::value_type complexe;
    complexe zero, y, s, t, x, zz, z;
    reel norm, tst1, tst2(0);

    SetComplexZero(zero);
    W.Fill(zero);
    int n = H.GetM();
    int ierr = 0, en, itn, its(0), enm1(0), lp1;
    
    Z.SetIdentity();
    
    for (int i = high-1; i >= low+1; i--)
      {
	if ((Q(i) != zero)&&(H(i,i-1) != zero))
	  {
	    norm = real(H(i,i-1)*conj(Q(i)));
	    for (int k = i+1; k <= high; k++)
	      Q(k) = H(k, i-1);
	    
	    for (int j = i; j <= high; j++)
	      {
		s = zero;
		for (int k = i; k <= high; k++)
		  s += conj(Q(k))*Z(k, j);
		
		s /= norm;
		
		for (int k = i; k <= high; k++)
		  Z(k, j) += s*Q(k);
	      }
	  }
      }

    //  .......... create real subdiagonal elements ..........
    int l = low + 1;
    
    for (int i = l; i <= high; i++)
      {
	int ll = min(i+1, high);
	if (imag(H(i,i-1)) != zero)
	  {	      
	    norm = abs(H(i, i-1));
	    y = H(i,i-1) / norm;
	    H(i,i-1) = norm;
	    
	    for (int j = i; j < n; j++)
	      H(i,j) = conj(y)*H(i,j);
	    
	    for (int j = 0; j <= ll; j++)
	      H(j,i) = y*H(j,i);
	    
	    for (int j = low; j <= high; j++)
	      Z(j, i) = y*Z(j,i);
	  }
      }
	
    
    //     .......... store roots isolated by cbal ..........
    for (int i = 0; i < n; i++)
      if ((i<low)|| (i>high))
	{
	  W(i) = H(i,i);
	}

    en = high;
    t = zero;
    itn = 30*n;
    bool init_iteration = true;
    //     .......... search for next eigenvalue ..........
    while (en >= low)
      {
	if (init_iteration)
	  {
	    its = 0;
	    enm1 = en - 1;
	  }
	
	//     .......... look for single small sub-diagonal element
	//                for l=en step -1 until low d0 -- ..........
	for (int ll = low; ll <= en; ll++)
	  {
	    l = en + low - ll;
	    if (l != low)
	      {
		tst1 = abs(real(H(l-1,l-1))) + abs(imag(H(l-1,l-1)))
		  + abs(real(H(l,l))) + abs(imag(H(l,l)));
		tst2 = tst1 + abs(real(H(l,l-1)));
		if (tst2==tst1)
		  break;
	      }
	    else
	      break;
	  }
	    
	
	//     .......... form shift ..........
	if (l==en)
	  {
	    H(en, en) += t;
	    W(en) = H(en,en);
	    en = enm1;
	    init_iteration = true;
	  }
	else
	  {
	    if (itn==0) 
	      return en;
	    
	    if ((its==10)||(its==20))
	      {
		//     .......... form exceptional shift ..........
		s = abs(real(H(en,enm1))) + abs(real(H(enm1,en-2)));
	      }
	    else
	      {
		s = H(en,en);
		x = H(enm1,en) * real(H(en,enm1));
		if (x != zero)
		  {
		    y = reel(0.5)*(H(enm1,enm1) - s);
		    z = y*y + x;
		    zz = sqrt(z);
		    tst1 = real(y*conj(zz));
		    if (tst1 < reel(0))
		      zz = -zz;
			
		    x = x/(y+zz);
		    s -= x;
		  }
	      }
		
	    for (int i = low; i <= en; i++)
	      H(i,i) -= s;
	    
	    t += s;
	    its ++;
	    itn --;
	    //     .......... reduce to triangle (rows) ..........
	    lp1 = l + 1;
	    
	    for (int i = lp1; i <= en; i++)
	      {
		s = real(H(i,i-1));
		norm = abs(complexe(abs(H(i-1, i-1)), real(s)));
		x = H(i-1,i-1) / norm;
		W(i-1) = x;
		H(i-1, i-1) = norm;
		H(i, i-1) = complexe(0, real(s)/norm);
		
		for (int j = i; j < n;j++)
		  {
		    y = H(i-1,j);
		    zz = H(i,j);
		    H(i-1, j) = conj(x)*y + zz*imag(H(i, i-1));
		    H(i, j) = x*zz - y*imag(H(i,i-1));
		  }
	      }
		
	    s = imag(H(en, en));
	    if (abs(s) != zero)
	      {
		norm = abs(H(en, en));
		s = H(en, en) / norm;
		H(en,en) = norm;
		
		for (int j = en+1; j < n; j++)
		  H(en, j) *= conj(s);
	      }
	    
	    init_iteration = false;
	    //     .......... inverse operation (columns) ..........
	    for (int j = lp1; j <= en; j++)
	      {
		x = W(j-1);
		
		for (int i = 0; i <= j; i++)
		  {
		    y = real(H(i,j-1));
		    zz = H(i,j);
		    if (i != j)
		      {
			y = H(i,j-1);
			H(i,j-1) = complexe(real(H(i, j-1)),
                                            imag(x*y) + imag(H(j, j-1)) * imag(zz));
		      }
		    
		    H(i,j-1) = complexe(real(x*y) + imag(H(j, j-1))*real(zz), imag(H(i, j-1)));
		    H(i,j) = conj(x)*zz - imag(H(j,j-1)) * y;
		  }
		
		for (int i = low; i <= high; i++)
		  {
		    y = Z(i, j-1);
		    zz = Z(i, j);
		    Z(i, j-1) = x*y + zz*imag(H(j, j-1));
		    Z(i, j) = conj(x)*zz - y*imag(H(j, j-1));
		  }
	      }
	    init_iteration = false;
		
	    if (abs(imag(s)) != zero)
	      {
		for (int i = 0; i <= en; i++)
		  {
		    y = H(i, en);
		    H(i, en) = s*y;
		  }
		
		for (int i = low; i <= high; i++)
		  Z(i, en) *= s;
	      }
	
	    init_iteration = false;
	  }
      }
    
    // .......... all roots found.  backsubstitute to find
    //                vectors of upper triangular form ..........
    norm = reel(0);
    
    reel tr;
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
	{
	  tr = abs(real(H(i, j))) + abs(imag(H(i, j)));
	  if (tr > norm)
	    norm = tr;
	}
    
    for (int nn = 1; nn < n; nn++)
      {
	en = n - nn;
	x = W(en);
	H(en, en) = 1.0;
	enm1 = en - 1;
	for (int i = en-1; i >= 0; i--)
	  {
	    zz = zero;
	    for (int j = i+1; j <= en; j++)
	      zz += H(i, j)*H(j, en);
	    
	    y = x - W(i);
	    reel yr = real(y), yi = imag(y);
	    if (y == zero)
	      {
		tst1 = norm;
		yr = tst1;
		while (tst2 > tst1)
		  {
		    yr = 0.01*yr;
		    tst2 = norm + yr;
		  }
	      }
	    y = complexe(yr, yi);
	    
	    H(i, en) = zz/y;
	    tr = abs(real(H(i, en))) + abs(imag(H(i, en)));
	    if (tr != zero)
	      {
		tst1 = tr;
		tst2 = tst1 + 1.0/tst1;
		if (tst2 <= tst1)
		  for (int j = i; j <= en; j++)
		    H(j, en) /= tr;
	      }
	  }
      }
	
    //.......... end backsubstitution ..........
    //     .......... vectors of isolated roots ..........
    for (int i = 0; i < n; i++)
      {
	if ((i < low)||(i > high))
	  for (int j = i; j < n; j++)
	    Z(i, j) = H(i, j);	
      }	

    //     .......... multiply by transformation matrix to give
    //                vectors of original full matrix.
    for (int j = n-1; j >= low; j--)
      {
	int m = min(j, high);
	for (int i = low; i <= high; i++)
	  {
	    zz = zero;
	    for (int k = low; k <= m; k++)
	      zz += Z(i,k)*H(k,j);
	    
	    
	    Z(i, j) = zz;
	  }
      }
    
    // normalizing eigenvectors
    for (int j = 0; j < n; j++)
      {
	norm = 0;
	for (int i = 0; i < n; i++)
	  norm += abs(Z(i,j)*conj(Z(i,j)));
	
	norm = reel(1)/sqrt(norm);
	for (int i = 0; i < n; i++)
	  Z(i,j) *= norm;
      }
   
    return ierr;
  }
  
  
  //! computation of eigenvalues of a general complex linear system
  /*!
    \param[inout] A original matrix, this matrix is modified during computation
    \param[out] L eigenvalues
   */
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<T>, Prop, Storage, Allocator1>& A,
		      Vector<complex<T>, VectFull, Allocator2>& L, LapackInfo& info)
  {
    complex<T> zero(T(0), T(0));
    int n = A.GetM();
    if (n <= 0)
      return;
    
    L.Reallocate(n); L.Fill(zero);
    Vector<complex<T> > Q(n); Q.Fill(zero);
    
    if (n == 1)
      {
	L.Reallocate(n);
	L(0) = A(0,0);
      }
    else
      {
	int is1 = 0, is2 = n-1;
	TransformMatrix_to_HessenbergForm(is1, is2, A, Q, T(0));
	ComputeEigenValuesQR_Complex(is1, is2, A, L, T(0));
      }
  }
  
  
  //! computation of eigenvalues and eigenvectors of a general complex matrix
  /*!
    \param[inout] A original matrix, this matrix is modified during computation
    \param[out] L eigenvalues
    \param[out] V eigenvectors stored as a matrix
    Each column of V is an eigenvector of matrix A
   */
  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<T>, General, Storage, Allocator1>& A,
				  Vector<complex<T>, Vect_Full, Allocator2>& L,
				  Matrix<complex<T>, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info)
  {
    complex<T> zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    int n = A.GetM();
    if (n <= 0)
      return;
    
    L.Reallocate(n); L.Fill(zero);
    Vector<complex<T> > Q(n); Q.Fill(zero);
    V.Reallocate(n, n); V.Fill(zero);

    if (n == 1)
      {
	L.Reallocate(n);
	L(0) = A(0,0);
	V(0,0) = one;
      }
    else
      {
	int is1 = 0, is2 = n-1;
	TransformMatrix_to_HessenbergForm(is1, is2, A, Q, T(0));
	ComputeEigenValuesQR_Complex(is1, is2, A, Q, L, V, T(0));	
      }
  }
  

  /////////////////////////
  // General real matrix //
  /////////////////////////
  
  
  //! computation of a general real matrix
  /*!
    \param[in] A original matrix, modified during computation
    \param[out] Lreal real part of eigenvalues
    \param[out] Limag imaginary part of eigenvalues
    \param[out] V eigenvectors stored as matrix
    Each column of V is an eigenvector.
    For complex pairs, the first column contains the real
    part and second column the imaginary part 
    (and you have the two eigenvectors v1 + i v2, v1 - i v2)
    See dgeev Lapack routine for more details
   */
  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Allocator3, class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<T, General, Storage, Allocator1>& A,
				  Vector<T, VectFull, Allocator2>& Lreal,
				  Vector<T, VectFull, Allocator3>& Limag,
				  Matrix<T, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info)
  {
    int n = A.GetM();
    if (n <= 0)
      return;
    
    T one(1), zero(0);
    
    if (n == 1)
      {
	Lreal.Reallocate(1); Lreal(0) = A(0,0);
	Limag.Reallocate(1); Limag(0) = zero;
	V.Reallocate(1, 1); V(0,0) = one;
	return;
      }

    T x, y; IVect permut(n);    
    Lreal.Reallocate(n); Lreal.Fill(zero);
    Limag.Reallocate(n); Limag.Fill(zero);
    V.Reallocate(n, n); V.SetIdentity();

    int low = 0, high = n-1;

    for (int m = low+1; m <= high-1; m++)
      {
	x = 0;
	int i = m;
	for (int j = m; j <= high; j++)
	  {
	    if (abs(A(j, m-1)) > abs(x))
	      {
		x = A(j, m-1);
		i = j;
	      }
	  }
	
	permut(m) = i;
	if (i != m)
	  {
	    // interchange rows and columns of a
	    for (int j = m-1; j < n; j++)
	      {
		y = A(i, j);
		A(i, j) = A(m, j);
		A(m, j) = y;
	      }
	    
	    for (int j = 0; j <= high; j++)
	      {
		y = A(j, i);
		A(j, i) = A(j, m);
		A(j, m) = y;
	      }
	    
	  }
	
	if (x != zero)
	  {
	    for (i = m+1; i <= high; i++)
	      {
		y = A(i, m-1);
		if (y != zero)
		  {
		    y /= x;
		    A(i, m-1) = y;
		    for (int j = m; j < n; j++)
		      A(i, j) -= y*A(m, j);
		    
		    for (int j = 0; j <= high; j++)
		      A(j, m) += y*A(j, i);
		  }
	      }		
	  }	
      }
    
    for (int mp = high-1; mp >= low+1; mp--)
      {
	for (int i = mp+1; i <= high; i++)
	  V(i, mp) = A(i, mp-1);
	
	int i = permut(mp);
	if (i != mp)
	  {
	    for (int j = mp; j <= high; j++)
	      {
		V(mp, j) = V(i, j);
		V(i, j) = zero;
	      }
	    V(i, mp) = one;
	  }	
      }
    
    // norm of the matrix
    T norm = 0.0;
    for (int i = 0; i < n; i++)
      for (int j = max(0,i-1); j < n; j++)
	norm += abs(A(i, j));

    int en = high, m;
    T t = zero, s(0), tst1, tst2, w, zz(0), r(0), p(0), q(0);
    int itn = 30*n, its = 0, na = 0, enm2 = 0;
    bool init_iteration = true;
    // search for next eigenvalues
    while (en >= low)
      {
	if (init_iteration)
	  {
	    its = 0;
	    na = en-1;
	    enm2 = na-1;
	  }

	int l = en; init_iteration = false;
	for (int ll = en; ll >= low; ll--)
	  {
	    l = ll;
	    if (l == low)
	      break;
	    
	    s = abs(A(ll-1, ll-1)) + abs(A(ll, ll));
	    if (s == zero)
	      s = norm;
	    
	    tst1 = s;
	    tst2 = tst1 + abs(A(ll, ll-1));
	    if (tst2 == tst1)
	      break;
	  }
	
	// form shift
	x = A(en, en);
	
	if (l != en)
	  {
	    y = A(na, na);
	    w = A(en, na)* A(na, en);
	    if (itn == 0)
	      return;
	    
	    if (l != na)
	      {
		if ((its == 10)||(its ==20))
		  {
		    // form exceptional shift
		    t += x;
		    for (int i = low; i <= en; i++)
		      A(i, i) -= x;
		    
		    s = abs(A(en, na)) + abs(A(na, enm2));
		    x = 0.75*s; y = x;
		    w = -0.4375*s*s;
		  }
		
		its ++; itn--;
		// look for two consecutive small sub-diagonal elements
		int mp2 = en;
		for (m = en-2; m >= l; m--)
		  {
		    mp2 = m+2;
		    zz = A(m, m);
		    r = x - zz;
		    s = y - zz;
		    p = (r * s - w) / A(m+1,m) + A(m,m+1);
		    q = A(m+1,m+1) - zz - r - s;
		    r = A(m+2,m+1);
		    s = abs(p) + abs(q) + abs(r);
		    p = p / s;
		    q = q / s;
		    r = r / s;
		    if (m == l)
		      break;
		    
		    tst1 = abs(p)*abs(A(m-1, m-1)) + abs(zz) + abs(A(m+1, m+1));
		    tst2 = tst1 + abs(A(m, m-1))*(abs(q) + abs(r));
		    if (tst1 == tst2)
		      break;
		  }
		
		for (int i = mp2; i <= en; i++)
		  {
		    A(i, i-2) = zero;
		    if (i != mp2)
		      A(i, i-3)  = zero;
		  }
		
		// double qr step involving rows l to en and columns m to na
		for (int k = m; k <= na; k++)
		  {
		    bool notlas = (k != na);
		    if (k != m)
		      {
			p = A(k, k-1);
			q = A(k+1, k-1);
			r = zero;
			if (notlas)
			  r = A(k+2, k-1);
			
			x = abs(p) + abs(q) + abs(r);
			if (x == zero)
			  continue;
			
			p /= x; q /= x; r /= x;
		      }
		    s = Sqrt(p*p+q*q+r*r);
		    if (p < 0)
		      s = -s;
		    
		    if (k != m)
		      A(k, k-1) = -s*x;
		    else
		      if (l != m)
			A(k, k-1) = -A(k, k-1);
		    
		    p += s; x = p/s; y = q/s;
		    zz = r/s;
		    q /= p; r /= p;
		    
		    if (!notlas)
		      {
			for (int j = k; j < n; j++)
			  {
			    p = A(k,j) + q * A(k+1,j);
			    A(k,j) = A(k,j) - p * x;
			    A(k+1,j) = A(k+1,j) - p * y;
			  }
			
			for (int i = 0; i <= min(en, k+3); i++)
			  {
			    p = x * A(i,k) + y * A(i,k+1);
			    A(i,k) = A(i,k) - p;
			    A(i,k+1) = A(i,k+1) - p * q;
			  }
			
			for (int i = low; i <= high; i++)
			  {
			    p = x * V(i,k) + y * V(i,k+1);
			    V(i,k) = V(i,k) - p;
			    V(i,k+1) = V(i,k+1) - p * q;
			  }
		      }
		    else
		      {
			// row modification
			for (int j = k; j < n; j++)
			  {
			    p = A(k,j) + q * A(k+1,j) + r * A(k+2,j);
			    A(k,j) = A(k,j) - p * x;
			    A(k+1,j) = A(k+1,j) - p * y;
			    A(k+2,j) = A(k+2,j) - p * zz;
			  }
			
			
			// column modification ..........
			for (int i = 0; i <= min(en, k+3); i++)
			  {
			    p = x * A(i,k) + y * A(i,k+1) + zz * A(i,k+2);
			    A(i,k) = A(i,k) - p;
			    A(i,k+1) = A(i,k+1) - p * q;
			    A(i,k+2) = A(i,k+2) - p * r;
			  }
			
			// accumulate transformations
			for (int i = 0; i <= high; i++)
			  {
			    p = x * V(i,k) + y * V(i,k+1) + zz * V(i,k+2);
			    V(i,k) = V(i,k) - p;
			    V(i,k+1) = V(i,k+1) - p * q;
			    V(i,k+2) = V(i,k+2) - p * r;
			  }

		      }
		  }
		
	    
	      }
	    else
	      {
		// two roots found
		p = (y - x) / 2;
		q = p * p + w;
		zz = Sqrt(abs(q));
		A(en,en) = x + t;
		x = A(en,en);
		A(na,na) = y + t;
		if (q >= zero)
		  {
		    // real pair
		    if (p < 0)
		      zz = p - zz;
		    else
		      zz = p + zz;
		    Lreal(na) = x + zz;
		    Lreal(en) = Lreal(na);
		    if (zz != zero)
		      Lreal(en) = x - w / zz;
		    
		    Limag(na) = zero;
		    Limag(en) = zero;
		    x = A(en,na);
		    s = abs(x) + abs(zz);
		    p = x / s;
		    q = zz / s;
		    r = Sqrt(p*p+q*q);
		    p = p / r;
		    q = q / r;
		    // row modification
		    for (int j = na; j < n; j++)
		      {
			zz = A(na,j);
			A(na,j) = q * zz + p * A(en,j);
			A(en,j) = q * A(en,j) - p * zz;
		      }
		    // column modification
		    for (int i = 0; i <= en; i++)
		      {
			zz = A(i,na);
			A(i,na) = q * zz + p * A(i,en);
			A(i,en) = q * A(i,en) - p * zz;
		      }
		      
		    // accumulate transformations
		    for (int i = low; i <= high; i++)
		      {
			zz = V(i,na);
			V(i,na) = q * zz + p * V(i,en);
			V(i,en) = q * V(i,en) - p * zz;
		      }
		  }
		else
		  {
		    // complex pair
		    Lreal(na) = x + p;
		    Lreal(en) = x + p;
		    Limag(na) = zz;
		    Limag(en) = -zz;
		  }
		
		en = enm2;
		init_iteration = true;
	      }
	  }
	else
	  {
	    // one root found
	    A(en,en) = x + t;
	    Lreal(en) = A(en,en);
	    Limag(en) = zero;
	    en = na;
	    init_iteration = true;
	  }
      }
    
    // all roots found, backsubstitue to find vectors of upper triangular form
    if (norm == zero)
      return;
    
    complex<T> zcplx; T ra, sa, vr, vi;
    for (en = n-1; en >= 0; en--)
      {
	p = Lreal(en);
	q = Limag(en);
	na = en - 1;
	if (q == zero)
	  {
	    // real vector
	    m = en;
	    A(en, en) = one;
	    if (na != -1)
	      {
		for (int i = en-1; i >= 0; i--)
		  {
		    w = A(i, i) - p;
		    r = zero;
		    for (int j = m; j <= en; j++)
		      r += A(i, j) * A(j, en);
		    
		    if (Limag(i) >= zero)
		      {
			m = i;
			if (Limag(i) != zero)
			  {
			    // solve real equations
			    x = A(i, i+1);
			    y = A(i+1, i);
			    q = (Lreal(i) - p) * (Lreal(i) - p) + Limag(i)*Limag(i);
			    t = (x * s - zz * r) / q;
			    A(i, en) = t;
			    if (abs(x) > abs(zz))
			      {
				A(i+1, en) = (-r - w * t) / x;
			      }
			    else
			      {
				A(i+1, en) = (-s - y * t) / zz;
			      }
			  }
			else
			  {
			    t = w;
			    if (t == zero)
			      {
				tst1 = norm;
				t = tst1; tst2 = 2*tst1 + one;
				while (tst2 > tst1)
				  {
				    t *= one/100;
				    tst2 = norm + t;
				  }
			      }
			    A(i, en) = -r/t;
			  }
			
			t = abs(A(i, en));
			if (t != zero)
			  {
			    tst1 = t;
			    tst2 = tst1 + one/tst1;
			    if (tst2 <= tst1)
			      for (int j = i; j <= en; j++)
				A(j,en) = A(j,en)/t;
			  }
		      }
		    else
		      {
			zz = w;
			s = r;
		      }
		  }
	      }
	  }
	else if (q < 0)
	  {
	    // complex vector
	    m = na;
	    // last vector component chosen imaginary so that
	    // eigenvector matrix is triangular
	    if (abs(A(en,na)) <= abs(A(na,en)))
	      {
		zcplx = complex<T>(zero, -A(na, en))/complex<T>(A(na, na)-p,q);
		A(na, na) = real(zcplx); A(na, en) = imag(zcplx);
	      }
	    else
	      {
		A(na, na) = q / A(en,na);
		A(na,en) = -(A(en,en) - p) / A(en,na);
	      }
	    
	    A(en, na) = zero;
	    A(en, en) = one;
	    enm2 = na - 1;
	    for (int i = en-2; i >= 0; i--)
	      {
		w = A(i, i) - p;
		ra = zero; sa = zero;
		for (int j = m; j <= en; j++)
		  {
		    ra += A(i, j) * A(j, na);
		    sa += A(i, j) * A(j, en);
		  }
		
		if (Limag(i) >= 0)
		  {
		    m = i;
		    if (Limag(i) != zero)
		      {
			x = A(i, i+1);
			y = A(i+1, i);
			vr = (Lreal(i) - p) * (Lreal(i) - p) + Limag(i) * Limag(i) - q * q;
			vi = (Lreal(i) - p) * 2 * q;
			if ((vr == zero) && (vi == zero))
			  {
			    tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz));
			    vr = tst1; tst2 = 2*tst1 + one;
			    while (tst2 > tst1)
			      {
				vr *= 0.01;
				tst2 = tst1 + vr;
			      }
			  }
			zcplx = complex<T>(x*r-zz*ra+q*sa, x*s-zz*sa-q*ra) / complex<T>(vr, vi);
			A(i, na) = real(zcplx); A(i, en) = imag(zcplx);
			if (abs(x) > abs(zz) + abs(q))
			  {
			    A(i+1,na) = (-ra - w * A(i,na) + q * A(i,en)) / x;
			    A(i+1,en) = (-sa - w * A(i,en) - q * A(i,na)) / x;
			  }
			else
			  {
			    zcplx = complex<T>(-r-y*A(i,na), -s-y*A(i,en))/complex<T>(zz,q);
			    A(i+1,na) = real(zcplx); A(i+1,en) = imag(zcplx);
			  }
		      }
		    else
		      {
			zcplx = complex<T>(-ra, -sa)/complex<T>(w, q);
			A(i, na) = real(zcplx); A(i, en) = imag(zcplx);
		      }
		    t = max(abs(A(i, na)), abs(A(i, en)));
		    if (t != zero)
		      {
			tst1 = t;
			tst2 = tst1 + one/tst1;
			if (tst2 <= tst1)
			  for (int j = i; j <= en; j++)
			    {
			      A(j, na) /= t;
			      A(j, en) /= t;
			    }
		      }
		  }
		else
		  {
		    zz = w; 
		    r = ra;
		    s = sa;
		  }
	      }
	  }
      }

    // end back substitution.
    // vectors of isolated roots
    for (int i = 0; i < n-1; i++)
      {
	if ((i < low) || (i > high))
	  for (int j = i; j < n; j++)
	    V(i, j) = A(i, j);
      }
    
    // multiply by transformation matrix to give
    // vectors of original full matrix.
    for (int j = n-1; j >= low; j--)
      {
	m = min(j, high);
	for (int i = low; i <= high; i++)
	  {
	    zz = zero;
	    for (int k = low; k <= m; k++)
	      zz += V(i, k)*A(k, j);
	    
	    V(i, j) = zz;
	  }
      }
    
    // normalizing eigenvectors
    m = 0;
    while (m < n)
      {
	if (Limag(m) == zero)
	  {
	    norm = 0;
	    for (int i = 0; i < n; i++)
	      norm += V(i,m)*V(i,m);
	
	    norm = one/Sqrt(norm);
	    for (int i = 0; i < n; i++)
	      V(i,m) *= norm;
	    
	    m++;
	  }
	else
	  {
	    norm = 0;
	    for (int i = 0; i < n; i++)
	      norm += V(i,m)*V(i,m) + V(i,m+1)*V(i,m+1);
	
	    norm = one/Sqrt(norm);
	    for (int i = 0; i < n; i++)
	      {
		V(i,m) *= norm;
		V(i,m+1) *= norm;
	      }
	    
	    m += 2;
	  }
      }

  }

  
  
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2,
	   class Allocator3>
  void GetEigenvalues(Matrix<T, Prop, Storage, Allocator1>& A,
                      Vector<T, VectFull, Allocator2>& Lreal,
                      Vector<T, VectFull, Allocator3>& Limag,
		      LapackInfo& info)
  {
    cout << "Not implemented " << endl;
    abort();
  }

  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<T, Prop, Storage, Allocator1>& A,
                      Vector<T, VectFull, Allocator2>& Lreal, LapackInfo& info)
  {
    cout << "Not implemented " << endl;
    abort();
  }

  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<T, Symmetric, Storage, Allocator1>& A,
				  Vector<T, VectFull, Allocator2>& L,
				  Matrix<T, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info)
  {
    cout << "not implemented" << endl;
    abort();
  }


  template<class T, class Storage, class Allocator1, class Allocator2,
	   class Prop2, class Storage2, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<T>, Symmetric, Storage, Allocator1>& A,
				  Vector<complex<T>, VectFull, Allocator2>& L,
				  Matrix<complex<T>, Prop2, Storage2, Allocator4>& V,
				  LapackInfo& info)
  {
    cout << "not implemented" << endl;
    abort();
  }
#endif
  
  
  //! Computes f(A) by using decomposition in eigenvalues of A
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrixReal(Matrix<T, General, Storage, Allocator1>& A, 
			     complex<T> (&f)(const complex<T>&))
  {
    int n = A.GetM();
    Vector<T> Lr(n), Li(n);
    Vector<complex<T> > scaleL(n);
    T zero(0);
    Matrix<T, General, Storage, Allocator1> Pr(n, n), invPr;
    Matrix<complex<T>, General, Storage> P, invP, Atmp;
    
    GetEigenvaluesEigenvectors(A, Lr, Li, Pr);

    // checking if all eigenvalues are real
    bool real_eig = true;
    for (int i = 0; i < n; i++)
      if (Li(i) != zero)
	real_eig = false;

    if (real_eig)
      {
	invPr = Pr;
	GetInverse(invPr);
	
	for (int i = 0; i < Lr.GetM(); i++)
	  Lr(i) = realpart(f(complex<T>(Lr(i), 0)));
	
	ScaleLeftMatrix(invPr, Lr);
	
	Mlt(Pr, invPr, A);
      }
    else
      {
	// transforming matrix Pr into the complex version
	A.Clear();
	P.Reallocate(n, n);
	int j = 0;
	while (j < n)
	  {
	    if (Li(j) == zero)
	      {
		for (int i = 0; i < n; i++)
		  P(i, j) = complex<T>(Pr(i, j), zero);
		
		j++;
	      }
	    else
	      {
		for (int i = 0; i < n; i++)
		  {
		    P(i, j) = complex<T>(Pr(i, j), Pr(i, j+1));
		    P(i, j+1) = complex<T>(Pr(i, j), -Pr(i, j+1));
		  }
		
		j += 2;
	      }
	  }    
	
	// computes in complex numbers P |L| P^{-1}
	Pr.Clear();
	invP = P;
	GetInverse(invP);
	
	for (int i = 0; i < Lr.GetM(); i++)
	  scaleL(i) = f(complex<T>(Lr(i), Li(i)));
	
	ScaleLeftMatrix(invP, scaleL);
	
	Atmp.Reallocate(n, n);
	Mlt(P, invP, Atmp);
	P.Clear(); invP.Clear();
	
	A.Reallocate(n, n);
	for (int i = 0; i < n; i++)
	  for (j = 0; j < n; j++)
	    A(i, j) = realpart(Atmp(i, j));
      }
  }


  //! Computes f(A) by using decomposition in eigenvalues of A
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrix(Matrix<complex<T>, General, Storage, Allocator1>& A,
			 complex<T> (&f)(const complex<T>&))
  {
    int n = A.GetM();
    Vector<complex<T> > L(n);
    Matrix<complex<T>, General, Storage, Allocator1> P(n, n), invP;
    
    GetEigenvaluesEigenvectors(A, L, P);
    
    invP = P;
    GetInverse(invP);
    
    for (int i = 0; i < n; i++)
      L(i) = f(L(i));
    
    ScaleLeftMatrix(invP, L);
    
    Mlt(P, invP, A);
  }


  //! Computes f(A) by using decomposition in eigenvalues
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrixReal(Matrix<T, Symmetric, Storage, Allocator1>& A, complex<T> (&f)(const complex<T>&))
  {
    int n = A.GetM();
    T one, zero; SetComplexOne(one);
    SetComplexZero(zero);
    Vector<T> d(n);
    Matrix<T, General, typename
      ClassStorageTypes<Storage>::DenseStorage> P(n, n), abs_d_mlt_p, Atmp;
    
    GetEigenvaluesEigenvectors(A, d, P);
    
    abs_d_mlt_p.Reallocate(n, n);
    Atmp.Reallocate(n, n);
    
    // eigenvectors orthonormals, the inverse of P is its transpose
    // we compute P f(D) P^{-1} 
    for (int i = 0; i < n; i++)
      {
	d(i) = realpart(f(complex<T>(d(i), 0)));
	for (int j = 0; j < n; j++)
	  abs_d_mlt_p(i,j) = P(j,i)*d(i);
      }
    
    MltAdd(one, SeldonNoTrans, P, SeldonNoTrans, abs_d_mlt_p, zero, Atmp);
    A.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
	A.Get(i, j) = Atmp(i, j);
    
  }


  //! Computes f(A) by using decomposition in eigenvalues of A
  template<class T, class Allocator1, class Storage>
  void GetFunctionMatrix(Matrix<complex<T>, Symmetric, Storage, Allocator1>& A,
			 complex<T> (&f)(const complex<T>&))
  {
    int n = A.GetM();
    Vector<complex<T> > L(n);
    Matrix<complex<T>, General, typename
      ClassStorageTypes<Storage>::DenseStorage> P(n, n), invP, Atmp;
    
    GetEigenvaluesEigenvectors(A, L, P);
    
    invP = P;
    GetInverse(invP);
    
    for (int i = 0; i < n; i++)
      L(i) = f(L(i));
    
    ScaleLeftMatrix(invP, L);

    Atmp.Reallocate(n, n);
    Mlt(P, invP, Atmp);

    A.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
	A.Get(i, j) = Atmp(i, j);    
  }


  //! Computes |A| by using decomposition in eigenvalues of A
  template<class T, class Prop, class Storage, class Allocator>
  void GetAbsoluteValue(Matrix<T, Prop, Storage, Allocator>& A, bool take_abs)
  {
    if (take_abs)
      GetFunctionMatrixReal(A, abs_cplx);
    else
      GetFunctionMatrixReal(A, abs_real);
  }
  

  //! Computes |A| by using decomposition in eigenvalues of A
  template<class T, class Prop, class Storage, class Allocator>
  void GetAbsoluteValue(Matrix<complex<T>, Prop, Storage, Allocator>& A, bool take_abs)
  {
    if (take_abs)
      GetFunctionMatrix(A, abs_cplx);
    else
      GetFunctionMatrix(A, abs_real);
  }
  

  //! Computes exp(A) by using decomposition in eigenvalues of A
  template<class T, class Prop, class Storage, class Allocator>
  void GetExponential(Matrix<T, Prop, Storage, Allocator>& A)
  {
    GetFunctionMatrixReal(A, exp);
  }


  //! Computes exp(A) by using decomposition in eigenvalues of A
  template<class T, class Prop, class Storage, class Allocator>
  void GetExponential(Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    GetFunctionMatrix(A, exp);
  }


  //! Computes sqrt(A) by using decomposition in eigenvalues of A
  template<class T, class Prop, class Storage, class Allocator>
  void GetSquareRoot(Matrix<T, Prop, Storage, Allocator>& A)
  {
    GetFunctionMatrixReal(A, sqrt);
  }


  //! Computes sqrt(A) by using decomposition in eigenvalues of A
  template<class T, class Prop, class Storage, class Allocator>
  void GetSquareRoot(Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    GetFunctionMatrix(A, sqrt);
  }

  //! replaces matrix A by its absolute value
  template<class T, int m>
  void GetAbsoluteValue(TinyMatrix<T, Symmetric, m, m>& A, bool take_abs)
  {
    Matrix<T, Symmetric, RowSymPacked> B(m, m);
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
	B(i, j) = A(i, j);
    
    GetAbsoluteValue(B, take_abs);

    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
	A(i, j) = B(i, j);    
  }
  
  
  //! replaces matrix A by its absolute value
  template<class T, int m>
  void GetAbsoluteValue(TinyMatrix<T, General, m, m>& A, bool take_abs)
  {
    Matrix<T, General, RowMajor> B(m, m);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
	B(i, j) = A(i, j);

    GetAbsoluteValue(B, take_abs);

    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
	A(i, j) = B(i, j);    
  }

}

#define SELDON_FILE_EIGENVALUE_CXX
#endif
