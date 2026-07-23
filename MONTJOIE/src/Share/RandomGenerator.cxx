#ifndef MONTJOIE_FILE_RANDOM_GENERATOR_CXX

namespace Montjoie
{
  /*********************
   * Random Generators *
   *********************/
  
  VirtualRandomGenerator::~VirtualRandomGenerator()
  {
  }


  //! default constructor
  VarRandomGenerator::VarRandomGenerator()
  {
    num_maximal_prime = 1000;
    
    // default random generator : rand()
    type_generator = RANDOM;

    var_random = NULL;    
  }
  
  
  //! destructor
  VarRandomGenerator::~VarRandomGenerator()
  {
    Clear();
  }
  
  
  //! releases memory used
  void VarRandomGenerator::Clear()
  {
    if (var_random != NULL)
      {
	var_random->Clear();
	var_random = NULL;
      }
    
    type_generator = RANDOM;
    prime_numbers.Clear();
    UsedPrimeNumber.Clear();
  }


  //! modifying an integer n so that n = nx*ny, where nx and ny are integers
  /*!
    ny is selected close to alpha nx
   */
  void VarRandomGenerator::RoundToSquare(int& n, const Real_wp& alpha,
					 int& nx, int& ny)
  { 
    // setting number of particles for the first divisions
    nx = toInteger(sqrt(n/alpha));
    ny = n/nx; n = nx*ny;
  }
  
  
  //! modifying an integer n so that n = nx*ny*nz, where nx, ny and nz are integers
  /*!
    ny is selected close to alpha nx
    nz is selected close to beta nx
   */
  void VarRandomGenerator::RoundToSquare(int& n, const Real_wp& alpha, const Real_wp& beta,
					 int& nx, int& ny, int& nz)
  {
    // setting number of particles for the first divisions
    nx = toInteger(pow(n/(alpha*beta), Real_wp(1)/3));
    ny = toInteger(alpha*nx); nz = toInteger(beta*nx); n = nx*ny*nz;
  }


  //! returns a random number between 0 and 1
  Real_wp VarRandomGenerator::GetRand()
  {
    switch (type_generator)
      {
      case RANDOM :
	return Real_wp(rand())/RAND_MAX;
      case BIT_REVERSED:
	cout << "not available for a single number" <<endl;
	abort();
	break;
      default :
	return var_random->GetRand();
      }
  }
  

  //! generation of random numbers with gsl function
  /*!
    \param[in] n how many random numbers are created ?
    \param[out] Rx list of random numbers
   */
  void VarRandomGenerator::GenerateRandomNumbers(int n, VectReal_wp& Rx)
  {    
    Rx.Reallocate(n);

    if (type_generator == BIT_REVERSED)
      GenerateBitReversedNumbers(n, Rx);
    else if (type_generator == RANDOM)
      {
	Rx.FillRand();
	Mlt(Real_wp(1)/RAND_MAX, Rx);
      }
    else
      var_random->GetRand(Rx);
  }


  //! obtention of a random permutation
  /*!
    the array RandomPermut will contain on output the values 
    0, 1, ..., N-1 sorted randomly
  */
  void VarRandomGenerator::GenerateRandomPermutation(int N, IVect& RandomPermut)
  {
    int num;
    RandomPermut.Reallocate(N); RandomPermut.Fill(-1);
    VectBool UsedNumber(N); UsedNumber.Fill(false);
    // generation of a first random set
    int nb_unaffected_int = 0;
    for (int i = 0; i < N; i++)      
      {
	num = toInteger(floor(GetRand()*N));
	if (!UsedNumber(num))
	  {
	    UsedNumber(num) = true;
	    RandomPermut(i) = num;
	  }
	else
	  nb_unaffected_int++;
      }
    
    // loop to affect all numbers
    while (nb_unaffected_int > 0)
      {
	int Np = nb_unaffected_int;
	IVect row_number_to_affect(nb_unaffected_int);
	IVect row_number_unused(nb_unaffected_int);	
	int nb1 = 0; int nb2 = 0;
	for (int i = 0; i < N; i++)
	  {
	    if (RandomPermut(i) < 0)
	      row_number_to_affect(nb1++) = i;
	    
	    if (!UsedNumber(i))
	      row_number_unused(nb2++) = i;
	  }
	
	nb_unaffected_int = 0;
	for (int i = 0; i < Np; i++)
	  {
	    nb1 = toInteger(floor(GetRand()*Np));
	    num = row_number_unused(nb1);
	    if (!UsedNumber(num))
	      {
		UsedNumber(num) = true;
		RandomPermut(row_number_to_affect(i)) = num;
	      }
	    else
	      nb_unaffected_int++;
	  }
      }
  }

  
  //! random permutation of a vector
  /*!
    \param[in,out] Rx the values of this vector are permuted randomly
  */
  void VarRandomGenerator::ApplyRandomPermutation(VectReal_wp& Rx)
  {
    // obtaining a random permutation
    IVect RandomPermut;
    GenerateRandomPermutation(Rx.GetM(), RandomPermut);
    
    VectReal_wp Ry = Rx;
    // now applying the permutation
    for (int i = 0; i < Rx.GetM(); i++)
      Rx(i) = Ry(RandomPermut(i));
  }

  
  //! selecting the type of quasi-random generator (uniform, bit-reversed or random)
  void VarRandomGenerator::SetRandomGenerator(int type)
  {
    // clearing previous random generator
    Clear();
    
    if (type == BIT_REVERSED)
      {
	FindPrimeNumbers(prime_numbers, num_maximal_prime);
	UsedPrimeNumber.Reallocate(prime_numbers.GetM());
	UsedPrimeNumber.Fill(false);
      }
    else if (type == RANDOM)
      srand(time(NULL));
    else
      {
#ifdef MONTJOIE_WITH_GSL
	if (var_random != NULL)
	  delete var_random;

	var_random = new VarRandomGeneratorGsl();
	var_random->Construct(type);
#else
	cout << "Recompile Montjoie with GSL" << endl;
	abort();
#endif
      }
  }
  

  //! selecting the type of quasi-random generator (uniform, bit-reversed or random)
  void VarRandomGenerator::SetInputRandomGenerator(const string& keyword)
  {
    int type = 0;
    if (!keyword.compare("BIT_REVERSED"))
      type = BIT_REVERSED;
    else if (!keyword.compare("RANDOM"))
      type = RANDOM;
    else if (keyword == "GSL_MT19937")
      type = GSL_MT19937;
    else if (keyword == "GSL_RANLXD1")
      type = GSL_RANLXD1;
    else if (keyword == "GSL_RANLXD2")
      type = GSL_RANLXD2;
    else if (keyword == "GSL_RANLUX")
      type = GSL_RANLUX;
    else if (keyword == "GSL_RANLUX389")
      type = GSL_RANLUX389;
    else if (keyword == "GSL_CMRG")
      type = GSL_CMRG;
    else if (keyword == "GSL_MRG")
      type = GSL_MRG;
    else if (keyword == "GSL_TAUS")
      type = GSL_TAUS;
    else if (keyword == "GSL_TAUS2")
      type = GSL_TAUS2;
    else if (keyword == "GSL_GFSR4")
      type = GSL_GFSR4;
    
    SetRandomGenerator(type);
  }

  
  //! generation of quasi-random numbers with the bit-reversed technique
  void VarRandomGenerator::GenerateBitReversedNumbers(int n, VectReal_wp& Rx)
  {    
    // finding a free prime number
    int base = 2;
    for (int i = 0; i < UsedPrimeNumber.GetM(); i++)
      if (!UsedPrimeNumber(i))
	{
	  UsedPrimeNumber(i) = true;
	  base = prime_numbers(i);
	  break;
	}
    
    Real_wp invBase = Real_wp(1)/base;
    
    int max_lvl = toInteger(ceil(log(n+1)/log(Real_wp(base))));
    IVect digits(max_lvl); digits.Fill(0);
    
    // quasi-random numbers with bit-reversed technique
    Rx.Reallocate(n);
    for (int i = 1; i <= n; i++) 
      {
	// decomposition of i in the base
	int n0 = i; int num = 0;
	while (n0 > 0)
	  {
	    digits(num) = n0%base;
	    n0 = n0/base;
	    num++;
	  }
	
	Real_wp x = 0;
	for (int n = 1; n <= num; n++)
	  x = invBase*(x + digits(num-n));
	
	Rx(i-1) = x;
      }
  }


#ifdef MONTJOIE_WITH_GSL
  //! default constructor
  VarRandomGeneratorGsl::VarRandomGeneratorGsl()
  {
    type_random_generator = NULL;
    instance_random_generator = NULL;
  }


  //! destructor
  VarRandomGeneratorGsl::~VarRandomGeneratorGsl()
  {
    Clear();
  }


  //! clears random generator
  void VarRandomGeneratorGsl::Clear()
  {
    if (instance_random_generator != NULL)
      {
	instance_random_generator = NULL;
	type_random_generator = NULL;
	gsl_rng_free(instance_random_generator);
      }
  }


  //! returns random numbers
  Real_wp VarRandomGeneratorGsl::GetRand()
  {
    return gsl_rng_uniform_pos(instance_random_generator);
  }


  //! fills Rx of random numbers
  void VarRandomGeneratorGsl::GetRand(VectReal_wp& Rx)
  {
    int n = Rx.GetM();
    for (int i = 0; i < n; i++) 
      Rx(i) = gsl_rng_uniform_pos(instance_random_generator);
  }


  //! Construct generator
  void VarRandomGeneratorGsl::Construct(int type)
  {
    gsl_rng_env_setup();
    if (type == VarRandomGenerator::GSL_MT19937)
      type_random_generator = gsl_rng_mt19937;
    else if (type == VarRandomGenerator::GSL_RANLXD1)
      type_random_generator = gsl_rng_ranlxd1;
    else if (type == VarRandomGenerator::GSL_RANLXD2)
      type_random_generator = gsl_rng_ranlxd2;
    else if (type == VarRandomGenerator::GSL_RANLUX)
      type_random_generator = gsl_rng_ranlux;
    else if (type == VarRandomGenerator::GSL_RANLUX389)
      type_random_generator = gsl_rng_ranlux389;
    else if (type == VarRandomGenerator::GSL_CMRG)
      type_random_generator = gsl_rng_cmrg;
    else if (type == VarRandomGenerator::GSL_MRG)
      type_random_generator = gsl_rng_mrg;
    else if (type == VarRandomGenerator::GSL_TAUS)
      type_random_generator = gsl_rng_taus;
    else if (type == VarRandomGenerator::GSL_TAUS2)
      type_random_generator = gsl_rng_ranlxd1;
    else if (type == VarRandomGenerator::GSL_GFSR4)
      type_random_generator = gsl_rng_gfsr4;
    
    instance_random_generator = gsl_rng_alloc(type_random_generator);
    
    gsl_rng_set(instance_random_generator, time(NULL));
  }
#endif

  
  /***********************************
   * Prime numbers and factorisation *
   ***********************************/
  
  
  //! computation of a bounded number of prime numbers
  /*!
    \param[in] max_num prime numbers are searched under the condition p < 2*max_num 
    \param[out] prime_numbers list of computed prime numbers
   */
  void FindPrimeNumbers(IVect& prime_numbers, int max_num)
  {
    if (max_num <= 1)
      {
        prime_numbers.Reallocate(1);
        prime_numbers(0) = 2;
        return;
      }
    
    // basic sieve of eratosthenes (even numbers are eliminated by construction)
    // first number of this array if 1
    Vector<bool> is_prime(max_num+2); is_prime.Fill(true); 
    // first prime : 2
    int prime = 2; // prime_numbers.Reallocate(max_num);
    prime_numbers.Reallocate(max_num+1);
    prime_numbers(0) = prime;
    // second prime : 3
    prime++; prime_numbers(1) = prime;
    int nb_primes = 2, index = 1;
    while (index <= max_num)
      {
	is_prime(index) = false;
	for (int i = 1; i <= max_num/prime; i++)
	  if ((prime*i+index) <= (max_num+1))
	      is_prime(prime*i+index) = false;
	
	// finding next prime number
	while ((index <= max_num)&&(is_prime(index) != true))
	  index++;
	
	if (index < max_num)
	  {
	    prime = 2*index+1;
	    prime_numbers(nb_primes) = prime;
	    nb_primes++;
	  }
      }
    
    prime_numbers.Resize(nb_primes);
  }
    

  //! Factorization of a number with prime numbers
  /*!
    \param[inout] n number to factorize
    \param[out] prime_decomp list of prime factors
    \param[in] prime_numbers prime numbers accepted in the decomposition
   */
  void FindPrimeDecomposition(int n, IVect& prime_decomp, const IVect& prime_numbers)
  {
    // loop over prime numbers
    int np = prime_numbers.GetM(); IVect pow(np); pow.Fill(0);
    int prod = 1; int nfactor = 0;
    for (int i = 0; i < np; i++)
      {
	int p = prime_numbers(i);
	if (n%p == 0)
	  {
	    pow(i) = 1; prod *= p; nfactor++;
	    int n0 = n/p;
	    while (n0%p == 0)
	      {
		pow(i)++;
		prod *= p; nfactor++;
		n0 = n0/p;
	      }
	  }
      }
    
    if (prod < n)
      prime_decomp.Reallocate(nfactor+1);
    else
      prime_decomp.Reallocate(nfactor);
    
    nfactor = 0;
    for (int i = 0; i < prime_numbers.GetM(); i++)
      for (int j = 0; j < pow(i); j++)
        prime_decomp(nfactor++) = prime_numbers(i);
    
    if (prod < n)
      prime_decomp(nfactor) = n/prod;
  }
  
  
  //! Factorization of a number with prime numbers
  /*!
    if n cannot be factorized with prime numbers of the array prime_numbers,
    n is increased to the nearest number where it is possible
    \param[inout] n number to factorize
    \param[out] prime_decomp list of prime factors
    \param[in] prime_numbers prime numbers accepted in the decomposition
   */
  void FindPrimeFactorization(int& n, IVect& prime_decomp, const IVect& prime_numbers)
  {
    // loop over prime numbers
    int np = prime_numbers.GetM(); IVect pow(np); pow.Fill(0);
    int prod = 1; int nfactor = 0;
    for (int i = 0; i < np; i++)
      {
	int p = prime_numbers(i);
	if (n%p == 0)
	  {
	    pow(i) = 1; prod *= p; nfactor++;
	    int n0 = n/p;
	    while (n0%p == 0)
	      {
		pow(i)++;
		prod *= p; nfactor++;
		n0 = n0/p;
	      }
	  }
      }
    
    if (prod != n)
      {
	// there is a prime number dividing n not in the list
	// we try with a greater n
	n++;
	FindPrimeFactorization(n, prime_decomp, prime_numbers);
      }
    else
      {
	prime_decomp.Reallocate(nfactor);
	nfactor = 0;
	for (int i = 0; i < prime_numbers.GetM(); i++)
	  for (int j = 0; j < pow(i); j++)
	    prime_decomp(nfactor++) = prime_numbers(i);
      }
    
  }

  
  //! Factorization of a number with prime numbers
  /*!
    if n cannot be factorized with prime numbers less or equal to pmax,
    n is increased to the nearest number where it is possible
    \param[inout] n number to factorize
    \param[out] prime_decomp list of prime factors
    \param[in] pmax maximal prime number considered
   */
  void FindPrimeFactorization(int& n, IVect& prime_decomp, int pmax)
  {
    IVect prime_numbers;
    FindPrimeNumbers(prime_numbers, (pmax+1)/2);
    
    FindPrimeFactorization(n, prime_decomp, prime_numbers);
  }


  //! Factorization of a number with prime numbers
  void FindPrimeDecomposition(int n, IVect& prime_decomp)
  {
    if (n <= 1)
      {
        prime_decomp.Reallocate(1);
        prime_decomp(0) = 1;
        return;
      }
    
    int pmax = toInteger(floor(sqrt(Real_wp(n))));
    IVect prime_numbers;
    FindPrimeNumbers(prime_numbers, (pmax+1)/2);
    
    FindPrimeDecomposition(n, prime_decomp, prime_numbers);
  }
  
  
  //! Finds two integers such that nx ny = n
  void FindTwoFactors(int n, int& nx, int& ny)
  {
    if (n <= 1)
      {
        nx = 1;
        ny = 1;
        return;
      }
    
    int pmax = toInteger(floor(sqrt(Real_wp(n))));
    IVect prime_numbers;
    FindPrimeNumbers(prime_numbers, (pmax+1)/2);
    
    IVect prime_decomp;
    FindPrimeDecomposition(n, prime_decomp, prime_numbers);
    nx = 1; ny = 1;
    int i = 0;
    while (ny <= pmax)
      {
        nx = ny;
        if (i < prime_decomp.GetM())
          {
            ny *= prime_decomp(i);
            i++;
          }
      }
    
    ny = n / nx;
    
  }
  
  
  //! finds the closest power of 2
  int FindClosestPow2(int N)
  {
    if (N <= 1)
      return N;

    int p = toInteger(ceil(log2(N)));
    int k = toInteger(pow(Real_wp(2.0), Real_wp(p)));
    if (k/2 >= N)
      return k/2;
    
    return k;
  }
  
  
  //! finds the closest 2^k 3^l to N
  int FindClosestPow23(int N)
  {
    if (N <= 1)
      return N;

    // brute force function
    Real_wp logN = log(Real_wp(N));
    int p2 = toInteger(ceil(logN/log(2.0)));
    int p3 = toInteger(ceil(logN/log(3.0)));
    IVect Pow2(p2+1);
    IVect Pow3(p3+1);
    Pow2(0) = 1; Pow3(0) = 1;
    for (int k = 1; k <= p2; k++)
      Pow2(k) = Pow2(k-1)*2;

    for (int k = 1; k <= p3; k++)
      Pow3(k) = Pow3(k-1)*3;

    int closest_number = 2, best_diff = N;
    // trying all combinations of 2^k 3^l until finding the closest one
    for (int k = 1; k <= p2; k++)
      for (int l = 0; l <= p3; l++)
	{
	  int n = Pow2(k)*Pow3(l);
	  if (n >= N)
	    {
	      if (n - N <= best_diff)
		{
		  closest_number = n;
		  best_diff = n-N;
		}
	    }
	}
    
    return closest_number;
  }

  
  //! finds the closest 2^k 3^l 5^m to N
  int FindClosestPow235(int N)
  {
    if (N <= 1)
      return N;

    // brute force strategy
    Real_wp logN = log(Real_wp(N));
    int p2 = toInteger(ceil(logN/log(2.0)));
    int p3 = toInteger(ceil(logN/log(3.0)));
    int p5 = toInteger(ceil(logN/log(5.0)));
    IVect Pow2(p2+1), Pow3(p3+1), Pow5(p5+1);
    Pow2(0) = 1; Pow3(0) = 1; Pow5(0) = 1;
    for (int k = 1; k <= p2; k++)
      Pow2(k) = Pow2(k-1)*2;

    for (int k = 1; k <= p3; k++)
      Pow3(k) = Pow3(k-1)*3;

    for (int k = 1; k <= p5; k++)
      Pow5(k) = Pow5(k-1)*5;

    int closest_number = 2, best_diff = N;
    // trying all combinations of 2^k 3^l 5^m until finding the closest one
    for (int k = 1; k <= p2; k++)
      for (int l = 0; l <= p3; l++)
	for (int m = 0; m <= p5; m++)
	  {
	    int n = Pow2(k)*Pow3(l)*Pow5(m);
	    if (n >= N)
	      {
		if (n - N <= best_diff)
		  {
		    closest_number = n;
		    best_diff = n-N;
		  }
	      }
	  }
    
    return closest_number;
  }


  //! finds the closest 2^k 3^l 5^m to N
  int FindClosestPow2357(int N)
  {
    if (N <= 1)
      return N;
    
    // brute force strategy
    Real_wp logN = log(Real_wp(N));
    int p2 = toInteger(ceil(logN/log(2.0)));
    int p3 = toInteger(ceil(logN/log(3.0)));
    int p5 = toInteger(ceil(logN/log(5.0)));
    int p7 = toInteger(ceil(logN/log(7.0)));
    IVect Pow2(p2+1), Pow3(p3+1), Pow5(p5+1), Pow7(p7+1);
    Pow2(0) = 1; Pow3(0) = 1; Pow5(0) = 1; Pow7(0) = 1;
    for (int k = 1; k <= p2; k++)
      Pow2(k) = Pow2(k-1)*2;

    for (int k = 1; k <= p3; k++)
      Pow3(k) = Pow3(k-1)*3;

    for (int k = 1; k <= p5; k++)
      Pow5(k) = Pow5(k-1)*5;

    for (int k = 1; k <= p7; k++)
      Pow7(k) = Pow7(k-1)*7;
    
    int closest_number = 2, best_diff = N;
    // trying all combinations of 2^k 3^l 5^m 7^n until finding the closest one
    for (int k = 1; k <= p2; k++)
      for (int l = 0; l <= p3; l++)
	for (int m = 0; m <= p5; m++)
	  for (int m7 = 0; m7 <= p7; m7++)
	    {
	      int n = Pow2(k)*Pow3(l)*Pow5(m)*Pow7(m7);
	      if (n >= N)
		{
		  if (n - N <= best_diff)
		    {
		      closest_number = n;
		      best_diff = n-N;
		    }
		}
	    }
    
    return closest_number;
  }

  
  //! extends the interval [i0, i1] such that the interval contains N elements
  /*!
    The subscripts i0 and i1 are extended within the interval <0, Nt-1>
    such that i1-i0+1 = N.
   */
  void ExtendIntegerInterval(int& i0, int& i1, int Nt, int N)
  {
    if (N > Nt)
      {
	i0 = 0;
	i1 = Nt-1;
	return;
      }
    
    int taille = i1-i0+1;
    if (taille >= N)
      return;
    
    int delta = N - taille;
    Real_wp coef0 = Real_wp(i0)/(Nt-1)*delta;
    Real_wp coef1 = Real_wp(Nt-1-i1)/(Nt-1)*delta;
    
    int p0 = toInteger(floor(coef0/(coef0+coef1)*delta));
    int p1 = toInteger(floor(coef1/(coef0+coef1)*delta));
    i0 -= p0;
    i1 += p1;
    
    if ((i1-i0+1) == N)
      return;
    else
      {
	if (i1-i0+1 == N-1)
	  {
	    if (i0 > 0)
	      {
		i0--;
		return;
	      }
	    else if (i1 < Nt-1)
	      {
		i1++;
		return;
	      }
	    else
	      {
		abort();
	      }
	  }
	else if (i1-i0+1 == N-2)
	  {
	    if (i0 > 0)
	      {
		i0--;
		if (i0 > 0)
		  {
		    i0--;
		    return;
		  }
		else if (i1 < Nt - 1)
		  {
		    i1++;
		    return;
		  }
		else
		  {
		    abort();
		  }
	      }
	    else if (i1 < Nt - 2)
	      {
		i1 += 2;
		return;
	      }
	    else
	      {
		abort();
	      }
	  }
      }
  }
  

  //! evaluation of a real number with his decomposition in continued fractions
  /*
    \param[in] decomp_coef decomposition of the real
    \param[out] coef resulting real number
  */
  template<class Integer, class T>
  void EvaluateContinuedFraction(const Vector<Integer>& decomp_coef, T& coef)
  {
    int n = decomp_coef.GetM();
    coef = T(decomp_coef(n-1));
    for (int i = (n-2); i >= 0; i--)
      coef = T(1)/coef + T(decomp_coef(i));
  }
  
  
  //! decomposition of a real number in continued fraction
  /*
    \param[in] coef real number to decompose
    \param[in] epsilon precision of the continued fraction
    \param[out] decomp_coef fractional decomposition
   */
  template<class Integer, class T>
  void DecomposeContinuedFraction(const T& coef, const T& epsilon,
				  Vector<Integer>& decomp_coef)
  {
    decomp_coef.Clear();
    Integer coef_int; T new_coef, coef_eval, coef_floor, residu(1);
    new_coef = coef;
    int nb_iter = 0;
    while ((residu > epsilon)&&(nb_iter < 10))
      {
	coef_floor = round(new_coef);
	coef_int = toInteger(coef_floor);
	decomp_coef.PushBack(coef_int); 
	new_coef = T(1)/(new_coef-coef_floor);
	EvaluateContinuedFraction(decomp_coef, coef_eval);
	residu = abs(coef-coef_eval);
	nb_iter++;
      }
    // DISP(coef); 
    // for (unsigned int i = 0; i < decomp_coef.size(); i++)
    // DISP(decomp_coef[i]);
  }
  
  
  //! computes the global fraction from the decomposition in continued fractions
  /*
    \param[in] decomp_coef decomposition in continued fractions
    \param[out] num resulting numerator
    \param[out] denom resulting denominator
   */
  template<class Integer>
  void GetNumeratorDenominator(const Vector<Integer>& decomp_coef, Integer& num, Integer& denom)
  {
    int n = decomp_coef.GetM();
    num = decomp_coef(n-1); denom = 1; Integer tmp;
    for (int i = (n-2); i >= 0; i--)
      {
	tmp = num;
	num = decomp_coef(i)*num+denom;
	denom = tmp;
      }
  }
  
}

#define MONTJOIE_FILE_RANDOM_GENERATOR_CXX
#endif  

