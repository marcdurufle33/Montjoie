#ifndef MONTJOIE_FILE_RANDOM_GENERATOR_HXX

namespace Montjoie
{

  //! base class for random generator
  class VirtualRandomGenerator
  {
  public:
    virtual ~VirtualRandomGenerator();

    virtual void Clear() = 0;

    virtual Real_wp GetRand() = 0;
    virtual void GetRand(VectReal_wp& Rx) = 0;

    virtual void Construct(int type) = 0;
    
  };
  
  
  //! class for generating random numbers
  class VarRandomGenerator
  {
  protected :
    // random generator
    //! maximum number used to generate the list of prime numbers
    int num_maximal_prime;
    //! list of prime numbers
    IVect prime_numbers;
    //! prime number already used ?
    VectBool UsedPrimeNumber;

    VirtualRandomGenerator* var_random;
    
    int type_generator;
    
  public:
    // type of random numbers
    enum{RANDOM, BIT_REVERSED, GSL_MT19937, GSL_RANLXD1, GSL_RANLXD2, GSL_RANLUX,
	 GSL_RANLUX389, GSL_CMRG, GSL_MRG, GSL_TAUS, GSL_TAUS2, GSL_GFSR4, UNIFORM};        
    
    VarRandomGenerator();
    ~VarRandomGenerator();
    
    void Clear();
    
    static void RoundToSquare(int& n, const Real_wp& alpha, int& nx, int& ny);
    static void RoundToSquare(int& n, const Real_wp& alpha,
			      const Real_wp& beta, int& nx, int& ny, int& nz);
    
    // generation of random, uniform or quasi-random numbers
    Real_wp GetRand();
    
    void GenerateRandomNumbers(int n, VectReal_wp& Rx);
    void ApplyRandomPermutation(VectReal_wp& Rx);
    void GenerateRandomPermutation(int n, IVect& num);
    
    void SetRandomGenerator(int);
    void SetInputRandomGenerator(const string& keyword);
    
    void GenerateBitReversedNumbers(int n, VectReal_wp& Rx);

  };

  
#ifdef MONTJOIE_WITH_GSL
  class VarRandomGeneratorGsl : public VirtualRandomGenerator
  {
  protected:
    
    //! GSL pointer for the generation of random numbers 
    const gsl_rng_type* type_random_generator;
    //! GSL pointer for the generation of random numbers
    gsl_rng* instance_random_generator;

  public:
    VarRandomGeneratorGsl();
    ~VarRandomGeneratorGsl();

    void Clear();

    Real_wp GetRand();
    void GetRand(VectReal_wp& Rx);

    void Construct(int type);
    
  };
#endif

  void FindPrimeNumbers(IVect& prime_numbers, int max_num);
  
  void FindPrimeDecomposition(int n, IVect& prime_decomp, const IVect& prime_numbers);
  void FindPrimeDecomposition(int n, IVect& prime_decomp);
  
  void FindPrimeFactorization(int& n, IVect& prime_decomp, const IVect& prime_numbers);
  void FindPrimeFactorization(int& n, IVect& prime_decomp, int pmax);
  
  void FindTwoFactors(int n, int& nx, int& ny);

  int FindClosestPow2(int N);
  int FindClosestPow23(int N);
  int FindClosestPow235(int N);
  int FindClosestPow2357(int N);
  void ExtendIntegerInterval(int& i0, int& i1, int Nt, int N);
  
  template<class Integer, class T>
  void EvaluateContinuedFraction(const Vector<Integer>& decomp_coef, T& coef);
  
  template<class Integer, class T>
  void DecomposeContinuedFraction(const T& coef, const T& epsilon,
				  Vector<Integer>& decomp_coef);

  template<class Integer>
  void GetNumeratorDenominator(const Vector<Integer>& decomp_coef, Integer& num, Integer& denom);
  
}

#define MONTJOIE_FILE_RANDOM_GENERATOR_HXX
#endif

  
