#ifndef MONTJOIE_FILE_MONTJOIE_TYPES_CXX

#ifndef SELDON_WITH_COMPILED_LIBRARY
double size_all_communications(0);
unsigned int Seed = 0;
bool is_seed_defined = false;
int rank_global_processor(0);
#endif

namespace Montjoie
{

#ifndef SELDON_WITH_COMPILED_LIBRARY
#ifdef MONTJOIE_WITH_MPFR
  //! Pi in working precision
  Real_wp pi_wp(0, MONTJOIE_DEFAULT_PRECISION) ;
  //! sqrt(-1) in working precision
  complex<Real_wp> Iwp(pi_wp, pi_wp);
  //! machine precision (1e-7 in simple, 2e-16 in double)
  Real_wp epsilon_machine(1e-32, MONTJOIE_DEFAULT_PRECISION);
#else

  //! pi variable at working precision
  Real_wp pi_wp(M_PI);

  //! complex number I
  complex<Real_wp> Iwp(0.0, 1.0);

  //! machine precision (1e-7 in simple, 2e-16 in double)
  Real_wp epsilon_machine(2e-16);
#endif
  
#ifdef SELDON_WITH_SCALAPACK
  BlacsHandle global_blacs_handle;
#endif
  
#endif

#ifndef SELDON_WITH_COMPILED_LIBRARY
  //! Initialization of different variables
  void InitMontjoie(int argc, char** argv, int print_level)
  {
    // initialization of MPI
#ifdef SELDON_WITH_MPI

#ifdef MONTJOIE_WITHOUT_THREAD
    MPI_Init(&argc, &argv);
#else
    int required = MPI_THREAD_MULTIPLE;
    int provided = -1;
    int rank;
    MPI_Init_thread(&argc, &argv, required, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ((rank == 0) && (print_level > 0))
      {
        switch (provided)
          {
          case MPI_THREAD_SINGLE:
            cout << "MPI_Init_thread level = MPI_THREAD_SINGLE" << endl;
            break;
          case MPI_THREAD_FUNNELED:
            cout << "MPI_Init_thread level = MPI_THREAD_FUNNELED" << endl;
            break;
          case MPI_THREAD_SERIALIZED:
            cout << "MPI_Init_thread level = MPI_THREAD_SERIALIZED" << endl;
            break;
          case MPI_THREAD_MULTIPLE:
            cout << "MPI_Init_thread level = MPI_THREAD_MULTIPLE" << endl;
            break;
          default:
            cout << "MPI_Init_thread level = ???" << endl;
          }
      }
#endif

    rank_global_processor = rank;

#endif

    // initial seed (to have something common to all processors and reproductible)
    srand(0);
    
    // setting precision for outputs and mpfr
    cout.precision(15);
    
#ifdef MONTJOIE_WITH_MPFR
    mpreal::set_default_prec(MONTJOIE_DEFAULT_PRECISION);
    mpfr_set_default_prec(MONTJOIE_DEFAULT_PRECISION);
    mpf_set_default_prec(MONTJOIE_DEFAULT_PRECISION);
#endif
    
    // constants
    InitPrecisionConstants();
    
#ifdef MONTJOIE_FILE_MESH_BASE_HXX
    // initialization of mesh constants
    Dimension2 dim; Real_wp x;
    InitStaticMeshData(dim, x);

#ifdef MONTJOIE_WITH_THREE_DIM
    Dimension3 dim3;
    InitStaticMeshData(dim3, x);
#endif
    
#ifdef MONTJOIE_FILE_VAR_FINITE_ELEMENT_HXX

#ifdef MONTJOIE_WITH_ONE_DIM
#ifdef MONTJOIE_WITH_NODAL_H1
    VarFiniteElementEnum<Dimension1, 1>::InitStaticData();
#endif
#endif

#ifdef MONTJOIE_WITH_TWO_DIM
#ifdef MONTJOIE_WITH_NODAL_H1
    VarFiniteElementEnum<Dimension2, 1>::InitStaticData();
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
    VarFiniteElementEnum<Dimension2, 2>::InitStaticData();
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
    VarFiniteElementEnum<Dimension2, 3>::InitStaticData();
#endif
#endif
    
#ifdef MONTJOIE_WITH_THREE_DIM

#ifdef MONTJOIE_WITH_NODAL_H1
    VarFiniteElementEnum<Dimension3, 1>::InitStaticData();
#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL
    VarFiniteElementEnum<Dimension3, 2>::InitStaticData();
#endif

#ifdef MONTJOIE_WITH_NODAL_HDIV
    VarFiniteElementEnum<Dimension3, 3>::InitStaticData();
#endif
    
#endif

#ifdef MONTJOIE_WITH_HCURL_AXI
    VarFiniteElementEnum<Dimension2, 4>::InitStaticData();
#endif
    
#endif
    
#endif

#ifdef _OPENMP
#ifdef MONTJOIE_FFTW_FFT
    fftw_init_threads();
#endif
#endif

#ifdef SELDON_WITH_SCALAPACK
    global_blacs_handle.Init();
#endif

#ifdef SELDON_WITH_SLEPC
    string help("Slepc interface");  
    SlepcInitialize(&argc, &argv, (char*)0, help.data());
#endif
    
  }
  
  
  int FinalizeMontjoie()
  {

#ifdef _OPENMP
#ifdef MONTJOIE_FFTW_FFT
    fftw_cleanup_threads();
#endif
#endif

#ifdef SELDON_WITH_SLEPC
    SlepcFinalize();
#endif
    
#ifdef SELDON_WITH_MPI
    MPI_Finalize();
#endif

    return 0;
  }
#endif

}

namespace Seldon
{
#ifndef SELDON_WITH_COMPILED_LIBRARY
  double Norm2(const complex<double>& c) { return abs(c); }
  double Norm2(const double& c) { return abs(c); }

#ifdef MONTJOIE_WITH_MULTIPLE
  Real_wp Norm2(const Complex_wp& c) { return abs(c); }
  Real_wp Norm2(const Real_wp& c) { return abs(c); }
#endif

#endif
}

#define MONTJOIE_FILE_MONTJOIE_TYPES_CXX
#endif
