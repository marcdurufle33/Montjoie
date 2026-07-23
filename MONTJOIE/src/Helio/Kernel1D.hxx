#ifndef MONTJOIE_FILE_KERNEL1D_HXX


namespace Montjoie
{
  class Kernel1D
  {
    
  private :
    int lMaxK;
    Vector<Real_wp> theta1, theta2, phi1, phi2;
    Real_wp omega;
    Real_wp Rsun; // solar radius
    Matrix<Complex_wp> G; // Green's function as a function of l and depth
    Matrix<Complex_wp> dG; // Radial derivative of the Green's function as a function of l and depth
    Matrix<Complex_wp> C; // cross covariance as a function of l and depth
    Matrix<Complex_wp> dC; // Radial derivative of the cross covariance as a function of l ande depth
    Vector<Real_wp> rho; // density
    Vector<Real_wp> c; // sound speed
    Vector<Real_wp> r; // radius where the kernel is computed
    Vector<Real_wp> frequencies;
    int indFreq; // frequency index for XS measurements
    Vector<Complex_wp> W; // weighting function from travel-time/amplitude
    Vector<Complex_wp> gamma; // damping as a function of frequency
    string typeOfFilter; // filter on the observations
    Real_wp centerFilter, widthFilter; // parameters of the filter
    string typeOfKernel;
    string typeOfObservable; // XS or tau
    string dirOut;
    int nbKernels;
    // Quantities averaged over frequencies for travel time computation
    Vector<Matrix<Complex_wp> > fllpr; 
    Vector<Matrix<Complex_wp> > gllpr;
    
  public :

    string kernelDimension;
    int lMaxG; // lMax for the Green's function
    int nr; // number of depths
    Vector<Matrix<Complex_wp> > kernels;

    Kernel1D();

    void getGreen(string input_file);

    void ReadInputFile(const string filename, const int indFreqMin);

    //    Vector<Matrix<Complex_wp> > mSumforKernel(const Vector<int> lKs, const int mK);

    Vector<Matrix<Complex_wp> > computeKernel(const Vector<int> lKs, const int mK, int rank_proc);

    void computeKernel3D_depth(const Vector<int> indz);
 
    void computeKernel3D();

    void compute(const int mK);

  };

}


#define MONTJOIE_FILE_KERNEL1D_HXX
#endif
