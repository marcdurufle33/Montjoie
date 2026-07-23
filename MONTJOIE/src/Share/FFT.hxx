#ifndef MONTJOIE_FILE_FFT_HXX

#ifdef MONTJOIE_WITH_GSL
#include <gsl/gsl_fft_complex.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#endif

#ifdef MONTJOIE_WITH_FFTW
#include<fftw3.h>
#endif

#ifdef SELDON_WITH_MKL
#include "mkl_dfti.h"
#endif

namespace Montjoie
{
  template<class T>
  class FftInterface;

  //! base class for complex fft
  template<class T>
  class VirtualFftObject
  {
    template<class T0>
    friend class FftInterface;
    
  protected:
    bool allocate_datasol;
    int nb_modes_x, nb_modes_y, nb_modes_z;
    Vector<T> data, xsol, phase_x, phase_y, phase_z;
    
  public:
    VirtualFftObject();
    virtual ~VirtualFftObject();

    int GetNbPoints() const;
    size_t GetMemorySize() const;

    virtual void SetNbThreads(int n);

    virtual void Init(int n);
    virtual void Init(int nx, int ny);
    virtual void Init(int nx, int ny, int nz);

    virtual void ApplyForward(Vector<T>& x) = 0;
    virtual void ApplyInverse(Vector<T>& x) = 0;
    
  };


  //! FFT using naive formulae (O(n^2) complexity)
  template<class T>
  class ManualFftObject : public VirtualFftObject<T>
  {
    void CheckRadixN(int N);
    void cfft(Vector<T>& x, T phi);
    
  public:
    void Init(int n);
    void Init(int nx, int ny);
    void Init(int nx, int ny, int nz);

    void ApplyForward(Vector<T>& x);
    void ApplyInverse(Vector<T>& x);
    
  };

  
#ifdef MONTJOIE_WITH_GSL
  //! interface with FFT proposed in GSL
  template<class T>
  class GslFftObject : public VirtualFftObject<T>
  {
  protected:
    gsl_fft_complex_wavetable * wavetable_x, *wavetable_y, *wavetable_z;
    gsl_fft_complex_workspace * workspace_x, *workspace_y, *workspace_z;

  public:
    GslFftObject();
    ~GslFftObject();
    
    void Init(int n);
    void Init(int nx, int ny);
    void Init(int nx, int ny, int nz);

    void ApplyForward(Vector<T>& x);
    void ApplyInverse(Vector<T>& x);
    
  };
#endif


#ifdef MONTJOIE_WITH_FFTW
  //! interface with FFT proposed in FFTW
  template<class T>
  class FftwFftObject : public VirtualFftObject<T>
  {
  protected:
    fftw_plan plan_forward;
    fftw_plan plan_backward;
    int nb_threads;

  public:
    FftwFftObject();
    ~FftwFftObject();

    void SetNbThreads(int n);

    void Init(int n);
    void Init(int nx, int ny);
    void Init(int nx, int ny, int nz);

    void ApplyForward(Vector<T>& x);
    void ApplyInverse(Vector<T>& x);
    
  };
#endif

#ifdef SELDON_WITH_MKL
  //! interface with FFT proposed in MKL
  template<class T>
  class MklFftObject : public VirtualFftObject<T>
  {
  protected:
    DFTI_DESCRIPTOR_HANDLE hand_mkl;

  public:
    ~MklFftObject();
    
    void Init(int n);
    void Init(int nx, int ny);
    void Init(int nx, int ny, int nz);

    void ApplyForward(Vector<T>& x);
    void ApplyInverse(Vector<T>& x);
    
  };
#endif

  
  //! object implementing Discrete Fourier Transform (interface with Gsl, Mkl or Fftw)
  template<class T>
  class FftInterface
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    
  protected :    
    VirtualFftObject<T>* var_fft;
    
  public :
    // available interfaces
    enum {MANUAL, FFT_MKL, FFTW, FFT_GSL};
    
    FftInterface();    
    ~FftInterface();

    void SelectFftAlgorithm(int type);
    
    int GetNbPoints() const;
    size_t GetMemorySize() const;
    void SetNbThreads(int);
    
    // Initialization (required before effectively computing a fft)
    void Init(int n);
    void Init(int nx, int ny);
    void Init(int nx, int ny, int nz);
    
    // overwrites x with its forward Discrete Fourier Transform
    void ApplyForward(Vector<T>& x);
    
    // overwrites x with its inverse Discrete Fourier Transform
    void ApplyInverse(Vector<T>& x);
    
    // computation of a single value of inverse Discrete Fourier Transform
    void ApplyForwardPoint(int k, const Vector<T>& x, T& u);
    void ApplyForwardPoint(int kx, int ky, const Vector<T>& x, T& u);
    void ApplyForwardPoint(int kx, int ky, int kz, const Vector<T>& x, T& u);

    void ApplyInversePoint(int k, const Vector<T>& x, T& u);
    void ApplyInversePoint(int kx, int ky, const Vector<T>& x, T& u);
    void ApplyInversePoint(int kx, int ky, int kz, const Vector<T>& x, T& u);

    T GetCoefficient(int ix, int nx);
    T GetCoefficient(int ix, int iy, int nx, int ny);
    T GetCoefficient(int ix, int iy, int iz, int nx, int ny, int nz);
    
    void GetCosSinAlpha(int n, Treal& cos_, Treal& sin_) const;
    
  };


  class FftRealInterface;

  
  //! base class for real fft
  template<class T>
  class VirtualFftRealObject
  {
    friend class FftRealInterface;

  protected:
    Vector<T> data, phase, phase_y;
    Vector<complex<T> > xsol;
    int nb_modes, nb_modes_y;

  public:
    VirtualFftRealObject();
    virtual ~VirtualFftRealObject();

    size_t GetMemorySize() const;
    virtual void SetNbThreads(int n);
    
    virtual void Init(int n);
    virtual void Init(int nx, int ny);
    
    virtual void ApplyForward(const Vector<T>&, Vector<complex<T> >&) = 0;
    virtual void ApplyInverse(const Vector<complex<T> >&, Vector<T>&) = 0;
    
  };


#ifdef MONTJOIE_WITH_GSL
  //! interface with FFT proposed in GSL
  template<class T>
  class GslFftRealObject : public VirtualFftRealObject<T>
  {
  protected:
    gsl_fft_real_wavetable * wavetable_real;
    gsl_fft_halfcomplex_wavetable * wavetable_cplx;
    gsl_fft_real_workspace * workspace_real;
    
  public:
    GslFftRealObject();
    ~GslFftRealObject();

    void Init(int n);
    void ApplyForward(const Vector<T>&, Vector<complex<T> >&);
    void ApplyInverse(const Vector<complex<T> >&, Vector<T>&);
    
  };
#endif


#ifdef MONTJOIE_WITH_FFTW
  //! interface with FFT proposed in FFTW
  template<class T>
  class FftwFftRealObject : public VirtualFftRealObject<T>
  {
  protected:
    fftw_plan plan_forward;
    fftw_plan plan_backward;
    int nb_threads;

  public:
    FftwFftRealObject();
    ~FftwFftRealObject();

    void SetNbThreads(int n);

    void Init(int n);
    void Init(int nx, int ny);
    
    void ApplyForward(const Vector<T>&, Vector<complex<T> >&);
    void ApplyInverse(const Vector<complex<T> >&, Vector<T>&);
    
  };
#endif

#ifdef SELDON_WITH_MKL
  //! interface with FFT proposed in MKL
  template<class T>
  class MklFftRealObject : public VirtualFftRealObject<T>
  {
  protected:
    DFTI_DESCRIPTOR_HANDLE hand_mkl;

  public:
    MklFftRealObject();
    ~MklFftRealObject();
    
    void Init(int n);
    void ApplyForward(const Vector<T>&, Vector<complex<T> >&);
    void ApplyInverse(const Vector<complex<T> >&, Vector<T>&);
    
  };
#endif


  //! object implementing Discrete Fourier Transform (interface with Gsl or Fftw)
  class FftRealInterface
  {
  protected :
    VirtualFftRealObject<Real_wp>* var_fft;
    
  public :
    // available interfaces
    enum {MANUAL, FFT_MKL, FFTW, FFT_GSL};

    FftRealInterface();    
    ~FftRealInterface();
    
    void SelectFftAlgorithm(int type);
    
    int GetNbPoints() const;
    size_t GetMemorySize() const;
    void SetNbThreads(int);
    
    // Initialization (required before effectively computing a fft)
    void Init(int n);
    void Init(int nx, int ny);
    
    // computes y = fft(x)
    void ApplyForward(const Vector<Real_wp>& x,
                      Vector<Complex_wp>& y);
    
    // computes y = ifft(x)
    void ApplyInverse(const Vector<Complex_wp>& x,
                      Vector<Real_wp>& y);

  };
  
}

#define MONTJOIE_FILE_FFT_HXX
#endif
