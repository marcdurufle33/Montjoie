#ifndef MONTJOIE_FILE_FFT_CXX

namespace Montjoie
{

  /********************
   * VirtualFftObject *
   ********************/
  
  
  //! default constructor
  template<class T>
  VirtualFftObject<T>::VirtualFftObject()
  {
    allocate_datasol = false;
    nb_modes_x = 0;
    nb_modes_y = 0;
    nb_modes_z = 0;
  }
  
  
  //! destructor
  template<class T>
  VirtualFftObject<T>::~VirtualFftObject()
  {
  }


  //! returns the number of points for the 1-D FFT
  template<class T>
  int VirtualFftObject<T>::GetNbPoints() const
  {
    return nb_modes_x;
  }


  //! returns the memory used by the object in bytes
  template<class T>
  size_t VirtualFftObject<T>::GetMemorySize() const
  {
    size_t taille = data.GetMemorySize() + xsol.GetMemorySize() + phase_x.GetMemorySize()
      + phase_y.GetMemorySize() + phase_z.GetMemorySize();
    
    taille += sizeof(*this);
    return taille;
  }


  //! modifies the number of threads
  template<class T>
  void VirtualFftObject<T>::SetNbThreads(int n)
  {    
  }
  

  //! Inits a 1-D FFT with n points
  template<class T>
  void VirtualFftObject<T>::Init(int n)
  {
    nb_modes_y = 0;
    nb_modes_z = 0;
    
    if (n <= 0)
      return;

    nb_modes_x = n;

    if (allocate_datasol)
      {
        data.Reallocate(n);
        xsol.Reallocate(n);
      }
    
    phase_x.Reallocate(n);
    
    SetComplexOne(phase_x(0));
    typename ClassComplexType<T>::Treal alpha = 2.0*pi_wp/n;
    T coef = exp(Iwp*alpha);
    for (int i = 0; i < n-1; i++)
      phase_x(i+1) = coef*phase_x(i);
  }


  //! Inits a 2-D FFT with (nx, ny) points
  template<class T>
  void VirtualFftObject<T>::Init(int nx, int ny)
  {
    nb_modes_z = 0;
    
    if ((nx <= 0) || (ny <= 0))
      return;

    nb_modes_x = nx;
    nb_modes_y = ny;
    if (allocate_datasol)
      {
        data.Reallocate(nx*ny);
        xsol.Reallocate(nx*ny);
      }
    
    phase_x.Reallocate(nx);
    phase_y.Reallocate(ny);
    
    SetComplexOne(phase_x(0));
    typename ClassComplexType<T>::Treal alpha = 2.0*pi_wp/nx;
    T coef = exp(Iwp*alpha);
    for (int i = 0; i < nx-1; i++)
      phase_x(i+1) = coef*phase_x(i);

    SetComplexOne(phase_y(0));
    alpha = 2.0*pi_wp/ny;
    coef = exp(Iwp*alpha);
    for (int i = 0; i < ny-1; i++)
      phase_y(i+1) = coef*phase_y(i);
    
  }


  //! Inits a 3-D FFT with (nx, ny, nz) points
  template<class T>
  void VirtualFftObject<T>::Init(int nx, int ny, int nz)
  {
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
      return;

    nb_modes_x = nx; 
    nb_modes_y = ny;
    nb_modes_z = nz;
    if (allocate_datasol)
      {
        data.Reallocate(nx*ny*nz);
        xsol.Reallocate(nx*ny*nz);
      }
    
    phase_x.Reallocate(nx);
    SetComplexOne(phase_x(0));
    typename ClassComplexType<T>::Treal alpha = 2.0*pi_wp/nx;
    T coef = exp(Iwp*alpha);
    for (int i = 0; i < nx-1; i++)
      phase_x(i+1) = coef*phase_x(i);

    phase_y.Reallocate(ny);
    SetComplexOne(phase_y(0));
    alpha = 2.0*pi_wp/ny;
    coef = exp(Iwp*alpha);
    for (int i = 0; i < ny-1; i++)
      phase_y(i+1) = coef*phase_y(i);

    phase_z.Reallocate(nz);
    SetComplexOne(phase_z(0));
    alpha = 2.0*pi_wp/nz;
    coef = exp(Iwp*alpha);
    for (int i = 0; i < nz-1; i++)
      phase_z(i+1) = coef*phase_z(i);
    
  }


  /*******************
   * ManualFftObject *
   *******************/
  

  //! Checks that N is a power of 2
  template<class T>
  void ManualFftObject<T>::CheckRadixN(int N)
  {
    int m = (int)log2(N);
    if (pow(2, m) != N)
      {
        cout << N << " is not a power of 2 " << endl;
        cout << "Link with an external Fft Library (Fftw, Mkl, etc)" << endl;
        abort();
      }
  }

  
  //! initializes 1-D FFT
  template<class T>
  void ManualFftObject<T>::Init(int Nx)
  {
    VirtualFftObject<T>::Init(Nx);
    CheckRadixN(Nx);
  }


  //! initializes 2-D FFT
  template<class T>
  void ManualFftObject<T>::Init(int Nx, int Ny)
  {
    VirtualFftObject<T>::Init(Nx, Ny);
    CheckRadixN(Nx); CheckRadixN(Ny);
  }


  //! initializes 3-D FFT
  template<class T>
  void ManualFftObject<T>::Init(int Nx, int Ny, int Nz)
  {
    VirtualFftObject<T>::Init(Nx, Ny, Nz);
    CheckRadixN(Nx); CheckRadixN(Ny); CheckRadixN(Nz);
  }
  

  //! Auxiliary function (1-D Fft)
  template<class T>
  void ManualFftObject<T>::cfft(Vector<T>& x, T phiT)
  {
    int N = x.GetM(), k = N, n;
    T tmp, t;
    while (k > 1)
      {
        n = k;
        k >>= 1;
        SetComplexOne(tmp);
        if (n == N)
          {
            for (int l = 0; l < k; l++)
              {
                int b = l + k;
                t = x(l) - x(b);
                x(l) += x(b);
                x(b) = t*tmp;
                tmp *= phiT;
              }                    
          }
        else
          {
            for (int l = 0; l < k; l++)
              {
                for (int a = l; a < N; a += n)
                  {
                    int b = a + k;
                    t = x(a) - x(b);
                    x(a) += x(b);
                    x(b) = t * tmp;
                  }
                
                tmp *= phiT;
              }
          }

        phiT = phiT * phiT;
      }

    // Decimate
    int m = (int)log2(N); unsigned int Nu(N);
    for (unsigned int a = 0; a < Nu; a++)
      {
        unsigned int b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a)
          {
            t = x(a);
            x(a) = x(b);
            x(b) = t;
          }
      }
  }
  
  
  //! Overwrites x by its Forward Discrete Fourier Transform
  template<class T>
  void ManualFftObject<T>::ApplyForward(Vector<T>& x)
  {
    if (this->nb_modes_z <= 1)
      {
        if (this->nb_modes_y <= 1)
          {
            // 1-D FFT
            if (this->nb_modes_x > 1)
              cfft(x, conjugate(this->phase_x(1)));
          }
        else
          {
            // 2-D FFT
            Vector<T> ux(this->nb_modes_x), uy(this->nb_modes_y);
            for (int i = 0; i < this->nb_modes_x; i++)
              {
                int offset = i*this->nb_modes_y;
                for (int j = 0; j < this->nb_modes_y; j++)
                  uy(j) = x(offset + j);
                
                if (this->nb_modes_y > 1)
                  cfft(uy, conjugate(this->phase_y(1)));
                
                for (int j = 0; j < this->nb_modes_y; j++)
                  x(offset + j) = uy(j);        
              }
            
            for (int j = 0; j < this->nb_modes_y; j++)
              {
                for (int i = 0; i < this->nb_modes_x; i++)
                  ux(i) = x(i*this->nb_modes_y + j);

                if (this->nb_modes_x > 1)
                  cfft(ux, conjugate(this->phase_x(1)));
                
                for (int i = 0; i < this->nb_modes_x; i++)
                  x(i*this->nb_modes_y + j) = ux(i);        
              }            
          }
      }
    else
      {
        // 3-D FFT
        Vector<T> ux(this->nb_modes_x), uy(this->nb_modes_y), uz(this->nb_modes_z);
        for (int ix = 0; ix < this->nb_modes_x; ix++)
          for (int iy = 0; iy < this->nb_modes_y; iy++)
            {
              int offset = this->nb_modes_z*(ix*this->nb_modes_y + iy);
              for (int j = 0; j < this->nb_modes_z; j++)
                uz(j) = x(offset + j);

              if (this->nb_modes_z > 1)
                cfft(uz, conjugate(this->phase_z(1)));

              for (int j = 0; j < this->nb_modes_z; j++)
                x(offset + j) = uz(j);
            }
        
        for (int ix = 0; ix < this->nb_modes_x; ix++)
          for (int iz = 0; iz < this->nb_modes_z; iz++)
            {
              int offset = ix*this->nb_modes_y;
              for (int j = 0; j < this->nb_modes_y; j++)
                uy(j) = x(this->nb_modes_z*(offset + j) + iz);

              if (this->nb_modes_y > 1)
                cfft(uy, conjugate(this->phase_y(1)));

              for (int j = 0; j < this->nb_modes_y; j++)
                x(this->nb_modes_z*(offset + j) + iz) = uy(j);
            }
        
        
        for (int iy = 0; iy < this->nb_modes_y; iy++)
          for (int iz = 0; iz < this->nb_modes_z; iz++)
            {
              for (int j = 0; j < this->nb_modes_x; j++)
                ux(j) = x(this->nb_modes_z*(j*this->nb_modes_y + iy) + iz);

              if (this->nb_modes_x > 1)
                cfft(ux, conjugate(this->phase_x(1)));
              
              for (int j = 0; j < this->nb_modes_x; j++)
                x(this->nb_modes_z*(j*this->nb_modes_y + iy) + iz) = ux(j);
            }
      }
  }


  //! Overwrites x by its Backward Discrete Fourier Transform
  template<class T>
  void ManualFftObject<T>::ApplyInverse(Vector<T>& x)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    if (this->nb_modes_z <= 1)
      {
        if (this->nb_modes_y <= 1)
          {
            // 1-D FFt
            if (this->nb_modes_x > 1)
              {
                cfft(x, this->phase_x(1));
                Mlt(Treal(1)/this->nb_modes_x, x);
              }
          }
        else
          {
            // 2-D FFT
            Vector<T> ux(this->nb_modes_x), uy(this->nb_modes_y);
            for (int i = 0; i < this->nb_modes_x; i++)
              {
                int offset = i*this->nb_modes_y;
                for (int j = 0; j < this->nb_modes_y; j++)
                  uy(j) = x(offset + j);
                
                if (this->nb_modes_y > 1)
                  cfft(uy, this->phase_y(1));
                
                for (int j = 0; j < this->nb_modes_y; j++)
                  x(offset + j) = uy(j);        
              }
            
            for (int j = 0; j < this->nb_modes_y; j++)
              {
                for (int i = 0; i < this->nb_modes_x; i++)
                  ux(i) = x(i*this->nb_modes_y + j);
                
                if (this->nb_modes_x > 1)
                  cfft(ux, this->phase_x(1));
                
                for (int i = 0; i < this->nb_modes_x; i++)
                  x(i*this->nb_modes_y + j) = ux(i);        
              }
            
            Mlt(Treal(1)/(this->nb_modes_x*this->nb_modes_y), x);
          }        
      }
    else
      {
        // 3-D FFT
        Vector<T> ux(this->nb_modes_x), uy(this->nb_modes_y), uz(this->nb_modes_z);
        for (int ix = 0; ix < this->nb_modes_x; ix++)
          for (int iy = 0; iy < this->nb_modes_y; iy++)
            {
              int offset = this->nb_modes_z*(ix*this->nb_modes_y + iy);
              for (int j = 0; j < this->nb_modes_z; j++)
                uz(j) = x(offset + j);

              if (this->nb_modes_z > 1)
                cfft(uz, this->phase_z(1));

              for (int j = 0; j < this->nb_modes_z; j++)
                x(offset + j) = uz(j);
            }
        
        for (int ix = 0; ix < this->nb_modes_x; ix++)
          for (int iz = 0; iz < this->nb_modes_z; iz++)
            {
              int offset = ix*this->nb_modes_y;
              for (int j = 0; j < this->nb_modes_y; j++)
                uy(j) = x(this->nb_modes_z*(offset + j) + iz);

              if (this->nb_modes_y > 1)
                cfft(uy, this->phase_y(1));

              for (int j = 0; j < this->nb_modes_y; j++)
                x(this->nb_modes_z*(offset + j) + iz) = uy(j);
            }
        
        
        for (int iy = 0; iy < this->nb_modes_y; iy++)
          for (int iz = 0; iz < this->nb_modes_z; iz++)
            {
              for (int j = 0; j < this->nb_modes_x; j++)
                ux(j) = x(this->nb_modes_z*(j*this->nb_modes_y + iy) + iz);

              if (this->nb_modes_x > 1)
                cfft(ux, this->phase_x(1));

              for (int j = 0; j < this->nb_modes_x; j++)
                x(this->nb_modes_z*(j*this->nb_modes_y + iy) + iz) = ux(j);
            }
        
        Mlt(Treal(1)/(this->nb_modes_x*this->nb_modes_y*this->nb_modes_z), x);
      }    
  }
  

  /****************
   * GslFftObject *
   ****************/

  
#ifdef MONTJOIE_WITH_GSL
  //! default constructor
  template<class T>
  GslFftObject<T>::GslFftObject()
  {
    wavetable_x = NULL;
    workspace_x = NULL;

    wavetable_y = NULL;
    workspace_y = NULL;

    wavetable_z = NULL;
    workspace_z = NULL;
  }


  //! destructor
  template<class T>
  GslFftObject<T>::~GslFftObject()
  {
    if (wavetable_x != NULL)
      {
	gsl_fft_complex_wavetable_free (wavetable_x);
	gsl_fft_complex_workspace_free (workspace_x);
	wavetable_x = NULL;
	workspace_x = NULL;
      }
    
    if (wavetable_y != NULL)
      {
	gsl_fft_complex_wavetable_free (wavetable_y);
	gsl_fft_complex_workspace_free (workspace_y);
	wavetable_y = NULL;
	workspace_y = NULL;
      }
    
    if (wavetable_z != NULL)
      {
	gsl_fft_complex_wavetable_free (wavetable_z);
	gsl_fft_complex_workspace_free (workspace_z);
	wavetable_z = NULL;
	workspace_z = NULL;
      }
  }

  
  //! Inits a 1-D FFT with n points
  template<class T>
  void GslFftObject<T>::Init(int n)
  {
    VirtualFftObject<T>::Init(n);
    if (n <= 0)
      return;
    
    wavetable_x = gsl_fft_complex_wavetable_alloc(n);
    workspace_x = gsl_fft_complex_workspace_alloc(n);
  }


  //! Inits a 2-D FFT with (nx, ny) points
  template<class T>
  void GslFftObject<T>::Init(int nx, int ny)
  {
    VirtualFftObject<T>::Init(nx, ny);
    if ((nx <= 0) || (ny <= 0))
      return;

    wavetable_x = gsl_fft_complex_wavetable_alloc(nx);
    workspace_x = gsl_fft_complex_workspace_alloc(nx);

    wavetable_y = gsl_fft_complex_wavetable_alloc(ny);
    workspace_y = gsl_fft_complex_workspace_alloc(ny);
  }


  //! Inits a 3-D FFT with (nx, ny, nz) points
  template<class T>
  void GslFftObject<T>::Init(int nx, int ny, int nz)
  {
    VirtualFftObject<T>::Init(nx, ny, nz);
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
      return;
    
    wavetable_x = gsl_fft_complex_wavetable_alloc(nx);
    workspace_x = gsl_fft_complex_workspace_alloc(nx);

    wavetable_y = gsl_fft_complex_wavetable_alloc(ny);
    workspace_y = gsl_fft_complex_workspace_alloc(ny);

    wavetable_z = gsl_fft_complex_wavetable_alloc(nz);
    workspace_z = gsl_fft_complex_workspace_alloc(nz);
  }


  //! Overwrites x by its Forward Discrete Fourier Transform
  template<class T>
  void GslFftObject<T>::ApplyForward(Vector<T>& x)
  {
    if (this->nb_modes_z == 0)
      {
        if (this->nb_modes_y == 0)
          {
            // 1-D FFT
            gsl_fft_complex_forward(reinterpret_cast<double*> (x.GetData()), 1,
                                    this->nb_modes_x, wavetable_x, workspace_x);
          }
        else
          {
            // 2-D FFT
            for (int i = 0; i < this->nb_modes_x; i++)
              gsl_fft_complex_forward(reinterpret_cast<double*> (&x(i*this->nb_modes_y)), 1,
                                      this->nb_modes_y, wavetable_y, workspace_y);
            
            for (int i = 0; i < this->nb_modes_y; i++)
              gsl_fft_complex_forward(reinterpret_cast<double*> (&x(i)), this->nb_modes_y,
                                      this->nb_modes_x, wavetable_x, workspace_x);
          }
      }
    else
      {
        // 3-D FFT
        for (int i = 0; i < this->nb_modes_x; i++)
          for (int j = 0; j < this->nb_modes_y; j++)
            gsl_fft_complex_forward(reinterpret_cast<double*>
                                    (&x(this->nb_modes_z*(i*this->nb_modes_y + j))), 1,
                                    this->nb_modes_z, wavetable_z, workspace_z);

        for (int i = 0; i < this->nb_modes_x; i++)
          for (int j = 0; j < this->nb_modes_z; j++)
            gsl_fft_complex_forward(reinterpret_cast<double*>
                                    (&x(this->nb_modes_z*this->nb_modes_y*i + j)), this->nb_modes_z,
                                    this->nb_modes_y, wavetable_y, workspace_y);
            
        for (int i = 0; i < this->nb_modes_y; i++)
          for (int j = 0; j < this->nb_modes_z; j++)
            gsl_fft_complex_forward(reinterpret_cast<double*>
                                    (&x(i*this->nb_modes_z + j)), this->nb_modes_y*this->nb_modes_z,
                                    this->nb_modes_x, wavetable_x, workspace_x);
      }
  }


  //! Overwrites x by its Backward Discrete Fourier Transform
  template<class T>
  void GslFftObject<T>::ApplyInverse(Vector<T>& x)
  {
    if (this->nb_modes_z == 0)
      {
        if (this->nb_modes_y == 0)
          {
            // 1-D FFT
            gsl_fft_complex_inverse(reinterpret_cast<double*> (x.GetData()), 1,
                                    this->nb_modes_x, wavetable_x, workspace_x);
          }
        else
          {
            // 2-D FFT
            for (int i = 0; i < this->nb_modes_x; i++)
              gsl_fft_complex_inverse(reinterpret_cast<double*> (&x(i*this->nb_modes_y)), 1,
                                      this->nb_modes_y, wavetable_y, workspace_y);
            
            for (int i = 0; i < this->nb_modes_y; i++)
              gsl_fft_complex_inverse(reinterpret_cast<double*> (&x(i)), this->nb_modes_y,
                                      this->nb_modes_x, wavetable_x, workspace_x);
          }
      }
    else
      {
        // 3-D FFT
        for (int i = 0; i < this->nb_modes_x; i++)
          for (int j = 0; j < this->nb_modes_y; j++)
            gsl_fft_complex_inverse(reinterpret_cast<double*>
                                    (&x(this->nb_modes_z*(i*this->nb_modes_y + j))), 1,
                                    this->nb_modes_z, wavetable_z, workspace_z);
        
        for (int i = 0; i < this->nb_modes_x; i++)
          for (int j = 0; j < this->nb_modes_z; j++)
            gsl_fft_complex_inverse(reinterpret_cast<double*>
                                    (&x(this->nb_modes_z*this->nb_modes_y*i + j)), this->nb_modes_z,
                                    this->nb_modes_y, wavetable_y, workspace_y);
        
        for (int i = 0; i < this->nb_modes_y; i++)
          for (int j = 0; j < this->nb_modes_z; j++)
            gsl_fft_complex_inverse(reinterpret_cast<double*>
                                    (&x(i*this->nb_modes_z + j)), this->nb_modes_y*this->nb_modes_z,
                                    this->nb_modes_x, wavetable_x, workspace_x);
      }
  }
#endif


  /****************
   * FftwFftObject *
   ****************/

  
#ifdef MONTJOIE_WITH_FFTW
  //! default constructor
  template<class T>
  FftwFftObject<T>::FftwFftObject()
  {
#ifdef _OPENMP
    nb_threads = omp_get_max_threads();
#else
    nb_threads = 1;
#endif
    this->allocate_datasol = true;
  }


  //! destructor
  template<class T>
  FftwFftObject<T>::~FftwFftObject()
  {
    if (this->nb_modes_x > 0)
      {
	fftw_destroy_plan(plan_forward);
        fftw_destroy_plan(plan_backward);

	this->nb_modes_x = 0;
      }
  }

  
  //! modifies the number of threads
  template<class T>
  void FftwFftObject<T>::SetNbThreads(int n)
  {
    nb_threads = n;
  }


  //! Inits a 1-D FFT with n points
  template<class T>
  void FftwFftObject<T>::Init(int n)
  {
    VirtualFftObject<T>::Init(n);
    if (n <= 0)
      return;

#ifdef _OPENMP
    fftw_plan_with_nthreads(nb_threads);
#endif
    
    fftw_complex* in = reinterpret_cast<fftw_complex*>(this->data.GetData());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(this->xsol.GetData());
    plan_forward = fftw_plan_dft_1d(n, in, out,
                                    FFTW_FORWARD, FFTW_ESTIMATE);
    
    plan_backward = fftw_plan_dft_1d(n, in, out,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);    
  }


  //! Inits a 2-D FFT with (nx, ny) points
  template<class T>
  void FftwFftObject<T>::Init(int nx, int ny)
  {
    VirtualFftObject<T>::Init(nx, ny);
    if ((nx <= 0) || (ny <= 0))
      return;

#ifdef _OPENMP
    fftw_plan_with_nthreads(nb_threads);
#endif
    
    fftw_complex* in = reinterpret_cast<fftw_complex*>(this->data.GetData());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(this->xsol.GetData());
    plan_forward = fftw_plan_dft_2d(nx, ny, in, out,
                                    FFTW_FORWARD, FFTW_ESTIMATE);

    plan_backward = fftw_plan_dft_2d(nx, ny, in, out,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
  }


  //! Inits a 3-D FFT with (nx, ny, nz) points
  template<class T>
  void FftwFftObject<T>::Init(int nx, int ny, int nz)
  {
    VirtualFftObject<T>::Init(nx, ny, nz);
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
      return;

#ifdef _OPENMP
    fftw_plan_with_nthreads(nb_threads);
#endif

    fftw_complex* in = reinterpret_cast<fftw_complex*>(this->data.GetData());
    fftw_complex* out = reinterpret_cast<fftw_complex*>(this->xsol.GetData());
    plan_forward = fftw_plan_dft_3d(nx, ny, nz, in, out,
                                    FFTW_FORWARD, FFTW_ESTIMATE);

    plan_backward = fftw_plan_dft_3d(nx, ny, nz, in, out,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    
  }


  //! Overwrites x by its Forward Discrete Fourier Transform
  template<class T>
  void FftwFftObject<T>::ApplyForward(Vector<T>& x)
  {
    for (int i = 0; i < this->data.GetM(); i++)
      this->data(i) = x(i);
            
    fftw_execute(plan_forward);
    
    for (int i = 0; i < this->xsol.GetM(); i++)
      x(i) = this->xsol(i);
  }


  //! Overwrites x by its Backward Discrete Fourier Transform
  template<class T>
  void FftwFftObject<T>::ApplyInverse(Vector<T>& x)
  {
    // 2-D FFT
    for (int i = 0; i < this->data.GetM(); i++)
      this->data(i) = x(i);
            
    fftw_execute(plan_backward);
    
    for (int i = 0; i < this->xsol.GetM(); i++)
      x(i) = this->xsol(i);
            
    if (this->nb_modes_z == 0)
      {
	if (this->nb_modes_y == 0)
	  Mlt(1.0/this->nb_modes_x, x);
	else
	  Mlt(1.0/(this->nb_modes_x*this->nb_modes_y), x);
      }
    else
        Mlt(1.0/(this->nb_modes_x*this->nb_modes_y*this->nb_modes_z), x);
  }
#endif


  /****************
   * MklFftObject *
   ****************/

  
#ifdef SELDON_WITH_MKL
  //! destructor
  template<class T>
  MklFftObject<T>::~MklFftObject()
  {
    if (this->nb_modes_x > 0)
      DftiFreeDescriptor(&hand_mkl);
  }


  //! Inits a 1-D FFT with n points
  template<class T>
  void MklFftObject<T>::Init(int n)
  {
    VirtualFftObject<T>::Init(n);
    if (n <= 0)
      return;

    DftiCreateDescriptor(&hand_mkl, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
    DftiSetValue(hand_mkl, DFTI_BACKWARD_SCALE, 1.0/n);
    //DftiSetValue(hand_mkl, DFTI_THREAD_LIMIT, nb_threads);
    DftiCommitDescriptor(hand_mkl);
  }


  //! Inits a 2-D FFT with (nx, ny) points
  template<class T>
  void MklFftObject<T>::Init(int nx, int ny)
  {
    VirtualFftObject<T>::Init(nx, ny);
    if ((nx <= 0) || (ny <= 0))
      return;
    
    MKL_LONG n[2] = {nx, ny};
    MKL_LONG strides[3] = {0, ny, 1};
    DftiCreateDescriptor(&hand_mkl, DFTI_DOUBLE, DFTI_COMPLEX, 2, n);
    DftiSetValue(hand_mkl, DFTI_BACKWARD_SCALE, 1.0/(nx*ny));
    DftiSetValue(hand_mkl, DFTI_INPUT_STRIDES, strides);
    //DftiSetValue(hand_mkl, DFTI_THREAD_LIMIT, nb_threads);
    DftiCommitDescriptor(hand_mkl);
  }


  //! Inits a 3-D FFT with (nx, ny, nz) points
  template<class T>
  void MklFftObject<T>::Init(int nx, int ny, int nz)
  {
    VirtualFftObject<T>::Init(nx, ny, nz);
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
      return;

    MKL_LONG n[3] = {nx, ny, nz};
    MKL_LONG strides[4] = {0, nz*ny, nz, 1};
    DftiCreateDescriptor(&hand_mkl, DFTI_DOUBLE, DFTI_COMPLEX, 3, n);
    DftiSetValue(hand_mkl, DFTI_BACKWARD_SCALE, 1.0/(nx*ny*nz));
    DftiSetValue(hand_mkl, DFTI_INPUT_STRIDES, strides);
    //DftiSetValue(hand_mkl, DFTI_THREAD_LIMIT, nb_threads);
    DftiCommitDescriptor(hand_mkl);

  }


  //! Overwrites x by its Forward Discrete Fourier Transform
  template<class T>
  void MklFftObject<T>::ApplyForward(Vector<T>& x)
  {
    DftiComputeForward(hand_mkl, x.GetDataVoid());
  }


  //! Overwrites x by its Backward Discrete Fourier Transform
  template<class T>
  void MklFftObject<T>::ApplyInverse(Vector<T>& x)
  {
    DftiComputeBackward(hand_mkl, x.GetDataVoid());
  }
#endif


  /****************
   * FftInterface *
   ****************/

  
  //! default constructor
  template<class T>
  FftInterface<T>::FftInterface()
  {
    var_fft = NULL;
    int type_fft = MANUAL;
    
    // using mkl if possible, then fftw, then gsl, then manual
#ifdef MONTJOIE_WITH_GSL
    type_fft = FFT_GSL;
#endif

#ifdef MONTJOIE_WITH_FFTW
    type_fft = FFTW;
#endif
    
#ifdef SELDON_WITH_MKL
    type_fft = FFT_MKL;
#endif

    SelectFftAlgorithm(type_fft);
  }

  
  //! destructor
  template<class T>
  FftInterface<T>::~FftInterface()
  {
    delete var_fft;
  }
  

  //! selects a different interface
  template<class T>
  void FftInterface<T>::SelectFftAlgorithm(int type)
  {
    if (var_fft != NULL)
      delete var_fft;
    
    switch (type)
      {
      case MANUAL:
	var_fft = new ManualFftObject<T>;
	break;
      case FFT_MKL:
	{
#ifdef SELDON_WITH_MKL
	  var_fft = new MklFftObject<T>;
#else
	  cout << "Recompile with SELDON_WITH_MKL" << endl;
	  abort();
#endif
	}
	break;
      case FFTW:
	{
#ifdef MONTJOIE_WITH_FFTW
	  var_fft = new FftwFftObject<T>;
#else
	  cout << "Recompile with MONTJOIE_WITH_FFTW" << endl;
	  abort();
#endif
	}
	break;
      case FFT_GSL:
	{
#ifdef MONTJOIE_WITH_GSL
	  var_fft = new GslFftObject<T>;
#else
	  cout << "Recompile with MONTJOIE_WITH_GSL" << endl;
	  abort();
#endif
	}
	break;
      default:
	cout << "Unknown type" << endl;
	abort();
      }
  }
  
  
  //! returns the number of points in the fft (1-D)
  template<class T>
  int FftInterface<T>::GetNbPoints() const
  {
    return var_fft->GetNbPoints();
  }


  //! returns the size used by the object in bytes
  template<class T>
  size_t FftInterface<T>::GetMemorySize() const
  {
    return var_fft->GetMemorySize();
  }
  
  
  //! modifies the number of threads
  template<class T>
  void FftInterface<T>::SetNbThreads(int n)
  {
    var_fft->SetNbThreads(n);
  }

  
  //! initialization for computation of 1-D Discrete Fourier transform
  template<class T>
  void FftInterface<T>::Init(int n)
  {
    var_fft->Init(n);
  }


  //! Initialization for computation of 2-D Discrete Fourier Transform
  template<class T>
  void FftInterface<T>::Init(int nx, int ny)
  {
    var_fft->Init(nx, ny);
  }
  
  
  //! Initialization for computation of 3-D Discrete Fourier Transform
  template<class T>
  void FftInterface<T>::Init(int nx, int ny, int nz)
  {
    var_fft->Init(nx, ny, nz);
  }
  

  //! Overwrites x by its Forward Discrete Fourier Transform
  template<class T>
  void FftInterface<T>::ApplyForward(Vector<T>& x)
  {
    var_fft->ApplyForward(x);
  }
  

  //! Overwrites x by its Backward Discrete Fourier Transform
  template<class T>
  void FftInterface<T>::ApplyInverse(Vector<T>& x)
  {
    var_fft->ApplyInverse(x);
  }
  

  template<class T>
  void FftInterface<T>::
  ApplyForwardPoint(int k, const Vector<T>& x, T& u)
  {
    int p = 0;
    SetComplexZero(u);
    for (int n = 0; n < var_fft->nb_modes_x; n++)
      {        
        u += conj(var_fft->phase_x(p))*x(n);
        p = (p+k)%var_fft->nb_modes_x;
      }
  }

  
  template<class T>
  void FftInterface<T>::
  ApplyForwardPoint(int ix, int iy, const Vector<T>& x, T& u)
  {
    SetComplexZero(u);
    int px = 0;
    for (int nx = 0; nx < var_fft->nb_modes_x; nx++)
      {   
        int py = 0;
        for (int ny = 0; ny < var_fft->nb_modes_y; ny++)
          {
            u += conj(var_fft->phase_x(px)*var_fft->phase_y(py))*x(nx*var_fft->nb_modes_y + ny);
            py = (py + iy)%var_fft->nb_modes_y;
          }
        
        px = (px + ix)%var_fft->nb_modes_x;
      }
  }


  template<class T>
  void FftInterface<T>::
  ApplyForwardPoint(int ix, int iy, int iz, const Vector<T>& x, T& u)
  {
    SetComplexZero(u);
    int px = 0;
    for (int nx = 0; nx < var_fft->nb_modes_x; nx++)
      {   
        int py = 0;
        for (int ny = 0; ny < var_fft->nb_modes_y; ny++)
          {
            int pz = 0;
            for (int nz = 0; nz < var_fft->nb_modes_z; nz++)
              {                
                u += conj(var_fft->phase_x(px)*var_fft->phase_y(py)*var_fft->phase_z(pz))
		  *x(var_fft->nb_modes_z*(nx*var_fft->nb_modes_y + ny)+nz);
		
                pz = (pz + iz)%var_fft->nb_modes_z;
              }
            py = (py + iy)%var_fft->nb_modes_y;
          }        
        px = (px + ix)%var_fft->nb_modes_x;
      }
  }


  template<class T>
  void FftInterface<T>::
  ApplyInversePoint(int k, const Vector<T>& x, T& u)
  {
    SetComplexZero(u);
    int p = 0;
    for (int n = 0; n < var_fft->nb_modes_x; n++)
      {        
        u += var_fft->phase_x(p)*x(n);
        p = (p+k)%var_fft->nb_modes_x;
      }
    
    u *= Treal(1)/var_fft->nb_modes_x;
  }


  template<class T>
  void FftInterface<T>::
  ApplyInversePoint(int ix, int iy, const Vector<T>& x, T& u)
  {
    SetComplexZero(u);
    int px = 0;
    for (int nx = 0; nx < var_fft->nb_modes_x; nx++)
      {   
        int py = 0;
        for (int ny = 0; ny < var_fft->nb_modes_y; ny++)
          {
            u += var_fft->phase_x(px)*var_fft->phase_y(py)*x(nx*var_fft->nb_modes_y + ny);
            py = (py + iy)%var_fft->nb_modes_y;
          }
        
        px = (px + ix)%var_fft->nb_modes_x;
      }
    
    u *= Treal(1)/(var_fft->nb_modes_x*var_fft->nb_modes_y);
  }


  template<class T>
  void FftInterface<T>::
  ApplyInversePoint(int ix, int iy, int iz, const Vector<T>& x, T& u)
  {
    SetComplexZero(u);
    int px = 0;
    for (int nx = 0; nx < var_fft->nb_modes_x; nx++)
      {   
        int py = 0;
        for (int ny = 0; ny < var_fft->nb_modes_y; ny++)
          {
            int pz = 0;
            for (int nz = 0; nz < var_fft->nb_modes_z; nz++)
              {                
                u += var_fft->phase_x(px)*var_fft->phase_y(py)*var_fft->phase_z(pz)
		  *x(var_fft->nb_modes_z*(nx*var_fft->nb_modes_y + ny)+nz);
                pz = (pz + iz)%var_fft->nb_modes_z;
              }
            py = (py + iy)%var_fft->nb_modes_y;
          }        
        px = (px + ix)%var_fft->nb_modes_x;
      }
    
    u *= Treal(1)/(var_fft->nb_modes_x*var_fft->nb_modes_y*var_fft->nb_modes_z);
  }
  

  template<class T>
  T FftInterface<T>::GetCoefficient(int ix, int nx)
  {
    T coef = var_fft->phase_x((ix*nx)%var_fft->nb_modes_x)/Treal(var_fft->nb_modes_x);    
    return coef;
  }
    

  template<class T>
  T FftInterface<T>
  ::GetCoefficient(int ix, int iy, int nx, int ny)
  {
    T coef = var_fft->phase_x((ix*nx)%var_fft->nb_modes_x)*var_fft->phase_y((iy*ny)%var_fft->nb_modes_y)
      /Treal(var_fft->nb_modes_x*var_fft->nb_modes_y);
    
    return coef;
  }


  template<class T>
  T FftInterface<T>
  ::GetCoefficient(int ix, int iy, int iz, int nx, int ny, int nz)
  {
    T coef = var_fft->phase_x((ix*nx)%var_fft->nb_modes_x)*var_fft->phase_y((iy*ny)%var_fft->nb_modes_y)
      *var_fft->phase_z((iz*nz)%var_fft->nb_modes_z)
      / Treal(var_fft->nb_modes_x*var_fft->nb_modes_y*var_fft->nb_modes_z);
    
    return coef;
  }
  
  
  template<class T>
  void FftInterface<T>::GetCosSinAlpha(int n, Treal& cos_, Treal& sin_) const
  {
    cos_ = real(var_fft->phase_x(n));
    sin_ = imag(var_fft->phase_x(n));
  }


  /************************
   * VirtualFftRealObject *
   ************************/


  template<class T>
  VirtualFftRealObject<T>::VirtualFftRealObject()
  {
    nb_modes = 0; nb_modes_y = 0;
  }


  template<class T>
  VirtualFftRealObject<T>::~VirtualFftRealObject()
  {
  }


  template<class T>
  void VirtualFftRealObject<T>::SetNbThreads(int)
  {
  }
  

  template<class T>
  size_t VirtualFftRealObject<T>::GetMemorySize() const
  {
    size_t taille = data.GetMemorySize() + phase.GetMemorySize()
      + xsol.GetMemorySize();
    
    taille += sizeof(*this);
    return taille;
  }

  
  template<class T>
  void VirtualFftRealObject<T>::Init(int n)
  {
    nb_modes = n; nb_modes_y = 1;
    
    data.Reallocate(n);
    xsol.Reallocate(n/2+1);
  }


  template<class T>
  void VirtualFftRealObject<T>::Init(int nx, int ny)
  {
    nb_modes = nx; nb_modes_y = ny;
    
    data.Reallocate(nx*ny);
    xsol.Reallocate(nx*(ny/2+1));
  }


  /********************
   * GslFftRealObject *
   ********************/
  
  
#ifdef MONTJOIE_WITH_GSL  
  template<class T>
  GslFftRealObject<T>::GslFftRealObject()
  {
    wavetable_real = NULL;
    wavetable_cplx = NULL;
    workspace_real = NULL;
  }

  
  template<class T>
  GslFftRealObject<T>::~GslFftRealObject()
  {
    if (wavetable_real != NULL)
      {
	gsl_fft_real_wavetable_free(wavetable_real);
	wavetable_real = NULL;
      }
    
    if (wavetable_cplx != NULL)
      {
	gsl_fft_halfcomplex_wavetable_free(wavetable_cplx);
	wavetable_cplx = NULL;
      }
    
    if (workspace_real != NULL)
      {
	gsl_fft_real_workspace_free(workspace_real);
	workspace_real = NULL;
      }
  }


  template<class T>
  void GslFftRealObject<T>::Init(int n)
  {
    VirtualFftRealObject<T>::Init(n);
    if (n <= 0)
      return;
    
    wavetable_real = gsl_fft_real_wavetable_alloc(n);
    workspace_real = gsl_fft_real_workspace_alloc(n);
    wavetable_cplx = gsl_fft_halfcomplex_wavetable_alloc(n);
  }
  

  template<class T>
  void GslFftRealObject<T>::ApplyForward(const Vector<T>& x, Vector<complex<T> >& y)
  {
    if (this->nb_modes < 2)
      {
        y(0) = x(0);
        return;
      }
    
    for (int i = 0; i < this->nb_modes/2; i++)
      y(i) = complex<T>(x(2*i), x(2*i+1));
    
    y(this->nb_modes/2) = 0;
    
    gsl_fft_real_transform(reinterpret_cast<double*> (y.GetData()), 1,
                           this->nb_modes, wavetable_real, workspace_real);
    
    Real_wp yr = imag(y(0)), yi, yr_next;
    y(0) = real(y(0));
    for (int i = 1; i < this->nb_modes/2; i++)
      {
        yr_next = imag(y(i));
        yi = real(y(i));
        y(i) = complex<T>(yr, yi);
        yr = yr_next;
      }
    
    y(this->nb_modes/2) = yr;
  }

  
  template<class T>
  void GslFftRealObject<T>::ApplyInverse(const Vector<complex<T> >& x, Vector<T>& y)
  {
    y(0) = real(x(0));
    for (int i = 1; i < this->nb_modes/2; i++)
      {
        y(2*i-1) = real(x(i));
        y(2*i) = imag(x(i));
      }
    
    y(this->nb_modes-1) = real(x(this->nb_modes/2));
    
    gsl_fft_halfcomplex_inverse(reinterpret_cast<double*> (y.GetData()), 1,
                                this->nb_modes, wavetable_cplx, workspace_real);
    
  }
#endif


  /*********************
   * FftwFftRealObject *
   *********************/
  
  
#ifdef MONTJOIE_WITH_FFTW
  template<class T>
  FftwFftRealObject<T>::FftwFftRealObject()
  {
#ifdef _OPENMP
    nb_threads = omp_get_max_threads();
#else
    nb_threads = 1;
#endif
  }

  
  template<class T>
  FftwFftRealObject<T>::~FftwFftRealObject()
  {
    if (this->nb_modes > 0)
      {
	fftw_destroy_plan(plan_forward);
        fftw_destroy_plan(plan_backward);
      }
  }


  template<class T>
  void FftwFftRealObject<T>::SetNbThreads(int n)
  {
    nb_threads = n;
  }

  
  template<class T>
  void FftwFftRealObject<T>::Init(int n)
  {
    VirtualFftRealObject<T>::Init(n);
    if (n <= 0)
      return;

#ifdef _OPENMP
    fftw_plan_with_nthreads(nb_threads);
#endif

    fftw_complex* out = reinterpret_cast<fftw_complex*>(this->xsol.GetData());
    Real_wp* in = reinterpret_cast<T*>(this->data.GetData());
    plan_forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    
    plan_backward = fftw_plan_dft_c2r_1d(n, out, in, FFTW_ESTIMATE);    
  }


  template<class T>
  void FftwFftRealObject<T>::Init(int nx, int ny)
  {
    VirtualFftRealObject<T>::Init(nx, ny);
    if ((nx <= 0) || (ny <= 0))
      return;

#ifdef _OPENMP
    fftw_plan_with_nthreads(nb_threads);
#endif

    fftw_complex* out = reinterpret_cast<fftw_complex*>(this->xsol.GetData());
    Real_wp* in = reinterpret_cast<T*>(this->data.GetData());
    plan_forward = fftw_plan_dft_r2c_2d(nx, ny, in, out, FFTW_ESTIMATE);
    
    plan_backward = fftw_plan_dft_c2r_2d(nx, ny, out, in, FFTW_ESTIMATE);    
  }
  

  template<class T>
  void FftwFftRealObject<T>::ApplyForward(const Vector<T>& x, Vector<complex<T> >& y)
  {
    cout << "c'est bien cela" << endl;
    for (int i = 0; i < this->data.GetM(); i++)
      this->data(i) = x(i);
        
    fftw_execute(plan_forward);
        
    for (int i = 0; i < this->xsol.GetM(); i++)
      y(i) = this->xsol(i);
  }

  
  template<class T>
  void FftwFftRealObject<T>::ApplyInverse(const Vector<complex<T> >& x, Vector<T>& y)
  {
    for (int i = 0; i < this->xsol.GetM(); i++)
      this->xsol(i) = x(i);
    
    fftw_execute(plan_backward);
    
    for (int i = 0; i < this->data.GetM(); i++)
      y(i) = this->data(i);
    
    Mlt(1.0/(this->nb_modes*this->nb_modes_y), y);
  }
#endif


  /********************
   * MklFftRealObject *
   ********************/

  
#ifdef SELDON_WITH_MKL
  template<class T>
  MklFftRealObject<T>::MklFftRealObject()
  {
  }

  
  template<class T>
  MklFftRealObject<T>::~MklFftRealObject()
  {
    if (this->nb_modes > 0)
      {
	DftiFreeDescriptor(&hand_mkl);
      }

  }


  template<class T>
  void MklFftRealObject<T>::Init(int n)
  {
    VirtualFftRealObject<T>::Init(n);
    if (n <= 0)
      return;

    DftiCreateDescriptor(&hand_mkl, DFTI_DOUBLE, DFTI_REAL, 1, n);
    DftiSetValue(hand_mkl, DFTI_BACKWARD_SCALE, 1.0/n);
    DftiSetValue(hand_mkl, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    DftiSetValue(hand_mkl, DFTI_CONJUGATE_EVEN_STORAGE,
                 DFTI_COMPLEX_COMPLEX);
    
    //DftiSetValue(hand_mkl, DFTI_THREAD_LIMIT, nb_threads);
    DftiCommitDescriptor(hand_mkl);
  }
  

  template<class T>
  void MklFftRealObject<T>::ApplyForward(const Vector<T>& x, Vector<complex<T> >& y)
  {
    DftiComputeForward(hand_mkl, x.GetData(), y.GetDataVoid());
  }

  
  template<class T>
  void MklFftRealObject<T>::ApplyInverse(const Vector<complex<T> >& x, Vector<T>& y)
  {
    DftiComputeBackward(hand_mkl, x.GetDataVoid(), y.GetData());
  }
#endif

  
  /********************
   * FftRealInterface *
   ********************/
  
  
  FftRealInterface::FftRealInterface()
  {
    var_fft = NULL;
    int type_fft = MANUAL;
    
    // using mkl if possible, then fftw, then gsl, then manual
#ifdef MONTJOIE_WITH_GSL
    type_fft = FFT_GSL;
#endif

#ifdef MONTJOIE_WITH_FFTW
    type_fft = FFTW;
#endif
    
#ifdef SELDON_WITH_MKL
    type_fft = FFT_MKL;
#endif

    SelectFftAlgorithm(type_fft);
  }

  
  FftRealInterface::~FftRealInterface()
  {
    delete var_fft;
  }

  
  void FftRealInterface::SelectFftAlgorithm(int type)
  {
    if (var_fft != NULL)
      delete var_fft;
    
    switch (type)
      {
      case FFT_MKL:
	{
#ifdef SELDON_WITH_MKL
	  var_fft = new MklFftRealObject<Real_wp>;
#else
	  cout << "Recompile with SELDON_WITH_MKL" << endl;
	  abort();
#endif
	}
	break;
      case FFTW:
	{
#ifdef MONTJOIE_WITH_FFTW
	  var_fft = new FftwFftRealObject<Real_wp>;
#else
	  cout << "Recompile with MONTJOIE_WITH_FFTW" << endl;
	  abort();
#endif
	}
	break;
      case FFT_GSL:
	{
#ifdef MONTJOIE_WITH_GSL
	  var_fft = new GslFftRealObject<Real_wp>;
#else
	  cout << "Recompile with MONTJOIE_WITH_GSL" << endl;
	  abort();
#endif
	}
	break;
      default:
	cout << "Unknown type" << endl;
	abort();
      }
  }
  
  
  //! returns the number of points of the fft (1-D)
  int FftRealInterface::GetNbPoints() const
  {
    return var_fft->nb_modes;
  }
  
  
  //! returns the size used by the object in bytes
  size_t FftRealInterface::GetMemorySize() const
  {
    return var_fft->GetMemorySize();
  }
  
    
  void FftRealInterface::SetNbThreads(int n)
  {
    var_fft->SetNbThreads(n);
  }
  
  
  void FftRealInterface::Init(int n)
  {
    var_fft->Init(n);
  }

  
  void FftRealInterface::Init(int nx, int ny)
  {
    var_fft->Init(nx, ny);
  }

  
  void FftRealInterface::ApplyForward(const Vector<Real_wp>& x,
                                      Vector<Complex_wp>& y)
  {
    var_fft->ApplyForward(x, y);
  }

  
  void FftRealInterface::ApplyInverse(const Vector<Complex_wp>& x,
                                      Vector<Real_wp>& y)
  {
    var_fft->ApplyInverse(x, y);
  }
  
}

#define MONTJOIE_FILE_FFT_CXX
#endif
