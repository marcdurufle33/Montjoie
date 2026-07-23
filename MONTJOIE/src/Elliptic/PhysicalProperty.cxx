#ifndef MONTJOIE_FILE_PHYSICAL_PROPERTY_CXX

#include "OneDimPhysicalIndex.cxx"

namespace Montjoie
{
  
  /*******************
   * RectangleCutOff *
   *******************/
  
  
  //! default constructor
  template<class Dimension, class T>
  RectangleCutOff<Dimension, T>::RectangleCutOff()
  {
    type_cut_off = NO_CUTOFF;    
  }
  
  
  //! Returns the size used by the object in bytes
  template<class Dimension, class T>
  size_t RectangleCutOff<Dimension, T>::GetMemorySize() const
  {
    return sizeof(*this);
  }

  
  //! Initialisation of bounds and coefficients
  template<class Dimension, class T>
  void RectangleCutOff<Dimension, T>::
  Init(const R_N& xmin_, const R_N& xmax_, const R_N& coef_)
  {
    envelope(0) = xmin_;
    envelope(1) = xmax_;
    coef_gaussian = coef_;
    bool coef_null = true;
    T zero; SetComplexZero(zero);
    for (int k = 0; k < Dimension::dim_N; k++)
      if (coef_gaussian(k) != zero)
        coef_null = false;
    
    if (coef_null)
      type_cut_off = NO_CUTOFF;
    else
      type_cut_off = GAUSSIAN_CUTOFF;
    
  }
  
  
  //! computation of val = val*tronc(x)
  template<class Dimension, class T> template<class T0>
  void RectangleCutOff<Dimension, T>::ApplyCutOff(const R_N& point, T0& val) const
  {
    if (type_cut_off == NO_CUTOFF)
      {
        // no cut-off
        // we don't modify
        return;
      }

    T arg, val1(1), val2(1);
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        if (abs(point(i)-envelope(0)(i)) > R_N::threshold)
          {
            arg = -coef_gaussian(i)/square(Real_wp(point(i)-envelope(0)(i)));
            val1 = exp(arg);
          }
        else
          {
            if (coef_gaussian(i) != T(0))
              val1 = 0;
            else
              val1 = 1.0;
          }
        
        if (abs(point(i)-envelope(1)(i)) > R_N::threshold)
          {
            arg = -coef_gaussian(i)/square(Real_wp(point(i)-envelope(1)(i)));
            val2 = exp(arg);
          }
        else
          {
            if (coef_gaussian(i) != T(0))
              val2 = 0;
            else
              val2 = 1.0;
          }
        
        val *= val1*val2;
      }
  }
  
  
  //! computation of val = val*tronc(x) and grad(val tronc(x))
  template<class Dimension, class T> template<class T0>
  void RectangleCutOff<Dimension, T>::
  ApplyCutOff(const R_N& point, T0& val, TinyVector<T0, Dimension::dim_N>& grad) const
  {
    if (type_cut_off == NO_CUTOFF)
      {
        // no cut-off
        // we don't modify
        return;
      }
    
    T arg, val1(1), dval1(0), val2(1), dval2(0), prod, dprod;
    for (int i = 0; i < Dimension::dim_N; i++)
      {
        if (abs(point(i)-envelope(0)(i)) > R_N::threshold)
          {
            arg = -coef_gaussian(i)/square(Real_wp(point(i)-envelope(0)(i)));
            val1 = exp(arg);
            dval1 = -2.0*arg/(point(i) - envelope(0)(i))*val1;
          }
        else
          {
            if (coef_gaussian(i) != T(0))
              val1 = 0;
            else
              val1 = 1.0;
            
            dval1 = 0;
          }
        
        if (abs(point(i)-envelope(1)(i)) > R_N::threshold)
          {
            arg = -coef_gaussian(i)/square(Real_wp(point(i)-envelope(1)(i)));
            val2 = exp(arg);
            dval2 = -2.0*arg/(point(i) - envelope(1)(i))*val2;
          }
        else
          {
            if (coef_gaussian(i) != T(0))
              val2 = 0;
            else
              val2 = 1.0;
            
            dval2 = 0;
          }
        
        prod = val1*val2;
        dprod = dval1*val2 + val1*dval2;

        grad *= prod;
        grad(i) += dprod*val;
        val *= prod;
      }        
  }
  
  
  /***************************
   * PhysicalSinusoidalMedia *
   ***************************/
  

  //! Returns the size used by the object in bytes
  template<class Dimension, class T>
  size_t PhysicalSinusoidalMedia<Dimension, T>::GetMemorySize() const
  {
    return sizeof(*this);
  }

    
  //! computes index and its gradient at a given point
  template<class Dimension, class T>
  void PhysicalSinusoidalMedia<Dimension, T>::
  ComputeIndexAtPoint(const R_N& point, T& coef, TinyVector<T, Dimension::dim_N>& grad) const
  {
    bool non_null_coef = PointInsideBoundingBox(point, this->envelope);
    
    T val3, dval3, zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    if (non_null_coef)
      {
        // evaluation of the sinus
        coef = T(1);
        grad.Fill(zero);
        
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            val3 = sin(point(i)*kwave(i));
            dval3 = kwave(i)*cos(point(i)*kwave(i));
            
            grad *= val3;
            grad(i) += coef*dval3;
            coef *= val3;
          }
        
        // applying cut-off function
        this->ApplyCutOff(point, coef, grad);
        
        // last operations
        coef = coef*amplitude + mu0;
        grad *= amplitude;
                    
        if (inverse)
          {
            coef = one/coef;
            grad *= -coef*coef;
          }
      }
    else
      {
        if (inverse)
          coef = one/mu0;
        else
          coef = mu0;
        
        grad.Fill(zero);
      }
  }
  
  
  //! computes index at a given point
  template<class Dimension, class T>
  void PhysicalSinusoidalMedia<Dimension, T>::
  ComputeIndexAtPoint(const R_N& point, T& coef) const
  {
    bool non_null_coef = PointInsideBoundingBox(point, this->envelope);
    
    T val3, one;
    SetComplexOne(one);
    if (non_null_coef)
      {
        // evaluation of the sinus
        coef = one;

        for (int i = 0; i < Dimension::dim_N; i++)
          {
            val3 = sin(point(i)*kwave(i));
            coef *= val3;
          }

        // applying cut-off function
        this->ApplyCutOff(point, coef);
        
        // last operations
        coef = coef*amplitude + mu0;
                            
        if (inverse)
          coef = one/coef;
      }
    else
      {
        if (inverse)
          coef = one/mu0;
        else
          coef = mu0;
      }
  }
  
  
  //! returns an evaluation of physical index on quadrature points
  /*!
    \param[in] mesh computational mesh
    \param[in] ElementRho array containing local element numbers
    among elements of reference equal to ref
    \param[in] ref reference
    \param[in] Points quadrature points (on each element)
    \param[in] compute_grad if true, gradient is also computed
    \param[out] rho evaluation of index on quadrature points
    \param[out] grad_rho gradient of index on quadrature points
  */
  template<class Dimension, class T>
  void PhysicalSinusoidalMedia<Dimension, T>::
  ComputeIndex(const Mesh<Dimension>& mesh, const IVect& ElementRho,
               int ref, Vector<Vector<R_N> >& Points,
               bool compute_grad, bool compute_hess, Vector<Vector<T> >& rho,
               Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
	       Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& hess_rho)
  {
    for (int ne = 0; ne < mesh.GetNbElt(); ne++)
      if (mesh.Element(ne).GetReference() == ref)
        {
          int i1 = ElementRho(ne);
	  if (compute_hess)
	    {
	      abort();
	    }
	  
          if (compute_grad)
            {
              for (int j = 0; j < Points(ne).GetM(); j++)              
                ComputeIndexAtPoint(Points(ne)(j), rho(i1)(j), grad_rho(i1)(j));
            }
          else
            {
              for (int j = 0; j < Points(ne).GetM(); j++)              
                ComputeIndexAtPoint(Points(ne)(j), rho(i1)(j));
            }
        }
  }
  
  
  /************************
   * PhysicalRegularMedia *
   ************************/
  
  
  //! default constructor
  template<class Dimension, class T>
  PhysicalRegularMedia<Dimension, T>::PhysicalRegularMedia()
    : RectangleCutOff<Dimension, Real_wp>()
  {
    nbx.Fill(10);
    SetComplexZero(amplitude);
    SetComplexOne(mu0);
    inverse = false;
    complex_index = false;
    periodic = USUAL;
  }
  

  //! Returns the size used by the object in bytes
  template<class Dimension, class T>
  size_t PhysicalRegularMedia<Dimension, T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    taille += xstep.GetMemorySize() + ystep.GetMemorySize() + zstep.GetMemorySize();
    taille += nu.GetMemorySize();
    return taille;
  }

    
  //! kind of constructor
  /*!
    \param[in] xmin_ extrema xmin, ymin, zmin of the regular grid
    \param[in] xmax_ extrema xmax, ymax, zmax of the regular grid
    \param[in] coef_ cut-off coefficients
    \param[in] offset index is equal to offset + ampli*h
    \param[in] ampli index is equal to offset + ampli*h
    \param[in] dimension 2 or 3
    \param[in] data_file name of the binary file containing values
    \param[in] complex_number if true the data file is containing complex values
    \param[in] double_prec if true, the data file is written in double precision
  */
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::
  Init(const R_N& xmin_, const R_N & xmax_, const R_N& coef_, const T& offset, 
       const T& ampli, int dim_N, const string& data_file,
       bool complex_number, bool double_prec)
  {
    RectangleCutOff<Dimension, Real_wp>::Init(xmin_, xmax_, coef_);
    
    // we read coefficients h_i,j in the data file
    ReadFile(data_file, dim_N, complex_number, double_prec, T(0));
    
    // checking if the number of points is greater than 5
    for (int i = 0; i < dim_N; i++)
      {
        if (nbx(i) < 5)
          {
            cout << "The number of points on each dimension should be greater"
                 << "or equal to 5" << endl;
            
            abort();
          }
      }
    
    amplitude = ampli;
    mu0 = offset;    
    // delta x, delta y and delta z
    for (int i = 0; i < Dimension::dim_N; i++)
      step_delta_x(i) = (xmax_(i) - xmin_(i))/(this->nbx(i) - 1);
    
    inverse =  false;
  }
  
  
  //! clearing datas
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::Clear()
  {
    nu.Clear();
    xstep.Clear();
    ystep.Clear();
    zstep.Clear();
  }
  
  
  //! physical index is read in a file
  /*!
    \param[in] file_name file name
    \param[in] dim_N dimension of the grid (2 or 3)
    \param[in] complex_number true if numbers stored in file are complex
    \param[in] double_prec true if numbers stored in file are written in double precision
    \param[in] cte dummy parameter
  */
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::
  ReadFile(const string& file_name, int dim_N,
           bool complex_number, bool double_prec, const Real_wp& cte)
  {
    complex_index = complex_number;
    ifstream file_in(file_name.data());
    if (!file_in.is_open())
      {
	cout << "File containing random media doesn't exist " << endl;
	cout << file_name << endl;
	abort();
      }
    
    nbx.Fill(0);
    // file is read in binary format, float or double precision
    if (dim_N == 2)
      {
	if (double_prec)
	  {
	    Matrix<double> mat_coef;
	    mat_coef.Read(file_in);
	    int m = mat_coef.GetM();
            int n = mat_coef.GetN();
	    nu.Reallocate(m*n);
	    for (int i = 0; i < m; i++)
	      for (int j = 0; j < n; j++)
		nu(i*n + j) = mat_coef(i, j);
	    
	    nbx(0) = m;
	    nbx(1) = n;
	  }
	else
	  {
	    Matrix<float> mat_coef;
	    mat_coef.Read(file_in);
	    int m = mat_coef.GetM();
	    int n = mat_coef.GetN();
	    nu.Reallocate(m*n);
	    for (int i = 0; i < m; i++)
	      for (int j = 0; j < n; j++)
		nu(i*n+j) = mat_coef(i, j);
	    
	    nbx(0) = m;
	    nbx(1) = n;
	  }
      }
    else
      {
	if (double_prec)
	  {
	    Array3D<double> mat_coef;
	    mat_coef.Read(file_in);
	    nbx(0) = mat_coef.GetLength1();
            nbx(1) = mat_coef.GetLength2();
            nbx(2) = mat_coef.GetLength3();
	    nu.Reallocate(nbx(0)*nbx(1)*nbx(2));
	    int nb = 0;
	    for (int i = 0; i < nbx(0); i++)
	      for (int j = 0; j < nbx(1); j++)
		for (int k = 0; k < nbx(2); k++)
		  nu(nb++) = mat_coef(i,j,k);
	    
	  }
	else
	  {
	    Array3D<float> mat_coef;
	    mat_coef.Read(file_in);
	    nbx(0) = mat_coef.GetLength1(); nbx(1) = mat_coef.GetLength2();
            nbx(2) = mat_coef.GetLength3();
	    nu.Reallocate(nbx(0)*nbx(1)*nbx(2));
	    int nb = 0;
	    for (int i = 0; i < nbx(0); i++)
	      for (int j = 0; j < nbx(1); j++)
		for (int k = 0; k < nbx(2); k++)
		  nu(nb++) = mat_coef(i,j,k);
	    
	  }
      }
    
    file_in.close();
  }
  
    
  //! physical index is read in a file
  /*!
    \param[in] file_name file name
    \param[in] dim_N dimension of the grid (2 or 3)
    \param[in] complex_number true if numbers stored in file are complex
    \param[in] double_prec true if numbers stored in file are written in double precision
    \param[in] cte dummy parameter
  */
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::
  ReadFile(const string& file_name, int dim_N,
           bool complex_number, bool double_prec, const Complex_wp& cte)
  {
    complex_index = complex_number;
    ifstream file_in(file_name.data());
    if (!file_in.is_open())
      {
	cout<<"File containing random media doesn't exist "<<endl;
	cout<<file_name<<endl;
	abort();
      }
    
    nbx.Fill(0);
    if (dim_N == 2)
      {
	if (double_prec)
	  {
	    if (complex_number)
	      {
		Matrix<complex<double> > mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetM(); nbx(1) = mat_coef.GetN();
		nu.Reallocate(nbx(0)*nbx(1));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    to_complex(mat_coef(i,j), nu(nb++));
		
	      }
	    else
	      {
		Matrix<double> mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetM(); nbx(1) = mat_coef.GetN();
		nu.Reallocate(nbx(0)*nbx(1));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    to_complex(mat_coef(i,j), nu(nb++));
		
	      }
	  }
	else
	  {
	    if (complex_number)
	      {
		Matrix<complex<float> > mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetM();
                nbx(1) = mat_coef.GetN();
		nu.Reallocate(nbx(0)*nbx(1));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    to_complex(mat_coef(i,j), nu(nb++));
		
	      }
	    else
	      {
		Matrix<float> mat_coef; 
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetM();
                nbx(1) = mat_coef.GetN();
		nu.Reallocate(nbx(0)*nbx(1));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    to_complex(mat_coef(i,j), nu(nb++));
                
	      }            
	  }
      }
    else
      {
	if (double_prec)
	  {
	    if (complex_number)
	      {
		Array3D<complex<double> > mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetLength1();
                nbx(1) = mat_coef.GetLength2();
		nbx(2) = mat_coef.GetLength3();
		nu.Reallocate(nbx(0)*nbx(1)*nbx(2));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    for (int k = 0; k < nbx(2); k++)
		      to_complex(mat_coef(i,j,k), nu(nb++));
	      }
	    else
	      {
		Array3D<double> mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetLength1();
                nbx(1) = mat_coef.GetLength2();
		nbx(2) = mat_coef.GetLength3();
		nu.Reallocate(nbx(0)*nbx(1)*nbx(2));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    for (int k = 0; k < nbx(2); k++)
		      to_complex(mat_coef(i,j,k), nu(nb++));
	      }
	  }
	else
	  {
	    if (complex_number)
	      {
		Array3D<complex<float> > mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetLength1();
                nbx(1) = mat_coef.GetLength2();
		nbx(2) = mat_coef.GetLength3();
		nu.Reallocate(nbx(0)*nbx(1)*nbx(2));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    for (int k = 0; k < nbx(2); k++)
		      to_complex(mat_coef(i,j,k), nu(nb++));
	      }
	    else
	      {
		Array3D<float> mat_coef;
		mat_coef.Read(file_in);
		nbx(0) = mat_coef.GetLength1();
                nbx(1) = mat_coef.GetLength2();
		nbx(2) = mat_coef.GetLength3();
		nu.Reallocate(nbx(0)*nbx(1)*nbx(2));
		int nb = 0;
		for (int i = 0; i < nbx(0); i++)
		  for (int j = 0; j < nbx(1); j++)
		    for (int k = 0; k < nbx(2); k++)
		      to_complex(mat_coef(i,j,k), nu(nb++));
	      }
	  }
      }
    
    file_in.close();
  }
  
  
  //! computation of coef = h(point)  and grad_coef = nabla h (point)
  /*!
    Interpolation of values h_i,j in order to obtain h at a given point    
  */
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::
  ComputeIndexAtPoint(const R_N& point2, T& coef, 
                      TinyVector<T, Dimension::dim_N>& grad_coef) const
  {
    R_N point = point2;
    T one; SetComplexOne(one);
    T zero; SetComplexZero(zero);
    coef = mu0;
    if (inverse)
      coef = one/mu0;
    
    grad_coef.Zero();
    
    T val1, dval1, val2, dval2, val3, dval3, val;
    bool non_null_coef = true;
    if ((periodic == PERIODICITY)||(periodic == QUASI_PERIODICITY))
      {
        bool center_cell = true;
        // point is moved to the cell [xmin, xmax] by using periodicity
        for (int k = 0; k < Dimension::dim_N; k++)
          {
            Real_wp Nf = (point(k)-this->envelope(0)(k))/(this->envelope(1)(k)-this->envelope(0)(k));
            int Nx = toInteger(floor(Nf));
            if (abs(Nf - round(Nf)) <= epsilon_machine)
              Nx = toInteger(round(Nf));
            
            point(k) -= Nx*(this->envelope(1)(k) - this->envelope(0)(k));
            if (Nx != 0)
              center_cell = false;
          }
        
        if (periodic == QUASI_PERIODICITY)
          if (center_cell)
            non_null_coef = false;
      }
    else
      non_null_coef = PointInsideBoundingBox(point, this->envelope);
    
    if (non_null_coef)
      {
        R_N pt_loc;
        TinyVector<int, Dimension::dim_N> imin, imax;
        
        // evaluating coef by fourth-order interpolation
        
        // finding the local coordinates pt_loc
        // and left-lower point in imin, right-upper in imax
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            int i1 = toInteger(floor((point(i)-this->envelope(0)(i))/step_delta_x(i)));
            i1 = max(0, i1); i1 = min(nbx(i), i1);
            Real_wp dx = (point(i)-this->envelope(0)(i))/step_delta_x(i) - i1;
                
            pt_loc(i) = GetIndexBounds(i1, imin(i), imax(i), nbx(i), dx);
          }
        
        // evaluating Lagrange basis functions
        TinyVector<Real_wp, 5> phix, phiy, dphix, dphiy;                    
        phix(0) = 32.0/3*(pt_loc(0)-0.25)*(pt_loc(0)-0.5)*(pt_loc(0)-0.75)*(pt_loc(0)-1.0);
        dphix(0) = 1.0/3*(-25.0 + pt_loc(0)*(140.0 + pt_loc(0)*(-240.0 + pt_loc(0)*128.0)));
        
        phix(1) = -128.0/3*pt_loc(0)*(pt_loc(0)-0.5)*(pt_loc(0)-0.75)*(pt_loc(0)-1.0);
        dphix(1) = 1.0/3*(48.0 + pt_loc(0)*(-416.0 + pt_loc(0)*(864.0 - pt_loc(0)*512.0)));
        
        phix(2) = 64.0*pt_loc(0)*(pt_loc(0)-0.25)*(pt_loc(0)-0.75)*(pt_loc(0)-1.0);
        dphix(2) = -12.0 + pt_loc(0)*(152.0 + pt_loc(0)*(-384.0 + pt_loc(0)*256.0));
        
        phix(3) = -128.0/3*pt_loc(0)*(pt_loc(0)-0.25)*(pt_loc(0)-0.5)*(pt_loc(0)-1.0);
        dphix(3) = 1.0/3*(16.0 + pt_loc(0)*(-224.0 + pt_loc(0)*(672.0 - 512.0*pt_loc(0))));
        
        phix(4) = 32.0/3*pt_loc(0)*(pt_loc(0)-0.25)*(pt_loc(0)-0.5)*(pt_loc(0)-0.75);
        dphix(4) = 1.0/3*(-3.0 + pt_loc(0)*(44.0 + pt_loc(0)*(-144.0 + pt_loc(0)*128.0)));
        
        phiy(0) = 32.0/3*(pt_loc(1)-0.25)*(pt_loc(1)-0.5)*(pt_loc(1)-0.75)*(pt_loc(1)-1.0);
        dphiy(0) = 1.0/3*(-25.0 + pt_loc(1)*(140.0 + pt_loc(1)*(-240.0 + pt_loc(1)*128.0)));
        
        phiy(1) = -128.0/3*pt_loc(1)*(pt_loc(1)-0.5)*(pt_loc(1)-0.75)*(pt_loc(1)-1.0);
        dphiy(1) = 1.0/3*(48.0 + pt_loc(1)*(-416.0 + pt_loc(1)*(864.0 - pt_loc(1)*512.0)));
        
        phiy(2) = 64.0*pt_loc(1)*(pt_loc(1)-0.25)*(pt_loc(1)-0.75)*(pt_loc(1)-1.0);
        dphiy(2) = -12.0 + pt_loc(1)*(152.0 + pt_loc(1)*(-384.0 + pt_loc(1)*256.0));
        
        phiy(3) = -128/3.0*pt_loc(1)*(pt_loc(1)-0.25)*(pt_loc(1)-0.5)*(pt_loc(1)-1.0);
        dphiy(3) = 1.0/3*(16.0 + pt_loc(1)*(-224.0 + pt_loc(1)*(672.0 - 512.0*pt_loc(1))));
        
        phiy(4) = 32/3.0*pt_loc(1)*(pt_loc(1)-0.25)*(pt_loc(1)-0.5)*(pt_loc(1)-0.75);
        dphiy(4) = 1.0/3*(-3.0 + pt_loc(1)*(44.0 + pt_loc(1)*(-144.0 + pt_loc(1)*128.0)));
        
        dphix *= Real_wp(0.25)/step_delta_x(0);
        dphiy *= Real_wp(0.25)/step_delta_x(1);
            
        coef = 0;
        grad_coef.Fill(zero);
        if (Dimension::dim_N == 2)
          {                  
            for (int i = imin(0); i <= imax(0); i++)
              for (int j = imin(1); j <= imax(1); j++)
                {
                  val1 = phix(i-imin(0));
                  val2 = phiy(j-imin(1));
                  dval1 = dphix(i-imin(0));
                  dval2 = dphiy(j-imin(1));
                  coef += val1*val2*nu(i*nbx(1) + j);
                  grad_coef(0) += dval1*val2*nu(i*nbx(1) + j);
                  grad_coef(1) += val1*dval2*nu(i*nbx(1) + j);
                }
          }
        else
          {
            TinyVector<Real_wp, 5> phiz, dphiz;
            phiz(0) = 32.0/3*(pt_loc(2)-0.25)*(pt_loc(2)-0.5)*(pt_loc(2)-0.75)*(pt_loc(2)-1.0);
            dphiz(0) = 1.0/3*(-25.0 + pt_loc(2)*(140.0 + pt_loc(2)*(-240.0 + pt_loc(2)*128.0)));
            
            phiz(1) = -128.0/3*pt_loc(2)*(pt_loc(2)-0.5)*(pt_loc(2)-0.75)*(pt_loc(2)-1.0);
            dphiz(1) = 1.0/3*(48.0 + pt_loc(2)*(-416.0 + pt_loc(2)*(864.0 - pt_loc(2)*512.0)));
            
            phiz(2) = 64.0*pt_loc(2)*(pt_loc(2)-0.25)*(pt_loc(2)-0.75)*(pt_loc(2)-1.0);
            dphiz(2) = -12.0 + pt_loc(2)*(152.0 + pt_loc(2)*(-384.0 + pt_loc(2)*256.0));
            
            phiz(3) = -128/3.0*pt_loc(2)*(pt_loc(2)-0.25)*(pt_loc(2)-0.5)*(pt_loc(2)-1.0);
            dphiz(3) = 1.0/3*(16.0 + pt_loc(2)*(-224.0 + pt_loc(2)*(672.0 - 512.0*pt_loc(2))));
            
            phiz(4) = 32/3.0*pt_loc(2)*(pt_loc(2)-0.25)*(pt_loc(2)-0.5)*(pt_loc(2)-0.75);
            dphiz(4) = 1.0/3*(-3.0 + pt_loc(2)*(44.0 + pt_loc(2)*(-144.0 + pt_loc(2)*128.0)));
            
            dphiz *= 4.0/step_delta_x(2);
            
            for (int i = imin(0); i <= imax(0); i++)
              for (int j = imin(1); j <= imax(1); j++)
                for (int k = imin(2); k <= imax(2); k++)
                  {
                    val1 = phix(i-imin(0));
                    val2 = phiy(j-imin(1));
                    val3 = phiz(k-imin(2));
                    dval1 = dphix(i-imin(0));
                    dval2 = dphiy(j-imin(1));
                    dval3 = dphiz(k-imin(2));
                    val = nu(nbx(2)*(i*nbx(1) + j) + k);
                    coef += val1*val2*val3*val;
                    grad_coef(0) += dval1*val2*val3*val;
                    grad_coef(1) += val1*dval2*val3*val;
                    grad_coef(2) += val1*val2*dval3*val;
                  }
            
          }
        
        // applying cut-off
        this->ApplyCutOff(point, coef, grad_coef);
        
        // last operations
        grad_coef *= amplitude;
        coef = mu0 + amplitude*coef;
        if (inverse)
          {
            coef = T(1)/coef;
            grad_coef *= -coef*coef;
          }
            
      }
  }

  
  //! computation of coef = h(point)  and grad_coef = nabla h (point)
  /*!
    Interpolation of values h_i,j in order to obtain h at a given point    
  */
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::
  ComputeIndexAtPoint(const R_N& point2, T& coef) const
  {
    R_N point = point2;
    T one; SetComplexOne(one);
    coef = mu0;
    if (inverse)
      coef = one/mu0;
    
    bool non_null_coef = true;
    if ((periodic == PERIODICITY)||(periodic == QUASI_PERIODICITY))
      {
        bool center_cell = true;
        // point is moved to the cell [xmin, xmax] by using periodicity
        for (int k = 0; k < Dimension::dim_N; k++)
          {
            Real_wp Nf = (point(k)-this->envelope(0)(k))/(this->envelope(1)(k)-this->envelope(0)(k));
            int Nx = toInteger(floor(Nf));
            if (abs(Nf - round(Nf)) <= epsilon_machine)
              Nx = toInteger(round(Nf));
            
            point(k) -= Nx*(this->envelope(1)(k) - this->envelope(0)(k));
            if (Nx != 0)
              center_cell = false;
          }
        
        if (periodic == QUASI_PERIODICITY)
          if (center_cell)
            non_null_coef = false;
      }
    else
      non_null_coef = PointInsideBoundingBox(point, this->envelope);
    
    if (non_null_coef)
      {
        R_N pt_loc;
        TinyVector<int, Dimension::dim_N> imin, imax;
        
        // evaluating coef by fourth-order interpolation
        
        // finding the local coordinates pt_loc
        // and left-lower point in imin, right-upper in imax
        for (int i = 0; i < Dimension::dim_N; i++)
          {
            int i1 = toInteger(floor((point(i)-this->envelope(0)(i))/step_delta_x(i)));
            i1 = max(0, i1); i1 = min(nbx(i), i1);
            Real_wp dx = (point(i)-this->envelope(0)(i))/step_delta_x(i) - i1;
            
            pt_loc(i) = GetIndexBounds(i1, imin(i), imax(i), nbx(i), dx);
          }
        
        // evaluating Lagrange basis functions
        TinyVector<Real_wp, 5> phix, phiy;
        phix(0) = 10.666666666666666
          *(pt_loc(0)-0.25)*(pt_loc(0)-0.5)*(pt_loc(0)-0.75)*(pt_loc(0)-1.0);
        phix(1) = -42.666666666666666*pt_loc(0)*(pt_loc(0)-0.5)*(pt_loc(0)-0.75)*(pt_loc(0)-1.0);
        phix(2) = 64.0*pt_loc(0)*(pt_loc(0)-0.25)*(pt_loc(0)-0.75)*(pt_loc(0)-1.0);
        phix(3) = -42.666666666666666*pt_loc(0)*(pt_loc(0)-0.25)*(pt_loc(0)-0.5)*(pt_loc(0)-1.0);
        phix(4) = 10.666666666666666*pt_loc(0)*(pt_loc(0)-0.25)*(pt_loc(0)-0.5)*(pt_loc(0)-0.75);
        
        phiy(0) = 10.666666666666666
          *(pt_loc(1)-0.25)*(pt_loc(1)-0.5)*(pt_loc(1)-0.75)*(pt_loc(1)-1.0);
        phiy(1) = -42.666666666666666*pt_loc(1)*(pt_loc(1)-0.5)*(pt_loc(1)-0.75)*(pt_loc(1)-1.0);
        phiy(2) = 64.0*pt_loc(1)*(pt_loc(1)-0.25)*(pt_loc(1)-0.75)*(pt_loc(1)-1.0);
        phiy(3) = -42.666666666666666*pt_loc(1)*(pt_loc(1)-0.25)*(pt_loc(1)-0.5)*(pt_loc(1)-1.0);
        phiy(4) = 10.666666666666666*pt_loc(1)*(pt_loc(1)-0.25)*(pt_loc(1)-0.5)*(pt_loc(1)-0.75);
        
        coef = 0;
        if (Dimension::dim_N == 2)
          {
            for (int i = imin(0); i <= imax(0); i++)
              for (int j = imin(1); j <= imax(1); j++)
                coef += phix(i-imin(0))*phiy(j-imin(1))*nu(i*nbx(1) + j);
          }
        else
          {
            TinyVector<Real_wp, 5> phiz;
            phiz(0) = 32.0/3*(pt_loc(2)-0.25)*(pt_loc(2)-0.5)*(pt_loc(2)-0.75)*(pt_loc(2)-1.0);
            phiz(1) = -128.0/3*pt_loc(2)*(pt_loc(2)-0.5)*(pt_loc(2)-0.75)*(pt_loc(2)-1.0);
            phiz(2) = 64.0*pt_loc(2)*(pt_loc(2)-0.25)*(pt_loc(2)-0.75)*(pt_loc(2)-1.0);
            phiz(3) = -128/3.0*pt_loc(2)*(pt_loc(2)-0.25)*(pt_loc(2)-0.5)*(pt_loc(2)-1.0);
            phiz(4) = 32/3.0*pt_loc(2)*(pt_loc(2)-0.25)*(pt_loc(2)-0.5)*(pt_loc(2)-0.75);
            
            for (int i = imin(0); i <= imax(0); i++)
              for (int j = imin(1); j <= imax(1); j++)
                for (int k = imin(2); k <= imax(2); k++)
                  coef += phix(i-imin(0))*phiy(j-imin(1))*phiz(k-imin(2))
                    *nu(nbx(2)*(i*nbx(1) + j) + k);
          }
        
        // applying cut-off
        this->ApplyCutOff(point, coef);
        
        // last operations
        coef = mu0 + amplitude*coef;
        if (inverse)
          coef = T(1)/coef;
        
      }
  }
  
  
  //! returns an evaluation of physical index on quadrature points
  /*!
    \param[in] mesh computational mesh
    \param[in] ElementRho array containing local element numbers
    among elements of reference equal to ref
    \param[in] ref reference
    \param[in] Points quadrature points (on each element)
    \param[in] compute_grad if true, gradient is also computed
    \param[out] rho evaluation of index on quadrature points
    \param[out] grad_rho gradient of index on quadrature points
  */
  template<class Dimension, class T>
  void PhysicalRegularMedia<Dimension, T>::
  ComputeIndex(const Mesh<Dimension>& mesh, const IVect& ElementRho,
               int ref, Vector<Vector<R_N> >& Points,
               bool compute_grad, bool compute_hess, Vector<Vector<T> >& rho,
               Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
	       Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& hess_rho)
  {
    for (int ne = 0; ne < mesh.GetNbElt(); ne++)
      if (mesh.Element(ne).GetReference() == ref)
        {
          int i1 = ElementRho(ne);
	  if (compute_hess)
	    {
	      abort();
	    }
	  
          if (compute_grad)
            {
              for (int j = 0; j < Points(ne).GetM(); j++)              
                ComputeIndexAtPoint(Points(ne)(j), rho(i1)(j), grad_rho(i1)(j));
            }
          else
            {
              for (int j = 0; j < Points(ne).GetM(); j++)              
                ComputeIndexAtPoint(Points(ne)(j), rho(i1)(j));
            }
        }
  }
  
    
  /*********************
   * PhysicalMeshMedia *
   *********************/
  

  //! Returns the size used by the object in bytes
  template<class Dimension, class T>
  size_t PhysicalMeshMedia<Dimension, T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + Seldon::GetMemorySize(nu)
      + Seldon::GetMemorySize(grad_nu) + Seldon::GetMemorySize(hessian_nu);
    
    return taille;
  }

  
  //! values of index are cleared
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>::Clear()
  {
    nu.Clear();
    grad_nu.Clear();
    hessian_nu.Clear();
  }
  
  
  //! physical index is read in a file
  /*!
    \param[in] name_mesh file name for mesh (file .mesh)
    \param[in] file_name file name for index (file .bb)
    \param[in] order order of approximation
    The mesh is a Q1 mesh, whereas the index can be an index
    computed by a high-order finite element method
  */
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>::ReadFile(const string & name_mesh,
                                                 const string& file_name, int order) 
  {
    mesh_file_name = name_mesh;
    mesh_order = order;
    ReadMeshData(nu, file_name);
  }


  //! kind of constructor
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>::Init(const string& file_mesh, const T& offset,
					     const T& amplitude,
					     int dim_N, const string& data_file, int order)
  {
    ReadFile(file_mesh, data_file, order);
    for (int i = 0; i < nu.GetM(); i++)
      {
        for (int j = 0; j < nu(i).GetM(); j++)
          nu(i)(j) = offset + amplitude*nu(i)(j);
      }
  }

  
  //! we keep only values on local elements (elements of the current processor)
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>::RemoveGlobalElements(const IVect& NumLoc)
  {
    Vector<Vector<T> > nu_glob = nu;
    bool grad = false, hess = false;
    if (grad_nu.GetM() > 0)
      grad = true;
    
    if (hessian_nu.GetM() > 0)
      hess = true;
    
    Vector<Vector<TinyVector<T, Dimension::dim_N> > > grad_nu_glob = grad_nu;
    Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >
      hessian_nu_glob = hessian_nu;
    
    int nb_elt = NumLoc.GetM();
    nu.Reallocate(nb_elt);
    if (grad)
      grad_nu.Reallocate(nb_elt);
    
    if (hess)
      hessian_nu.Reallocate(nb_elt);
    
    for (int i = 0; i < nb_elt; i++)
      {
        int ne = NumLoc(i);
        nu(i) = nu_glob(ne);
        if (grad)
	  grad_nu(i) = grad_nu_glob(ne);
	
	if (hess)
	  hessian_nu(i) = hessian_nu_glob(ne);
      }
  }
  
  
  //! checking if the mesh is coherent with the index read nu
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>::CheckMesh(const Mesh<Dimension>& mesh)
  {
    if (nu.GetM() != mesh.GetNbElt())
      {
        cout << "The mesh contains " << mesh.GetNbElt() << " elements whereas "
             << "the index is written with " << nu.GetM() << " elements" << endl;
        abort();
      }
    
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        int nb_nodes = mesh.GetNbPointsNodalElt(i);
        if (nu(i).GetM() != nb_nodes)
          {
            cout << "The number of nodes differ on element " << i << endl;
            cout << "Number of nodes for nu = " << nu(i).GetM() << endl; 
            cout << "Number of nodes for the mesh = " << nb_nodes << endl; 
            abort();
          }
      }
  }

  
  //! computation of gradient of the index
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>::ComputeGradient(const Mesh<Dimension>& mesh,
							bool hessian)
  {
    typedef typename Dimension::R_N R_N;
    Vector<R_N> s;
    SetPoints<Dimension> PointsElem;
    SetMatrices<Dimension> MatricesElem;
    grad_nu.Reallocate(nu.GetM());
    if (hessian)
      hessian_nu.Reallocate(nu.GetM());
    
    Vector<T> nu_loc;
    Vector<TinyVector<T, Dimension::dim_N> > hessian_loc;
    // main loop over elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        // computation of points and jacobian matrices
        mesh.GetVerticesElement(i, s);
        mesh.FjElemNodal(s, PointsElem, mesh, i);
        mesh.DFjElemNodal(s, PointsElem, MatricesElem, mesh, i);
        
	// allocating arrays
        int nb_nodes = nu(i).GetM();
        grad_nu(i).Reallocate(nb_nodes);
	if (hessian)
	  hessian_nu(i).Reallocate(nb_nodes);
	
        // effective computation of gradient
        mesh.ComputeNodalGradient(MatricesElem, nu(i), grad_nu(i), i);
	
	// computing hessian if needed
	if (hessian)
	  {
	    nu_loc.Reallocate(nb_nodes);
	    for (int k = 0; k < Dimension::dim_N; k++)
	      {
		for (int j = 0; j < nb_nodes; j++)
		  nu_loc(j) = grad_nu(i)(j)(k);
		
		mesh.ComputeNodalGradient(MatricesElem, nu_loc, hessian_loc, i);
		
		for (int j = 0; j < nb_nodes; j++)
		  for (int m = 0; m < Dimension::dim_N; m++)
		    hessian_nu(i)(j)(k, m) = hessian_loc(j)(m);
	      }
	  }
      }
  }

  
  //! evaluates the index on a point of element i
  template<class Dimension, class T>
  void PhysicalMeshMedia<Dimension, T>
  ::ComputeIndexAtPoint(int i, const VectReal_wp& phi, T& coef) const
  {
    if (phi.GetM() != nu(i).GetM())
      {
        cout << "Use the same order for the index and the discretization" << endl;
        abort();
      }
    
    SetComplexZero(coef);
    for (int j = 0; j < phi.GetM(); j++)
      coef += phi(j)*nu(i)(j);
  }
  
  
  /**********************
   * RadialVaryingMedia *
   **********************/
  

  //! Returns the size used by the object in bytes
  template<class Dimension, class T>
  size_t RadialVaryingMedia<Dimension, T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this);
    if (index_radial != NULL)
      taille += index_radial->GetMemorySize();
    
    return taille;
  }


  //! inits the 1-D index with cubic spline
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>
  ::InitSpline(const string& data_file, const T& offset, const T& amplitude, T& cte)
  {
    Clear();
    
    SplineInterpolationIndex<T>* spline;
    spline = new SplineInterpolationIndex<T>();

    // values are read on a file
    Matrix<T> data;
    data.ReadText(data_file);

    // first column : abscissae x_i
    // second column : values of the index at at points x_i
    int N = data.GetM();
    Vector<Real_wp> x(N);
    Vector<T> y(N);
    for (int i = 0; i < N; i++)
      {
        x(i) = realpart(data(i, 0));
        y(i) = offset + amplitude*data(i, 1);
      }
    
    spline->Init(x, y, cte, false);
    index_radial = spline;
  }
  

  //! returns true if the index is complex
  template<class Dimension, class T>
  bool RadialVaryingMedia<Dimension, T>::IsComplex() const
  {
    return complex;
  }
  
  
  //! tells to return the inverse
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>::GetInverse()
  {
    inverse = !inverse;
  }
  
  
  //! multiplies values by a given coefficient
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>::Mlt(const T& coef)
  {
    if (index_radial != NULL)
      index_radial->Mlt(realpart(coef));
  }
  
  
  //! releases memory used by the arrays
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>::Clear()
  {
    if (index_radial != NULL)
      {
	delete index_radial;
	index_radial = NULL;
      }
  }
  
  
  //! evaluates at a given point
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>::ComputeIndexAtPoint(const R_N& point, T& coef) const
  {
    Real_wp r = Norm2(point);
    coef = index_radial->Evaluate(r);
    if (inverse)
      coef = 1.0/coef;
  }
  
  
  //! evaluates index and its gradient at a given point
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>
  ::ComputeIndexAtPoint(const R_N& point, T& coef, 
                        TinyVector<T, Dimension::dim_N>& grad_coef,
                        TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>& hess_coef,
                        bool hessian) const
  {
    T y, dy, d2y, zero;

    SetComplexZero(zero);
    Real_wp r = Norm2(point);
    if (hessian)
      index_radial->EvaluateSecondDerivative(r, y, dy, d2y);
    else
      index_radial->EvaluateDerivative(r, y, dy);
    
    if (inverse) 
      {
        if (abs(y) > epsilon_machine)
          {
            // g = 1/y
            T g, dg, d2g;
            g = Real_wp(1) / y;
            dg = -dy * g*g;
            d2g = g*g*(-d2y + Real_wp(2)*dy*dy*g);

            // y, dy d2y are substituted by g, dg, d2g
            y = g; dy = dg; d2y = d2g;
          }
      }
    
    Real_wp invR = Real_wp(1)/r;
    TinyVector<Real_wp, Dimension::dim_N> gradR;
    gradR = invR*point;

    coef = y;
    grad_coef = dy*gradR;
    
    if (hessian)
      {
        // hessian of r : I/r - x_i x_j / r^3
        TinyMatrix<Real_wp, Symmetric, Dimension::dim_N, Dimension::dim_N> hess_r, grad_rirj;
        GetNormalProjector(point, grad_rirj);
        grad_rirj *= invR*invR*invR;
        
        hess_r.SetDiagonal(invR);
        hess_r -= grad_rirj;
        
        // hessian of f : f''(r) dr/dx_i dr/dx_j + f'(r) hess_r
        GetNormalProjector(gradR, grad_rirj);
        hess_coef = d2y * grad_rirj + dy*hess_r;        
      }    
  }


  //! evalues index and its gradient for a set of points
  template<class Dimension, class T>
  void RadialVaryingMedia<Dimension, T>
  ::ComputeIndex(const Mesh<Dimension>& mesh, const IVect& ElementRho, int ref,
                 Vector<Vector<R_N> >& Points, bool compute_grad, bool compute_hess, Vector<Vector<T> >& rho,
                 Vector<Vector<TinyVector<T, Dimension::dim_N> > >& grad_rho,
		 Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >& hess_rho)
  {
    T rho0; TinyVector<T, Dimension::dim_N> grad_rho0;
    TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> hess_rho0;
    for (int ne = 0; ne < mesh.GetNbElt(); ne++)
      if (mesh.Element(ne).GetReference() == ref)
        {
          int i1 = ElementRho(ne);
          for (int j = 0; j < Points(ne).GetM(); j++)
            {

              if (compute_grad || compute_hess)
                ComputeIndexAtPoint(Points(ne)(j), rho0, grad_rho0, hess_rho0, compute_hess);
              else
                ComputeIndexAtPoint(Points(ne)(j), rho0);
              
              rho(i1)(j) = rho0;
              if (compute_grad)
                grad_rho(i1)(j) = grad_rho0;
              
              if (compute_hess)
                hess_rho(i1)(j) = hess_rho0;
            }
        }
  }
  
  
  /************************
   * PhysicalVaryingMedia *
   ************************/
  

  //! true if the index is complex
  template<class Dimension, class T>
  bool PhysicalVaryingMedia<Dimension, T>::IsComplex() const
  {
    switch (type)
      {
      case CONSTANT :
        {
          if (abs(imagpart(cte_coef)) > ScalarPhysicalIndice<Dimension, T>::threshold)
            return true;
        }
        break;
      case SINUSOIDE :
        return index_sinus.IsComplex();
      case REGULAR_GRID :
        return index_regular.IsComplex();
      case MESH :
        return index_mesh.IsComplex();
      case RADIAL :
        return index_radial.IsComplex();
      case USER :
        for (int i = 0; i < eval_coef.GetM(); i++)
          for (int j = 0; j < eval_coef(i).GetM(); j++)
            if (abs(imagpart(eval_coef(i)(j)))
                > ScalarPhysicalIndice<Dimension, T>::threshold)
              return true;
        
        break;
      }
    
    return false;
  }
  
  
  //! user asks to inverse physical index
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::GetInverse()
  {
    T one; SetComplexOne(one);
    if (type == CONSTANT)
      {
	if (cte_coef != T(0))
	  cte_coef = T(1)/cte_coef;
      }
    else
      {
        if (type == SINUSOIDE)
          index_sinus.GetInverse();
        if (type == REGULAR_GRID)
          index_regular.GetInverse();
        if (type == RADIAL)
          index_radial.GetInverse();
        if (type == MESH)
          index_mesh.GetInverse();
        
        // inverting discrete values        
        for (int i = 0; i < eval_coef.GetM(); i++)
          for (int j = 0; j < eval_coef(i).GetM(); j++)
            eval_coef(i)(j) = one/eval_coef(i)(j);
          
        // modification of gradient too
        for (int i = 0; i < grad_coef.GetM(); i++)
          for (int j = 0; j < grad_coef(i).GetM(); j++)
            grad_coef(i)(j) *= -square(eval_coef(i)(j));
        
        if (hessian_coef.GetM() > 0)
          {
            cout << "hessian of inverse not implemented" << endl;
            abort();
          }        
      }
  }
  
  
  //! multiplication of the physical index by a coefficient
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::Mlt(const T& coef)
  {
    switch (type)
      {
      case CONSTANT : cte_coef *= coef; break;
      case SINUSOIDE : index_sinus.Mlt(coef); break;
      case REGULAR_GRID : index_regular.Mlt(coef); break;
      case RADIAL : index_radial.Mlt(coef); break;
      case MESH : index_mesh.Mlt(coef); break;
      }
    
    // multiplying discrete values        
    for (int i = 0; i < eval_coef.GetM(); i++)
      for (int j = 0; j < eval_coef(i).GetM(); j++)
        eval_coef(i)(j) *= coef;
          
    for (int i = 0; i < grad_coef.GetM(); i++)
      for (int j = 0; j < grad_coef(i).GetM(); j++)
        grad_coef(i)(j) *= coef;

    for (int i = 0; i < hessian_coef.GetM(); i++)
      for (int j = 0; j < hessian_coef(i).GetM(); j++)
        hessian_coef(i)(j) *= coef;
    
  }
  
  
  //! sets the index to zero
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::Zero()
  {
    Clear();
    
    type = CONSTANT;
    SetComplexZero(cte_coef);
  }
  
  
  //! clears index
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::Clear()
  {
    eval_coef.Clear();
    grad_coef.Clear();
    hessian_coef.Clear();
    index_regular.Clear();
    index_mesh.Clear();
    index_radial.Clear();
  }


  //! clears gradient and hessian
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::ClearHessianGradient()
  {
    grad_coef.Clear();
    hessian_coef.Clear();
  }
  
  
  //! evaluates the index on element i and global point pt_glob
  template<class Dimension, class T>
  T PhysicalVaryingMedia<Dimension, T>
  ::GetIndexAtPoint(int i, const VectReal_wp& phi, const R_N& pt_glob) const
  {
    T val;
    switch (type)
      {
      case CONSTANT :
        val = cte_coef;
        break;
      case SINUSOIDE :
        index_sinus.ComputeIndexAtPoint(pt_glob, val);
        break;
      case REGULAR_GRID :
        index_regular.ComputeIndexAtPoint(pt_glob, val);
        break;
      case MESH :
        index_mesh.ComputeIndexAtPoint(i, phi, val);
        break;
      case RADIAL :
        index_radial.ComputeIndexAtPoint(pt_glob, val);
        break;
      default :
        cout << "not implemented" << endl;
        abort();
      }
    
    return val;
  }
  
  
  //! sets physical index to a constant
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::InitConstant(const T& constante)
  {
    type = CONSTANT;
    cte_coef = constante;
  }
  
    
  //! sets physical index to a sinusoidal field
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::
  InitSinus(const R_N& xmin, const R_N& xmax, const R_N& coef,
	    const T& offset, const T& amplitude, const R_N& kwave)
  {
    type = SINUSOIDE;
    index_sinus.Init(xmin, xmax, coef, offset, amplitude, kwave);
  }
  
  
  //! sets physical index to a field defined on a regular grid
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::
  InitRandom(const R_N& xmin, const R_N& xmax, const R_N& coef,
	     const T& offset, const T& amplitude,
	     int dim_N, const string& data_file,
	     bool complex_number, bool double_prec)
  {
    type = REGULAR_GRID;
    index_regular.Init(xmin, xmax, coef, offset, amplitude, dim_N,
		       data_file, complex_number, double_prec);
  }
  
  
  //! sets physical index to a field defined on a mesh
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::
  InitMesh(const string& mesh_file, const T& offset, const T& amplitude,
	   int dim_N, const string& data_file, int order)
  {
    type = MESH;
    index_mesh.Init(mesh_file, offset, amplitude, dim_N, data_file, order);
  }
  
  
  //! sets physical index to a field defined radially
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>
  ::InitRadial(const string& data_file, const T& offset,
               const T& amplitude, T& cte)
  {
    type = RADIAL;
    
    // cubic splines are used
    index_radial.InitSpline(data_file, offset, amplitude, cte_coef);
    cte = cte_coef;
  }


  //! sets physical index to a field defined by the user
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>
  ::InitUser(const T& offset, const T& amplitude, const T& cte)
  {
    type = USER;
    offset_coef = offset;
    amplitude_coef = amplitude;
    cte_coef = cte;
  }
  
  
  //! Returns the memory used by the object in bytes
  template<class Dimension, class T>
  size_t PhysicalVaryingMedia<Dimension, T>::GetMemorySize() const
  {
    size_t taille = sizeof(T)*15;
    taille += Seldon::GetMemorySize(eval_coef) + Seldon::GetMemorySize(grad_coef)
      + Seldon::GetMemorySize(hessian_coef) + index_regular.GetMemorySize()
      + index_sinus.GetMemorySize() + index_mesh.GetMemorySize()
      + index_radial.GetMemorySize();
    
    return taille;
  }


  //! allocation of array eval_coef
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::Reallocate(int N, bool compute_grad, bool compute_hess)
  {
    if (type != CONSTANT)
      {
        eval_coef.Reallocate(N);
        if (compute_grad)
          grad_coef.Reallocate(N);
	
	if (compute_hess)
	  hessian_coef.Reallocate(N);
      }
  }
  
  
  //! allocation of array eval_coef
  template<class Dimension, class T>
  void PhysicalVaryingMedia<Dimension, T>::
  ReallocateOnElement(int i, int N, bool compute_grad, bool compute_hess)
  {
    if (type != CONSTANT)
      {
        eval_coef(i).Reallocate(N);
        if (compute_grad)
          grad_coef(i).Reallocate(N);
	
	if (compute_hess)
	  hessian_coef(i).Reallocate(N);
      }
  }
  
  
  /************************
   * ScalarPhysicalIndice *
   ************************/


  //! Returns the memory used by the object in bytes
  template<class Dimension, class T>
  size_t ScalarPhysicalIndice<Dimension, T>::GetMemorySize() const
  {
    size_t taille = sizeof(T) + sizeof(int) + fct_rho.GetMemorySize();
    return taille;
  }  
  
  
  //! returns the number of varying coefficients (0 or 1)
  template<class Dimension, class T>
  int ScalarPhysicalIndice<Dimension, T>::GetNbVaryingMedia() const
  {
    if (type_media == VARYING)
      return 1;
    
    return 0;
  }
  
  
  //! provides variables coefficients in array rho
  template<class Dimension, class T>
  void ScalarPhysicalIndice<Dimension, T>::
  GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho,
                  IVect& num_component)
  {
    if (type_media == VARYING)
      {
        rho(nb) = &fct_rho;
        num_component(nb) = 0;
        nb++;
      }
  }
  
  
  //! clearing variable coefficients
  template<class Dimension, class T>
  void ScalarPhysicalIndice<Dimension, T>::Clear()
  {
    fct_rho.Clear();
  }
  

  /************************
   * VectorPhysicalIndice *
   ************************/


  //! Returns the memory used by the object in bytes
  template<class Dimension, int m, class T>
  size_t VectorPhysicalIndice<Dimension, m, T>::GetMemorySize() const
  {
    size_t taille = m*sizeof(T) + sizeof(int);
    for (int i = 0; i < m; i++)
      taille += fct_rho(i).GetMemorySize();
    
    return taille;
  }  
  
  
  //! returns true if a component of the constant index is complex
  template<class Dimension, int m, class T>
  bool VectorPhysicalIndice<Dimension, m, T>::IsComplex() const
  {
    if (IsVarying())
      {
	for (int i = 0; i < m; i++)
	  if (fct_rho(i).IsComplex())
	    return true;
	
	return false;
      }
    
    for (int i = 0; i < m; i++)
      if (abs(imagpart(cte_rho(i))) > threshold)
	return true;
    
    return false;
  }
  
  
  //! returns true if a component of the physical index is varying
  template<class Dimension, int m, class T>
  bool VectorPhysicalIndice<Dimension, m, T>::IsVarying() const
  {
    if (type_media == CONSTANT)
      return false;
    
    for (int i = 0; i < m; i++)
      if (fct_rho(i).IsVarying())
	return true;
    
    return false;
  }
  
  
  //! returns true if the index is null (all components equal to 0)
  template<class Dimension, int m, class T>
  bool VectorPhysicalIndice<Dimension, m, T>::IsZero() const
  {
    if (type_media == CONSTANT)
      return cte_rho.IsZero();
    
    for (int i = 0; i < m; i++)
      if (!fct_rho(i).IsZero())
	return false;
    
    return true;
  }
  
  
  //! sets all components of the vector to a given type of varying media
  template<class Dimension, int m, class T>
  void VectorPhysicalIndice<Dimension, m, T>::SetMediaType(int type_)
  {
    for (int i = 0; i < m; i++)
      fct_rho(i).SetMediaType(type_);
    
    if (type_ > 0)
      type_media = VARYING;
    else
      type_media = CONSTANT;
  }
  
  
  //! sets a component of the vector to a given type of varying media
  template<class Dimension, int m, class T>
  void VectorPhysicalIndice<Dimension, m, T>::SetMediaType(int i, int type_)
  {
    fct_rho(i).SetMediaType(type_);
    
    if (type_ > 0)
      type_media = VARYING;
    else
      {
        type_media = CONSTANT;
        for (int j = 0; j < m; j++)
          if (fct_rho(i).GetMediaType() > 0)
            type_media = VARYING;
      }
  }
  
    
  //! returns the number of varying coefficients (between 0 and m)
  template<class Dimension, int m, class T>
  int VectorPhysicalIndice<Dimension, m, T>::GetNbVaryingMedia() const
  {
    int nb = 0;
    if (type_media == VARYING)
      for (int i = 0; i < m; i++)
        if (fct_rho(i).IsVarying())
          nb++;
    
    return nb;
  }
  
  
  //! provides pointers on varying coefficients and associated component numbers
  template<class Dimension, int m, class T>
  void VectorPhysicalIndice<Dimension, m, T>::
  GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho,
                  IVect& num_component)
  {
    if (type_media == VARYING)
      for (int i = 0; i < m; i++)
        if (fct_rho(i).IsVarying())
          {
            rho(nb) = &fct_rho(i);
            num_component(nb) = i;
            nb++;
          }
  }
  
  
  //! clears varying coefficients
  template<class Dimension, int m, class T>
  void VectorPhysicalIndice<Dimension, m, T>::Clear()
  {
    for (int i = 0; i < m; i++)
      fct_rho(i).Clear();
  }


  //! provides coefficient and its gradient
  template<class Dimension, int m, class T> 
  void VectorPhysicalIndice<Dimension, m, T>::
  GetCoefGradient(const VarPhysicalProblem& var, int i, int j, TinyVector<T, m>& coef,
		  TinyMatrix<T, General, m, Dimension::dim_N>& grad) const
  {
    if (type_media != CONSTANT)
      {
	int i1 = var.ElementRho(i);
	grad.Zero();
	for (int k = 0; k < m; k++)
	  {
	    if (fct_rho(k).IsVarying())
	      {
		coef(k) = fct_rho(k).GetCoefficient(i1, j);
		SetRow(fct_rho(k).GetCoefGradient(i1, j), k, grad);
	      }
	    else
	      coef(k) = fct_rho(k).GetConstant();
	  }
      }
    else
      {
	coef = cte_rho;
	grad.Zero();
      }
  }
  
  
  //! retrieves coefficient, gradient and its hessian
  template<class Dimension, int m, class T> 
  void VectorPhysicalIndice<Dimension, m, T>
  ::GetCoefHessian(const VarPhysicalProblem& var, int i, int j,
		   TinyVector<T, m>& coef, TinyMatrix<T, General, m, Dimension::dim_N>& grad,
		   TinyVector<TinyMatrix<T, Symmetric,
		   Dimension::dim_N, Dimension::dim_N>, m>& hess)
  {
    if (type_media != CONSTANT)
      {
	int i1 = var.ElementRho(i);
	grad.Zero(); 
	for (int k = 0; k < m; k++)
	  {
	    hess(k).Zero();
	    if (fct_rho(k).IsVarying())
	      {
		coef(k) = fct_rho(k).GetCoefficient(i1, j);
		SetRow(fct_rho(k).GetCoefGradient(i1, j), k, grad);
		hess(k) = fct_rho(k).GetCoefHessian(i1, j);
	      }
	    else
	      coef(k) = fct_rho(k).GetConstant();
	  }
      }
    else
      {
	coef = cte_rho;
	grad.Zero();
	hess.Zero();
      }    
  }
  
  
  //! reads the physical index in a data file
  /*!
    \param[in] var given problem
    \param[in] nb subscript to access array parameters
    \param[in] parameters values retrieved from data file
    \param[in] keyword keyword of the line of data file
   */
  template<class Dimension, int m, class T>
  void VectorPhysicalIndice<Dimension, m, T>::
  SetInputData(int& nb, const VectString& parameters,
	       const string& keyword)
  {
    type_media = CONSTANT;
    for (int k = 0; k < m; k++)
      {
	this->SetInputVaryingMedia(nb, fct_rho(k), cte_rho(k), parameters);
	if (fct_rho(k).IsVarying())
	  type_media = VARYING;
      }
  }
    
  
  /************************
   * TensorPhysicalIndice *
   ************************/
  

  //! Returns the memory used by the object in bytes
  template<class Dimension, int m, class T>
  size_t TensorPhysicalIndice<Dimension, m, T>::GetMemorySize() const
  {
    size_t taille = cte_mu.GetSize()*sizeof(T) + 2*sizeof(int);
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        taille += fct_mu(i, j).GetMemorySize();
    
    return taille;
  }  
  
  
  //! sets all components of tensor to a given varying media
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::SetMediaType(int type_)
  {
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
	fct_mu(i, j).SetMediaType(type_);
    
    if (type_ > 0)
      type_media = VARYING;
    else
      type_media = CONSTANT;
  }
  
  
  //! sets a single component of tensor to a given varying media
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::
  SetMediaType(int i, int j, int type_)
  {
    fct_mu(i, j).SetMediaType(type_);
    
    type_media = CONSTANT;
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        if (fct_mu(i, j).GetMediaType() > 0)
          type_media = VARYING;
  }
  

  //! computes v = C u where C is the current tensor index
  template<class Dimension, int m, class T>
  template<int p, class T0, class T1>
  void TensorPhysicalIndice<Dimension, m, T>::
  MltMatrix(const VarPhysicalProblem& var, int i, int j,
            const TinyVector<T0, p>& u, TinyVector<T1, p>& v) const
  {
    if (type_media == CONSTANT)
      Seldon::Mlt(cte_mu, u, v);
    else
      {
        int i1 = var.ElementRho(i);
        if (type_anisotropy == ISOTROPE)
          {
            v = u;
            v *= fct_mu(0, 0).GetCoefficient(i1, j);
          }
        else if (type_anisotropy == ORTHOTROPE)
          {
            v = u;
            for (int k = 0; k < m; k++)
              v(k) *= fct_mu(k, k).GetCoefficient(i1, j);
          }
        else
          {
            v.Fill(0);
            for (int p1 = 0; p1 < m; p1++)
	      for (int q = 0; q < m; q++)
                v(p1) += fct_mu(p1, q).GetCoefficient(i1, j)*u(q);
          }
      }
  }

  
  //! computes v = C u where C is the current tensor index
  template<class Dimension, int m, class T> template<class T0, class T1>
  void TensorPhysicalIndice<Dimension, m, T>::
  MltMatrix(const VarPhysicalProblem& var, int i, int j,
            const TinyVector<T0, 2>& u, TinyVector<T1, 2>& v) const
  {
    // 2-D case
    if (type_media == CONSTANT)
      Seldon::Mlt(cte_mu, u, v);
    else
      {
        int i1 = var.ElementRho(i);
        if (type_anisotropy == ISOTROPE)
          {
            v = u;
            v *= fct_mu(0, 0).GetCoefficient(i1, j);
          }
        else if (type_anisotropy == ORTHOTROPE)
          {
            v = u;
	    v(0) *= fct_mu(0, 0).GetCoefficient(i1, j);
	    v(1) *= fct_mu(1, 1).GetCoefficient(i1, j);
          }
        else
          {
	    v(0) = fct_mu(0, 0).GetCoefficient(i1, j)*u(0)
	      + fct_mu(0, 1).GetCoefficient(i1, j)*u(1);

	    v(1) = fct_mu(0, 1).GetCoefficient(i1, j)*u(0)
	      + fct_mu(1, 1).GetCoefficient(i1, j)*u(1);
          }
      }
  }


  //! computes v = C u where C is the current tensor index
  template<class Dimension, int m, class T> template<class T0, class T1>
  void TensorPhysicalIndice<Dimension, m, T>::
  MltMatrix(const VarPhysicalProblem& var, int i, int j,
            const TinyVector<T0, 3>& u, TinyVector<T1, 3>& v) const
  {
    // 3-D case
    if (type_media == CONSTANT)
      Seldon::Mlt(cte_mu, u, v);
    else
      {
        int i1 = var.ElementRho(i);
        if (type_anisotropy == ISOTROPE)
          {
            v = u;
            v *= fct_mu(0, 0).GetCoefficient(i1, j);
          }
        else if (type_anisotropy == ORTHOTROPE)
          {
            v = u;
	    v(0) *= fct_mu(0, 0).GetCoefficient(i1, j);
	    v(1) *= fct_mu(1, 1).GetCoefficient(i1, j);
	    v(2) *= fct_mu(2, 2).GetCoefficient(i1, j);
          }
        else
          {
	    v(0) = fct_mu(0, 0).GetCoefficient(i1, j)*u(0)
	      + fct_mu(0, 1).GetCoefficient(i1, j)*u(1)
	      + fct_mu(0, 2).GetCoefficient(i1, j)*u(2);

	    v(1) = fct_mu(1, 0).GetCoefficient(i1, j)*u(0)
	      + fct_mu(1, 1).GetCoefficient(i1, j)*u(1)
	      + fct_mu(1, 2).GetCoefficient(i1, j)*u(2);

	    v(2) = fct_mu(2, 0).GetCoefficient(i1, j)*u(0)
	      + fct_mu(2, 1).GetCoefficient(i1, j)*u(1)
	      + fct_mu(2, 2).GetCoefficient(i1, j)*u(2);
          }
      }
  }
      
        
  //! returns the number of varying coefficients (between 0 and m*m)
  template<class Dimension, int m, class T>
  int TensorPhysicalIndice<Dimension, m, T>::GetNbVaryingMedia() const
  {
    int nb = 0;
    if (type_media == VARYING)
      for (int i = 0; i < m; i++)
        for (int j = i; j < m; j++)
          if (fct_mu(i, j).IsVarying())
            nb++;
    
    return nb;
  }
  
  
  //! provides pointers on varying coefficients and associated component numbers
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::
  GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho,
                  IVect& num_component)
  {
    if (type_media == VARYING)
      {
        int num = 0;
        for (int i = 0; i < m; i++)
          for (int j = i; j < m; j++)
            {
              if (fct_mu(i, j).IsVarying())
                {
                  rho(nb) = &fct_mu(i, j);
                  num_component(nb) = num;
                  nb++;
                }
              num++;
            }
      }
  }
  

  //! sets the physical index to alpha I
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::SetDiagonal(const T& c)
  {
    type_anisotropy = ISOTROPE;
    type_media = CONSTANT;
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        fct_mu(i, j).SetMediaType(fct_mu(i, j).CONSTANT);
    
    cte_mu.SetDiagonal(c);
  }
  
  
  //! sets physical index to a constant matrix
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::
  SetConstant(const TinyMatrix<T, Symmetric, m, m>& coef)
  {
    type_anisotropy = ISOTROPE;
    type_media = CONSTANT;
    cte_mu = coef;
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        fct_mu(i, j).SetMediaType(fct_mu(i, j).CONSTANT);
    
    for (int i = 0; i < m; i++)
      for (int j = i+1; j < m; j++)
	if (abs(coef(i,j)) > threshold)
	  type_anisotropy = ANISOTROPE;
    
    if (type_anisotropy != ANISOTROPE)
      {
	for (int i = 1; i < m; i++)
	  if (abs(coef(i,i)-coef(0,0)) > threshold)
	    type_anisotropy = ORTHOTROPE;
      }
  }
  
    
  //! clears varying coefficients
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::Clear()
  {
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
	fct_mu(i, j).Clear();
  }
  

  //! sets physical index to null matrix
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::Zero()
  {
    type_anisotropy = ISOTROPE;
    type_media = CONSTANT;
    cte_mu.Fill(0);
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        fct_mu(i, j).SetMediaType(fct_mu(i, j).CONSTANT);
    
  }
  
  
  //! matrix is inversed
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::GetInverse()
  {
    // DISP(cte_mu);
    TinyMatrix<T, Symmetric, m, m> mu;
    if (!cte_mu.IsZero())
      {
        Seldon::GetInverse(cte_mu, mu);
        cte_mu = mu;
      }
    
    if (type_media == VARYING)
      {
	if (type_anisotropy == ISOTROPE)
	  fct_mu(0,0).GetInverse();
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    for (int k = 0; k < fct_mu.GetM(); k++)
	      fct_mu(k, k).GetInverse();
	  }
	else
	  {
	    // anisotropic case, only coefficients can be inverted
	    
	    // finding a variable coefficient
	    int nb_elt = 0, i0 = -1, j0 = -1;
            for (int i = 0; i < m; i++)
              for (int j = i; j < m; j++)
                {
		  if (fct_mu(i, j).GetValue().GetM() > 0)
		    {
		      nb_elt = fct_mu(i, j).GetValue().GetM();
		      i0 = i; j0 = j;
		    }
		}
	    
	    // constant coefficients are transformed into variable ones
	    Vector<Vector<T> >& mu00 = fct_mu(i0, j0).GetValue();
	    for (int i = 0; i < m; i++)
              for (int j = i; j < m; j++)
                {
		  if (fct_mu(i, j).GetMediaType() == fct_mu(i, j).CONSTANT)
		    {
		      T cte = fct_mu(i, j).GetConstant();
		      fct_mu(i, j).SetMediaType(fct_mu(i, j).USER);
		      fct_mu(i, j).Reallocate(nb_elt, false, false);
		      Vector<Vector<T> >& mu_ij = fct_mu(i, j).GetValue();		      
		      for (int k = 0; k < nb_elt; k++)
			{
			  fct_mu(i, j).ReallocateOnElement(k, mu00(k).GetM(), false, false);
			  mu_ij(k).Fill(cte);
			}
		    }
		  else
		    {
		      fct_mu(i, j).SetMediaType(fct_mu(i, j).USER);		  
		      fct_mu(i, j).ClearHessianGradient();
		    }
                }
            
	    // achieving the inversion of the index on quadrature points
            for (int iquad = 0; iquad < mu00.GetM(); iquad++)
              for (int jloc = 0; jloc < mu00(iquad).GetM(); jloc++)
                {
                  for (int i = 0; i < m; i++)
                    for (int j = i; j < m; j++)
                      mu(i, j) = fct_mu(i, j).GetCoefficient(iquad, jloc);
                  
                  Seldon::GetInverse(mu);
                  for (int i = 0; i < m; i++)
                    for (int j = i; j < m; j++)
                      fct_mu(i, j).SetCoefficient(iquad, jloc, mu(i, j));
                }
	  }
      }
    // DISP(cte_mu);
  }
  
  
  //! multiplication of physical index by a coefficient
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::Mlt(const T& coef)
  {
    // DISP(cte_mu);
    Seldon::Mlt(coef, cte_mu);
    
    if (type_media == VARYING)
      {
	if (type_anisotropy == ISOTROPE)
	  fct_mu(0,0).Mlt(coef);
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    for (int k = 0; k < m; k++)
	      fct_mu(k,k).Mlt(coef);
	  }
	else
	  {
            for (int i = 0; i < m; i++)
              for (int j = i; j < m; j++)
                fct_mu(i, j).Mlt(coef);
	  }
      }
    // DISP(cte_mu);
  }
  
  
  //! returns true if an entry of the constant index is complex
  template<class Dimension, int m, class T>
  bool TensorPhysicalIndice<Dimension, m, T>::IsComplex() const
  {
    if (type_media == CONSTANT)
      {
	for (int i = 0; i < m; i++)
	  for (int j = i; j < m; j++)
	    if (abs(imagpart(cte_mu(i,j))) > threshold)
	      return true;
      }
    else
      {
	if (type_anisotropy == ISOTROPE)
	  return fct_mu(0,0).IsComplex();
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    for (int i = 0; i < m; i++)
	      if (fct_mu(i,i).IsComplex())
		return true;
	  }
	else
	  {
	    for (int i = 0; i < m; i++)
	      for (int j = i; j < m; j++)
		if ( fct_mu(i,j).IsComplex())
		  return true;
	  }
      }
    
    return false;
  }
  
  
  //! returns true if an entry of the physical index is varying
  template<class Dimension, int m, class T>
  bool TensorPhysicalIndice<Dimension, m, T>::IsVarying() const
  {
    if (type_media == CONSTANT)
      return false;
    else
      {
	if (type_anisotropy == ISOTROPE)
	  return fct_mu(0,0).IsVarying();
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    for (int i = 0; i < m; i++)
	      if (fct_mu(i,i).IsVarying())
		return true;
	  }
	else
	  {
	    for (int i = 0; i < m; i++)
	      for (int j = i; j < m; j++)
		if ( fct_mu(i,j).IsVarying())
		  return true;
	  }
      }
    
    return false;
  }    
  
  
  //! returns true if an entry of the physical index is varying
  template<class Dimension, int m, class T>
  bool TensorPhysicalIndice<Dimension, m, T>::IsZero() const
  {
    if (type_media == CONSTANT)
      return cte_mu.IsZero();
    else
      {
	if (type_anisotropy == ISOTROPE)
	  return fct_mu(0,0).IsZero();
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    for (int i = 0; i < m; i++)
	      if (!fct_mu(i,i).IsZero())
		return false;
	  }
	else
	  {
	    for (int i = 0; i < m; i++)
	      for (int j = i; j < m; j++)
		if ( !fct_mu(i,j).IsZero())
		  return false;
	  }
      }
    
    return true;
  }    
  
  
  //! computes v = C u where C is the current tensor index
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::
  MltMatrixPoint(int i1, const VectReal_wp& phi,
		 const typename Dimension::R_N& pt_glob,
		 const TinyVector<T, m>& u, TinyVector<T, m>& v) const
  {
    if (type_media == CONSTANT)
      Seldon::Mlt(cte_mu, u, v);
    else
      {
        if (type_anisotropy == ISOTROPE)
          {
            v = u;
            v *= fct_mu(0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          }
        else if (type_anisotropy == ORTHOTROPE)
          {
            v = u;
            for (int k = 0; k < m; k++)
              v(k) *= fct_mu(k, k).GetIndexAtPoint(i1, phi, pt_glob);
          }
        else
          {
            v.Fill(0);
            for (int p = 0; p < m; p++)
	      for (int q = p; q < m; q++)
                {
		  T coef = fct_mu(p, q).GetIndexAtPoint(i1, phi, pt_glob);
		  v(p) += coef*u(q);
		  if (q != p)
		    v(q) += coef*u(p);
		}
          }
      }
  }

  
  //! returns interpolation of tensor C at point pt_glob
  template<class Dimension, int m, class T> TinyMatrix<T, Symmetric, m, m>
  TensorPhysicalIndice<Dimension, m, T>::GetMatrixPoint(int i1, const VectReal_wp& phi, 
							const typename Dimension::R_N& pt_glob) const
  {
    if (type_media == CONSTANT)
      return cte_mu;
    else
      {
	TinyMatrix<T, Symmetric, m, m> C;
        if (type_anisotropy == ISOTROPE)
	  {
	    C.SetDiagonal(fct_mu(0, 0).GetIndexAtPoint(i1, phi, pt_glob));
	  }
        else if (type_anisotropy == ORTHOTROPE)
          {
            for (int k = 0; k < m; k++)
	      C(k, k) = fct_mu(k, k).GetIndexAtPoint(i1, phi, pt_glob);
	  }
	else
	  {
	    for (int p = 0; p < m; p++)
	      for (int q = p; q < m; q++)
		C(p, q) = fct_mu(p, q).GetIndexAtPoint(i1, phi, pt_glob);
	  }
	
	return C;
      }
  }

    
  //! reads the physical index in a data file
  /*!
    \param[in] var given problem
    \param[in] nb subscript to access array parameters
    \param[in] parameters values retrieved from data file
    \param[in] keyword keyword of the line of data file
   */
  template<class Dimension, int m, class T>
  void TensorPhysicalIndice<Dimension, m, T>::
  SetInputData(int& nb,
	       const VectString& parameters, const string& keyword)
  {
    T zero; SetComplexZero(zero);
 
    if (!keyword.compare("ISOTROPE"))
      {
	type_anisotropy = ISOTROPE; cte_mu.Fill(zero);
	this->SetInputVaryingMedia(nb, fct_mu(0,0), cte_mu(0,0), parameters);
	type_media = CONSTANT;
	if (fct_mu(0,0).IsVarying())
	  type_media = VARYING;
	
	for (int i = 1; i < cte_mu.GetM(); i++)
	  cte_mu(i,i) = cte_mu(0,0);
	
      }
    else if (!keyword.compare("ORTHOTROPE"))
      {
	type_anisotropy = ORTHOTROPE; cte_mu.Fill(zero);
	type_media = CONSTANT;
	for (int i = 0; i < cte_mu.GetM(); i++)
	  {
	    this->SetInputVaryingMedia(nb, fct_mu(i,i), cte_mu(i,i), parameters);
	    if (fct_mu(i,i).IsVarying())
	      type_media = VARYING;
	  }
      }
    else if (!keyword.compare("ANISOTROPE"))
      {
	type_anisotropy = ANISOTROPE; cte_mu.Fill(zero);
	type_media = CONSTANT;
	for (int i = 0; i < cte_mu.GetM(); i++)
	  for (int j = i; j < cte_mu.GetM(); j++)
	    {
	      this->SetInputVaryingMedia(nb, fct_mu(i,j), cte_mu(i,j), parameters);
	      if (fct_mu(i,j).IsVarying())
		type_media = VARYING;
	    }
      }
  }
  

  //! returns physical index of quadrature point (i,j)
  template<class Dimension, int m, class T> 
  TinyMatrix<T, Symmetric, m, m>
  TensorPhysicalIndice<Dimension, m, T>::GetCoefficient(const VarPhysicalProblem& var,
                                                        int i, int j) const
  {
    if (type_media == CONSTANT)
      return cte_mu;
    else
      {
	int i1 = var.ElementRho(i);
        TinyMatrix<T, Symmetric, m, m> mu;
	if (type_anisotropy == ISOTROPE)
	  {
	    T coef = fct_mu(0,0).GetCoefficient(i1, j);
	    mu.SetDiagonal(coef);
	  }
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    for (int k = 0; k < m; k++)
	      mu(k, k) = fct_mu(k,k).GetCoefficient(i1, j);
	  }
	else
	  {
	    for (int p = 0; p < m; p++)
	      for (int q = p; q < m; q++)
		mu(p, q) = fct_mu(p, q).GetCoefficient(i1, j);
	  }
        
	return mu;
      }
    
    return cte_mu;
  }
  
  
  //! returns coefficient mu(i, j) and its gradient
  template<class Dimension, int m, class T> 
  void TensorPhysicalIndice<Dimension, m, T>::
  GetCoefGradient(const VarPhysicalProblem& var, int i, int j,
                  TinyMatrix<T, Symmetric, m, m>& coef,
		  TinyArray3D<T, m, m, Dimension::dim_N>& grad) const
  {
    if (type_media != CONSTANT)
      {
	abort();
      }
    else
      {
	coef = cte_mu;
	grad.Zero();
      }
  }
  
  
  //! returns entry (m,n) of physical index of quadrature point (i,j)
  template<class Dimension, int m, class T> 
  T TensorPhysicalIndice<Dimension, m, T>::
  GetCoefficient(const VarPhysicalProblem& var, int i, int j, int p, int q) const
  {
    if (type_media == CONSTANT)
      return cte_mu(p, q);
    else
      {
	int i1 = var.ElementRho(i);
	if (type_anisotropy == ISOTROPE)
	  return fct_mu(0,0).GetCoefficient(i1,j);
	else if (type_anisotropy == ORTHOTROPE)
	  {
	    if (p == q)
	      return fct_mu(p, p).GetCoefficient(i1, j);
	    else
	      return T(0);
	  }
	else
	  return fct_mu(p, q).GetCoefficient(i1, j);
      }
    
    return cte_mu(p, q);
  }

  
  //! displays informations about class TensorPhysicalIndice
  template<class Dimension, int m, class T>
  ostream& operator <<(ostream& out, const TensorPhysicalIndice<Dimension, m, T>& var)
  {
    if (var.type_media == var.CONSTANT)
      out<<var.cte_mu;
    else
      out<<"varying media";
    
    return out;
  }
  
  
  /*************************
   * ElasticPhysicalIndice *
   *************************/
  

  //! Returns the memory used by the object in bytes
  template<class Dimension, int m, class T>
  size_t ElasticPhysicalIndice<Dimension, m, T>::GetMemorySize() const
  {
    size_t taille = cte_C.GetSize()*sizeof(T) + 2*sizeof(int);
    for (int i = 0; i < nb_components; i++)
      for (int j = i; j < nb_components; j++)
        taille += fct_C(i, j).GetMemorySize();
    
    return taille;
  }  
  
    
  //! returns the number of varying scalar indexes
  template<class Dimension, int m, class T>
  int ElasticPhysicalIndice<Dimension, m, T>::GetNbVaryingMedia() const
  {
    int nb = 0;
    if (type_media == VARYING)
      for (int i = 0; i < nb_components; i++)
        for (int j = i; j < nb_components; j++)
          if (fct_C(i, j).IsVarying())
            nb++;
    
    return nb;
  }
  
  
  //! fills the arrays containing varying indices
  template<class Dimension, int m, class T>
  void ElasticPhysicalIndice<Dimension, m, T>::
  GetVaryingMedia(int& nb, Vector<PhysicalVaryingMedia<Dimension, T>* >& rho, IVect& num_component)
  {
    if (type_media == VARYING)
      {
        int num = 0;
        for (int i = 0; i < nb_components; i++)
          for (int j = i; j < nb_components; j++)
            {
              if (fct_C(i, j).IsVarying())
                {
                  rho(nb) = &fct_C(i, j);
                  num_component(nb) = num;
                  nb++;
                }
              
              num++;
            }
      }
  }
  

  //! returns physical index of quadrature point (num_elem,jloc)
  template<class Dimension, int m, class T> template<class T0>
  void ElasticPhysicalIndice<Dimension, m, T>::
  GetCoefficient(const VarPhysicalProblem& var, int num_elem, int jloc,
		 TinyMatrix<TinyMatrix<T0, General, m, m>, General, m, m>& mat) const
  {
    FillZero(mat);
    
    if (type_media == VARYING)
      {	
        int i1 = var.ElementRho(num_elem);
	switch (type_anisotropy)
          {
          case ISOTROPE :
            {
              T lambda = fct_C(0,0).GetCoefficient(i1, jloc);
              T mu = fct_C(1,1).GetCoefficient(i1, jloc);
              for (int i = 0; i < m; i++)
                {
                  mat(i, i)(i, i) = lambda + T(2)*mu;
                  for (int j = 0; j < m; j++)
                    if (j != i)
                      {
                        mat(i, j)(i, j) = lambda;
                        mat(i, i)(j, j) = mu;
                        mat(i, j)(j, i) = mu;
                      }
                }
            }
            break;
          case ORTHOTROPE :
            {
              for (int i = 0; i < m; i++)
                for (int j = i; j < m; j++)
                  {
                    mat(i, j)(i, j) = fct_C(i, i, j, j).GetCoefficient(i1, jloc);
                    mat(j, i)(j, i) = mat(i, j)(i, j);
                  }
              
              for (int i = 0; i < m; i++)
		for (int j = (i+1); j < m; j++)
                  {
                    mat(i, i)(j, j) = fct_C(i, j, i, j).GetCoefficient(i1, jloc);
                    mat(j, j)(i, i) = mat(i, i)(j, j);
                    mat(j, i)(i, j) = mat(i, i)(j, j);
                    mat(i, j)(j, i) = mat(i, i)(j, j);
                  }
              }
            break;
          case ANISOTROPE :
            {
              for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                  for (int k = 0; k < m; k++)
                    for (int l = 0; l < m; l++)
                      mat(i, k)(j, l) = fct_C(i, j, k, l).GetCoefficient(i1, jloc);
            }
            break;
          }
      }
    else
      {                                         
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < m; j++)
	    for (int k = 0; k < m; k++)
	      for (int l = 0; l < m; l++)
		mat(i, k)(j, l) = cte_C(i, j, k, l);	
      }
  }
   

  //! sets directly tensor C
  template<class Dimension, int m, class T>
  void ElasticPhysicalIndice<Dimension, m, T>::
  SetConstant(const TinySymmetricTensor<T, m>& C)
  {
    type_media = CONSTANT;
    cte_C = C;
    
    // detecting type of anisotropy
    T zero; SetComplexZero(zero);
    type_anisotropy = ORTHOTROPE;
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        for (int k = 0; k < m; k++)
          for (int l = k+1; l < m; l++)
            if ((i != k) || (j != l))
              if (cte_C(i, j, k, l) != zero)
                type_anisotropy = ANISOTROPE;
  }
  
  
  //! constant tensor is inversed
  template<class Dimension, int m, class T>
  void ElasticPhysicalIndice<Dimension, m, T>::GetInverse()
  {    
    cte_C.GetInverse();
    
    if (type_media == VARYING)
      {
        switch (type_anisotropy)
          {
          case ISOTROPE:
            {
              Vector<Vector<T> >& Lambda = fct_C(0, 0).GetValue();
              Vector<Vector<T> >& Mu = fct_C(1, 1).GetValue();
              T L, mu;
              for (int i = 0; i < Lambda.GetM(); i++)
                for (int j = 0; j < Lambda(i).GetM(); j++)
                  {
                    L = Lambda(i)(j); mu = Mu(i)(j);
                    Lambda(i)(j) = -L/(4.0*mu*(L+mu));
                    Mu(i)(j) = 1.0/(4.0*mu);
                  }
              
              fct_C(0, 0).SetMediaType(fct_C(0, 0).USER);
              fct_C(0, 0).ClearHessianGradient();
              fct_C(1, 1).SetMediaType(fct_C(1, 1).USER);
              fct_C(1, 1).ClearHessianGradient();
            }
            break;
          case ORTHOTROPE:
          case ANISOTROPE :
            {
              for (int i = 0; i < nb_components; i++)
                for (int j = i; j < nb_components; j++)
                  {
                    fct_C(i, j).ClearHessianGradient();
                    fct_C(i, j).SetMediaType(fct_C(0, 0).USER);
                  }
              
              TinySymmetricTensor<T, m> C;
              Vector<Vector<T> >& coef = fct_C(0, 0).GetValue();
              for (int iquad = 0; iquad < coef.GetM(); iquad++)
                for (int jloc = 0; jloc < coef(iquad).GetM(); jloc++)
                  {
                    for (int i = 0; i < nb_components; i++)
                      for (int j = i; j < nb_components; j++)
                        C(i, j) = fct_C(i, j).GetCoefficient(iquad, jloc);
                    
                    C.GetInverse();
                    
                    for (int i = 0; i < nb_components; i++)
                      for (int j = i; j < nb_components; j++)
                        fct_C(i, j).SetCoefficient(iquad, jloc, C(i, j));
                  }
            }
            break;
          }
      }
  }
  

  //! multiplication of physical index by a coefficient
  template<class Dimension, int m, class T>
  void ElasticPhysicalIndice<Dimension, m, T>::Mlt(const T& coef)
  {
    // DISP(cte_mu);
    Seldon::Mlt(coef, cte_C);
    
    if (type_media == VARYING)
      {
        TinyMatrix<int, Symmetric, m, m> indice_doublet;
        indice_doublet.Fill();
        
	if (type_anisotropy == ISOTROPE)
          {
            fct_C(0,0).Mlt(coef);
            fct_C(1,1).Mlt(coef);
          }
	else if (type_anisotropy == ORTHOTROPE)
	  {
            for (int i = 0; i < m; i++)
              for (int j = i; j < m; j++)
                {
                  int num1 = indice_doublet(i,i), num2 = indice_doublet(j,j);
                  fct_C(num1, num2).Mlt(coef);
                }
            
            for (int i = 0; i < m; i++)
              for (int j = (i+1); j < m; j++)
                {
                  int num = indice_doublet(i,j);
                  fct_C(num, num).Mlt(coef);
                }
	  }
	else
	  {
            for (int i = 0; i < nb_components; i++)
              for (int j = i; j < nb_components; j++)
                fct_C(i, j).Mlt(coef);
	  }
      }
  }

    
  template<class Dimension, int m, class T> template<class Vector1>
  void ElasticPhysicalIndice<Dimension, m, T>::
  MltMatrixVar(int i1, int jloc, const TinyVector<int, 2>& xt,
               const TinyVector<Vector1, 2>& dU, TinyVector<Vector1, 2>& dV) const
  {
    switch (type_anisotropy)
      {
      case ISOTROPE :
        {
          T lambda = fct_C(0, 0).GetCoefficient(i1, jloc);
          T mu = fct_C(1, 1).GetCoefficient(i1, jloc);
          dV(0)(0) = (lambda+2.0*mu)*dU(0)(0) + lambda*dU(1)(1);
          dV(1)(1) = (lambda+2.0*mu)*dU(1)(1) + lambda*dU(0)(0);
          dV(0)(1) = mu*(dU(0)(1) + dU(1)(0));
          dV(1)(0) = dV(0)(1);
          break;
        }
      case ORTHOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          dV(0)(0) = c0000*dU(0)(0) + c0011*dU(1)(1);
          dV(1)(1) = c0011*dU(0)(0) + c1111*dU(1)(1);
          dV(0)(1) = c0101*(dU(0)(1)+dU(1)(0));
          dV(1)(0) = dV(0)(1);
        }
        break;
      case ANISOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c0001 = fct_C(0, 0, 0, 1).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c1110 = fct_C(1, 1, 1, 0).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          dV(0)(0) = c0000*dU(0)(0) + c0011*dU(1)(1) + c0001*(dU(0)(1)+dU(1)(0));
          dV(1)(1) = c0011*dU(0)(0) + c1111*dU(1)(1) + c1110*(dU(0)(1)+dU(1)(0));
          dV(0)(1) = c0001*dU(0)(0) + c1110*dU(1)(1) + c0101*(dU(0)(1)+dU(1)(0));
          dV(1)(0) = dV(0)(1);
        }
        break;
      }
  }


  template<class Dimension, int m, class T> template<class Vector1>
  void ElasticPhysicalIndice<Dimension, m, T>::
  MltMatrixVar(int i1, int jloc, const TinyVector<int, 3>& xt,
               const TinyVector<Vector1, 3>& dU, TinyVector<Vector1, 3>& dV) const
  {
    switch (type_anisotropy)
      {
      case ISOTROPE :
        {
          T lambda = fct_C(0, 0).GetCoefficient(i1, jloc);
          T mu = fct_C(1, 1).GetCoefficient(i1, jloc);
          dV(0)(0) = (lambda+2.0*mu)*dU(0)(0) + lambda*dU(1)(1) + lambda*dU(2)(2);
          dV(1)(1) = (lambda+2.0*mu)*dU(1)(1) + lambda*dU(0)(0) + lambda*dU(2)(2);
          dV(2)(2) = (lambda+2.0*mu)*dU(2)(2) + lambda*dU(0)(0) + lambda*dU(1)(1);
          dV(0)(1) = mu*(dU(0)(1) + dU(1)(0));
          dV(1)(0) = dV(0)(1);
          dV(0)(2) = mu*(dU(0)(2) + dU(2)(0));
          dV(2)(0) = dV(0)(2);
          dV(1)(2) = mu*(dU(1)(2) + dU(2)(1));
          dV(2)(1) = dV(1)(2);
        }
        break;
      case ORTHOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c0022 = fct_C(0, 0, 2, 2).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c1122 = fct_C(1, 1, 2, 2).GetCoefficient(i1, jloc);
          T c2222 = fct_C(2, 2, 2, 2).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          T c0202 = fct_C(0, 2, 0, 2).GetCoefficient(i1, jloc);
          T c1212 = fct_C(1, 2, 1, 2).GetCoefficient(i1, jloc);
          dV(0)(0) = c0000*dU(0)(0) + c0011*dU(1)(1) + c0022*dU(2)(2);
          dV(1)(1) = c0011*dU(0)(0) + c1111*dU(1)(1) + c1122*dU(2)(2);
          dV(2)(2) = c0022*dU(0)(0) + c1122*dU(1)(1) + c2222*dU(2)(2);
          dV(0)(1) = c0101*(dU(0)(1)+dU(1)(0));
          dV(1)(0) = dV(0)(1);
          dV(0)(2) = c0202*(dU(0)(2)+dU(2)(0));
          dV(2)(0) = dV(0)(2);
          dV(1)(2) = c1212*(dU(1)(2)+dU(2)(1));
          dV(2)(1) = dV(1)(2);
        }
        break;
      case ANISOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c0022 = fct_C(0, 0, 2, 2).GetCoefficient(i1, jloc);
          T c0001 = fct_C(0, 0, 0, 1).GetCoefficient(i1, jloc);
          T c0002 = fct_C(0, 0, 0, 2).GetCoefficient(i1, jloc);
          T c0012 = fct_C(0, 0, 1, 2).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c1122 = fct_C(1, 1, 2, 2).GetCoefficient(i1, jloc);
          T c1101 = fct_C(1, 1, 0, 1).GetCoefficient(i1, jloc);
          T c1102 = fct_C(1, 1, 0, 2).GetCoefficient(i1, jloc);
          T c1112 = fct_C(1, 1, 1, 2).GetCoefficient(i1, jloc);
          T c2222 = fct_C(2, 2, 2, 2).GetCoefficient(i1, jloc);
          T c2201 = fct_C(2, 2, 0, 1).GetCoefficient(i1, jloc);
          T c2202 = fct_C(2, 2, 0, 2).GetCoefficient(i1, jloc);
          T c2212 = fct_C(2, 2, 1, 2).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          T c0102 = fct_C(0, 1, 0, 2).GetCoefficient(i1, jloc);
          T c0112 = fct_C(0, 1, 1, 2).GetCoefficient(i1, jloc);
          T c0202 = fct_C(0, 2, 0, 2).GetCoefficient(i1, jloc);
          T c0212 = fct_C(0, 2, 1, 2).GetCoefficient(i1, jloc);
          T c1212 = fct_C(1, 2, 1, 2).GetCoefficient(i1, jloc);
          dV(0)(0) = c0000*dU(0)(0) + c0011*dU(1)(1) + c0022*dU(2)(2) + 
            c0001*(dU(0)(1)+dU(1)(0)) + c0002*(dU(0)(2)+dU(2)(0)) + c0012*(dU(1)(2)+dU(2)(1));

          dV(1)(1) = c0011*dU(0)(0) + c1111*dU(1)(1) + c1122*dU(2)(2) + 
            c1101*(dU(0)(1)+dU(1)(0)) + c1102*(dU(0)(2)+dU(2)(0)) + c1112*(dU(1)(2)+dU(2)(1));

          dV(2)(2) = c0022*dU(0)(0) + c1122*dU(1)(1) + c2222*dU(2)(2) + 
            c2201*(dU(0)(1)+dU(1)(0)) + c2202*(dU(0)(2)+dU(2)(0)) + c2212*(dU(1)(2)+dU(2)(1));

          dV(0)(1) = c0001*dU(0)(0) + c1101*dU(1)(1) + c2201*dU(2)(2) + 
            c0101*(dU(0)(1)+dU(1)(0)) + c0102*(dU(0)(2)+dU(2)(0)) + c0112*(dU(1)(2)+dU(2)(1));
          dV(1)(0) = dV(0)(1);
          
          dV(0)(2) = c0002*dU(0)(0) + c1102*dU(1)(1) + c2202*dU(2)(2) + 
            c0102*(dU(0)(1)+dU(1)(0)) + c0202*(dU(0)(2)+dU(2)(0)) + c0212*(dU(1)(2)+dU(2)(1));
          dV(2)(0) = dV(0)(2);
          
          dV(1)(2) = c0012*dU(0)(0) + c1112*dU(1)(1) + c2212*dU(2)(2) + 
            c0112*(dU(0)(1)+dU(1)(0)) + c0212*(dU(0)(2)+dU(2)(0)) + c1212*(dU(1)(2)+dU(2)(1));
          dV(2)(1) = dV(1)(2);
        }
        break;
      }
  }


  template<class Dimension, int m, class T> template<class Vector1>
  void ElasticPhysicalIndice<Dimension, m, T>::
  MltVectorVar(int i1, int jloc, const TinyVector<int, 2>& xt,
               const Vector1& dU, Vector1& dV) const
  {
    switch (type_anisotropy)
      {
      case ISOTROPE:
        {
          T lambda = fct_C(0, 0).GetCoefficient(i1, jloc);
          T mu = fct_C(1, 1).GetCoefficient(i1, jloc);
          dV(0) = (lambda+2.0*mu)*dU(0) + lambda*dU(3);
          dV(3) = (lambda+2.0*mu)*dU(3) + lambda*dU(0);
          dV(2) = mu*(dU(1) + dU(2));
          dV(1) = dV(2);
        }
        break;
      case ORTHOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          dV(0) = c0000*dU(0) + c0011*dU(3);
          dV(3) = c0011*dU(0) + c1111*dU(3);
          dV(1) = c0101*(dU(1)+dU(2));
          dV(2) = dV(1);
        }
        break;
      case ANISOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c0001 = fct_C(0, 0, 0, 1).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c1110 = fct_C(1, 1, 1, 0).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          dV(0) = c0000*dU(0) + c0011*dU(3) + c0001*(dU(1)+dU(2));
          dV(3) = c0011*dU(0) + c1111*dU(3) + c1110*(dU(1)+dU(2));
          dV(1) = c0001*dU(0) + c1110*dU(3) + c0101*(dU(1)+dU(2));
          dV(2) = dV(1);
        }
        break;
      }
  }


  template<class Dimension, int m, class T> template<class Vector1>
  void ElasticPhysicalIndice<Dimension, m, T>::
  MltVectorVar(int i1, int jloc, const TinyVector<int, 3>& xt,
               const Vector1& dU, Vector1& dV) const
  {
    switch (type_anisotropy)
      {
      case ISOTROPE:
        {
          T lambda = fct_C(0, 0).GetCoefficient(i1, jloc);
          T mu = fct_C(1, 1).GetCoefficient(i1, jloc);
          dV(0) = (lambda+2.0*mu)*dU(0) + lambda*dU(4) + lambda*dU(8);
          dV(4) = (lambda+2.0*mu)*dU(4) + lambda*dU(0) + lambda*dU(8);
          dV(8) = (lambda+2.0*mu)*dU(8) + lambda*dU(0) + lambda*dU(4);
          dV(1) = mu*(dU(1) + dU(3));
          dV(3) = dV(1);
          dV(2) = mu*(dU(2) + dU(6));
          dV(6) = dV(2);
          dV(5) = mu*(dU(5) + dU(7));
          dV(7) = dV(5);
        }
        break;
      case ORTHOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c0022 = fct_C(0, 0, 2, 2).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c1122 = fct_C(1, 1, 2, 2).GetCoefficient(i1, jloc);
          T c2222 = fct_C(2, 2, 2, 2).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          T c0202 = fct_C(0, 2, 0, 2).GetCoefficient(i1, jloc);
          T c1212 = fct_C(1, 2, 1, 2).GetCoefficient(i1, jloc);
          dV(0) = c0000*dU(0) + c0011*dU(4) + c0022*dU(8);
          dV(4) = c0011*dU(0) + c1111*dU(4) + c1122*dU(8);
          dV(8) = c0022*dU(0) + c1122*dU(4) + c2222*dU(8);
          dV(1) = c0101*(dU(1)+dU(3));
          dV(3) = dV(1);
          dV(2) = c0202*(dU(2)+dU(6));
          dV(6) = dV(2);
          dV(5) = c1212*(dU(5)+dU(7));
          dV(7) = dV(5);
        }
        break;
      case ANISOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetCoefficient(i1, jloc);
          T c0011 = fct_C(0, 0, 1, 1).GetCoefficient(i1, jloc);
          T c0022 = fct_C(0, 0, 2, 2).GetCoefficient(i1, jloc);
          T c0001 = fct_C(0, 0, 0, 1).GetCoefficient(i1, jloc);
          T c0002 = fct_C(0, 0, 0, 2).GetCoefficient(i1, jloc);
          T c0012 = fct_C(0, 0, 1, 2).GetCoefficient(i1, jloc);
          T c1111 = fct_C(1, 1, 1, 1).GetCoefficient(i1, jloc);
          T c1122 = fct_C(1, 1, 2, 2).GetCoefficient(i1, jloc);
          T c1101 = fct_C(1, 1, 0, 1).GetCoefficient(i1, jloc);
          T c1102 = fct_C(1, 1, 0, 2).GetCoefficient(i1, jloc);
          T c1112 = fct_C(1, 1, 1, 2).GetCoefficient(i1, jloc);
          T c2222 = fct_C(2, 2, 2, 2).GetCoefficient(i1, jloc);
          T c2201 = fct_C(2, 2, 0, 1).GetCoefficient(i1, jloc);
          T c2202 = fct_C(2, 2, 0, 2).GetCoefficient(i1, jloc);
          T c2212 = fct_C(2, 2, 1, 2).GetCoefficient(i1, jloc);
          T c0101 = fct_C(0, 1, 0, 1).GetCoefficient(i1, jloc);
          T c0102 = fct_C(0, 1, 0, 2).GetCoefficient(i1, jloc);
          T c0112 = fct_C(0, 1, 1, 2).GetCoefficient(i1, jloc);
          T c0202 = fct_C(0, 2, 0, 2).GetCoefficient(i1, jloc);
          T c0212 = fct_C(0, 2, 1, 2).GetCoefficient(i1, jloc);
          T c1212 = fct_C(1, 2, 1, 2).GetCoefficient(i1, jloc);
          dV(0) = c0000*dU(0) + c0011*dU(4) + c0022*dU(8) + 
            c0001*(dU(1)+dU(3)) + c0002*(dU(2)+dU(6)) + c0012*(dU(5)+dU(7));

          dV(4) = c0011*dU(0) + c1111*dU(4) + c1122*dU(8) + 
            c1101*(dU(1)+dU(3)) + c1102*(dU(2)+dU(6)) + c1112*(dU(5)+dU(7));

          dV(8) = c0022*dU(0) + c1122*dU(4) + c2222*dU(8) + 
            c2201*(dU(1)+dU(3)) + c2202*(dU(2)+dU(6)) + c2212*(dU(5)+dU(7));

          dV(1) = c0001*dU(0) + c1101*dU(4) + c2201*dU(8) + 
            c0101*(dU(1)+dU(3)) + c0102*(dU(2)+dU(6)) + c0112*(dU(5)+dU(7));
          dV(3) = dV(1);
          
          dV(2) = c0002*dU(0) + c1102*dU(4) + c2202*dU(8) + 
            c0102*(dU(1)+dU(3)) + c0202*(dU(2)+dU(6)) + c0212*(dU(5)+dU(7));
          dV(6) = dV(2);
          
          dV(5) = c0012*dU(0) + c1112*dU(4) + c2212*dU(8) + 
            c0112*(dU(1)+dU(3)) + c0212*(dU(2)+dU(6)) + c1212*(dU(5)+dU(7));
          dV(7) = dV(5);
        }
        break;
      }
  }


  template<class Dimension, int m, class T> template<class Vector1>
  void ElasticPhysicalIndice<Dimension, m, T>::
  MltVectorPoint(int i1, const VectReal_wp& phi,
                 const R2& pt_glob, const Vector1& dU, Vector1& dV) const
  {
    switch (type_anisotropy)
      {
      case ISOTROPE:
        {
          T lambda = fct_C(0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T mu = fct_C(1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          dV(0) = (lambda+2.0*mu)*dU(0) + lambda*dU(3);
          dV(3) = (lambda+2.0*mu)*dU(3) + lambda*dU(0);
          dV(2) = mu*(dU(1) + dU(2));
          dV(1) = dV(2);        
        }
        break;
      case ORTHOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T c0011 = fct_C(0, 0, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c1111 = fct_C(1, 1, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0101 = fct_C(0, 1, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          dV(0) = c0000*dU(0) + c0011*dU(3);
          dV(3) = c0011*dU(0) + c1111*dU(3);
          dV(1) = c0101*(dU(1)+dU(2));
          dV(2) = dV(1);
        }
        break;
      case ANISOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T c0011 = fct_C(0, 0, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0001 = fct_C(0, 0, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c1111 = fct_C(1, 1, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c1110 = fct_C(1, 1, 1, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T c0101 = fct_C(0, 1, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          dV(0) = c0000*dU(0) + c0011*dU(3) + c0001*(dU(1)+dU(2));
          dV(3) = c0011*dU(0) + c1111*dU(3) + c1110*(dU(1)+dU(2));
          dV(1) = c0001*dU(0) + c1110*dU(3) + c0101*(dU(1)+dU(2));
          dV(2) = dV(1);
        }
        break;
      }
  }


  template<class Dimension, int m, class T> template<class Vector1>
  void ElasticPhysicalIndice<Dimension, m, T>::
  MltVectorPoint(int i1, const VectReal_wp& phi,
                 const R3& pt_glob, const Vector1& dU, Vector1& dV) const
  {
    switch (type_anisotropy)
      {
      case ISOTROPE:
        {
          T lambda = fct_C(0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T mu = fct_C(1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          dV(0) = (lambda+2.0*mu)*dU(0) + lambda*(dU(4) + dU(8));
          dV(4) = (lambda+2.0*mu)*dU(4) + lambda*(dU(0) + dU(8));
          dV(8) = (lambda+2.0*mu)*dU(8) + lambda*(dU(0) + dU(4));
          dV(1) = mu*(dU(1) + dU(3));
          dV(3) = dV(1);        
          dV(2) = mu*(dU(2) + dU(6));
          dV(6) = dV(2);        
          dV(5) = mu*(dU(5) + dU(7));
          dV(7) = dV(5);        
        }
        break;
      case ORTHOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T c0011 = fct_C(0, 0, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0022 = fct_C(0, 0, 2, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c1111 = fct_C(1, 1, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c1122 = fct_C(1, 1, 2, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c2222 = fct_C(2, 2, 2, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0101 = fct_C(0, 1, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0202 = fct_C(0, 2, 0, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c1212 = fct_C(1, 2, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          dV(0) = c0000*dU(0) + c0011*dU(4) + c0022*dU(8);
          dV(4) = c0011*dU(0) + c1111*dU(4) + c1122*dU(8);
          dV(8) = c0022*dU(0) + c1122*dU(4) + c2222*dU(8);
          dV(1) = c0101*(dU(1)+dU(3));
          dV(3) = dV(1);
          dV(2) = c0202*(dU(2)+dU(6));
          dV(6) = dV(2);
          dV(5) = c1212*(dU(5)+dU(7));
          dV(7) = dV(5);
        }
        break;
      case ANISOTROPE :
        {
          T c0000 = fct_C(0, 0, 0, 0).GetIndexAtPoint(i1, phi, pt_glob);
          T c0011 = fct_C(0, 0, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0022 = fct_C(0, 0, 2, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0001 = fct_C(0, 0, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0002 = fct_C(0, 0, 0, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0012 = fct_C(0, 0, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c1111 = fct_C(1, 1, 1, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c1122 = fct_C(1, 1, 2, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c1101 = fct_C(1, 1, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c1102 = fct_C(1, 1, 0, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c1112 = fct_C(1, 1, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c2222 = fct_C(2, 2, 2, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c2201 = fct_C(2, 2, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c2202 = fct_C(2, 2, 0, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c2212 = fct_C(2, 2, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0101 = fct_C(0, 1, 0, 1).GetIndexAtPoint(i1, phi, pt_glob);
          T c0102 = fct_C(0, 1, 0, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0112 = fct_C(0, 1, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0202 = fct_C(0, 2, 0, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c0212 = fct_C(0, 2, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          T c1212 = fct_C(1, 2, 1, 2).GetIndexAtPoint(i1, phi, pt_glob);
          dV(0) = c0000*dU(0) + c0011*dU(4) + c0022*dU(8) + 
            c0001*(dU(1)+dU(3)) + c0002*(dU(2)+dU(6)) + c0012*(dU(5)+dU(7));

          dV(4) = c0011*dU(0) + c1111*dU(4) + c1122*dU(8) + 
            c1101*(dU(1)+dU(3)) + c1102*(dU(2)+dU(6)) + c1112*(dU(5)+dU(7));

          dV(8) = c0022*dU(0) + c1122*dU(4) + c2222*dU(8) + 
            c2201*(dU(1)+dU(3)) + c2202*(dU(2)+dU(6)) + c2212*(dU(5)+dU(7));

          dV(1) = c0001*dU(0) + c1101*dU(4) + c2201*dU(8) + 
            c0101*(dU(1)+dU(3)) + c0102*(dU(2)+dU(6)) + c0112*(dU(5)+dU(7));
          dV(3) = dV(1);
          
          dV(2) = c0002*dU(0) + c1102*dU(4) + c2202*dU(8) + 
            c0102*(dU(1)+dU(3)) + c0202*(dU(2)+dU(6)) + c0212*(dU(5)+dU(7));
          dV(6) = dV(2);
          
          dV(5) = c0012*dU(0) + c1112*dU(4) + c2212*dU(8) + 
            c0112*(dU(1)+dU(3)) + c0212*(dU(2)+dU(6)) + c1212*(dU(5)+dU(7));
          dV(7) = dV(5);
        }
        break;
      }
  }
  

  //! returns true if index is complex
  template<class Dimension, int m, class T>
  bool ElasticPhysicalIndice<Dimension, m, T>::IsComplex() const
  {
    TinyMatrix<int, Symmetric, m, m> indice_doublet;
    indice_doublet.Fill();
    if (type_media == CONSTANT)
      {
        for (int i = 0; i < nb_components; i++)
	  for (int j = i; j < nb_components; j++)
	    if (abs(imagpart(cte_C(i,j))) > 1e-10)
	      return true;
        
        return false;
      }
    else
      {
        if (type_anisotropy == ISOTROPE)
          {
            if (fct_C(0, 0).IsComplex())
              return true;
            
            if (fct_C(1, 1).IsComplex())
              return true;
            
            return false;
          }
        else if (type_anisotropy == ORTHOTROPE)
          {
            for (int i = 0; i < m; i++)
              for (int j = i; j < m; j++)
                {
                  int num1 = indice_doublet(i,i), num2 = indice_doublet(j,j);
                  if (fct_C(num1, num2).IsComplex())
                    return true;
                }
            
            for (int i = 0; i < m; i++)
              for (int j = (i+1); j < m; j++)
                {
                  int num = indice_doublet(i,j);
                  if (fct_C(num, num).IsComplex())
                    return true;                  
                }
            
            return false;
          }
        else if (type_anisotropy == ANISOTROPE)
          {
            for (int i = 0; i < nb_components; i++)
              for (int j = i; j < nb_components; j++)
                if (fct_C(i, j).IsComplex())
                  return true;
            
            return false;
          }
      }
    
    return false;
  }
  
      
  //! clears varying indices
  template<class Dimension, int m, class T>
  void ElasticPhysicalIndice<Dimension, m, T>::Clear()
  {
    for (int i = 0; i < nb_components; i++)
      for (int j = i; j < nb_components; j++)
	fct_C(i, j).Clear();
  }
  
  
  //! reads the physical index in a data file
  /*!
    \param[in] var given problem
    \param[in] nb subscript to access array parameters
    \param[in] parameters values retrieved from data file
    \param[in] keyword keyword of the line of data file
   */
  template<class Dimension, int m, class T>
  void ElasticPhysicalIndice<Dimension, m, T>::
  SetInputData(int& nb, const VectString& parameters,
	       const string& keyword)
  {
    T zero; SetComplexZero(zero);
    if (!keyword.compare("ISOTROPE"))
      {
	type_anisotropy = ISOTROPE; T lambda, mu;
	this->SetInputVaryingMedia(nb, fct_C(0,0), lambda, parameters);
	this->SetInputVaryingMedia(nb, fct_C(1,1), mu, parameters);
	cte_C.FillIsotrope(lambda, mu);
	if ((fct_C(0,0).IsVarying())||(fct_C(1,1).IsVarying()))
	  type_media = VARYING;
	
      }
    else if (!keyword.compare("YOUNG_POISSON"))
      {
	type_anisotropy = ISOTROPE; T lambda, mu, E, nu;
	this->SetInputVaryingMedia(nb, fct_C(0,0), E, parameters);
	this->SetInputVaryingMedia(nb, fct_C(1,1), nu, parameters);
	if (m == 2)
          {
            // relations between Lame coefficients, and E, nu in 2-D
            lambda = E*nu/(1.0 - nu*nu);
            mu = 0.5*E/(1.0+nu);
          }
        else
          {
            // relations in 3-D
            lambda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu));
            mu = 0.5*E/(1.0+nu);
          }
        
	cte_C.FillIsotrope(lambda, mu);
	if ((fct_C(0,0).IsVarying())||(fct_C(1,1).IsVarying()))
	  {
            cout << "variable case not implemented " <<endl;            
            type_media = VARYING;
            abort();
          }	
      }
    else if (!keyword.compare("ORTHOTROPE"))
      {
        for (int i = 0; i < nb_components; i++)
          for (int j = i; j < nb_components; j++)
            fct_C(i, j).Zero();
        
	type_anisotropy = ORTHOTROPE; cte_C.Fill(zero);
	type_media = CONSTANT;
	for (int i = 0; i < m; i++)
	  for (int j = i; j < m; j++)
	    {
	      this->SetInputVaryingMedia(nb, fct_C(i, i, j, j), cte_C(i, i, j, j),
                                         parameters);
	      if (fct_C(i, i, j, j).IsVarying())
                type_media = VARYING;               
	    }
	
	for (int i = 0; i < m; i++)
	  for (int j = (i+1); j < m; j++)
	    {
	      this->SetInputVaryingMedia(nb, fct_C(i, j, i, j),
                                         cte_C(i, j, i, j), parameters);
	      if (fct_C(i, j, i, j).IsVarying())
                type_media = VARYING;
	    }
      }
    else if (!keyword.compare("ORTHOTROPE_POISSON"))
      {
	type_anisotropy = ORTHOTROPE;
        cte_C.Fill(zero);
	type_media = CONSTANT;
        
        if ((m != 2) && (m != 3))
          {
            cout << "only implemented in 2-D and 3-D" << endl;
            abort();
          }
        
        // first, Young coefficients
        T Ex, Ey, Ez;
        this->SetInputVaryingMedia(nb, fct_C(0, 0), Ex, parameters);
        if (fct_C(0, 0).IsVarying())
          {
            cout << "variable case not implemented " <<endl;
            abort();
          }
        
        this->SetInputVaryingMedia(nb, fct_C(0, 0), Ey, parameters);
        if (fct_C(0, 0).IsVarying())
          {
            cout << "variable case not implemented " <<endl;
            abort();
          }
        
        if (m == 3)
          {
            this->SetInputVaryingMedia(nb, fct_C(0, 0), Ez, parameters);
            if (fct_C(0, 0).IsVarying())
              {
                cout << "variable case not implemented " <<endl;
                abort();
              }
          }

        // then Poisson coefficients
        T nu_xy, nu_xz, nu_yz, nu_yx, nu_zx, nu_zy;
        this->SetInputVaryingMedia(nb, fct_C(0, 0), nu_xy, parameters);
        if (fct_C(0, 0).IsVarying())
          {
            cout << "variable case not implemented " <<endl;
            abort();
          }
        
        if (m == 3)
          {
            this->SetInputVaryingMedia(nb, fct_C(0, 0), nu_xz, parameters);
            if (fct_C(0, 0).IsVarying())
              {
                cout << "variable case not implemented " <<endl;
                abort();
              }
            
            this->SetInputVaryingMedia(nb, fct_C(0, 0), nu_yz, parameters);
            if (fct_C(0, 0).IsVarying())
              {
                cout << "variable case not implemented " <<endl;
                abort();
              }
          }
        
        nu_yx = nu_xy*Ey/Ex;
        if (m == 3)
          {
            nu_zx = nu_xz*Ez/Ex;
            nu_zy = nu_yz*Ez/Ey;
          }
        
        // finally shear modulus
        T Gxy, Gxz, Gyz;
        this->SetInputVaryingMedia(nb, fct_C(0, 0), Gxy, parameters);
        if (fct_C(0, 0).IsVarying())
          {
            cout << "variable case not implemented " <<endl;
            abort();
          }
        
        if (m == 3)
          {
            this->SetInputVaryingMedia(nb, fct_C(0, 0), Gxz, parameters);
            if (fct_C(0, 0).IsVarying())
              {
                cout << "variable case not implemented " <<endl;
                abort();
              }
            
            this->SetInputVaryingMedia(nb, fct_C(0, 0), Gyz, parameters);
            if (fct_C(0, 0).IsVarying())
              {
                cout << "variable case not implemented " <<endl;
                abort();
              }
          }

        Real_wp theta(0.5*pi_wp), phi(0);
        if (parameters.GetM() > nb)
          {
            if (parameters(nb++) == "ANGLE")
              {
                if (m == 3)
                  {
		    if (parameters.GetM() <= nb+1)
		      {
			cout << "In SetInputData of ElasticPhysicalIndice" << endl;
			cout << "we need two angles in 3-D" << endl;
			cout << "Current parameters are : " << endl << parameters << endl;
			abort();
		      }

                    theta = to_num<Real_wp>(parameters(nb++))*pi_wp/180.0;
                    phi = to_num<Real_wp>(parameters(nb++))*pi_wp/180.0;
                  }
                else
		  {
		    if (parameters.GetM() <= nb)
		      {
			cout << "In SetInputData of ElasticPhysicalIndice" << endl;
			cout << "we need one angle in 2-D" << endl;
			cout << "Current parameters are : " << endl << parameters << endl;
			abort();
		      }
		    
		    phi = to_num<Real_wp>(parameters(nb++))*pi_wp/180.0;
		  }
                
              }
          }
        
        TinyMatrix<Real_wp, General, m, m> matQ;
        matQ.SetIdentity();
        if (m == 2)
          {
            matQ(0, 0) = cos(phi); matQ(0, 1) = -sin(phi);
            matQ(1, 0) = sin(phi); matQ(1, 1) = matQ(0, 0);
          }
        else
          {
            matQ(0, 0) = sin(theta)*cos(phi);
            matQ(1, 0) = sin(theta)*sin(phi);
            matQ(2, 0) = cos(theta);
            
            matQ(0, 1) = -sin(phi);
            matQ(1, 1) = cos(phi);
            matQ(2, 1) = 0;
            
            matQ(0, 2) = -cos(theta)*cos(phi);
            matQ(1, 2) = -cos(theta)*sin(phi);
            matQ(2, 2) = sin(theta);            
          }
        
        if (m == 3)
          {
            T Delta = (1.0 - nu_xy*nu_yx - nu_zx*nu_xz
                       - nu_yz*nu_zy - 2.0*nu_xy*nu_yz*nu_zx)/(Ex*Ey*Ez);
            
            cte_C(0, 0, 0, 0) = (1.0 - nu_yz*nu_zy)/(Ey*Ez*Delta);
            cte_C(1, 1, 1, 1) = (1.0 - nu_zx*nu_xz)/(Ex*Ez*Delta);
            cte_C(2, 2, 2, 2) = (1.0 - nu_xy*nu_yx)/(Ex*Ey*Delta);
            
            cte_C(0, 0, 1, 1) = (nu_yx + nu_zx*nu_yz)/(Ey*Ez*Delta);            
            cte_C(0, 0, 2, 2) = (nu_zx + nu_yx*nu_zy)/(Ey*Ez*Delta);
            cte_C(1, 1, 2, 2) = (nu_zy + nu_zx*nu_xy)/(Ez*Ex*Delta);
            
            cte_C(0, 1, 0, 1) = Gxy;
            cte_C(0, 2, 0, 2) = Gxz;
            cte_C(1, 2, 1, 2) = Gyz;
          }
        else
          {
            T Delta = 1.0 - nu_xy * nu_yx;
            
            cte_C(0, 0, 0, 0) = Ex/Delta;
            cte_C(1, 1, 1, 1) = Ey/Delta;

            cte_C(0, 0, 1, 1) = nu_yx*Ex/Delta;
            cte_C(0, 1, 0, 1) = Gxy;
          }
        
        if ((theta != Real_wp(0.5*pi_wp)) || (phi != Real_wp(0)))
          cte_C.ApplyRotation(matQ);
        
        type_anisotropy = ANISOTROPE;
      }
    else if (!keyword.compare("ANISOTROPE"))
      {
	type_anisotropy = ANISOTROPE; cte_C.Fill(zero);
	type_media = CONSTANT;
	IVect perm(nb_components);
	if (m == 2)
	  {
	    perm(0) = 0; perm(1) = 2; perm(2) = 1;
	  }
	else
	  {
	    perm(0) = 0; perm(1) = 3; perm(2) = 5;
	    perm(3) = 1; perm(4) = 2; perm(5) = 4;
	  }
	
	for (int i = 0; i < nb_components; i++)
	  for (int j = i; j < nb_components; j++)
	    {
	      this->SetInputVaryingMedia(nb, fct_C(perm(i), perm(j)),
                                         cte_C(perm(i), perm(j)), parameters);
	      if (fct_C(perm(i),perm(j)).IsVarying())
		type_media = VARYING;
	    }
      }
  }
  
  
  //! displays informations about class ElasticPhysicalIndice
  template<class Dimension, int m, class T>
  ostream& operator <<(ostream& out, const ElasticPhysicalIndice<Dimension, m, T>& var)
  {
    if (var.type_media == var.CONSTANT)
      out<<var.cte_C;
    else
      out<<"varying media";
    
    return out;
  }
  
  
  //! reads a physical index in a data file
  /*!
    \param[in] var given problem
    \param[in] nb subscript to access array parameters
    \param[in] parameters values retrieved from data file
    \param[out] var_media varying field
    \param[out] cte constant field
   */  
  template<class Dimension, class T>
  void GenericPhysicalIndice<Dimension, T>::
  SetInputVaryingMedia(int& nb, PhysicalVaryingMedia<Dimension, T>& var_media,
		       T& cte, const VectString& parameters)
  {
    T one; SetComplexOne(one);
    cte = T(0);
    
    if (nb >= parameters.GetM())
      return;
    
    if (!parameters(nb).compare("SINUS"))
      {
	nb++;
	R_N xmin, xmax, coef, kwave;
        T offset(0), amplitude(0);
	if (nb < (parameters.GetM()-3))
	  {
	    xmin(0) = to_num<Real_wp>(parameters(nb++));
            xmax(0) = to_num<Real_wp>(parameters(nb++));
	    xmin(1) = to_num<Real_wp>(parameters(nb++));
            xmax(1) = to_num<Real_wp>(parameters(nb++));
	  }
	
	if (Dimension::dim_N == 3)
	  {
	    if (nb < (parameters.GetM()-1))
	      {
		xmin(2) = to_num<Real_wp>(parameters(nb++));
		xmax(2) = to_num<Real_wp>(parameters(nb++));
	      }
	  }
	
	if (nb < (parameters.GetM()-1))
	  {
	    coef(0) = to_num<Real_wp>(parameters(nb++));
	    coef(1) = to_num<Real_wp>(parameters(nb++));
	  }
	
	if (Dimension::dim_N == 3)
	  if (nb < parameters.GetM())
	    coef(2) = to_num<Real_wp>(parameters(nb++));
	
	if (nb < (parameters.GetM()-3))
	  {
	    offset = to_num<T>(parameters(nb++));
	    amplitude = to_num<T>(parameters(nb++));
	
	    kwave(0) = 2.0*pi_wp/to_num<Real_wp>(parameters(nb++));
	    kwave(1) = 2.0*pi_wp/to_num<Real_wp>(parameters(nb++));
	  }
	
	if (Dimension::dim_N == 3)
	  if (nb < parameters.GetM())
	    kwave(2) = 2.0*pi_wp/to_num<Real_wp>(parameters(nb++));
	
	if (nb < parameters.GetM())
	  cte = to_num<Real_wp>(parameters(nb++));
	
	var_media.InitSinus(xmin, xmax, coef, offset, amplitude, kwave);
      }
    else if ((!parameters(nb).compare("MESH")) || (!parameters(nb).compare("SAME_MESH")))
      {
        if (!parameters(nb).compare("SAME_MESH"))
          var_media.SetSameMesh();
        
	nb++;
	T offset(0), amplitude(0); string data_file, mesh_file;

	if (parameters.GetM() <= nb+4)
	  {
	    cout << "In SetInputVaryingMedia of GenericPhysicalIndice" << endl;
	    cout << "We need at least five parameters : "
                 << "mesh_file offset amplitude data_file order" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	mesh_file = parameters(nb++);
	
	offset = to_num<T>(parameters(nb++));
	amplitude = to_num<T>(parameters(nb++));
	
	data_file = parameters(nb++);
	
	int order = to_num<int>(parameters(nb++));
	
	if (nb < parameters.GetM())
	  cte = to_num<Real_wp>(parameters(nb++));
	
	var_media.InitMesh(mesh_file, offset, amplitude, Dimension::dim_N, data_file, order);        
      }
    else if ((!parameters(nb).compare("RANDOM"))||(!parameters(nb).compare("PERIODIC"))
             ||(!parameters(nb).compare("QUASI_PERIODIC")))
      {
	if (!parameters(nb).compare("PERIODIC"))
	  var_media.index_regular.SetPeriodic();

	if (!parameters(nb).compare("QUASI_PERIODIC"))
	  var_media.index_regular.SetQuasiPeriodic();
	
	nb++;
	R_N xmin, xmax, coef, kwave;  T offset(0), amplitude(0); string data_file;
	if (nb < (parameters.GetM()-3))
	  {
	    xmin(0) = to_num<Real_wp>(parameters(nb++));
            xmax(0) = to_num<Real_wp>(parameters(nb++));
	    xmin(1) = to_num<Real_wp>(parameters(nb++)); 
            xmax(1) = to_num<Real_wp>(parameters(nb++));
	  }		
	
	if (Dimension::dim_N == 3)
	  {
	    if (nb < (parameters.GetM()-1))
	      {
		xmin(2) = to_num<Real_wp>(parameters(nb++));
		xmax(2) = to_num<Real_wp>(parameters(nb++));
	      }
	  }
	
	if (nb < (parameters.GetM()-1))
	  {
	    coef(0) = to_num<Real_wp>(parameters(nb++));
	    coef(1) = to_num<Real_wp>(parameters(nb++));
	  }
	
	if (Dimension::dim_N == 3)
	  if (nb < parameters.GetM())
	    coef(2) = to_num<Real_wp>(parameters(nb++));
	
	if (nb < (parameters.GetM()-2))
	  {
	    offset = to_num<T>(parameters(nb++));
	    amplitude = to_num<T>(parameters(nb++));
	    
	    data_file = parameters(nb++);
            cout << "Data file read: " << data_file << endl;
	  }
	
	string type_number("DOUBLE");
	if (nb < parameters.GetM())
	  type_number = parameters(nb++);
	
	bool complex_number = false, double_prec = false;
	if (!type_number.compare("DOUBLE"))
	  double_prec = true;
	else if (!type_number.compare("COMPLEX"))
	  complex_number = true;
	else if (!type_number.compare("COMPLEX_DBLE"))
	  {
	    complex_number = true;
	    double_prec = true;
	  }
		
	if (nb < parameters.GetM())
	  cte = to_num<Real_wp>(parameters(nb++));
		
	var_media.InitRandom(xmin, xmax, coef, offset, amplitude,
                             Dimension::dim_N, data_file, complex_number, double_prec);
      }
    else if (!parameters(nb).compare("RADIAL"))
      {
        nb++;
        string data_file = parameters(nb++);
	T offset = to_num<T>(parameters(nb++));
	T amplitude = to_num<T>(parameters(nb++));
        
        var_media.InitRadial(data_file, offset, amplitude, cte);
      }
    else if (!parameters(nb).compare("USER"))
      {
        nb++;
        
	T offset = to_num<T>(parameters(nb++));
	T amplitude = to_num<T>(parameters(nb++));
        cte = to_num<Real_wp>(parameters(nb++));
        
        var_media.InitUser(offset, amplitude, cte);
      }
    else if (!parameters(nb).compare("VELOCITY"))
      {
	nb++;
	if (nb < parameters.GetM())
	  {
	    cte = to_num<T>(parameters(nb++));
	    cte = one/(cte*cte);
	    var_media.InitConstant(cte);
	  }
      }
    else if (!parameters(nb).compare("SQUARE"))
      {
	nb++;
	if (nb < parameters.GetM())
	  {
	    cte = to_num<T>(parameters(nb++));
	    cte = cte*cte;
	    var_media.InitConstant(cte);
	  }
      }
    else if (!parameters(nb).compare("INVERSE"))
      {
	nb++;
	if (nb < parameters.GetM())
	  {
	    cte = to_num<T>(parameters(nb++));
	    cte = one/cte;
	    var_media.InitConstant(cte);
	  }
      }
    else if (!parameters(nb).compare("ANGLE"))
      {
	nb++;
	if (nb < parameters.GetM())
	  {
	    cte = to_num<T>(parameters(nb++));
	    cte = cte*pi_wp/Real_wp(180);
	    var_media.InitConstant(cte);
	  }
      }
    else
      {
	if (nb < parameters.GetM())
	  {
	    cte = to_num<T>(parameters(nb++));
	    var_media.InitConstant(cte);
	  }
      }
    
  }
  
}

#define MONTJOIE_FILE_PHYSICAL_PROPERTY_CXX
#endif

