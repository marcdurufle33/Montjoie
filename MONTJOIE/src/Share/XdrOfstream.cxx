#ifndef MONTJOIE_FILE_XDR_OFSTREAM_CXX

/*
 * This class is only alive when iostream is broken
 */
#ifdef BROKEN_IOSTREAM



// the template for reals
template <typename T>
OFStream& OFStream::scientific (const sizetype w,
				const T r)
{
  assert (w < 30);
  char buf[30];  
  char format[8];
  // form the format for r
  sprintf (format, "%%%de", w);
  // form string as desired
  sprintf (buf, format, r);
  *this << buf;
  return *this;
}



// full specialization in case of complex numbers
#if defined(USE_COMPLEX_NUMBERS) 

template <>
OFStream& OFStream::scientific (const sizetype w,
				const Complex r)
{
  assert (w < 30);
  char buf[60];  
  char format[16];
  // form the format for r
  sprintf (format, "%%%de %%%de", w, w);
  // form string as desired
  sprintf (buf, format, r.real(), r.imag());
  *this << buf;
  return *this;
}

#endif // if defined(USE_COMPLEX_NUMBERS) 



//--------------------------------------------------------------
// Explicit instantiations for reals
template OFStream& OFStream::scientific (const sizetype w,
					 const double r);

template OFStream& OFStream::scientific (const sizetype w,
					 const float r);



#endif // ifdef BROKEN_IOSTREAM

#define MONTJOIE_FILE_XDR_OFSTREAM_CXX
#endif
