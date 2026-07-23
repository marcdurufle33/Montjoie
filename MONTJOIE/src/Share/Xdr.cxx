#ifndef MONTJOIE_FILE_XDR_CXX

// C/C++ includes
#include <iostream>
#include <string.h>

// Local includes
#include "Xdr.hxx"

//-------------------------------------------------------------
// Xdr class implementation
Xdr::Xdr (const std::string& name, const XdrMODE m) :
  mode(m),
#ifdef HAVE_XDR
  xdrs(NULL),
  fp(NULL),
#endif
  comm_len(xdr_MAX_STRING_LENGTH)
{
  this->open(name);
}



Xdr::~Xdr()
{
  close();
}



void Xdr::open (const std::string& name)
{
  if (name == "")
    return;

  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	fp = fopen(name.c_str(), (mode == ENCODE) ? "w" : "r");
	assert (fp);
	xdrs = new XDR;
	xdrstdio_create (xdrs, fp, (mode == ENCODE) ? XDR_ENCODE : XDR_DECODE);
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();
	
#endif
	return;

      }

    case READ:
      {
	in.open(name.c_str(), std::ios::in);
	assert (in.good());
	return;
      }

    case WRITE:
      {
	out.open(name.c_str(), std::ios::out);
	assert (out.good());
	return;
      }
      
    default:
      abort();
    }  
}



void Xdr::close ()
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	if (xdrs)
	  {
	    xdr_destroy (xdrs);
	    delete xdrs;
	    xdrs = NULL;
	  }
	
	if (fp)
	  {
	    fflush(fp);
	    fclose(fp);
	    fp = NULL;
	  }
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();
	
#endif
	return;
      }
      
    case READ:
      {
	if (in.is_open()) 
	  in.close();      
	return;
      }

    case WRITE:
      {
	if (out.is_open()) 
	  out.close();      
	return;
      }

    default:
      abort();
    }
}



bool Xdr::is_open() const
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	if (fp)
	  if (xdrs)
	    return true;

	return false;

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

	return false;	

#endif

      }
      
    case READ:
      {
	return in.good();
      }

    case WRITE:
      {
	return out.good();
      }

    default:
      abort();
    }

  return false;
}



void Xdr::data (int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_int(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());
	
	out << a << "\t " << comment << '\n';

	return;
      }

    default:
      abort();
    }
}



void Xdr::data (unsigned int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (this->is_open());

	xdr_u_int(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      abort();
    }
}



void Xdr::data (short int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_short(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      abort();
    }
}



void Xdr::data (unsigned short int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_u_short(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      abort();
    }
}



void Xdr::data (float& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_float(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      abort();
    }
}



void Xdr::data (double& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_double(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      abort();
    }
}



#ifdef USE_COMPLEX_NUMBERS

void Xdr::data (std::complex<double>& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());
	double a_real = a.real(), a_imag = a.imag();
	xdr_double(xdrs, &a_real); 
	xdr_double(xdrs, &a_imag);
	a = std::complex<double>(a_real, a_imag);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	double _r, _i;
	in >> _r;
	in  >> _i;
	a = std::complex<double>(_r,_i);
        in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a.real() << "\t " 
	    << a.imag() << "\t "
	    << comment << '\n';
	
	return;
      }

    default:
      abort();
    }
}

#endif // USE_COMPLEX_NUMBERS


template<class Allocator>
void Xdr::data (Vector<int, Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();

	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(int),
		   (xdrproc_t) xdr_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();

	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(int),
		   (xdrproc_t) xdr_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (in.good());
	    in >> v(i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (out.good());
	    out << v(i) << " ";
	  }
	
	out << "\t " << comment << '\n';
	
	return;	
      }

    default:
      abort();
    }
}


template<class Allocator>
void Xdr::data (Vector<unsigned int, Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();
	
	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(unsigned int),
		   (xdrproc_t) xdr_u_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();
	
	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned int),
		   (xdrproc_t) xdr_u_int);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif	
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length = v.GetM();

	for (unsigned int i=0; i<length; i++)
	  {
	    assert (in.good());
	    in >> v(i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (out.good());
	    out << v(i) << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}


template<class Allocator>
void Xdr::data (Vector<short int, Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();
	
	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(short int),
		   (xdrproc_t) xdr_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();

	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(short int),
		   (xdrproc_t) xdr_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length = v.GetM();
	
	for (unsigned int i=0; i<length; i++)
	  {
	    assert (in.good());
	    in >> v(i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (out.good());
	    out << v(i) << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}


template<class Allocator>
void Xdr::data (Vector<unsigned short int,Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();
	
	xdr_vector(xdrs, 
		   (char*) &v.GetData(),
		   length,
		   sizeof(unsigned short int),
		   (xdrproc_t) xdr_u_short);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();

	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(unsigned short int),
		   (xdrproc_t) xdr_u_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (in.good());
	    in >> v(i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (out.good());
	    out << v(i) << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}


template<class Allocator>
void Xdr::data (Vector<float, Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();

	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(float),
		   (xdrproc_t) xdr_float);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.GetM();
	
	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(float),
		   (xdrproc_t) xdr_float);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (in.good());
	    in >> v(i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());
	
	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (out.good());
	    OFSRealscientific(out,12,v(i)) << " ";
 	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}


template<class Allocator>
void Xdr::data (Vector<double, Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (this->is_open());

	unsigned int length = v.GetM();

	xdr_vector(xdrs, 
		   (char*) v.GetData(),
		   length,
		   sizeof(double),
		   (xdrproc_t) xdr_double);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (this->is_open());

	unsigned int length = v.GetM();
	
	// Note: GCC 3.4.1 will crash in debug mode here if length
	// is zero and you attempt to access the zeroth index of v.
	if (length > 0)
	  xdr_vector(xdrs, 
		     (char*) v.GetData(),
		     length,
		     sizeof(double),
		     (xdrproc_t) xdr_double);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	// If you were expecting to read in a vector at this
	// point, it's not going to happen if length == 0!
	// assert (length != 0);
	
	for (int i = 0; i < v.GetM(); i++)
	  {
	    assert (in.good());
	    in >> v(i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	for (int i = 0; i < v.GetM(); i++)
	  {
	    assert (out.good());
	    OFSRealscientific(out,12,v(i)) << " ";
 	  }



	out << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}




#ifdef USE_COMPLEX_NUMBERS

template<class Allocator>
void Xdr::data (Vector< std::complex<double>, Vect_Full, Allocator>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());
	
	for (int i = 0; i < v.GetM(); i++)
	    data(v(i), "");
	
// Alternative code
// 	/*
// 	 * save complex values as two std::vectors<double>.
// 	 * Using just one buffer increases time for copying,
// 	 * but reduces memory consumption
// 	 */
// 	std::vector<double> buf;
// 	buf.resize(length);

// 	// real part
// 	std::vector< std::complex<double> >::iterator c_iter   = v.begin();
// 	std::vector<double>::iterator                 buf_iter = buf.begin();
// 	for (; c_iter != v.end(); ++c_iter)
// 	{
// 	  *buf_iter = c_iter->real();
// 	  ++buf_iter;
// 	}
// 	data(buf, "");

// 	// imaginary part
// 	c_iter   = v.begin();
// 	buf_iter = buf.begin();
// 	for (; c_iter != v.end(); ++c_iter)
// 	{
// 	  *buf_iter = c_iter->real();
// 	  ++buf_iter;
// 	}
// 	data(buf, "");

// 	buf.clear();


// did not work...
// 	// with null pointer, let XDR dynamically allocate?
// 	xdr_vector(xdrs, 
// 		   (char*) &v[0],
// 		   length,
// 		   sizeof(std::complex<double>),
// 		   (xdrproc_t) 0);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	for (int i = 0; i < v.GetM(); i++)
	    data(v(i), "");

// alternative code
// 	/*
// 	 * load complex values as two std::vector<double>
// 	 * one after the other, store them in two buffers,
// 	 * since we have @e no chance to get the real and complex
// 	 * part one after the other into the std::complex
// 	 * (apart from messing with += or so, which i don't want to)
// 	 */
// 	std::vector<double> real_buf, imag_buf;
// 	real_buf.resize(length);
// 	imag_buf.resize(length);

// 	// get real & imaginary part
// 	data(real_buf, "");
// 	data(imag_buf, "");

// 	// copy into vector
// 	std::vector< std::complex<double> >::iterator c_iter   = v.begin();
// 	std::vector<double>::iterator                 real_buf_iter = real_buf.begin();
// 	std::vector<double>::iterator                 imag_buf_iter = imag_buf.begin();

// 	for (; c_iter != v.end(); ++c_iter)
// 	{
// 	  *c_iter = std::complex<double>(*real_buf_iter, *imag_buf_iter);
// 	  ++real_buf_iter;
// 	  ++imag_buf_iter;
// 	}

// 	// clear up
// 	real_buf.clear();
// 	imag_buf.clear();


// did not work...
// 	xdr_vector(xdrs, 
// 		   (char*) &v[0],
// 		   length,
// 		   sizeof(std::complex<double>),
// 		   (xdrproc_t) 0);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	for (unsigned int i=0; i<v.GetM(); i++)
	  {
	    assert (in.good());
	
	    double _r, _i;
	    in >> _r;
	    in  >> _i;
	    v(i) = std::complex<double>(_r,_i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());
 	
	for (unsigned int i=0; i<v.GetM(); i++)
 	  {
 	    assert (out.good());
	    OFSNumberscientific(out,12,v(i)) << " ";
 	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}

#endif // ifdef USE_COMPLEX_NUMBERS




void Xdr::data (std::string& s, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	{
	  char* sptr = new char[s.size()+1];

	  for (unsigned int c=0; c<s.size(); c++)
	    sptr[c] = s[c];
	
	  sptr[s.size()] = '\0';
	  
	  xdr_string(xdrs,
		     &sptr,
		     strlen(sptr));

	  delete [] sptr;
	}

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	{
	  char* sptr = new char[xdr_MAX_STRING_LENGTH];
	  
	  xdr_string(xdrs,
		     &sptr,
		     xdr_MAX_STRING_LENGTH);

	  s.resize(strlen(sptr));

	  for (unsigned int c=0; c<s.size(); c++)
	    s[c] = sptr[c];
	  
	  delete [] sptr;  
	}

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	abort();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	in.getline(comm, comm_len);

//#ifndef BROKEN_IOSTREAM
//	s.clear();
//#else
	s = "";
//#endif

	for (unsigned int c=0; c<strlen(comm); c++)
	  {
	    if (comm[c] == '\t') 
	      break;
	    
	    s.push_back(comm[c]);
	  }

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	out << s << "\t " << comment << '\n';

	return;	
      }

    default:
      abort();
    }
}



#undef xdr_REAL

#define MONTJOIE_FILE_XDR_CXX
#endif
