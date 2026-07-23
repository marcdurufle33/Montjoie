#ifndef MONTJOIE_FILE_XDR_HXX

#include <assert.h>

#define HAVE_XDR
#define USE_COMPLEX_NUMBERS

// Local includes
namespace libMeshEnums {
  
  /**
   * Defines an \p enum for read/write mode in Xdr format.
   * \p READ, \p WRITE perform reading and writing in ASCII format,
   * and \p DECODE, \p ENCODE do the same in binary format.
   */
  enum XdrMODE
    {
      UNKNOWN = -1, ENCODE=0, DECODE, WRITE, READ
    };
}

using namespace libMeshEnums;

#include "lib_xdr/Xdr_Ofstream.cxx"

#ifndef BROKEN_IOSTREAM
 /*
  * ---------------------------------------------------------------------------------
  * Everything for a clean iostream
  */

# include <iomanip>

 /*
  * Outputs \p Real \p d  with width
  * \p v in scientific format to stream \p o.
  */
# define OFSRealscientific(o,v,d)       (o) << std::setw(v) << std::scientific << (d)

 /*
  * Outputs \p Number \p d, (note that \p Number
  * may be either real or complex) with width
  * \p v in scientific format to stream \p o.
  */
# if defined(USE_COMPLEX_NUMBERS) 
#  define OFSNumberscientific(o,v,d) (o) << std::setw(v) << std::scientific << (d).real() << " " \
                                          << std::setw(v) << std::scientific << (d).imag()
# else
#  define OFSNumberscientific(o,v,d)    (o) << std::setw(v) << std::scientific << (d)
# endif

 /*
  * class alias
  */
# define OFStream                    std::ofstream


#else
 /*
  * ---------------------------------------------------------------------------------
  * Everything for broken iostream
  */

# include <stdio.h>

 /*
  * Outputs \p Real \p d with width 
  * \p v in scientific format to stream \p o.
  */
# define OFSRealscientific(o,v,d)       (o).scientific( (v), (d) )

 /*
  * Outputs \p Number \p d, (note that \p Number
  * may be either real or complex) with width
  * \p v in scientific format to stream \p o.
  */
# define OFSNumberscientific(o,v,d)     (o).scientific( (v), (d) )

//  /*
//   * class alias
//   */
// # define OFSOFStream                    OFStream



 /**
  * This class provides a compatibility class for broken
  * features in the \p std::ofstream of the older \p GCC
  * versions.  Other compilers do not see this class.
  */

 // ------------------------------------------------------------
 // OFStream class definition
 class OFStream : public std::ofstream
 {
 public:

  /**
   * Default constructor.
   */
  OFStream () {};

  /**
   * Default destructor.
   */
  ~OFStream () {};

  /**
   * convenient typedef
   */
  typedef std::string::size_type sizetype;

  /**
   * @returns a \p OFStream, where \p r
   * was directed in scientific style with
   * size \p w.
   */
  template <typename T>
  OFStream& scientific (const sizetype w,
			const T r);

 };


 // ------------------------------------------------------------
 // OFStream inline methods
 

#endif // ifndef ... else ... BROKEN_IOSTREAM

// C++ includes
#ifdef HAVE_XDR
#  include <rpc/rpc.h>
#endif

#include <string>
#ifdef USE_COMPLEX_NUMBERS
# include <complex>
#endif

#define xdr_MAX_STRING_LENGTH 256

#ifndef SINGLE_PRECISION
#define xdr_REAL xdr_double
#else
#define xdr_REAL xdr_float
#endif



//--------------------------------------------------------------
// Xdr class definition

/**
 * This class implements a C++ interface to the XDR 
 * (eXternal Data Representation) format.  XDR is useful for
 * creating platform-independent binary files.  This class was
 * created to handle equation system output as a replacement for
 * XdrIO since that is somewhat limited.
 */

class Xdr
{
  
public:

  /**
   * Constructor.  Takes the filename and the mode.
   * Valid modes are ENCODE, DECODE, READ, and WRITE.
   */
  Xdr (const std::string& name="", const libMeshEnums::XdrMODE m=UNKNOWN);

  /**
   * Destructor.  Closes the file if it is open.
   */
  ~Xdr ();
  
  /**
   * Opens the file.
   */ 
  void open (const std::string& name);

  /**
   * Closes the file if it is open.
   */ 
  void close();

  /**
   * Returns true if the Xdr file is open, false
   * if it is closed.
   */ 
  bool is_open() const;

  /**
   * Returns true if the file is opened in a reading
   * state, false otherwise.
   */
  bool reading() const { return ((mode == DECODE) || (mode == READ)); }

  /**
   * Returns true if the file is opened in a writing
   * state, false otherwise.
   */
  bool writing() const { return ((mode == ENCODE) || (mode == WRITE)); }

  /**
   * Returns the mode used to access the file.  Valid modes
   * are ENCODE, DECODE, READ, or WRITE.
   */
  XdrMODE access_mode () const { return mode; }

  // Data access methods

  /**
   * Inputs or outputs a single integer.
   */
  void data(int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (int a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (int& a) { assert (reading()); data(a); return *this; }
  
  /**
   * Inputs or outputs a single unsigned integer.
   */
  void data(unsigned int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (unsigned int a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (unsigned int& a) { assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single short integer.
   */
  void data(short int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (short int a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (short int& a) { assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single unsigned short integer.
   */
  void data(unsigned short int& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (unsigned short int a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (unsigned short int& a) { assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single float.
   */
  void data(float& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (float a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (float& a) { assert (reading()); data(a); return *this; }

  /**
   * Inputs or outputs a single double.
   */
  void data(double& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (double a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (double& a) { assert (reading()); data(a); return *this; }


#ifdef USE_COMPLEX_NUMBERS

  /**
   * Inputs or outputs a single complex<double>.
   */
  void data(std::complex<double>& a, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (std::complex<double> a) { assert (writing()); data(a); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::complex<double>& a) { assert (reading()); data(a); return *this; }

#endif


  /**
   * Inputs or outputs a vector of integers.
   */
  template<class Allocator>
  void data(Vector<int, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<int, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<int, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of unsigned integers.
   */
  template<class Allocator>
  void data(Vector<unsigned int, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<unsigned int, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<unsigned int, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of short integers.
   */
  template<class Allocator>
  void data(Vector<short int, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<short int, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<short int, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of unsigned short integers.
   */
  template<class Allocator>
  void data(Vector<unsigned short int, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<unsigned short int, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<unsigned short int, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of floats.
   */
  template<class Allocator>
  void data(Vector<float, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<float, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<float, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }

  /**
   * Inputs or outputs a vector of doubles.
   */
  template<class Allocator>
  void data(Vector<double, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<double, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<double, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }


#ifdef USE_COMPLEX_NUMBERS

  /**
   * Inputs or outputs a vector of complex<double>.
   */
  template<class Allocator>
  void data(Vector< std::complex<double>, Vect_Full, Allocator>& v, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  template<class Allocator>
  Xdr& operator << (Vector<complex<double>, Vect_Full, Allocator>& v)
  { assert (writing()); data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  template<class Allocator>
  Xdr& operator >> (Vector<complex<double>, Vect_Full, Allocator>& v)
  { assert (reading()); data(v); return *this; }

#endif


  /**
   * Inputs or outputs a single string.
   */
  void data(std::string& s, const char* comment="");

  /**
   * Same, but provides an \p ostream like interface.
   */
  Xdr& operator << (const std::string& s)
  { assert (writing()); string v = s; data(v); return *this; }

  /**
   * Same, but provides an \p istream like interface.
   */
  Xdr& operator >> (std::string& v) { assert (reading()); data(v); return *this; }


private:

  /**
   * The mode used for accessing the file.
   */ 
  const XdrMODE mode;

#ifdef HAVE_XDR
  
  /**
   * Pointer to the standard @p xdr
   * struct.  See the standard
   * header file rpc/rpc.h
   * for more information.
   */
  XDR* xdrs;

  /**
   * File pointer.
   */
  FILE* fp;
  
#endif

  /**
   * The output file stream.
   * Use the customized class to enable
   * features also for compilers with broken
   * iostream
   */
  OFStream out;

  /**
   * The input file stream.
   */
  std::ifstream in;

  /**
   * A buffer to put comment strings into.
   */
  const int comm_len;
  char comm[xdr_MAX_STRING_LENGTH];  
};

#include <fstream>
#include <string>
#include <vector>



/*
 * Some compilers, at least HP \p aCC do not even
 * accept empty classes derived from \p std::ostringstream.
 * Therefore, resort to preprocessor definitions. 
 */


#define MONTJOIE_FILE_XDR_HXX
#endif
