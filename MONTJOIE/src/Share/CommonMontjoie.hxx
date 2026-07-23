#ifndef MONTJOIE_FILE_COMMON_HXX

namespace std
{
  template<class T>
  T square(const T& x);
}

namespace Seldon
{
  /* functions in CommonInline.cxx */

  void Sort (int& a, int& b);
  
  template<class T, int m>
  void Sort (TinyVector<T, m>& a, TinyVector<T, m>& b);

  void Sort (int& a, int& b, int& c);
  void Sort(int& i, int& j, int& k, int& l);
  
  template<class T, class Storage1, class Allocator1>
  void Check_And_ReallocateIfNecessary(Vector<T, Storage1, Allocator1>& X, int n);

  template<class T>
  void SetIncidentAngle(const T& omega, TinyVector<T, 2>& kwave,
			const T& teta, const T& phi);
  
  template<class T>
  void SetIncidentAngle(const T& omega, TinyVector<T, 3>& kwave,
			const T& teta, const T& phi);
  
  template<class T>
  void PrintNoBrackets(ostream& out, const Vector<T>& u);
  
  template<class T>
  void PrintNoBrackets(ostream& out, const T& u);
  
  template<class T, int m>
  void PrintNoBrackets(ostream& out, const TinyVector<T, m>& u);

  template<class T>
  size_t GetMemorySize(const Vector<Vector<T> >&);

  template<class T>
  size_t GetMemorySize(const Vector<Vector<T, VectSparse> >& x);

  template<class T, class Prop, class Storage>
  size_t GetMemorySize(const Vector<Matrix<T, Prop, Storage> >&);

  template<class T>
  size_t GetMemorySize(const Array2D<Matrix<T> >&);

  template<class T>
  size_t GetMemorySize(const Vector<Vector<Vector<T> > >&);

  template<class T, class Prop, class Storage>
  size_t GetMemorySize(const Vector<Vector<Matrix<T, Prop, Storage> > >&);

  template<class T, int q>
  size_t GetMemorySize(const Vector<TinyVector<Vector<T>, q> >&);

  template<class T, int q>
  size_t GetMemorySize(const TinyVector<Vector<T>, q>&);

  /* functions in Common.cxx */

  size_t GetMemorySize(const string&);
  size_t GetMemorySize(const Vector<string>&);
  
  void MakeLink(string name1, string name2);
  void CopyFile(string name1, string name2);
 
  void StringTokenize(const string& chaine, Vector<string>& param, const string& delim="");
  void StringTrim(string & s);
  void DeleteSpaceAtExtremityOfString(string & s);
  string NumberToString(int number, int nb_char = 4);
  
  template<class T>
  void Copy(list<T>& x, Vector<T>& y);
  
  template<class T0, class T1, class Allocator0, class Allocator1>
  void CopyVector(const Vector<T0, VectFull, Allocator0>& x,
		  Vector<T1, VectSparse, Allocator1>& y);
 
  template<class T>
  void Linspace(const T& a, const T& b, int N, Vector<T>& x);

  string GetHumanReadableMemory(size_t taille);
  
#ifdef SELDON_WITH_MPI
  void MPI_Bcast_string(string& s, int root, const MPI_Comm& comm);

#ifdef MONTJOIE_WITH_MPFR
  void FillVectorMPI(mpfr::mpreal* x, int n, Vector<int64_t>& xtmp);
  void FillVectorMPI(Vector<int64_t>& xtmp, int n, mpfr::mpreal* x);
  
  void FillVectorMPI(complex<mpfr::mpreal>* x, int n, Vector<int64_t>& xtmp);
  void FillVectorMPI(Vector<int64_t>& xtmp, int n, complex<mpfr::mpreal>* x);

  void ReduceVectorMPI(const MPI_Comm& comm, const MPI_Op& op,
                       Vector<int64_t>& xtmp, int n, mpfr::mpreal* x);

  void ReduceVectorMPI(const MPI_Comm& comm, const MPI_Op& op,
                       Vector<int64_t>& xtmp, int n, complex<mpfr::mpreal>* x);
  
  MPI_Request MpiIsend(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIsend(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIsend(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIsend(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIrecv(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIrecv(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIrecv(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  MPI_Request MpiIrecv(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  void MpiCompleteIrecv(mpfr::mpreal* x, Vector<int64_t>& xtmp, int n);
  void MpiCompleteIrecv(Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp, int n);
  void MpiCompleteIrecv(complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp, int n);
  void MpiCompleteIrecv(Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp, int n);

  void MpiSsend(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);

  void MpiSsend(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);

  void MpiSsend(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);

  void MpiSsend(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);

  void MpiSend(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);

  void MpiSend(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);

  void MpiSend(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);

  void MpiSend(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);

  void MpiGather(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                 mpfr::mpreal* y, int n, int proc);

  void MpiGather(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                 Vector<mpfr::mpreal>& y, int n, int proc);

  void MpiGather(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                 complex<mpfr::mpreal>* y, int n, int proc);

  void MpiGather(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                 Vector<complex<mpfr::mpreal> >& y, int n, int proc);

  void MpiReduce(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                 mpfr::mpreal* y, int n, const MPI_Op& op, int proc);
  
  void MpiReduce(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                 Vector<mpfr::mpreal>& y, int n, const MPI_Op& op, int proc);

  void MpiReduce(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                 complex<mpfr::mpreal>* y, int n, const MPI_Op& op, int proc);

  void MpiReduce(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                 Vector<complex<mpfr::mpreal> >& y, int n, const MPI_Op& op, int proc);
  
  void MpiRecv(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status);

  void MpiRecv(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status);

  void MpiRecv(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status);

  void MpiRecv(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status);

  void MpiBcast(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp, int n, int proc);
  void MpiBcast(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp, int n, int proc);

  void MpiBcast(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp, int n, int proc);
  void MpiBcast(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp, int n, int proc);
  
  void MpiAllreduce(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                    Vector<mpfr::mpreal>& y, int n, const MPI_Op& op);

  void MpiAllreduce(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                    Vector<complex<mpfr::mpreal> >& y, int n, const MPI_Op& op);

  void MpiAllreduce(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                    mpfr::mpreal* y, int n, const MPI_Op& op);

  void MpiAllreduce(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                    complex<mpfr::mpreal>* y, int n, const MPI_Op& op);
#endif
  
#endif 


  void PrintFactorizationFailed(int type, int ierr);
  
}

namespace Montjoie
{
  Real_wp GetMinimumDotProd(const R2& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax);
  
  Real_wp GetMinimumDotProd(const R3& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax);

  Real_wp GetMaximumDotProd(const R2& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax);
  
  Real_wp GetMaximumDotProd(const R3& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax);

  void CartesianToPolar(const Real_wp& x, const Real_wp& y, Real_wp& r, Real_wp& theta);
  
  void CartesianToSpherical(const Real_wp& x, const Real_wp& y, const Real_wp& z,
                            Real_wp& r, Real_wp& teta, Real_wp& phi,
                            Real_wp& cos_teta, Real_wp& sin_teta);

  void CartesianToSpherical(const Real_wp& x, const Real_wp& y, const Real_wp& z,
                            Real_wp& r, Real_wp& teta, Real_wp& phi);

  void SphericalToCartesian(const Real_wp& r, const Real_wp& teta, const Real_wp& phi,
                            Real_wp& x, Real_wp& y, Real_wp& z);

  template<class T>
  void ComputePadeCoefficientsSqrt(const T& alpha, int n, complex<T>& CoefC0,
				   Vector<complex<T> >& CoefAl, Vector<complex<T> >& CoefBl);
  
} // namespace Montjoie

#define MONTJOIE_FILE_COMMON_HXX
#endif
