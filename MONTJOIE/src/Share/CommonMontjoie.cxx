#ifndef MONTJOIE_FILE_COMMON_CXX

#include<list>
#include "Share/CommonMontjoie.hxx"
#include "Share/CommonMontjoieInline.cxx"

namespace Seldon
{
  
  size_t GetMemorySize(const string& s)
  {
    return s.size() + sizeof(s);
  }


  size_t GetMemorySize(const Vector<string>& x)
  {
    size_t taille = x.GetMemorySize();
    for (int i = 0; i < x.GetM(); i++)
      taille += x(i).size();
    
    return taille;
  }
  

  //! creation of a symbolic link 
  /*!
    \param[in] name1 file name of destination
    \param[in] name2 file name of source
   */
  void MakeLink(string name2, string name1)
  {
    if (name1.compare(name2))
      {
	int success = system("pwd > sortie");
	if (success == 0)
	  {
	    string path;
	    ifstream file_in("sortie"); file_in >> path; file_in.close();
	    string command = "rm -f " + name1 + "; ln -s " + path + "/" + name2 +" " + name1;
	    // DISP(command);
	    success = system(command.data());
	    if (success != 0)
	      {
		// nothing
	      }
	  }
      }
  }
  
  
  //! copies a file into another one
  void CopyFile(string name1, string name2)
  {
    if (name1.compare(name2))
      {
	string command = "cp " + name1 + " " + name2;
	int success = system(command.data());
        if (success != 0)
          {
            cerr << "Error: unable to copy " << name2 << " from " << name1
                 << endl << "Abort..." << endl;
            abort();
          }
      }
  }
  
    
  //! Splits a string into several pieces (=tokens) by taking into account parenthesis
  void StringTokenize(const string& chaine, Vector<string>& param, const string& delim)
  {
    list<string> param_;
    unsigned pos = 0;
    int level = 0;
    for (unsigned i = 0; i < chaine.size(); i++)
      {
	if (chaine[i] == '(')
	  level++;

	if (chaine[i] == ')')
	  level--;

	if ((delim.find(chaine[i]) != string::npos) && (level == 0))
	  {
	    // a delimiting character has been found => new token
	    if (i > pos)
	      param_.push_back(chaine.substr(pos, i-pos));
	    
	    pos = i+1;
	  }
      }
    
    if (chaine.size() > pos)
      param_.push_back(chaine.substr(pos, chaine.size()-pos));
    
    Copy(param_, param);
  }
  
  
  //! Removes space and tab characters of a string
  void StringTrim(string & s)
  {
    string s2(s);
    unsigned nb = 0;
    for (unsigned int i = 0; i < s2.size(); i++)
      if ((s2[i] != ' ') && (s2[i] != '\t'))
	nb++;
    
    s.resize(nb);
    nb = 0;
    for (unsigned int i = 0; i < s2.size(); i++)
      if ((s2[i] != ' ') && (s2[i] != '\t'))
	s[nb++] = s2[i];    
  }
  
  
  //! removes space characters at the two extremities of a string
  void DeleteSpaceAtExtremityOfString(string & s)
  {
    if (s.size() <= 0)
      return;
    
    // gets the initial and final position
    unsigned pos_begin = 0;
    for (unsigned int i = 0; i < s.size(); i++)
      {
	if ((s[i] == ' ') || (s[i] == '\t'))
	  pos_begin = i+1;
	else
	  break;
      }

    unsigned pos_end = s.size();
    for (unsigned int i = s.size(); i != 0; i--)
      {
	if ((s[i-1] == ' ') || (s[i-1] == '\t'))
	  pos_end = i-1;
	else
	  break;
      }
    
    if (pos_end <= pos_begin)
      s.clear();
    else
      s = s.substr(pos_begin, pos_end-pos_begin);
  }
  
  
  //! converts a number to a string with fixed format xxxx
  string NumberToString(int number, int nb_char)
  {
    if (number < 0)
      {
	cout << "The number must be positive" << endl;
	abort();
      }
    
    string res(nb_char, '0');
    int p = number, r;
    int i = nb_char;
    while (i > 0)
      {
	i--;
	r = p%10; p = p/10;
	res[i] = 48+r; // 48 -> ascii code for 0
      }
    
    if (p != 0)
      {
	cout << "The number is larger than the number of characters" << endl;
	abort();
      }
    
    return res;
  }

  
  //! conversion from a list to a vector
  template<class T>
  void Copy(list<T>& x, Vector<T>& y)
  {
    y.Reallocate(x.size());
    typename list<T>::iterator it; int nb = 0;
    for (it = x.begin(); it != x.end(); ++it)
      y(nb++) = *it;
    
  }
  
  
  //! conversion from a dense vector to a sparse vector
  template<class T0, class T1, class Allocator0, class Allocator1>
  void CopyVector(const Vector<T0, VectFull, Allocator0>& x,
		  Vector<T1, VectSparse, Allocator1>& y)
  {
    T0 zero; SetComplexZero(zero);
    // counting length of sparse vector
    int size_vec = 0;
    for (int i = 0; i < x.GetM(); i++)
      if (x(i) != zero)
	size_vec++;
    
    // filling y
    y.Clear();
    y.ReallocateVector(size_vec);
    size_vec = 0;
    for (int i = 0; i < x.GetM(); i++)
      if (x(i) != zero)
	{
	  y.Index(size_vec) = i;
	  y.Value(size_vec) = x(i);
	  size_vec++;
	}
  }  


  
  //! Linspace Matlab function
  template<class T>
  void Linspace(const T& a, const T& b, int N, Vector<T>& x)
  {
    if (N <= 0)
      {
        x.Clear();
        return;
      }
    
    x.Reallocate(N);
    if (N == 1)
      {
        x(0) = a;
        return;
      }
    
    T dx = (b-a)/(N-1);
    for (int i = 0; i < N-1; i++)
      x(i) = a + dx*i;
    
    x(N-1) = b;
  }


  //! returns the size n in a readable string
  string GetHumanReadableMemory(size_t taille)
  {
    string size_chaine;
    if (taille < 1024)
      size_chaine = to_str(taille) + " bytes";
    else if (taille < 1024*1024)
      {
        int n = (100*taille)/1024;
        size_chaine = to_str(double(n)/100) + " KiB";
      }
    else if (taille < 1024*1024*1024)
      {
        size_t n = (100*taille)/size_t(1024*1024);
        size_chaine = to_str(double(n)/100) + " MiB";
      }
    else if (taille/(1024*1024) < 1024*1024)
      {
        size_t n = (100*taille)/size_t(1024*1024*1024);
        size_chaine = to_str(double(n)/100) + " GiB";
      }
    else
      {
        size_t n = (100*taille)/int64_t(1024*1024*int64_t(1024*1024));
        size_chaine = to_str(double(n)/100) + " TiB";
      }
    
    return size_chaine;
  }
  
    
#ifdef SELDON_WITH_MPI
  //! broadcasting a string to all the processors of communicator
  void MPI_Bcast_string(string& s, int root, const MPI_Comm& comm)
  {
    int n = s.size();
    MPI_Bcast(&n, 1, MPI_INTEGER, root, comm);
    
    Vector<char> tab(n);
    int rank; MPI_Comm_rank(comm, &rank);
    if (rank == root)
      {
        for (int i = 0; i < n; i++)
          tab(i) = s[i];
      }
    else
      s.resize(n);
    
    MPI_Bcast(tab.GetData(), n, MPI_CHAR, root, comm);
    
    if (rank != root)
      for (int i = 0; i < n; i++)
        s[i] = tab(i);
  }

#ifdef MONTJOIE_WITH_MPFR
  void FillVectorMPI(mpfr::mpreal* x, int n, Vector<int64_t>& xtmp)
  {
    int p = (MONTJOIE_SIZE_MPFR+2)*n;
    xtmp.Reallocate(p);
    int nb = 0;
    mpfr_ptr xptr;
    for (int i = 0; i < n; i++)
      {
        xptr = x[i].mpfr_ptr();
        xtmp(nb++) = xptr->_mpfr_sign;
        xtmp(nb++) = xptr->_mpfr_exp;
        for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
          xtmp(nb++) = xptr->_mpfr_d[j];
      }
  }

  void FillVectorMPI(Vector<int64_t>& xtmp, int n, mpfr::mpreal* x)
  {
    int nb = 0;
    mpfr_ptr xptr;
    for (int i = 0; i < n; i++)
      {
        xptr = x[i].mpfr_ptr();
        xptr->_mpfr_sign = xtmp(nb++); 
        xptr->_mpfr_exp = xtmp(nb++);
        for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
          xptr->_mpfr_d[j] = xtmp(nb++);
      }
  }

  void FillVectorMPI(complex<mpfr::mpreal>* x, int n, Vector<int64_t>& xtmp)
  {
    int p = (MONTJOIE_SIZE_MPFR+2)*n*2;
    xtmp.Reallocate(p);
    int nb = 0;
    mpfr_ptr xptr;
    for (int i = 0; i < n; i++)
      {	
        xptr = reinterpret_cast<mpfr::mpreal*>(&x[i])[0].mpfr_ptr();
        xtmp(nb++) = xptr->_mpfr_sign;
        xtmp(nb++) = xptr->_mpfr_exp;
        for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
          xtmp(nb++) = xptr->_mpfr_d[j];

	xptr = reinterpret_cast<mpfr::mpreal*>(&x[i])[1].mpfr_ptr();
        xtmp(nb++) = xptr->_mpfr_sign;
        xtmp(nb++) = xptr->_mpfr_exp;
        for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
          xtmp(nb++) = xptr->_mpfr_d[j];
      }
  }

  void FillVectorMPI(Vector<int64_t>& xtmp, int n, complex<mpfr::mpreal>* x)
  {
    int nb = 0;
    mpfr_ptr xptr; mpfr::mpreal xr, xi;
    for (int i = 0; i < n; i++)
      {
        xptr = xr.mpfr_ptr();
        xptr->_mpfr_sign = xtmp(nb++); 
        xptr->_mpfr_exp = xtmp(nb++);
        for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
          xptr->_mpfr_d[j] = xtmp(nb++);

        xptr = xi.mpfr_ptr();
        xptr->_mpfr_sign = xtmp(nb++); 
        xptr->_mpfr_exp = xtmp(nb++);
        for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
          xptr->_mpfr_d[j] = xtmp(nb++);
        
        x[i] = complex<mpfr::mpreal>(xr, xi);
      }
  }

  void ReduceVectorMPI(const MPI_Comm& comm, const MPI_Op& op,
                       Vector<int64_t>& xtmp, int n, mpfr::mpreal* x)
  {
    int taille = (MONTJOIE_SIZE_MPFR+2)*n;
    mpfr_ptr xptr;
    mpfr::mpreal result, xval;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    for (int i = 0; i < n; i++)
      {
        for (int p = 0; p < nb_proc; p++)
          {
            int offset = taille*p + i*(MONTJOIE_SIZE_MPFR+2);
            xptr = xval.mpfr_ptr();
            xptr->_mpfr_sign = xtmp(offset); 
            xptr->_mpfr_exp = xtmp(offset + 1);
            for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
              xptr->_mpfr_d[j] = xtmp(offset + j + 2);
            
            if (p == 0)
              result = xval;
            else
              {
                if (op == MPI_MIN)
                  result = min(result, xval);
                else if (op == MPI_MAX)
                  result = max(result, xval);
                else if (op == MPI_SUM)
                  result += xval;
                else
                  {
                    cout << "not implemented" << endl;
                    abort();
                  }
              }
          }
        
        x[i] = result;
      }
  }

  void ReduceVectorMPI(const MPI_Comm& comm, const MPI_Op& op,
                       Vector<int64_t>& xtmp, int n, complex<mpfr::mpreal>* x)
  {
    int taille = (MONTJOIE_SIZE_MPFR+2)*n*2;
    mpfr_ptr xptr;
    mpfr::mpreal xr, xi; complex<mpfr::mpreal> result;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    for (int i = 0; i < n; i++)
      {
        for (int p = 0; p < nb_proc; p++)
          {
            int offset = taille*p + 2*i*(MONTJOIE_SIZE_MPFR+2);
            xptr = xr.mpfr_ptr();
            xptr->_mpfr_sign = xtmp(offset); 
            xptr->_mpfr_exp = xtmp(offset + 1);
            for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
              xptr->_mpfr_d[j] = xtmp(offset + j + 2);
            
            offset += MONTJOIE_SIZE_MPFR+2;
            xptr = xi.mpfr_ptr();
            xptr->_mpfr_sign = xtmp(offset); 
            xptr->_mpfr_exp = xtmp(offset + 1);
            for (int j = 0; j < MONTJOIE_SIZE_MPFR; j++)
              xptr->_mpfr_d[j] = xtmp(offset + j + 2);
            
            if (p == 0)
              result = complex<mpfr::mpreal>(xr, xi);
            else
              {
                if (op == MPI_SUM)
                  result += complex<mpfr::mpreal>(xr, xi);
                else
                  {
                    cout << "not implemented" << endl;
                    abort();
                  }
              }
          }
        
        x[i] = result;
      }
  }

  MPI_Request MpiIsend(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Request request;
    MPI_Isend(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, tag, comm, &request);
    return request;
  }

  MPI_Request MpiIsend(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    return MpiIsend(comm, x.GetData(), xtmp, n, proc, tag);
  }

  MPI_Request MpiIsend(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Request request;
    MPI_Isend(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, tag, comm, &request);
    return request;
  }

  MPI_Request MpiIsend(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    return MpiIsend(comm, x.GetData(), xtmp, n, proc, tag);
  }

  MPI_Request MpiIrecv(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    int p = (MONTJOIE_SIZE_MPFR+2)*n;
    xtmp.Reallocate(p);
    MPI_Request request;
    MPI_Irecv(xtmp.GetData(), p, MPI_INTEGER8, proc, tag, comm, &request);
    return request;
  }

  MPI_Request MpiIrecv(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    return MpiIrecv(comm, x.GetData(), xtmp, n, proc, tag);
  }

  MPI_Request MpiIrecv(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    int p = (MONTJOIE_SIZE_MPFR+2)*2*n;
    xtmp.Reallocate(p);
    MPI_Request request;
    MPI_Irecv(xtmp.GetData(), p, MPI_INTEGER8, proc, tag, comm, &request);
    return request;
  }

  MPI_Request MpiIrecv(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    return MpiIrecv(comm, x.GetData(), xtmp, n, proc, tag);
  }

  void MpiCompleteIrecv(mpfr::mpreal* x, Vector<int64_t>& xtmp, int n)
  {
    FillVectorMPI(xtmp, n, x);
  }

  void MpiCompleteIrecv(Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp, int n)
  {
    FillVectorMPI(xtmp, n, x.GetData());
  }

  void MpiCompleteIrecv(complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp, int n)
  {
    FillVectorMPI(xtmp, n, x);    
  }

  void MpiCompleteIrecv(Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp, int n)
  {
    FillVectorMPI(xtmp, n, x.GetData());    
  }

  void MpiSsend(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Send(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, tag, comm);
  }

  void MpiSsend(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    MpiSsend(comm, x.GetData(), xtmp, n, proc, tag);
  }

  void MpiSsend(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Send(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, tag, comm);
  }

  void MpiSsend(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    MpiSsend(comm, x.GetData(), xtmp, n, proc, tag);
  }

  void MpiSend(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Send(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, tag, comm);
  }

  void MpiSend(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag)
  {
    MpiSend(comm, x.GetData(), xtmp, n, proc, tag);
  }

  void MpiSend(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Send(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, tag, comm);
  }

  void MpiSend(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag)
  {
    MpiSend(comm, x.GetData(), xtmp, n, proc, tag);
  }

  void MpiGather(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                 mpfr::mpreal* y, int n, int proc)
  {
    Vector<int64_t> ytmp;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int p = (MONTJOIE_SIZE_MPFR+2)*n*nb_proc;
    if (rank_proc == proc)
      ytmp.Reallocate(p);
    
    FillVectorMPI(x, n, xtmp);
    MPI_Gather(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8,
               ytmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, comm);
    
    if (rank_proc == proc)
      FillVectorMPI(ytmp, n*nb_proc, y);
  }

  void MpiGather(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                 Vector<mpfr::mpreal>& y, int n, int proc)
  {
    MpiGather(comm, x.GetData(), xtmp, y.GetData(), n, proc);
  }

  void MpiGather(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                 complex<mpfr::mpreal>* y, int n, int proc)
  {
    Vector<int64_t> ytmp;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int p = (MONTJOIE_SIZE_MPFR+2)*n*2*nb_proc;
    if (rank_proc == proc)
      ytmp.Reallocate(p);
    
    FillVectorMPI(x, n, xtmp);
    MPI_Gather(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8,
               ytmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, comm);
    
    if (rank_proc == proc)
      FillVectorMPI(ytmp, n*nb_proc, y);
  }

  void MpiGather(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                 Vector<complex<mpfr::mpreal> >& y, int n, int proc)
  {
    MpiGather(comm, x.GetData(), xtmp, y.GetData(), n, proc);
  }

  void MpiReduce(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                 mpfr::mpreal* y, int n, const MPI_Op& op, int proc)
  {
    Vector<int64_t> ytmp;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int p = (MONTJOIE_SIZE_MPFR+2)*n*nb_proc;
    if (rank_proc == proc)
      ytmp.Reallocate(p);
    
    // regrouping all values in processor proc
    FillVectorMPI(x, n, xtmp);
    MPI_Gather(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8,
               ytmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, comm);
    
    // then performing reduction phase
    if (rank_proc == proc)
      ReduceVectorMPI(comm, op, ytmp, n, y);
  }

  void MpiReduce(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                 Vector<mpfr::mpreal>& y, int n, const MPI_Op& op, int proc)
  {
    MpiReduce(comm, x.GetData(), xtmp, y.GetData(), n, op, proc);
  }

  void MpiReduce(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                 complex<mpfr::mpreal>* y, int n, const MPI_Op& op, int proc)
  {
    Vector<int64_t> ytmp;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    int p = (MONTJOIE_SIZE_MPFR+2)*2*n*nb_proc;
    if (rank_proc == proc)
      ytmp.Reallocate(p);
    
    // regrouping all values in processor proc
    FillVectorMPI(x, n, xtmp);
    MPI_Gather(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8,
               ytmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, comm);
    
    // then performing reduction phase
    if (rank_proc == proc)
      ReduceVectorMPI(comm, op, ytmp, n, y);
  }

  void MpiReduce(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                 Vector<complex<mpfr::mpreal> >& y, int n, const MPI_Op& op, int proc)
  {
    MpiReduce(comm, x.GetData(), xtmp, y.GetData(), n, op, proc);
  }
  
  void MpiRecv(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status)
  {
    int p = (MONTJOIE_SIZE_MPFR+2)*n;
    xtmp.Reallocate(p);
    MPI_Recv(xtmp.GetData(), p, MPI_INTEGER8, proc, tag, comm, &status);
    FillVectorMPI(xtmp, n, x);
  }

  void MpiRecv(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status)
  {
    MpiRecv(comm, x.GetData(), xtmp, n, proc, tag, status);
  }

  void MpiRecv(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status)
  {
    int p = (MONTJOIE_SIZE_MPFR+2)*2*n;
    xtmp.Reallocate(p);
    MPI_Recv(xtmp.GetData(), p, MPI_INTEGER8, proc, tag, comm, &status);
    FillVectorMPI(xtmp, n, x);
  }

  void MpiRecv(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI_Status& status)
  {
    MpiRecv(comm, x.GetData(), xtmp, n, proc, tag, status);
  }

  void MpiBcast(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp, int n, int proc)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Bcast(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, comm);
    FillVectorMPI(xtmp, n, x);
  }

  void MpiBcast(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp, int n, int proc)
  {
    MpiBcast(comm, x.GetData(), xtmp, n, proc);
  }

  void MpiBcast(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp, int n, int proc)
  {
    FillVectorMPI(x, n, xtmp);
    MPI_Bcast(xtmp.GetData(), xtmp.GetM(), MPI_INTEGER8, proc, comm);
    FillVectorMPI(xtmp, n, x);
  }

  void MpiBcast(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp, int n, int proc)
  {
    MpiBcast(comm, x.GetData(), xtmp, n, proc);
  }

  void MpiAllreduce(const MPI_Comm& comm, mpfr::mpreal* x, Vector<int64_t>& xtmp,
                    mpfr::mpreal* y, int n, const MPI_Op& op)
  {
    // reducing the vector on processor, then broadcastin to other processors
    MpiReduce(comm, x, xtmp, y, n, op, 0);
    MpiBcast(comm, y, xtmp, n, 0);
  }

  void MpiAllreduce(const MPI_Comm& comm, Vector<mpfr::mpreal>& x, Vector<int64_t>& xtmp,
                    Vector<mpfr::mpreal>& y, int n, const MPI_Op& op)
  {
    MpiAllreduce(comm, x.GetData(), xtmp, y.GetData(), n, op);
  }

  void MpiAllreduce(const MPI_Comm& comm, complex<mpfr::mpreal>* x, Vector<int64_t>& xtmp,
                    complex<mpfr::mpreal>* y, int n, const MPI_Op& op)
  {
    // reducing the vector on processor, then broadcastin to other processors
    MpiReduce(comm, x, xtmp, y, n, op, 0);
    MpiBcast(comm, y, xtmp, n, 0);
  }

  void MpiAllreduce(const MPI_Comm& comm, Vector<complex<mpfr::mpreal> >& x, Vector<int64_t>& xtmp,
                    Vector<complex<mpfr::mpreal> >& y, int n, const MPI_Op& op)
  {
    MpiAllreduce(comm, x.GetData(), xtmp, y.GetData(), n, op);
  }
#endif

#endif


  //! displays error messages of direct solver
  /*!
    \param[in] ierr error code returned by SparseDirectSolver
   */
  void PrintFactorizationFailed(int type, int ierr)
  {
    if (type == SparseDirectSolver<double>::FACTO_OK)
      return;
    
    if (type == SparseDirectSolver<double>::STRUCTURALLY_SINGULAR_MATRIX)
      cout<<"Structurally singular matrix"<<endl;
    else if (type == SparseDirectSolver<double>::NUMERICALLY_SINGULAR_MATRIX)
      cout<<"Numerically singular matrix"<<endl;
    else if (type == SparseDirectSolver<double>::OUT_OF_MEMORY)
      {
	cout<<"Matrix too large, we are unable to factorize this matrix. \n";
	cout<<"Run the simulation on a machine with more memory or on more machines "<<endl;
      }
    else if (type == SparseDirectSolver<double>::INVALID_ARGUMENT)
      cout << "Invalid argument" << endl;
    else if (type == SparseDirectSolver<double>::INCORRECT_NUMBER_OF_ROWS)
      cout << "Number of rows of the matrix incorrect. Empty matrix ? " << endl;
    else if (type == SparseDirectSolver<double>::MATRIX_INDICES_INCORRECT)
      cout << "Row or column indices of the matrix incorrect. Check your matrix" << endl;
    else if (type == SparseDirectSolver<double>::INVALID_PERMUTATION)
      cout << "Invalidation permutation array" << endl;
    else if (type == SparseDirectSolver<double>::ORDERING_FAILED)
      cout << "Computation of the ordering failed" << endl;
    else if (type == SparseDirectSolver<double>::INTERNAL_ERROR)
      cout << "Internal error, check the documentation of the direct solver" << endl;
    
    cout<<"Factorization failed"<<endl;
    cout<<"The error code of the direct solver is "<<ierr<<endl;
    
    abort();
  }
    
}

namespace Montjoie
{ 
  
  //! returns minimum value of k.x on box [xmin, xmax] x [ymin, ymax]
  Real_wp GetMinimumDotProd(const R2& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax)
  {
    Real_wp res = k(0)*xmin + k(1)*ymin;
    res = min(res, k(0)*xmin + k(1)*ymax);
    res = min(res, k(0)*xmax + k(1)*ymin);
    res = min(res, k(0)*xmax + k(1)*ymax);
    return res;
  }


  //! returns minimum value of k.x on box [xmin, xmax] x [ymin, ymax] x [zmin, zmax]
  Real_wp GetMinimumDotProd(const R3& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax)
  {
    Real_wp res = k(0)*xmin + k(1)*ymin + k(2)*zmin;
    res = min(res, k(0)*xmin + k(1)*ymin + k(2)*zmax);
    res = min(res, k(0)*xmin + k(1)*ymax + k(2)*zmax);
    res = min(res, k(0)*xmin + k(1)*ymax + k(2)*zmin);
    res = min(res, k(0)*xmax + k(1)*ymin + k(2)*zmin);
    res = min(res, k(0)*xmax + k(1)*ymin + k(2)*zmax);
    res = min(res, k(0)*xmax + k(1)*ymax + k(2)*zmax);
    res = min(res, k(0)*xmax + k(1)*ymax + k(2)*zmin);
    return res;
  }


  //! returns maximum value of k.x on box [xmin, xmax] x [ymin, ymax]
  Real_wp GetMaximumDotProd(const R2& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax)
  {
    Real_wp res = k(0)*xmin + k(1)*ymin;
    res = max(res, k(0)*xmin + k(1)*ymax);
    res = max(res, k(0)*xmax + k(1)*ymin);
    res = max(res, k(0)*xmax + k(1)*ymax);
    return res;
  }


  //! returns minimum value of k.x on box [xmin, xmax] x [ymin, ymax] x [zmin, zmax]
  Real_wp GetMaximumDotProd(const R3& k, const Real_wp& xmin, const Real_wp& xmax,
                            const Real_wp& ymin, const Real_wp& ymax,
                            const Real_wp& zmin, const Real_wp& zmax)
  {
    Real_wp res = k(0)*xmin + k(1)*ymin + k(2)*zmin;
    res = max(res, k(0)*xmin + k(1)*ymin + k(2)*zmax);
    res = max(res, k(0)*xmin + k(1)*ymax + k(2)*zmax);
    res = max(res, k(0)*xmin + k(1)*ymax + k(2)*zmin);
    res = max(res, k(0)*xmax + k(1)*ymin + k(2)*zmin);
    res = max(res, k(0)*xmax + k(1)*ymin + k(2)*zmax);
    res = max(res, k(0)*xmax + k(1)*ymax + k(2)*zmax);
    res = max(res, k(0)*xmax + k(1)*ymax + k(2)*zmin);
    return res;
  }


  //! conversion from cartesian coordinates (x, y) to polar coordinates (r, theta)
  void CartesianToPolar(const Real_wp& x, const Real_wp& y, Real_wp& r, Real_wp& theta)
  {
    r = sqrt(x*x + y*y);
    theta = Real_wp(0);
    if (r > 0)
      {
        Real_wp arg = x/r;
        if (arg >= 1.0)
          theta = 0.0;
        else if (arg <= -1.0)
          theta = pi_wp;
        else
          {
            theta = acos(x/r);            
            if (y < 0)
              theta = -theta;
          }
      }
  }

  
  //! conversion from cartesian coordinates (x, y, z) to spherical coordinates(r, theta, phi)
  /*!
    x = r sin(teta) cos(phi)
    y = r sin(teta) sin(phi)
    z = r cos(teta)
    r = sqrt( x^2 + y^2 + z^2 )
    phi = arctan( y/x )     with phi between -pi/2 and pi/2
    teta = arccos(z/r)   with negative sign if sin(teta) is negative
   */
  void CartesianToSpherical(const Real_wp& x, const Real_wp& y, const Real_wp& z,
                            Real_wp& r, Real_wp& teta, Real_wp& phi,
                            Real_wp& cos_teta, Real_wp& sin_teta)
  {
    r = sqrt(x*x + y*y + z*z);
    if (r == 0)
      {
	teta = Real_wp(0); phi = Real_wp(0);
	sin_teta = Real_wp(0); cos_teta = Real_wp(1);
	return;
      }
    
    if (abs(x) >= epsilon_machine)
      {
        phi = atan(y/x);
        sin_teta = x/(r*cos(phi));
      }
    else
      {
        phi = pi_wp/2;
        sin_teta = y/r;
      }
    
    cos_teta = z/r;
    teta = acos(cos_teta);
    if (sin_teta < -epsilon_machine)
      teta = -teta;
    else if (abs(sin_teta) <= epsilon_machine)
      {
        sin_teta = epsilon_machine;
        if (teta > pi_wp/2)
          teta = pi_wp-epsilon_machine;
        else
          teta = epsilon_machine;
      }

  }


  //! conversion from cartesian coordinates (x, y, z) to spherical coordinates(r, theta, phi)
  /*!
    x = r sin(teta) cos(phi)
    y = r sin(teta) sin(phi)
    z = r cos(teta)
    r = sqrt( x^2 + y^2 + z^2 )
    phi = arctan( y/x )     with phi between -pi/2 and pi/2
    teta = arccos(z/r)   with negative sign if sin(teta) is negative
   */
  void CartesianToSpherical(const Real_wp& x, const Real_wp& y, const Real_wp& z,
                            Real_wp& r, Real_wp& teta, Real_wp& phi)
  {
    Real_wp cos_teta, sin_teta;
    CartesianToSpherical(x, y, z, r, teta, phi, cos_teta, sin_teta);
  }
  

  //! conversion from spherical coordinates (r, theta, phi) to cartesian coordinates (x, y, z)
  void SphericalToCartesian(const Real_wp& r, const Real_wp& teta, const Real_wp& phi,
                            Real_wp& x, Real_wp& y, Real_wp& z)
  {
    Real_wp sin_teta = sin(teta);
    x = r*sin_teta*cos(phi);
    y = r*sin_teta*sin(phi);
    z = r*cos(teta);
  }

  
  //! computes coefficients of Pade expansion of square root
  template<class T>
  void ComputePadeCoefficientsSqrt(const T& alpha, int n, complex<T>& CoefC0,
				   Vector<complex<T> >& CoefAl, Vector<complex<T> >& CoefBl)
  {
    // usual Pade coefficients a_l, b_l
    T teta = pi_wp/(2*n+1);
    Vector<T> coef_al(n), coef_bl(n);
    for (int l = 0; l < n; l++)
      {
	T sin_sq = square(sin((l+1)*teta));
	coef_bl(l) = 1.0 - sin_sq;
	coef_al(l) = 2*sin_sq/(2*n+1);
      }
    
    complex<T> exp_ialpha = exp(Iwp*alpha);
    complex<T> exp_ialpha05 = exp(0.5*Iwp*alpha);
    
    // coefficient C_0 = e(i alpha/2) R_N ( e^{-i \alpha} -1)
    complex<T> arg = conj(exp_ialpha) - 1.0;
    complex<T> Rn(1, 0);
    for (int l = 0; l < n; l++)
      Rn += coef_al(l)*arg / (1.0 + coef_bl(l)*arg);
    
    CoefC0 = exp_ialpha05*Rn;
    
    // coefficient A_l = exp(-i alpha/2) a_l / (1 + b_l (exp(-i alpha) - 1))^2
    // coefficient B_l = exp(-i alpha) b_l / (1 + b_l (exp(-i alpha) - 1))
    CoefAl.Reallocate(n);
    CoefBl.Reallocate(n);
    for (int l = 0; l < n; l++)
      {
	CoefAl(l) = conj(exp_ialpha05)*coef_al(l) / square(1.0 + coef_bl(l)*arg);
	CoefBl(l) = conj(exp_ialpha)*coef_bl(l) / (1.0 + coef_bl(l)*arg);
      }
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_COMMON_CXX
#endif
