#ifndef MONTJOIE_FILE_PARAMETER_OUTPUT_REPRISE_CXX

namespace Montjoie
{

  ParameterOutputReprise::ParameterOutputReprise()
  {
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#else
    rank_proc = 0;
#endif
    master_proc = rank_proc;
    nb_procs_slave = 1;
  }
  

  void ParameterOutputReprise
  ::RegroupWritingOnMaster(int nb_proc_per_node)
  {
    master_proc = rank_proc - rank_proc%nb_proc_per_node;
    nb_procs_slave = nb_proc_per_node;
  }
  
  
  void ParameterOutputReprise::StartReading(const string& input_file)
  {
    if (master_proc == rank_proc)
      {
        file_in.open(input_file.data());
        
        if (!file_in.is_open())
          {
            cout << "Unable to open file " << input_file << endl;
            abort();
          }
      }
  }
  
  
  void ParameterOutputReprise::Read(VectReal_wp& x, bool master_only)
  {
    if (master_proc == rank_proc)
      {
        x.Read(file_in);
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              VectReal_wp y; Vector<int64_t> ytmp;
              y.Read(file_in);
              
              int n = y.GetM();
              MPI_Ssend(&n, 1, MPI_INTEGER,
                        master_proc+i, 62, MPI_COMM_WORLD);
              
              if (n > 0)
                MpiSsend(MPI_COMM_WORLD, y, ytmp, y.GetM(), master_proc+i, 63);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        MPI_Status status;
        int n; Vector<int64_t> xtmp;
        MPI_Recv(&n, 1, MPI_INTEGER,
                 master_proc, 62, MPI_COMM_WORLD, &status);
        
        x.Reallocate(n);
        if (n > 0)
          MpiRecv(MPI_COMM_WORLD, x, xtmp, n, master_proc, 63, status);
#endif
      }
  }
  
  
  void ParameterOutputReprise::Read(Vector<int>& x, bool master_only)
  {
    if (master_proc == rank_proc)
      {
        x.Read(file_in);
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              Vector<int> y;
              y.Read(file_in);
              
              int n = y.GetM();
              MPI_Ssend(&n, 1, MPI_INTEGER,
                        master_proc+i, 64, MPI_COMM_WORLD);
              
              if (n > 0)
                MPI_Ssend(y.GetData(), y.GetM(), MPI_INTEGER,
                          master_proc+i, 65, MPI_COMM_WORLD);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        MPI_Status status;
        int n;
        MPI_Recv(&n, 1, MPI_INTEGER,
                 master_proc, 64, MPI_COMM_WORLD, &status);
        
        x.Reallocate(n);
        if (n > 0)
          MPI_Recv(x.GetData(), n, MPI_INTEGER,
                   master_proc, 65, MPI_COMM_WORLD, &status);
#endif
      }
  }
  
  
  template<int p>
  void ParameterOutputReprise::Read(Vector<TinyVector<Real_wp, p> >& x, bool master_only)
  {
    if (master_proc == rank_proc)
      {
        x.Read(file_in);
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              Vector<TinyVector<Real_wp, p> > y;
              y.Read(file_in);
              
              int n = y.GetM()*p;
              MPI_Ssend(&n, 1, MPI_INTEGER,
                        master_proc+i, 66, MPI_COMM_WORLD);
              
              if (n > 0)
                {
                  Vector<int64_t> ytmp;
                  MpiSsend(MPI_COMM_WORLD, reinterpret_cast<Real_wp*>(y.GetData()), ytmp,
                           n, master_proc+i, 67);
                }
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        MPI_Status status;
        int n; Vector<int64_t> xtmp;
        MPI_Recv(&n, 1, MPI_INTEGER,
                 master_proc, 66, MPI_COMM_WORLD, &status);
        
        x.Reallocate(n/p);
        if (n > 0)
          MpiRecv(MPI_COMM_WORLD, reinterpret_cast<Real_wp*>(x.GetData()),
                  xtmp, n, master_proc, 67, status);
#endif
      }
  }
  
  
  void ParameterOutputReprise::Read(Matrix<Real_wp>& A, bool master_only)
  {
    if (master_proc == rank_proc)
      {
        A.Read(file_in);
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              Matrix<Real_wp> B; Vector<int64_t> xtmp;
              B.Read(file_in);
              
              IVect sizeB(2);
              int m = B.GetM(), n = B.GetN();
              sizeB(0) = m; sizeB(1) = n;
              MPI_Ssend(sizeB.GetData(), 2, MPI_INTEGER,
                        master_proc+i, 68, MPI_COMM_WORLD);
              
              if (m*n > 0)
                MpiSsend(MPI_COMM_WORLD, B.GetData(), xtmp, m*n, master_proc+i, 69);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        MPI_Status status;
        IVect sizeB(2); Vector<int64_t> xtmp;
        MPI_Recv(sizeB.GetData(), 2, MPI_INTEGER,
                 master_proc, 68, MPI_COMM_WORLD, &status);
        
        int m = sizeB(0), n = sizeB(1);
        A.Reallocate(m, n);
        if (m*n > 0)
          MpiRecv(MPI_COMM_WORLD, A.GetData(), xtmp, m*n, master_proc, 69, status);
#endif
      }
  }
  
  
  void ParameterOutputReprise::Read(int& n, bool master_only)
  {
    if (master_proc == rank_proc)
      {
        file_in.read(reinterpret_cast<char*>(&n), sizeof(int));
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              int m;
              file_in.read(reinterpret_cast<char*>(&m), sizeof(int));
              
              MPI_Ssend(&m, 1, MPI_INTEGER,
                        master_proc+i, 70, MPI_COMM_WORLD);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        MPI_Status status;
        MPI_Recv(&n, 1, MPI_INTEGER,
                 master_proc, 70, MPI_COMM_WORLD, &status);
#endif
      }

  }
  
  
  void ParameterOutputReprise::FinalizeReading()
  {
    if (master_proc == rank_proc)
      file_in.close();
  }
  
  
  void ParameterOutputReprise::StartWriting(const string& output_file)
  {
    if (master_proc == rank_proc)
      {
        file_out.open(output_file.data());
        
        if (!file_out.is_open())
          {
            cout << "Unable to open file " << output_file << endl;
            abort();
          }
      }
  }
  
  
  void ParameterOutputReprise::Write(const VectReal_wp& x, bool master_only)
  {
    if (master_proc == rank_proc)
      {        
        x.Write(file_out);
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              MPI_Status status;
              int n; VectReal_wp y; Vector<int64_t> xtmp;
              MPI_Recv(&n, 1, MPI_INTEGER,
                       master_proc+i, 71, MPI_COMM_WORLD, &status);
              
              y.Reallocate(n);
              if (n > 0)
                MpiRecv(MPI_COMM_WORLD, y, xtmp, n, master_proc+i, 72, status);
              
              y.Write(file_out);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        int n = x.GetM(); Vector<int64_t> xtmp;
        MPI_Ssend(&n, 1, MPI_INTEGER,
                  master_proc, 71, MPI_COMM_WORLD);
        
        if (n > 0)
          MpiSsend(MPI_COMM_WORLD, x.GetData(), xtmp, n, master_proc, 72);
#endif
      }    
  }
  
  
  void ParameterOutputReprise::Write(const Vector<int>& x, bool master_only)
  {
    if (master_proc == rank_proc)
      {        
        x.Write(file_out);

#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              MPI_Status status;
              int n; VectReal_wp y;
              MPI_Recv(&n, 1, MPI_INTEGER,
                       master_proc+i, 73, MPI_COMM_WORLD, &status);
              
              y.Reallocate(n);
              if (n > 0)
                MPI_Recv(y.GetData(), n, MPI_INTEGER,
                         master_proc+i, 74, MPI_COMM_WORLD, &status);
              
              y.Write(file_out);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        int n = x.GetM();
        MPI_Ssend(&n, 1, MPI_INTEGER,
                  master_proc, 73, MPI_COMM_WORLD);
        
        if (n > 0)
          MPI_Ssend(x.GetData(), n, MPI_INTEGER,
                    master_proc, 74, MPI_COMM_WORLD);
#endif
      }
  }
  
  
  template<int p>
  void ParameterOutputReprise::Write(const Vector<TinyVector<Real_wp, p> >& x,
                                     bool master_only)
  {
    if (master_proc == rank_proc)
      {        
        x.Write(file_out);
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              MPI_Status status;
              int n; Vector<TinyVector<Real_wp, p> > y;
              Vector<int64_t> xtmp;
              MPI_Recv(&n, 1, MPI_INTEGER,
                       master_proc+i, 75, MPI_COMM_WORLD, &status);
              
              y.Reallocate(n/p);
              if (n > 0)
                MpiRecv(MPI_COMM_WORLD, reinterpret_cast<Real_wp*>(y.GetData()), 
                        xtmp, n, master_proc+i, 76, status);
              
              y.Write(file_out);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        int n = x.GetM()*p; Vector<int64_t> xtmp;
        MPI_Ssend(&n, 1, MPI_INTEGER,
                  master_proc, 75, MPI_COMM_WORLD);
        
        if (n > 0)
          MpiSsend(MPI_COMM_WORLD, reinterpret_cast<Real_wp*>(x.GetData()),
                   xtmp, n, master_proc, 76);
#endif
      }    
  }
  
  
  void ParameterOutputReprise::Write(const Matrix<Real_wp>& A, bool master_only)
  {
    if (master_proc == rank_proc)
      {        
        A.Write(file_out);
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              MPI_Status status;
              IVect sizeB(2); Vector<int64_t> xtmp;
              int m, n; Matrix<Real_wp> B;
              MPI_Recv(sizeB.GetData(), 2, MPI_INTEGER,
                       master_proc+i, 37, MPI_COMM_WORLD, &status);
              
              m = sizeB(0); n = sizeB(1);
              B.Reallocate(m, n);
              if (m*n > 0)
                MpiRecv(MPI_COMM_WORLD, B.GetData(), xtmp, m*n,
                        master_proc+i, 38, status);
              
              B.Write(file_out);
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        IVect sizeA(2); Vector<int64_t> xtmp;
        int m = A.GetM(), n = A.GetN();
        sizeA(0) = m; sizeA(1) = n;        
        MPI_Ssend(sizeA.GetData(), 2, MPI_INTEGER,
                  master_proc, 37, MPI_COMM_WORLD);
        
        if (m*n > 0)
          MpiSsend(MPI_COMM_WORLD, A.GetData(), 
                   xtmp, m*n, master_proc, 38);
#endif
      }    
  }
  
  
  void ParameterOutputReprise::Write(int n, bool master_only)
  {
    if (master_proc == rank_proc)
      {        
        file_out.write(reinterpret_cast<char*>(&n), sizeof(int));
        
#ifdef SELDON_WITH_MPI
        if (!master_only)
          for (int i = 1; i < nb_procs_slave; i++)
            {
              int m; MPI_Status status;
              MPI_Recv(&m, 1, MPI_INTEGER,
                       master_proc+i, 39, MPI_COMM_WORLD, &status);
              
              file_out.write(reinterpret_cast<char*>(&m), sizeof(int));
            }
#endif
      }
    else
      {
#ifdef SELDON_WITH_MPI
        MPI_Ssend(&n, 1, MPI_INTEGER,
                  master_proc, 39, MPI_COMM_WORLD);
#endif
      }    
        
  }
  
  
  void ParameterOutputReprise::FinalizeWriting()
  {
    if (master_proc == rank_proc)
      file_out.close();
  }
  
}

#define MONTJOIE_FILE_PARAMETER_OUTPUT_REPRISE_CXX
#endif
