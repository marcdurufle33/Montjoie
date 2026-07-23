#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

Real_wp threshold;
int rank_processor, nb_processors;

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  int type = rand()%3;
  //int type = 2;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int m)
{
  x.Reallocate(m);
  for (int i = 0; i < m; i++)
    GetRand(x(i));  
}


void GetParallelRowNumbers(int m, int nb_proc, Vector<IVect>& list_proc,
			   IVect& GlobalRowNumbers, IVect& OverlapRowNumbers,
			   IVect& OverlapProcNumbers, IVect& ProcSharingRows,
			   Vector<IVect>& SharingRowNumbers,
			   const MPI_Comm& comm)
{
  // for each dof determining processors that share this dof
  list_proc.Reallocate(m);
  for (int i = 0; i < m; i++)
    {
      int nb = rand()%nb_proc+1;
      
      IVect num(nb);
      for (int j = 0; j < nb; j++)
        num(j) = rand()%nb_processors;
      
      RemoveDuplicate(num);
      nb = num.GetM();
      list_proc(i) = num;
    }
  
  // constructing GlobalRowNumbers
  int nodl = 0, nodl_overlap = 0;
  IVect Glob_to_local(m); Glob_to_local.Fill(-1);
  Vector<int> OverlappedGlobal(m); OverlappedGlobal.Fill(-1);
  IVect NbSharedDofPerProc(nb_processors); NbSharedDofPerProc.Fill(0);
  for (int i = 0; i < list_proc.GetM(); i++)
    for (int j = 0; j < list_proc(i).GetM(); j++)
      if (list_proc(i)(j) == rank_processor)
        {
          Glob_to_local(i) = nodl;
          if (j > 0)
            {
              OverlappedGlobal(i) = list_proc(i)(0);
              nodl_overlap++;
            }
          
          for (int k = 0; k < list_proc(i).GetM(); k++)
            if (j != k)
              NbSharedDofPerProc(list_proc(i)(k))++;
          
          nodl++;
        }
  
  GlobalRowNumbers.Reallocate(nodl);
  for (int i = 0; i < m; i++)
    if (Glob_to_local(i) >= 0)
      GlobalRowNumbers(Glob_to_local(i)) = i;
  
  //DISP(nodl); DISP(GlobalRowNumbers);
  
  // then OverlapRowNumbers/OverlapProcNumbers
  OverlapRowNumbers.Reallocate(nodl_overlap);
  OverlapProcNumbers.Reallocate(nodl_overlap);
  nodl_overlap = 0;
  for (int i = 0; i < m; i++)
    if (OverlappedGlobal(i) >= 0)
      {
        OverlapRowNumbers(nodl_overlap) = Glob_to_local(i);
        OverlapProcNumbers(nodl_overlap) = OverlappedGlobal(i);
        nodl_overlap++;
      }
  
  //DISP(OverlapRowNumbers); DISP(OverlapProcNumbers);
  // and finally ProcSharingRows, SharingRowNumbers
  nb_proc = 0;
  IVect IndexProc(nb_processors); IndexProc.Fill(-1);
  for (int p = 0; p < nb_processors; p++)
    if (NbSharedDofPerProc(p) > 0)
      {
        IndexProc(p) = nb_proc;
        nb_proc++;
      }
  
  ProcSharingRows.Reallocate(nb_proc);
  SharingRowNumbers.Reallocate(nb_proc);
  for (int p = 0; p < nb_processors; p++)
    if (IndexProc(p) >= 0)
      {
        ProcSharingRows(IndexProc(p)) = p;
        SharingRowNumbers(IndexProc(p)).Reallocate(NbSharedDofPerProc(p));
        SharingRowNumbers(IndexProc(p)).Fill(-1);
      }
  
  NbSharedDofPerProc.Fill(0);
  for (int i = 0; i < list_proc.GetM(); i++)
    for (int j = 0; j < list_proc(i).GetM(); j++)
      if (list_proc(i)(j) == rank_processor)
        {
          for (int k = 0; k < list_proc(i).GetM(); k++)
            if (j != k)
              {
                int cpt = NbSharedDofPerProc(list_proc(i)(k));
                int proc = IndexProc(list_proc(i)(k));
                SharingRowNumbers(proc)(cpt) = Glob_to_local(i);
                NbSharedDofPerProc(list_proc(i)(k))++;
              }
        }
  
  for (int i = 0; i < nb_proc; i++)
    {
      //DISP(i); 
      //DISP(ProcSharingRows(i));
      //DISP(SharingRowNumbers(i));
    }

}


template<class T>
void CheckDistributedVector(Vector<T>&, MPI_Comm comm, string suffix)
{
  int m = 50, nb_proc = min(nb_processors, 4);
  Vector<T> xref, yref;
  GenerateRandomVector(xref, m);
  GenerateRandomVector(yref, m);

  IVect GlobalRowNumbers, OverlapRowNumbers, OverlapProcNumbers, ProcSharingRows;
  Vector<IVect> SharingRowNumbers, list_proc;

  GetParallelRowNumbers(m, nb_proc, list_proc, GlobalRowNumbers, OverlapRowNumbers,
			OverlapProcNumbers, ProcSharingRows, SharingRowNumbers, comm);
  
  int n = GlobalRowNumbers.GetM();
  Vector<T> x(n), y(n);
  for (int i = 0; i < n; i++)
    {
      x(i) = xref(GlobalRowNumbers(i));
      y(i) = yref(GlobalRowNumbers(i));
    }

  // testing SetData
  DistributedVector<T> xd(OverlapRowNumbers, comm);
  DistributedVector<T> yd(OverlapRowNumbers, comm);
  xd.SetData(x); yd.SetData(y);

  // testing GetNbOverlap/GetOverlapRow
  if (xd.GetNbOverlap() != OverlapRowNumbers.GetM())
    {
      cout << "GetNbOverlap incorrect" << endl;
      abort();
    }

  for (int i = 0; i < xd.GetNbOverlap(); i++)
    if (xd.GetOverlapRow(i) != OverlapRowNumbers(i))
      {
	cout << "GetOverlapRow incorrect" << endl;
	abort();
      }

  // testing GetCommunicator
  if (xd.GetCommunicator() != comm)
    {
      cout << "GetCommunicator incorrect" << endl;
      abort();
    }

  // testing DotProd
  T scal_ref = DotProd(xref, yref);
  T scal = DotProd(xd, yd);

  if ((abs(scal-scal_ref) > threshold) || isnan(abs(scal-scal_ref)))
    {
      cout << "DotProd incorrect" << endl;
      abort();
    }

  // testing DotProdConj
  scal_ref = DotProdConj(xref, yref);
  scal = DotProdConj(xd, yd);

  if ((abs(scal-scal_ref) > threshold) || isnan(abs(scal-scal_ref)))
    {
      cout << "DotProdConj incorrect" << endl;
      abort();
    }

  // testing Norm2
  scal_ref = Norm2(xref);
  scal = Norm2(xd);

  if ((abs(scal-scal_ref) > threshold) || isnan(abs(scal-scal_ref)))
    {
      cout << "Norm2 incorrect" << endl;
      abort();
    }

  // testing AssembleVectorMin
  Vector<int> n_ref(m), proc_ref(m);
  Vector<int> num(n), num_proc(n);
  num.Fill(-1); num_proc.Fill(-1);
  int nb = 0;
  for (int i = 0; i < m; i++)
    {
      int p = 0, j = 0, pmin = RAND_MAX, jmin = RAND_MAX;
      for (int k = 0; k < list_proc(i).GetM(); k++)
	{
	  int proc = list_proc(i)(k);
	  p = rand()%nb_processors;
	  j = rand();
	  if (p < pmin)
	    {
	      pmin = p;
	      jmin = j;
	    }
	  else if (p == pmin)
	    jmin = min(jmin, j);
	  
	  if (proc == rank_processor)
	    {
	      num(nb) = j;
	      num_proc(nb) = p;
	      nb++;
	    }
	}

      n_ref(i) = jmin;
      proc_ref(i) = pmin;
    }
  
  AssembleVectorMin(num, num_proc, ProcSharingRows, SharingRowNumbers,
		    comm, n, 1, 12);

  for (int i = 0; i < n; i++)
    if ( (num(i) != n_ref(GlobalRowNumbers(i))) ||
	 (num_proc(i) != proc_ref(GlobalRowNumbers(i))) )
      {
	DISP(i); DISP(num(i)); DISP(n_ref(GlobalRowNumbers(i)));
	DISP(num_proc(i)); DISP(proc_ref(GlobalRowNumbers(i)));
	cout << "AssembleVectorMin incorrect" << endl;
	abort();
      }

  // testing AssembleVector
  nb = 0;
  for (int i = 0; i < m; i++)
    {
      T val, val_ref; SetComplexZero(val_ref);
      for (int k = 0; k < list_proc(i).GetM(); k++)
	{
	  int proc = list_proc(i)(k);
	  GetRand(val);
	  val_ref += val;
	  if (rank_processor == proc)
	    {
	      xd(nb) = val;
	      nb++;
	    }
	}

      xref(i) = val_ref;
    }

  AssembleVector(xd, MPI_SUM, ProcSharingRows, SharingRowNumbers,
		 comm, n, 1, 12);

  for (int i = 0; i < n; i++)
    if ( ( abs(xd(i) - xref(GlobalRowNumbers(i))) > threshold) ||
	 isnan(abs(xd(i) - xref(GlobalRowNumbers(i)))))
      {
	cout << "AssembleVector incorrect" << endl;
	abort();
      }
  

  // testing ExchangeVector
  nb = 0;
  for (int i = 0; i < m; i++)
    {
      Vector<T> val_(list_proc(i).GetM());
      T val_ref; SetComplexZero(val_ref);
      for (int k = 0; k < list_proc(i).GetM(); k++)
	{
	  int proc = list_proc(i)(k);
	  GetRand(val_(k));
	  if (rank_processor != proc)
	    val_ref = val_(k);
	}

      if (list_proc(i).GetM() == 1)
	val_ref = val_(0);

      for (int k = 0; k < list_proc(i).GetM(); k++)
	{
	  int proc = list_proc(i)(k);
	  if (rank_processor == proc)
	    {
	      xd(nb) = val_(k);
	      yd(nb) = val_ref;
	      nb++;
	    }
	}
    }

  ExchangeVector(xd, ProcSharingRows, SharingRowNumbers,
		 comm, n, 1, 12);

  for (int i = 0; i < n; i++)
    if ( ( abs(xd(i) - yd(i)) > threshold) ||
	 isnan(abs(xd(i) - yd(i))))
      {
	cout << "ExchangeVector incorrect" << endl;
	abort();
      }

  // testing ExchangeRelaxVector
  for (int p = 0; p < nb_processors; p++)
    {
      Real_wp omega;
      GetRand(omega);

      nb = 0;
      for (int i = 0; i < m; i++)
	{
	  T val;
	  T val_ref; GetRand(val_ref);

	  bool p_present = false;
	  for (int k = 0; k < list_proc(i).GetM(); k++)
	    if (list_proc(i)(k) == p)
	      p_present = true;
	  
	  //DISP(i); DISP(list_proc(i)); DISP(nb);
	  if (list_proc(i).GetM() == 1)
	    {
	      if (list_proc(i)(0) == rank_processor)
		{
		  xd(nb) = val_ref;
		  yd(nb) = val_ref;
		  nb++;
		}
	    }
	  else
	    {	      
	      for (int k = 0; k < list_proc(i).GetM(); k++)
		{
		  int proc = list_proc(i)(k);
		  GetRand(val);
		  if (rank_processor == proc)
		    {
		      if (p != proc)
			{
			  if (p_present)
			    {
			      xd(nb) = val;
			      yd(nb) = (Real_wp(1)-omega)*val + omega*val_ref;
			    }
			  else
			    {
			      xd(nb) = val;
			      yd(nb) = val;
			    }
			}
		      else
			{
			  xd(nb) = val_ref;
			  yd(nb) = val_ref;
			}
		      
		      nb++;
		    }
		}
	    }
	}

      ExchangeRelaxVector(xd, omega, p, ProcSharingRows, SharingRowNumbers,
			  comm, n, 1, 12);

      for (int i = 0; i < n; i++)
	if ( ( abs(xd(i) - yd(i)) > threshold) ||
	     isnan(abs(xd(i) - yd(i))))
	  {
	    DISP(i); DISP(xd(i)); DISP(yd(i));
	    cout << "ExchangeRelaxVector incorrect" << endl;
	    abort();
	  }
    }
}


void SplitWorldCommunicator(MPI_Comm& comm, string& suffix)
{
  if (false)
    {
      // no split
      comm = MPI_COMM_WORLD;
      suffix = string("");
    }
  else
    {
      // splitting in two communicators
      int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
      if (nb_proc < 4)
	{
	  cout << "This unitary test must be launched on four processors at least" << endl;
	  abort();
	}
      
      int color = 0;
      int key = rank;
      suffix = "_C0";
      if (rank >= nb_proc/2)
	{
	  color = 1;
	  key -= nb_proc/2;
	  suffix = "_C1";
	}
      
      MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
    }
  
  MPI_Comm_rank(comm, &rank_processor);
  MPI_Comm_size(comm, &nb_processors);
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  MPI_Comm comm; string suffix;
  SplitWorldCommunicator(comm, suffix);

  threshold = 1e4*epsilon_machine;
  cout.precision(15);

  {
    Vector<Real_wp> V;
    CheckDistributedVector(V, comm, suffix);
  }

  {
    Vector<Complex_wp> V;
    CheckDistributedVector(V, comm, suffix);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank_processor == 0)
    cout << "All tests passed successfully" << endl;

  return FinalizeMontjoie();
}
