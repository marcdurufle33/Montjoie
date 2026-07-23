// file to be included => do not compile this file

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

void GenerateRandomPermutation(int n, IVect& permut)
{
  Vector<bool> NumUsed(n);
  NumUsed.Fill(false);
  permut.Reallocate(n);
  permut.Fill(-1);
  int nb = 0;
  // premiere iteration
  for (int i = 0; i < n; i++)
    {
      int i2 = rand()%n;
      if (!NumUsed(i2))
        {
          NumUsed(i2) = true;
          permut(i) = i2;
          nb++;
        }
    }
  
  while (nb < n)
    {
      // on recupere les numeros non-selectionnes
      IVect non_selec(n-nb);
      int k = 0;
      for (int i = 0; i < n; i++)
        if (!NumUsed(i))
          non_selec(k++) = i;
      
      // iteration suivante
      for (int i = 0; i < n; i++)
        if (permut(i) == -1)
          {
            int i2 = rand()%(n-nb);
            int j = non_selec(i2);
            if (!NumUsed(j))
              {
                NumUsed(j) = true;
                permut(i) = j;
                nb++;
              }
          }
    }
}

void PickRandomInteger(int nb, IVect& num)
{
  if (nb == num.GetM())
    return;
  
  IVect permut;
  GenerateRandomPermutation(num.GetM(), permut);
  
  IVect old_num(num);
  num.Reallocate(nb);
  for (int i = 0; i < nb; i++)
    num(i) = old_num(permut(i));
}

template<class T>
void PickRandomPartition(const T& x, int n, Vector<T>& eval)
{
  eval.Reallocate(n);
  double zero, sum; SetComplexZero(zero);
  sum = zero; Vector<double> coef(n);
  for (int i = 0; i < n; i++)
    {
      GetRand(coef(i));
      sum += coef(i);
    }
  
  Mlt(1.0/sum, coef);
  for (int i = 0; i < n; i++)
    eval(i) = coef(i)*x;
}

template<class MatrixSparse, class T>
void AddInteraction(MatrixSparse& A, int i, int j, const T& x,
                    const IVect& Glob_to_local, int proc_row, int proc_col)
{
  int iloc = Glob_to_local(i);
  int jloc = Glob_to_local(j);
  bool sym = IsSymmetricMatrix(A);
  if (i == j)
    sym = false;
  
  if (iloc >= 0)
    {
      if (jloc >= 0)
        {
          A.AddInteraction(iloc, jloc, x);
          if (sym)
            A.AddInteraction(jloc, iloc, x);
        }
      else
        {
          A.AddDistantInteraction(iloc, j, proc_col, x); 
          if (sym)
            A.AddRowDistantInteraction(j, iloc, proc_col, x); 
        }
    }
  else
    {
      A.AddRowDistantInteraction(i, jloc, proc_row, x);
      if (sym)
        A.AddDistantInteraction(jloc, i, proc_row, x);
    }
}

template<class MatrixSeq, class MatrixPar>
void DistributeMatrixProcessor(const MatrixSeq& Aref, MatrixPar& A, int rank_processor,
                               const Vector<IVect>& list_proc, const IVect& Glob_to_local)
{
  int n = A.GetM();
  A.Clear(); A.Reallocate(n, n);
  typedef typename MatrixSeq::entry_type T;
  T zero; SetComplexZero(zero);
  for (int i = 0; i < Aref.GetM(); i++)
    {
      int jlow = 0;
      if (IsSymmetricMatrix(Aref))
        jlow = i;
      
      for (int j = jlow; j < Aref.GetM(); j++)
        if (Aref(i, j) != zero)
          {
            IVect num(list_proc(i).GetM() + list_proc(j).GetM());
            for (int k = 0; k < list_proc(i).GetM(); k++)
              num(k) = list_proc(i)(k);
            
            for (int k = 0; k < list_proc(j).GetM(); k++)
              num(list_proc(i).GetM() + k) = list_proc(j)(k);          
            
            RemoveDuplicate(num);
            
            int nb = rand()%num.GetM() + 1;
            PickRandomInteger(nb, num);
            
            IVect num_row(nb), num_col(nb);
            for (int k = 0; k < nb; k++)
              {
                num_row(k) = rand()%list_proc(i).GetM();
                num_col(k) = rand()%list_proc(j).GetM();
              }
            
            Vector<T> eval;
            PickRandomPartition(Aref(i, j), nb, eval);
            
            for (int k = 0; k < nb; k++)
            if (rank_processor == num(k))
              AddInteraction(A, i, j, eval(k), Glob_to_local,
                             list_proc(i)(num_row(k)), list_proc(j)(num_col(k)));
          }
    }
}

template<class T, class Prop, class Storage>
void ScatterMatrix(Vector<DistributedMatrix<T, Prop, Storage> >& Ai,
                   DistributedMatrixIntegerArray& info,
                   int nb_processors, int rank_processor, MPI_Comm comm)
{
  int m = Ai(0).GetM();
  
  // for each dof determining processors that share this dof
  Vector<IVect> list_proc(m);
  for (int i = 0; i < m; i++)
    {
      int nb = rand()%nb_processors+1;
      //int nb = 1;
      
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
  
  IVect& GlobalRowNumbers = info.GlobalRowNumbers;
  GlobalRowNumbers.Reallocate(nodl);
  for (int i = 0; i < m; i++)
    if (Glob_to_local(i) >= 0)
      GlobalRowNumbers(Glob_to_local(i)) = i;

  // then OverlapRowNumbers/OverlapProcNumbers
  IVect& OverlapRowNumbers = info.OverlapRowNumbers;
  OverlapRowNumbers.Reallocate(nodl_overlap);
  IVect& OverlapProcNumbers = info.OverlapProcNumbers;
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
  int nodl_scalar = nodl; // nb_u = 1
  int nb_proc = 0;
  IVect IndexProc(nb_processors); IndexProc.Fill(-1);
  for (int p = 0; p < nb_processors; p++)
    if (NbSharedDofPerProc(p) > 0)
      {
        IndexProc(p) = nb_proc;
        nb_proc++;
      }
  
  IVect& ProcSharingRows = info.ProcSharingRows;
  ProcSharingRows.Reallocate(nb_proc);
  Vector<IVect>& SharingRowNumbers = info.SharingRowNumbers;
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

  // distributing matrices to all the processors
  for (int i = 0; i < Ai.GetM(); i++)
    {
      DistributedMatrix<T, Prop, Storage> Aref;
      Aref = Ai(i);
      Ai(i).Reallocate(nodl, nodl);
      Ai(i).Init(m, &GlobalRowNumbers, &OverlapRowNumbers, &OverlapProcNumbers,
                 nodl_scalar, 1, &ProcSharingRows, &SharingRowNumbers, comm);
      
      DistributeMatrixProcessor(Aref, Ai(i), rank_processor, list_proc, Glob_to_local);
    }
}
