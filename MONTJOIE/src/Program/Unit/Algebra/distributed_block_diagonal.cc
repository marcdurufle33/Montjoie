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

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRand(x);
      A.Set(i, j, x);
      if (IsSymmetricMatrix(A))
        A.Set(j, i, x);
    }
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
  Real_wp zero, sum; SetComplexZero(zero);
  sum = zero; VectReal_wp coef(n);
  for (int i = 0; i < n; i++)
    {
      GetRand(coef(i));
      sum += coef(i);
    }
  
  Mlt(Real_wp(1)/sum, coef);
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

template<class MatrixSparse, class MatrixSparse2>
void CheckMatrixMlt(const MatrixSparse& A,
                    const MatrixSparse2& Aref,
                    const string& fct_name, bool all_col = false)
{
  typedef typename MatrixSparse::entry_type T;
  int m = A.GetGlobalM();
  Vector<T> Xref, Yref;
  GenerateRandomVector(Xref, m);
  GenerateRandomVector(Yref, m);

  T alpha, beta;
  GetRand(alpha);   GetRand(beta);
  
  Vector<T> X(A.GetM()), Y(A.GetM());
  const IVect& global = A.GetGlobalRowNumber();
  for (int i = 0; i < global.GetM(); i++)
    {
      X(i) = Xref(global(i));
      Y(i) = Yref(global(i));
    }
  
  if (all_col)
    {
      Vector<T> Ones(A.GetM()), Ah_Ones(A.GetM());
      for (int j = 0; j < m; j++)
        {
          Ones.Fill(0);
          for (int k = 0; k < global.GetM(); k++)
            if (global(k) == j)
              SetComplexOne(Ones(k));
          
          Mlt(A, Ones, Ah_Ones);
          
          for (int i = 0; i < global.GetM(); i++)
            if ((abs(Ah_Ones(i) - Aref(global(i), j)) > threshold) || isnan(abs(Ah_Ones(i) - Aref(global(i), j))))
              {
                cout << fct_name << " incorrect" << endl;
                DISP(global(i)); DISP(j); DISP(Ah_Ones(i));
                DISP(Aref(global(i), j));
                abort();
              }
        }
    }

  MltAdd(alpha, A, X, beta, Y);
  MltAdd(alpha, Aref, Xref, beta, Yref);
  
  for (int i = 0; i < global.GetM(); i++)
    if ((abs(Yref(global(i)) - Y(i)) > threshold) || isnan(abs(Yref(global(i)) - Y(i))))
      {
        cout << fct_name << " incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(global(i)));
        abort();
      }

  MltAdd(alpha, SeldonTrans, A, X, beta, Y);
  MltAdd(alpha, SeldonTrans, Aref, Xref, beta, Yref);
  
  for (int i = 0; i < global.GetM(); i++)
    if ((abs(Yref(global(i)) - Y(i)) > threshold) || isnan(abs(Yref(global(i)) - Y(i))))
      {
        cout << fct_name << " (transpose) incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(global(i)));
        abort();
      }

  MltAdd(alpha, SeldonConjTrans, A, X, beta, Y);
  MltAdd(alpha, SeldonConjTrans, Aref, Xref, beta, Yref);
  
  for (int i = 0; i < global.GetM(); i++)
    if ((abs(Yref(global(i)) - Y(i)) > threshold) || isnan(abs(Yref(global(i)) - Y(i))))
      {
        cout << fct_name << " (conjugate transpose) incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(global(i)));
        DISP(Y);
        abort();
      }
  
}

template<class MatrixSeq, class MatrixPar>
void DistributeMatrixProcessor(const MatrixSeq& Aref, MatrixPar& A,
                               const Vector<IVect>& list_proc, const IVect& Glob_to_local)
{
  A.Fill(0);
  //A.SwitchToGlobalNumbers();
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
            
            //DISP(i); DISP(j); DISP(Aref(i, j)); DISP(eval); DISP(num_row); DISP(num_col);
            for (int k = 0; k < nb; k++)
              if (rank_processor == num(k))
                AddInteraction(A, i, j, eval(k), Glob_to_local,
                               list_proc(i)(num_row(k)), list_proc(j)(num_col(k)));
          }
    }
}

template<class T, class Storage, class Allocator, class T1, class Allocator1>
void ScaleLeftMatrix(Matrix<T, Symmetric, Storage, Allocator>&, const Vector<T1, VectFull, Allocator1>&)
{}

template<class T, class Storage, class Allocator, class T1, class Allocator1>
void ScaleRightMatrix(Matrix<T, Symmetric, Storage, Allocator>&, const Vector<T1, VectFull, Allocator1>&)
{}

template<class T, class Prop, class Storage, class Allocator,
         class T2, class Prop2, class Storage2, class Allocator2>
void CheckDistributedMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& mat,
                            DistributedMatrix<T2, Prop2, Storage2, Allocator2>& mat2,
			    MPI_Comm comm, string suffix)
{
  Matrix<T, Prop, Storage, Allocator> Aref;
  typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type T0;
  
  // random pattern
  int m = 250, size_max = 20; IVect permut;
  GenerateRandomPermutation(m, permut);
  
  IVect size_block;
  int nb = 0;
  while (nb < m)
    {
      int k = rand()%size_max + 1;
      if (nb+k > m)
        k = m-nb;
      
      size_block.PushBack(k);
      nb += k;
    }

  Vector<IVect> pattern;
  pattern.Reallocate(size_block.GetM());
  nb = 0;
  for (int i = 0; i < size_block.GetM(); i++)
    {
      pattern(i).Reallocate(size_block(i));
      for (int j = 0; j < size_block(i); j++)
        {
          pattern(i)(j) = permut(nb);
          nb++;
        }
      
      Sort(pattern(i));
    }
  
  //DISP(size_block); DISP(pattern);
  Aref.SetPattern(pattern);
  
  for (int i = 0; i < Aref.GetNbBlocks(); i++)
    {
      int size_block = Aref.GetBlockSize(i);
      for (int j = 0; j < size_block; j++)
        for (int k = 0; k < size_block; k++)
          GetRand(Aref.Value(i, j, k));
    }
  
  // for each dof determining processors that share this dof
  Vector<IVect> list_proc(m);
  for (int i = 0; i < m; i++)
    {
      nb = rand()%nb_processors+1;
      //nb = 1;
      
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
  
  IVect GlobalRowNumbers(nodl);
  for (int i = 0; i < m; i++)
    if (Glob_to_local(i) >= 0)
      GlobalRowNumbers(Glob_to_local(i)) = i;
  
  //DISP(nodl); DISP(GlobalRowNumbers);
  
  // then OverlapRowNumbers/OverlapProcNumbers
  IVect OverlapRowNumbers(nodl_overlap), OverlapProcNumbers(nodl_overlap);
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
  
  IVect ProcSharingRows(nb_proc);
  Vector<IVect> SharingRowNumbers(nb_proc);
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
  
  // distributing Aref to all the processors
  DistributedMatrix<T, Prop, Storage, Allocator> A;
  A.Reallocate(nodl, nodl);
  A.Init(m, &GlobalRowNumbers, &OverlapRowNumbers, &OverlapProcNumbers,
         nodl_scalar, 1, &ProcSharingRows, &SharingRowNumbers, comm);

  
  int nblock = 0;
  Vector<bool> block_included(pattern.GetM());
  for (int i = 0; i < pattern.GetM(); i++)
    {
      bool block_taken = false;
      for (int j = 0; j < pattern(i).GetM(); j++)
        if (Glob_to_local(pattern(i)(j)) >= 0)
          block_taken = true;
      
      block_included(i) = block_taken;
      if (block_taken)
        nblock++;
    }
  
  Vector<IVect> pattern_loc(nblock), proc_loc(nblock);
  nblock = 0;
  for (int i = 0; i < pattern.GetM(); i++)
    if (block_included(i))
      {    
        proc_loc(nblock).Reallocate(pattern(i).GetM());
        pattern_loc(nblock).Reallocate(pattern(i).GetM());      
        for (int j = 0; j < pattern(i).GetM(); j++)
          {
            if (Glob_to_local(pattern(i)(j)) >= 0)
              {
                proc_loc(nblock)(j) = rank_processor;
                pattern_loc(nblock)(j) = Glob_to_local(pattern(i)(j));
              }
            else
              {
                int p = list_proc(pattern(i)(j)).GetM()-1;
                proc_loc(nblock)(j) = list_proc(pattern(i)(j))(p);
                pattern_loc(nblock)(j) = pattern(i)(j);
              }
          }      
        
        nblock++;
      }  
  
  for (int i = 0; i < pattern_loc.GetM(); i++)
    {
      //DISP(i);
      //DISP(pattern_loc(i)); DISP(proc_loc(i));
    }
  
  A.SetPattern(pattern_loc, proc_loc);
  
  DistributeMatrixProcessor(Aref, A, list_proc, Glob_to_local);
  
  if (rank_processor == 0)
    Aref.WriteText(string("Dmat_ref" + suffix + ".dat"));
  
  A.PrepareMltAdd();
  A.WriteText(string("Dmat"+suffix+".dat"));

  // testing Mlt function
  CheckMatrixMlt(A, Aref, "Mlt");
  
  // Mlt with overwrite
  Vector<T> Yref, Y(nodl);
  GenerateRandomVector(Yref, m);
  for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
    Y(i) = Yref(GlobalRowNumbers(i));
  
  Mlt(A, Y);
  Mlt(Aref, Yref);
  for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
    if ((abs(Yref(GlobalRowNumbers(i)) - Y(i)) > threshold) || isnan(abs(Yref(GlobalRowNumbers(i)) - Y(i))))
      {
        cout << " Mlt overwrite incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(GlobalRowNumbers(i)));
        abort();
      }
  
  // testing conversion to a sparse matrix
  DistributedMatrix<T2, Prop2, Storage2, Allocator2> Asparse, Bsparse; 
  Copy(A, Asparse);
  
  CheckMatrixMlt(Asparse, Aref, "Copy");
  //Asparse.WriteText("mat_sp.dat");

  ConvertToSparse(A, Bsparse);
  CheckMatrixMlt(Bsparse, Aref, "ConvertToSparse");
  
  // testing conversion in the other way
  DistributedMatrix<T, Prop, Storage, Allocator> B;
  
  Copy(Bsparse, B);
  //A.WriteText("mat_modif.dat");
  //B.WriteText("mat_copy.dat");
  CheckMatrixMlt(B, Aref, "Copy");

  // testing function Add
  Matrix<bool> keep_interac(m, m);
  keep_interac.Fill(false);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        if (rand()%3 == 1)
          keep_interac(i, j) = true;
      }

  if (IsSymmetricMatrix(Aref))
    for (int i = 0; i < m; i++)
      for (int j = i+1; j < m; j++)
        keep_interac(j, i) = keep_interac(i, j);
  
  for (int i = 0; i < B.GetNbBlocks(); i++)
    {
      int size_block = B.GetBlockSize(i);
      for (int j = 0; j < size_block; j++)
        for (int k = 0; k < size_block; k++)
          {
            int irow = B.IndexGlobal(i, j);
            int jcol = B.IndexGlobal(i, k);
            if (!keep_interac(irow, jcol))
              SetComplexZero(B.Value(i, j, k));
          }
    }
  
  Copy(B, Bsparse);
  
  //B.WriteText("matBzero.dat");
  Bsparse.RemoveSmallEntry(threshold);

  //Bsparse.WriteText("matBzero_sp.dat");
  Copy(Bsparse, B);
  
  // B.WriteText("matBzero_bis.dat");
  for (int i = 0; i < B.GetNbBlocks(); i++)
    {
      int size_block = B.GetBlockSize(i);
      //DISP(i);
      for (int j = 0; j < size_block; j++)
        {
          //DISP(j); DISP(B.IndexGlobal(i, j));
        }
    }
  
  Matrix<T, Prop, Storage, Allocator> Bref;
  Bref.SetPattern(pattern);

  for (int i = 0; i < Bref.GetNbBlocks(); i++)
    {
      int size_block = Bref.GetBlockSize(i);
      for (int j = 0; j < size_block; j++)
        for (int k = 0; k < size_block; k++)
          {
            int irow = Bref.Index(i, j);
            int jcol = Bref.Index(i, k);
            if (keep_interac(irow, jcol))
              GetRand(Bref.Value(i, j, k));
          }
    }
  
  DistributeMatrixProcessor(Bref, B, list_proc, Glob_to_local);
  
  T alpha; GetRand(alpha);
  Add(alpha, B, A);
  Add(alpha, Bref, Aref);
  
  CheckMatrixMlt(A, Aref, "Add");  

  Matrix<T, Prop, Storage, Allocator> A0_ref(Aref);
  DistributedMatrix<T, Prop, Storage, Allocator> A0(A);
  
  // testing GetInverse
  //A.Assemble();
  GetInverse(A);
  GetInverse(Aref);
  
  Mlt(A, Y, false);
  Mlt(Aref, Yref);
  for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
    if ((abs(Yref(GlobalRowNumbers(i)) - Y(i)) > 1e4*threshold) || isnan(abs(Yref(GlobalRowNumbers(i)) - Y(i))))
      {
        cout << " GetInverse incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(GlobalRowNumbers(i)));
        abort();
      }

  // testing EraseRow
  int nb_row = 10; IVect num_row, num_row_loc;
  for (int i = 0; i < nb_row; i++)
    {
      int j = rand()%m;
      num_row.PushBack(j);

      int jloc = Glob_to_local(j);
      if (jloc >= 0)
	num_row_loc.PushBack(jloc);
    }

  RemoveDuplicate(num_row);
  RemoveDuplicate(num_row_loc);

  Matrix<T, Prop, Storage, Allocator> B0_ref(A0_ref);
  DistributedMatrix<T, Prop, Storage, Allocator> B0(A0);
  
  //if (rank_processor == 0)
  //A0_ref.WriteText("Bref.dat");

  //A0.WriteText("B.dat");

  EraseRow(num_row, A0_ref);
  EraseRow(num_row_loc, A0);

  //if (rank_processor == 0)
  //A0_ref.WriteText("Aref.dat");

  A0.PrepareMltAdd();
  //A0.WriteText("A.dat");

  CheckMatrixMlt(A0, A0_ref, "EraseRow", true);

  EraseCol(num_row, B0_ref);
  EraseCol(num_row_loc, B0);

  B0.PrepareMltAdd();
  CheckMatrixMlt(B0, B0_ref, "EraseCol", true);  

}


template<class T, class Prop, class Storage, class Allocator>
void CheckCholesky(DistributedMatrix<T, Prop, Storage, Allocator>& mat, MPI_Comm comm, string suffix)
{
  Matrix<T, Prop, Storage, Allocator> Aref, Lref;
  typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type T0;
  
  // random pattern
  int m = 50; IVect permut;
  GenerateRandomPermutation(m, permut);
  
  IVect size_block;
  int nb = 0;
  while (nb < m)
    {
      int k = rand()%10 + 1;
      if (nb+k > m)
        k = m-nb;
      
      size_block.PushBack(k);
      nb += k;
    }

  Vector<IVect> pattern;
  pattern.Reallocate(size_block.GetM());
  nb = 0;
  for (int i = 0; i < size_block.GetM(); i++)
    {
      pattern(i).Reallocate(size_block(i));
      for (int j = 0; j < size_block(i); j++)
        {
          pattern(i)(j) = permut(nb);
          nb++;
        }
      
      Sort(pattern(i));
    }
  
  //DISP(size_block); DISP(pattern);
  Lref.SetPattern(pattern);
  Aref.SetPattern(pattern);
  
  for (int i = 0; i < Lref.GetNbBlocks(); i++)
    {
      int size_block = Lref.GetBlockSize(i);
      for (int j = 0; j < size_block; j++)
        for (int k = 0; k < size_block; k++)
          GetRand(Lref.Value(i, j, k));
    }
  
  // Aref = Lref Lref to ensure positivity
  Mlt(Lref, Lref, Aref);
  
  // for each dof determining processors that share this dof
  Vector<IVect> list_proc(m);
  for (int i = 0; i < m; i++)
    {
      nb = rand()%nb_processors+1;
      //nb = 1;
      
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
  
  IVect GlobalRowNumbers(nodl);
  for (int i = 0; i < m; i++)
    if (Glob_to_local(i) >= 0)
      GlobalRowNumbers(Glob_to_local(i)) = i;
  
  //DISP(nodl); DISP(GlobalRowNumbers);
  
  // then OverlapRowNumbers/OverlapProcNumbers
  IVect OverlapRowNumbers(nodl_overlap), OverlapProcNumbers(nodl_overlap);
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
  
  IVect ProcSharingRows(nb_proc);
  Vector<IVect> SharingRowNumbers(nb_proc);
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
  
  // distributing Aref to all the processors
  DistributedMatrix<T, Prop, Storage, Allocator> A;
  A.Reallocate(nodl, nodl);
  A.Init(m, &GlobalRowNumbers, &OverlapRowNumbers, &OverlapProcNumbers,
         nodl_scalar, 1, &ProcSharingRows, &SharingRowNumbers, comm);
  
  int nblock = 0;
  Vector<bool> block_included(pattern.GetM());
  for (int i = 0; i < pattern.GetM(); i++)
    {
      bool block_taken = false;
      for (int j = 0; j < pattern(i).GetM(); j++)
        if (Glob_to_local(pattern(i)(j)) >= 0)
          block_taken = true;
      
      block_included(i) = block_taken;
      if (block_taken)
        nblock++;
    }
  
  Vector<IVect> pattern_loc(nblock), proc_loc(nblock);
  nblock = 0;
  for (int i = 0; i < pattern.GetM(); i++)
    if (block_included(i))
      {    
        proc_loc(nblock).Reallocate(pattern(i).GetM());
        pattern_loc(nblock).Reallocate(pattern(i).GetM());      
        for (int j = 0; j < pattern(i).GetM(); j++)
          {
            if (Glob_to_local(pattern(i)(j)) >= 0)
              {
                proc_loc(nblock)(j) = rank_processor;
                pattern_loc(nblock)(j) = Glob_to_local(pattern(i)(j));
              }
            else
              {
                int p = list_proc(pattern(i)(j)).GetM()-1;
                proc_loc(nblock)(j) = list_proc(pattern(i)(j))(p);
                pattern_loc(nblock)(j) = pattern(i)(j);
              }
          }      
        
        nblock++;
      }  
  
  for (int i = 0; i < pattern_loc.GetM(); i++)
    {
      //DISP(i);
      //DISP(pattern_loc(i)); DISP(proc_loc(i));
    }
  
  A.SetPattern(pattern_loc, proc_loc);
  
  DistributeMatrixProcessor(Aref, A, list_proc, Glob_to_local);
  
  if (rank_processor == 0)
    Aref.WriteText(string("mat_ref"+suffix+".dat"));
  
  A.PrepareMltAdd();
  A.WriteText(string("mat"+suffix+".dat"));
  
  // testing Mlt function
  CheckMatrixMlt(A, Aref, "Mlt");
  
  // testing GetCholesky
  GetCholesky(A);
  GetCholesky(Aref);
  
  if (rank_processor == 0)
    Aref.WriteText(string("chol_ref"+suffix+".dat"));
  
  A.WriteText(string("chol"+suffix+".dat"));
  
  Vector<T> Yref, Y(nodl);
  GenerateRandomVector(Yref, m);
  for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
    Y(i) = Yref(GlobalRowNumbers(i));

  MltCholesky(SeldonTrans, A, Y);
  MltCholesky(SeldonTrans, Aref, Yref);
  MltCholesky(SeldonNoTrans, A, Y);
  MltCholesky(SeldonNoTrans, Aref, Yref);
  
  for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
    if ((abs(Yref(GlobalRowNumbers(i)) - Y(i)) > 10*threshold) || isnan(abs(Yref(GlobalRowNumbers(i)) - Y(i))))
      {
        cout << " MltCholesky incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(GlobalRowNumbers(i)));
        abort();
      }

  SolveCholesky(SeldonNoTrans, A, Y);
  SolveCholesky(SeldonNoTrans, Aref, Yref);
  SolveCholesky(SeldonTrans, A, Y);
  SolveCholesky(SeldonTrans, Aref, Yref);

  for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
    if ((abs(Yref(GlobalRowNumbers(i)) - Y(i)) > 10*threshold) || isnan(abs(Yref(GlobalRowNumbers(i)) - Y(i))))
      {
        cout << " SolveCholesky incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(GlobalRowNumbers(i)));
        abort();
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
  
  {
    DistributedMatrix<Real_wp, General, BlockDiagRow> A;
    DistributedMatrix<Real_wp, General, ArrayRowSparse> B;
    CheckDistributedMatrix(A, B, comm, suffix);    
  }

  {
    DistributedMatrix<Complex_wp, General, BlockDiagRow> A;
    DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse> B;
    CheckDistributedMatrix(A, B, comm, suffix);        
  }

  {
    DistributedMatrix<Real_wp, Symmetric, BlockDiagRowSym> A;
    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> B;
    CheckDistributedMatrix(A, B, comm, suffix);    
    CheckCholesky(A, comm, suffix);
  }

  {
    DistributedMatrix<Complex_wp, Symmetric, BlockDiagRowSym> A;
    DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    CheckDistributedMatrix(A, B, comm, suffix);    
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank_processor == 0)
    cout << "All tests passed successfully" << endl;

  return FinalizeMontjoie();
}
