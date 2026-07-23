//#define MONTJOIE_WITH_ONE_DIM
//#define MONTJOIE_WITH_TWO_DIM
//#define MONTJOIE_WITH_THREE_DIM
//#define MONTJOIE_WITH_NODAL_H1
//#define MONTJOIE_WITH_HP_H1
//#define MONTJOIE_WITH_NODAL_DG
//#define MONTJOIE_WITH_ORTHO_DG
//#define MONTJOIE_WITH_LEGENDRE_DG
//#define MONTJOIE_WITH_HP_HCURL
//#define MONTJOIE_WITH_HCURL_AXI
//#define MONTJOIE_WITH_NODAL_HCURL
//#define MONTJOIE_WITH_HP_HDIV
//#define MONTJOIE_WITH_NODAL_HDIV

//#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"

#include "Algebra/MontjoieAlgebra.hxx"

//#include "vector/Vector2.cxx"
//#include "vector/Vector3.cxx"

using namespace Montjoie;

void CheckMemorySize()
{
  int test_input;
  
  cout << "Waiting... " << endl; cin >> test_input;

  /*****************************************
   * structures in folder array and vector *
   *****************************************/
  
  /* Array<double, 5> x;
  
  x.Reallocate(10, 20, 40, 20, 30);
  x.Fill();
  
  cout << "Memory taken by Array x = " << GetHumanReadableMemory(x.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Array4D<complex<double> > x4D;
  
  x4D.Reallocate(21, 32, 43, 55);
  x4D.Fill();
  
  cout << "Memory taken by Array4D = " << GetHumanReadableMemory(x4D.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Array3D<complex<double> > x3D;
  
  x3D.Reallocate(108, 64, 78);
  x3D.Fill();
  
  cout << "Memory taken by Array3D = " << GetHumanReadableMemory(x3D.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<int64_t> v;
  
  v.Reallocate(2000000);
  v.Fill();
  
  cout << "Memory taken by Vector = " << GetHumanReadableMemory(v.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  */

  Vector<Vector<double, VectSparse> > vSp;
  
  vSp.Reallocate(48000);
  for (int i = 0; i < vSp.GetM(); i++)
    {
      vSp(i).Reallocate(100);
      for (int j = 0; j < vSp(i).GetM(); j++)
	{
          vSp(i).Index(j) = rand()%1000;
          vSp(i).Value(j) = rand();
        }
    }
  
  cout << "Memory taken by VectSparse = " << GetHumanReadableMemory(GetMemorySize(vSp)) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  for (int i = 0; i < vSp.GetM(); i++)
    {
      vSp(i).Reallocate(80);
      for (int j = 0; j < vSp(i).GetM(); j++)
	{
          vSp(i).Index(j) = rand()%1000;
          vSp(i).Value(j) = rand();
        }
    }

  cout << "Memory taken by VectSparse = " << GetHumanReadableMemory(GetMemorySize(vSp)) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<double, VectSparse> vSp2;
  
  vSp2.Reallocate(2000000);
  for (int k = 0; k < vSp2.GetM(); k++)
    {
      vSp2.Index(k) = rand()%1000;
      vSp2.Value(k) = rand();
    }

  cout << "Memory taken by VectSparse = " << GetHumanReadableMemory(GetMemorySize(vSp)+vSp2.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  /*
  Vector2<complex<double> > v2;
  
  v2.Reallocate(100000);
  for (int i = 0; i < v2.GetLength(); i++)
    {
      v2.Reallocate(i, 10);
      v2(i).FillRand();
    }

  cout << "Memory taken by Vector2 = " << GetHumanReadableMemory(v2.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector3<complex<double> > v3;
  
  v3.Reallocate(500);
  for (int i = 0; i < v3.GetLength(); i++)
    {
      v3.Reallocate(i, 200);
      for (int j = 0; j < v3.GetLength(i); j++)
        {
          v3(i, j).Reallocate(20);
          v3(i, j).FillRand();
        }
    }
  
  cout << "Memory taken by Vector3 = " << GetHumanReadableMemory(v3.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */

  /*******************************
   * structures in folder matrix *
   *******************************/
  
  /*
  Matrix<complex<double>, Hermitian, RowHermPacked> Mhp;
  
  Mhp.Reallocate(2034, 2034);
  Mhp.Fill();
  
  cout << "Memory taken by Matrix_HermPacked = " << GetHumanReadableMemory(Mhp.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Hermitian, RowHerm> Mhe;
  
  Mhe.Reallocate(2034, 2034);
  Mhe.Fill();
  
  cout << "Memory taken by Matrix_Hermitian = " << GetHumanReadableMemory(Mhe.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Matrix<complex<double>, General, RowMajor> Mge;
  
  Mge.Reallocate(1228, 1228);
  Mge.Fill();
  
  cout << "Memory taken by Matrix_Pointers = " << GetHumanReadableMemory(Mge.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, RowSymPacked> Msp;
  
  Msp.Reallocate(1228, 1228);
  Msp.Fill();
  
  cout << "Memory taken by Matrix_SymPacked = " << GetHumanReadableMemory(Msp.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, RowSym> Msy;
  
  Msy.Reallocate(1228, 1228);
  Msy.Fill();
  
  cout << "Memory taken by Matrix_Symmetric = " << GetHumanReadableMemory(Msy.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Matrix<complex<double>, General, RowLoTriangPacked> Mtp;
  
  Mtp.Reallocate(1228, 1228);
  Mtp.Fill();
  
  cout << "Memory taken by Matrix_TriangPacked = " << GetHumanReadableMemory(Mtp.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<double, Symmetric, RowLoTriang> Mtr;
  
  Mtr.Reallocate(1228, 1228);
  Mtr.Fill();
  
  cout << "Memory taken by Matrix_Triangular = " << GetHumanReadableMemory(Mtr.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  */

  /**************************************
   * structures in folder matrix_sparse *
   **************************************/
  
  /* Matrix<complex<double>, General, RowSparse> A1;
  
  A1.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_Sparse = " << GetHumanReadableMemory(A1.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, RowSymSparse> A2;
  
  A2.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_SymSparse = " << GetHumanReadableMemory(A2.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Matrix<complex<double>, General, ArrayRowSparse> A3;
  
  A3.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_ArraySparse = " << GetHumanReadableMemory(A3.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, ArrayRowSymSparse> A4;
  
  A4.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_ArraySparse = " << GetHumanReadableMemory(A4.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, General, RowComplexSparse> Ac1;
  
  Ac1.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_ComplexSparse = " << GetHumanReadableMemory(Ac1.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, RowSymComplexSparse> Ac2;
  
  Ac2.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_SymComplexSparse = " << GetHumanReadableMemory(Ac2.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Matrix<complex<double>, General, ArrayRowComplexSparse> Ac3;
  
  Ac3.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_ArrayComplexSparse = " << GetHumanReadableMemory(Ac3.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, ArrayRowSymComplexSparse> Ac4;
  
  Ac4.ReadText("mat0.dat");
  
  cout << "Memory taken by Matrix_ArrayComplexSparse = " << GetHumanReadableMemory(Ac4.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Matrix<double, General, BandedCol> B;
  
  B.Reallocate(100000, 100000, 10, 20);
  B.Fill(1.0);
  
  cout << "Memory taken by Matrix_Band = " << GetHumanReadableMemory(B.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<double, General, ArrowCol> Ba;
  
  Ba.Reallocate(100000, 100000, 10, 20, 5, 8);
  Ba.Fill(1.0);
  
  cout << "Memory taken by Matrix_Arrow = " << GetHumanReadableMemory(Ba.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  DistributedMatrix<complex<double>, General, ArrayRowSparse> C;
  
  Copy(A1, C);
  for (int k = 0; k < 2000000; k++)
    {
      int iloc = rand()%C.GetM();
      int iglob = rand();
      complex<double> val = rand();
      C.AddDistantInteraction(iloc, iglob, rand()%4, val);
      val = rand();
      C.AddRowDistantInteraction(iglob, iloc, rand()%4, val);
    }
  
  cout << "Memory taken by DistributedMatrix = " << GetHumanReadableMemory(C.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  */
  
  /******************************************
   * structures in folder interfaces/direct *
   ******************************************/

  /*
#ifdef SELDON_WITH_CHOLMOD
  MatrixCholmod mat_chol;
  
  Matrix<double, Symmetric, RowSymSparse> A;
  A.ReadText("mass0.dat");
  
  mat_chol.ShowFullHistory();
  mat_chol.FactorizeMatrix(A);
  
  cout << "Memory taken by Choldmod = " << GetHumanReadableMemory(mat_chol.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif

#ifdef SELDON_WITH_MUMPS
  MatrixMumps<complex<double> > mat_mumps;
  
  Matrix<complex<double>, Symmetric, RowSymSparse> B;
  B.ReadText("mat0.dat");
  
  mat_mumps.HideMessages();
  mat_mumps.FactorizeMatrix(B);
  
  cout << "Memory taken by Mumps = " << GetHumanReadableMemory(mat_mumps.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif

#ifdef SELDON_WITH_PARDISO
  MatrixPardiso<complex<double> > mat_pardiso;
  
  B.ReadText("mat0.dat");
  
  mat_pardiso.HideMessages();
  mat_pardiso.FactorizeMatrix(B);
  
  cout << "Memory taken by Pardiso = " << GetHumanReadableMemory(mat_pardiso.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif

#ifdef SELDON_WITH_PASTIX
  MatrixPastix<complex<double> > mat_pastix;
  
  B.ReadText("mat0.dat");
  
  mat_pastix.HideMessages();
  mat_pastix.FactorizeMatrix(B);
  
  cout << "Memory taken by Pastix = " << GetHumanReadableMemory(mat_pastix.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif

#ifdef SELDON_WITH_SUPERLU
  MatrixPastix<complex<double> > mat_superlu;
  
  B.ReadText("mat0.dat");
  
  mat_superlu.HideMessages();
  mat_superlu.FactorizeMatrix(B);
  
  cout << "Memory taken by SuperLU = " << GetHumanReadableMemory(mat_superlu.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif

#ifdef SELDON_WITH_UMFPACK
  MatrixUmfPack<complex<double> > mat_umfpack;
  
  B.ReadText("mat0.dat");
  
  mat_umfpack.HideMessages();
  mat_umfpack.FactorizeMatrix(B);
  
  cout << "Memory taken by UmfPack = " << GetHumanReadableMemory(mat_umfpack.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif

#ifdef SELDON_WITH_WSMP
  MatrixWsmp<complex<double> > mat_wsmp;
  
  B.ReadText("mat0.dat");
  
  mat_wsmp.HideMessages();
  mat_wsmp.FactorizeMatrix(B);
  
  cout << "Memory taken by Wsmp = " << GetHumanReadableMemory(mat_wsmp.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
#endif
  */
  
  /********************************
   * structures in folder Algebra *
   ********************************/
  
  /* Matrix<double, General, RowSkyLine> As;
  
  As.Reallocate(20000, 20000);
  for (int i = 0; i < As.GetM(); i++)
    {
      int size_row = rand()%50+1;
      int idx_begin = max(0, i - rand()%size_row);
      As.ReallocateRow(i, size_row, idx_begin);
    }
  
  cout << "Memory taken by GeneralSkyLineMatrix = " << GetHumanReadableMemory(As.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<double, Symmetric, SymColSkyLine> Bs;
  Matrix<double, Symmetric, RowSymSparse> mass;
  
  mass.ReadText("massSmall0.dat");
  CopyMatrix(mass, Bs);
  mass.Clear();
  
  cout << "Memory taken by SymmetricSkyLineMatrix = " << GetHumanReadableMemory(Bs.GetMemorySize()) << endl;
  
  Matrix<complex<double>, General, BlockDiagRow> Bh;
  
  int N = 20000;
  Vector<int> NumBlock(N);
  int nb_block = 2000;
  IVect NbDofsPerBlock(nb_block);
  NbDofsPerBlock.Zero();
  for (int i = 0; i < N; i++)
    {
      NumBlock(i) = rand()%nb_block;
      NbDofsPerBlock(NumBlock(i))++;
    }

  int nb_actual_block = 0;
  for (int i = 0; i < nb_block; i++)
    if (NbDofsPerBlock(i) > 0)
      nb_actual_block++;
  
  Vector<IVect> NumDofs(nb_actual_block);
  nb_actual_block = 0;
  IVect IndexBlock(nb_block);
  for (int i = 0; i < nb_block; i++)
    if (NbDofsPerBlock(i) > 0)
      {
        IndexBlock(i) = nb_actual_block;
        NumDofs(nb_actual_block).Reallocate(NbDofsPerBlock(i));
        nb_actual_block++;
      }

  NbDofsPerBlock.Zero();
  for (int i = 0; i < N; i++)
    {
      int n = IndexBlock(NumBlock(i));
      NumDofs(n)(NbDofsPerBlock(n)) = i;
      NbDofsPerBlock(n)++;
    }
  
  IndexBlock.Clear(); NbDofsPerBlock.Clear(); NumBlock.Clear();
  cout << "Waiting... " << endl; cin >> test_input;
  
  Bh.SetPattern(NumDofs);  
  Bh.FillRand();
  
  cout << "Memory taken by Matrix_BlockDiagonal = " << GetHumanReadableMemory(Bh.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Matrix<complex<double>, Symmetric, BlockDiagRowSym> Bh_sym;
  
  Bh_sym.SetPattern(NumDofs);  
  Bh_sym.FillRand();
  
  cout << "Memory taken by MatrixSym_BlockDiagonal = " << GetHumanReadableMemory(Bh_sym.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  TinyBandMatrix<complex<double>, 5> C;
  
  C.Reallocate(1000000, 1000000);
  C.FillRand();
  
  cout << "Memory taken by TinyBandMatrix = " << GetHumanReadableMemory(C.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  TinyArrowMatrix<complex<double>, 6, 3> Ca;
  
  Ca.Reallocate(1000000, 1000000);
  Ca.FillRand();
  
  cout << "Memory taken by TinyArrowMatrix = " << GetHumanReadableMemory(Ca.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  */

  /************************************
   * structures in folder Computation *
   ************************************/

  /* MatrixVectorProductLevel level;
  Vector<IVect>& test = level.GetLevelArray();
  
  test.Reallocate(1000000);
  for (int i = 0; i < test.GetM(); i++)
    {
      test(i).Reallocate(30);
      test(i).Fill();
    }

  cout << "Memory taken by MatrixVectorProductLevel = " << GetHumanReadableMemory(level.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  ExtrapolVariablesProductFEM<complex<double>, LaplaceEquation<Dimension2> > var;
  
  var.extrapolU(0).Reallocate(500000); var.extrapolU(0).Fill();
  var.extrapolDun(0).Reallocate(600000); var.extrapolDun(0).Fill();
  var.extrapolCgradU(0).Reallocate(800000); var.extrapolCgradU(0).Fill();
  
  var.Un_quad.Reallocate(20000);
  for (int i = 0; i < var.Un_quad.GetM(); i++)
    {
      var.Un_quad(i)(0).Reallocate(30);
      var.Un_quad(i)(0).Fill();
    }

  var.Uloc(0).Reallocate(200000); var.Uloc(0).Fill();
  var.Uquad(0).Reallocate(300000); var.Uquad(0).Fill();
  var.Uneighbor(0).Reallocate(400000); var.Uneighbor(0).Fill();
  var.DunNeighbor(0).Reallocate(600000); var.DunNeighbor(0).Fill();
  var.CgradUNeighbor(0).Reallocate(500000); var.CgradUNeighbor(0).Fill();

  cout << "Memory taken by ExtrapolVariablesProductFEM = " << GetHumanReadableMemory(var.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  EllipticProblem<LaplaceEquation<Dimension2> > var_harmonic;
  FemMatrixFreeClass<complex<double>, LaplaceEquation<Dimension2> > mat(var_harmonic);
  
  mat.GetExtrapolVariables() = var;  
  mat.matCSR_iterative_unsym.ReadText("mat0.dat");
  
  cout << "Memory taken by FemMatrixFreeClass = " << GetHumanReadableMemory(mat.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */
  
  /*******************************************
   * structures in folder FiniteElement/Edge *
   *******************************************/
  /* 
  Vector<EdgeGauss> elt(1000);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%10+2;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by EdgeGauss = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<EdgeLobatto> eltLob(5000);
  size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%10+2;
      eltLob(i).ConstructFiniteElement(r);
      size += eltLob(i).GetMemorySize();
    }

  cout << "Memory taken by EdgeLobatto = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<EdgeHierarchic> eltHp(5000);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by EdgeHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  */

  /***********************************************
   * structures in folder FiniteElement/Triangle *
   ***********************************************/

  // problem with size in TriangleGeomReference ?
  /* Vector<TriangleGeomReference> elt(200);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%10+2;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleGeomReference = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleClassical> eltClassical(100);
  size = 0;
  for (int i = 0; i < eltClassical.GetM(); i++)
    {
      int r = rand()%10+2;
      eltClassical(i).ConstructFiniteElement(r);
      size += eltClassical(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleClassical = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleLobatto> eltLobatto(100);
  size = 0;
  for (int i = 0; i < eltLobatto.GetM(); i++)
    {
      int r = rand()%10+2;
      eltLobatto(i).ConstructFiniteElement(r);
      size += eltLobatto(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleLobatto = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleHierarchic> eltHp(100);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleDgOrtho> eltOrtho(100);
  size = 0;
  for (int i = 0; i < eltOrtho.GetM(); i++)
    {
      int r = rand()%10+2;
      eltOrtho(i).ConstructFiniteElement(r);
      size += eltOrtho(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleDgOrtho = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<TriangleHcurlFirstFamily> eltHcurl(100);
  size = 0;
  for (int i = 0; i < eltHcurl.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurl(i).ConstructFiniteElement(r);
      size += eltHcurl(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHcurlFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleHcurlHierarchic> eltHcurlHp(100);
  size = 0;
  for (int i = 0; i < eltHcurlHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlHp(i).ConstructFiniteElement(r);
      size += eltHcurlHp(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHcurlHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<TriangleHcurlOptimalHpFirstFamily> eltHcurlOpt(100);
  size = 0;
  for (int i = 0; i < eltHcurlOpt.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlOpt(i).ConstructFiniteElement(r);
      size += eltHcurlOpt(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHcurlOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleHcurlSecondFamily> eltHcurlSec(100);
  size = 0;
  for (int i = 0; i < eltHcurlSec.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlSec(i).ConstructFiniteElement(r);
      size += eltHcurlSec(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHcurlSecondFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleHdivFirstFamily> eltHdiv(100);
  size = 0;
  for (int i = 0; i < eltHdiv.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHdiv(i).ConstructFiniteElement(r);
      size += eltHdiv(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHdivFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TriangleHdivOptimalHpFirstFamily> eltHdivOpt(100);
  size = 0;
  for (int i = 0; i < eltHdivOpt.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHdivOpt(i).ConstructFiniteElement(r);
      size += eltHdivOpt(i).GetMemorySize();
    }

  cout << "Memory taken by TriangleHdivOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;*/

  /***********************************************
   * structures in folder FiniteElement/Quadrangle *
   ***********************************************/

  /* Vector<QuadrangleGeomReference> elt(200);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%10+2;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleGeomReference = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleDgGauss> eltDgGauss(10);
  size = 0;
  for (int i = 0; i < eltDgGauss.GetM(); i++)
    {
      int r = rand()%10+2;
      eltDgGauss(i).ConstructFiniteElement(r);
      size += eltDgGauss(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleDgGauss = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleDgOrtho> eltOrtho(10);
  size = 0;
  for (int i = 0; i < eltOrtho.GetM(); i++)
    {
      int r = rand()%10+2;
      eltOrtho(i).ConstructFiniteElement(r);
      size += eltOrtho(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleDgOrtho = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleGauss> eltGauss(10);
  size = 0;
  for (int i = 0; i < eltGauss.GetM(); i++)
    {
      int r = rand()%10+2;
      eltGauss(i).ConstructFiniteElement(r);
      size += eltGauss(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleGauss = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlAxi> eltHcurlAxi(5);
  size = 0;
  for (int i = 0; i < eltHcurlAxi.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlAxi(i).ConstructFiniteElement(r);
      size += eltHcurlAxi(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlAxi = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlFirstFamily> eltHcurl(5);
  size = 0;
  for (int i = 0; i < eltHcurl.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurl(i).ConstructFiniteElement(r);
      size += eltHcurl(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlGaussFirstFamily> eltHcurlGauss(5);
  size = 0;
  for (int i = 0; i < eltHcurlGauss.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlGauss(i).ConstructFiniteElement(r);
      size += eltHcurlGauss(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlGaussFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlHpFirstFamily> eltHcurlHp(5);
  size = 0;
  for (int i = 0; i < eltHcurlHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlHp(i).ConstructFiniteElement(r);
      size += eltHcurlHp(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlLobatto> eltHcurlLob(5);
  size = 0;
  for (int i = 0; i < eltHcurlLob.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlLob(i).ConstructFiniteElement(r);
      size += eltHcurlLob(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlLobatto = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlOptimalFirstFamily> eltHcurlOpt(5);
  size = 0;
  for (int i = 0; i < eltHcurlOpt.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlOpt(i).ConstructFiniteElement(r);
      size += eltHcurlOpt(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlOptimalHpFirstFamily> eltHcurlOptHp(5);
  size = 0;
  for (int i = 0; i < eltHcurlOptHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHcurlOptHp(i).ConstructFiniteElement(r);
      size += eltHcurlOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHcurlOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;


  Vector<QuadrangleHdivFirstFamily> eltHdiv(5);
  size = 0;
  for (int i = 0; i < eltHdiv.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHdiv(i).ConstructFiniteElement(r);
      size += eltHdiv(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHdivFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHdivHpFirstFamily> eltHdivHp(5);
  size = 0;
  for (int i = 0; i < eltHdivHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHdivHp(i).ConstructFiniteElement(r);
      size += eltHdivHp(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHdivHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHcurlOptimalFirstFamily> eltHdivOpt(5);
  size = 0;
  for (int i = 0; i < eltHdivOpt.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHdivOpt(i).ConstructFiniteElement(r);
      size += eltHdivOpt(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHdivOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHdivOptimalHpFirstFamily> eltHdivOptHp(5);
  size = 0;
  for (int i = 0; i < eltHdivOptHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHdivOptHp(i).ConstructFiniteElement(r);
      size += eltHdivOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHdivOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleHierarchic> eltHp(100);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%10+2;
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleLobatto> eltLob(100);
  size = 0;
  for (int i = 0; i < eltLob.GetM(); i++)
    {
      int r = rand()%10+2;
      eltLob(i).ConstructFiniteElement(r);
      size += eltLob(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleLobatto = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<QuadrangleRadau> eltRadau(100);
  size = 0;
  for (int i = 0; i < eltRadau.GetM(); i++)
    {
      int r = rand()%10+2;
      eltRadau(i).SetRadauInterpolationPoints();
      eltRadau(i).ConstructFiniteElement(r);
      size += eltRadau(i).GetMemorySize();
    }

  cout << "Memory taken by QuadrangleRadau = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  */

  /**************************************************
   * structures in folder FiniteElement/Tetrahedron *
   **************************************************/
  
  // problem of size in TetrahedronGeomReference ?
  /* Vector<TetrahedronGeomReference> elt(20);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%7+1;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronGeomReference = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TetrahedronClassical> eltClassical(10);
  size = 0;
  for (int i = 0; i < eltClassical.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltClassical(i).ConstructFiniteElement(r);
      size += eltClassical(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronClassical = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TetrahedronHierarchic> eltHp(10);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<TetrahedronDgOrtho> eltOrtho(10);
  size = 0;
  for (int i = 0; i < eltOrtho.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltOrtho(i).ConstructFiniteElement(r);
      size += eltOrtho(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronDgOrtho = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<TetrahedronHcurlFirstFamily> eltHcurl(10);
  size = 0;
  for (int i = 0; i < eltHcurl.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHcurl(i).ConstructFiniteElement(r);
      size += eltHcurl(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronHcurlFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<TetrahedronHcurlOptimalHpFirstFamily> eltHcurlHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHcurlHp(i).ConstructFiniteElement(r);
      size += eltHcurlHp(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronHcurlOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<TetrahedronHdivFirstFamily> eltHdiv(10);
  size = 0;
  for (int i = 0; i < eltHdiv.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdiv(i).ConstructFiniteElement(r);
      size += eltHdiv(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronHdivFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<TetrahedronHdivOptimalHpFirstFamily> eltHdivHp(10);
  size = 0;
  for (int i = 0; i < eltHdivHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdivHp(i).ConstructFiniteElement(r);
      size += eltHdivHp(i).GetMemorySize();
    }

  cout << "Memory taken by TetrahedronHdivOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */

  /**************************************************
   * structures in folder FiniteElement/Pyramid *
   **************************************************/
  
  // problem of size in PyramidGeomReference ?
  /* Vector<PyramidGeomReference> elt(20);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%7+1;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidGeomReference = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidClassical> eltClassical(10);
  size = 0;
  for (int i = 0; i < eltClassical.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltClassical(i).ConstructFiniteElement(r);
      size += eltClassical(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidClassical = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidHierarchic> eltHp(10);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<PyramidDgOrtho> eltOrtho(10);
  size = 0;
  for (int i = 0; i < eltOrtho.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltOrtho(i).ConstructFiniteElement(r);
      size += eltOrtho(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidDgOrtho = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidDgLegendre> eltLeg(10);
  size = 0;
  for (int i = 0; i < eltLeg.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltLeg(i).ConstructFiniteElement(r);
      size += eltLeg(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidDgLegendre = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidHcurlFirstFamily> eltHcurl(10);
  size = 0;
  for (int i = 0; i < eltHcurl.GetM(); i++)
    {
      int r = rand()%4+1; DISP(r);
      eltHcurl(i).ConstructFiniteElement(r);
      size += eltHcurl(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHcurlFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidHcurlHpFirstFamily> eltHcurlHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHcurlHp(i).ConstructFiniteElement(r);
      size += eltHcurlHp(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHcurlHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidHcurlOptimalFirstFamily> eltHcurlOpt(10);
  size = 0;
  for (int i = 0; i < eltHcurlOpt.GetM(); i++)
    {
      int r = rand()%4+1; DISP(r);
      eltHcurlOpt(i).ConstructFiniteElement(r);
      size += eltHcurlOpt(i).GetMemorySize();
    }


  cout << "Memory taken by PyramidHcurlOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<PyramidHcurlOptimalHpFirstFamily> eltHcurlOptHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlOptHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHcurlOptHp(i).ConstructFiniteElement(r);
      size += eltHcurlOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHcurlOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<PyramidHdivFirstFamily> eltHdiv(10);
  size = 0;
  for (int i = 0; i < eltHdiv.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      //eltHdiv(i).ConstructFiniteElement(r);
      //size += eltHdiv(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHdivFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidHdivHpFirstFamily> eltHdivHp(10);
  size = 0;
  for (int i = 0; i < eltHdivHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdivHp(i).ConstructFiniteElement(r);
      size += eltHdivHp(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHdivHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<PyramidHdivOptimalFirstFamily> eltHdivOpt(10);
  size = 0;
  for (int i = 0; i < eltHdivOpt.GetM(); i++)
    {
      int r = rand()%4+1; DISP(r);
      eltHdivOpt(i).ConstructFiniteElement(r);
      size += eltHdivOpt(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHdivOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<PyramidHdivOptimalHpFirstFamily> eltHdivOptHp(10);
  size = 0;
  for (int i = 0; i < eltHdivOptHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdivOptHp(i).ConstructFiniteElement(r);
      size += eltHdivOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by PyramidHdivOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */


  /**************************************************
   * structures in folder FiniteElement/Wedge *
   **************************************************/
  
  // problem of size in WedgeGeomReference ?
  /* Vector<WedgeGeomReference> elt(20);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%7+1;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeGeomReference = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeClassical> eltClassical(10);
  size = 0;
  for (int i = 0; i < eltClassical.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltClassical(i).ConstructFiniteElement(r);
      size += eltClassical(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeClassical = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeHierarchic> eltHp(10);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<WedgeDgOrtho> eltOrtho(10);
  size = 0;
  for (int i = 0; i < eltOrtho.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltOrtho(i).ConstructFiniteElement(r);
      size += eltOrtho(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeDgOrtho = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeDgLegendre> eltLeg(10);
  size = 0;
  for (int i = 0; i < eltLeg.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltLeg(i).ConstructFiniteElement(r);
      size += eltLeg(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeDgLegendre = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeHcurlFirstFamily> eltHcurl(10);
  size = 0;
  for (int i = 0; i < eltHcurl.GetM(); i++)
    {
      int r = rand()%4+1; DISP(r);
      eltHcurl(i).ConstructFiniteElement(r);
      size += eltHcurl(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHcurlFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeHcurlHpFirstFamily> eltHcurlHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHcurlHp(i).ConstructFiniteElement(r);
      size += eltHcurlHp(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHcurlHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeHcurlOptimalFirstFamily> eltHcurlOpt(10);
  size = 0;
  for (int i = 0; i < eltHcurlOpt.GetM(); i++)
    {
      int r = rand()%4+1; DISP(r);
      eltHcurlOpt(i).ConstructFiniteElement(r);
      size += eltHcurlOpt(i).GetMemorySize();
    }


  cout << "Memory taken by WedgeHcurlOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<WedgeHcurlOptimalHpFirstFamily> eltHcurlOptHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlOptHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHcurlOptHp(i).ConstructFiniteElement(r);
      size += eltHcurlOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHcurlOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<WedgeHdivFirstFamily> eltHdiv(10);
  size = 0;
  for (int i = 0; i < eltHdiv.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdiv(i).ConstructFiniteElement(r);
      size += eltHdiv(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHdivFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeHdivHpFirstFamily> eltHdivHp(10);
  size = 0;
  for (int i = 0; i < eltHdivHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdivHp(i).ConstructFiniteElement(r);
      size += eltHdivHp(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHdivHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<WedgeHdivOptimalFirstFamily> eltHdivOpt(10);
  size = 0;
  for (int i = 0; i < eltHdivOpt.GetM(); i++)
    {
      int r = rand()%4+1; DISP(r);
      eltHdivOpt(i).ConstructFiniteElement(r);
      size += eltHdivOpt(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHdivOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<WedgeHdivOptimalHpFirstFamily> eltHdivOptHp(10);
  size = 0;
  for (int i = 0; i < eltHdivOptHp.GetM(); i++)
    {
      int r = rand()%7+1; DISP(r);
      eltHdivOptHp(i).ConstructFiniteElement(r);
      size += eltHdivOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by WedgeHdivOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */


  /**************************************************
   * structures in folder FiniteElement/Hexahedron *
   **************************************************/
  
  /* Vector<HexahedronGeomReference> elt(20);
  int64_t size = 0;
  for (int i = 0; i < elt.GetM(); i++)
    {
      int r = rand()%10+1;
      elt(i).ConstructFiniteElement(r);
      size += elt(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronGeomReference = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronGauss> eltClassical(10);
  size = 0;
  for (int i = 0; i < eltClassical.GetM(); i++)
    {
      int r = rand()%10+1; DISP(r);
      eltClassical(i).ConstructFiniteElement(r);
      size += eltClassical(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronGauss = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronDgGauss> eltDgGauss(10);
  size = 0;
  for (int i = 0; i < eltDgGauss.GetM(); i++)
    {
      int r = rand()%10+1; DISP(r);
      eltDgGauss(i).ConstructFiniteElement(r);
      size += eltDgGauss(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronDgGauss = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronLobatto> eltLob(10);
  size = 0;
  for (int i = 0; i < eltLob.GetM(); i++)
    {
      int r = rand()%10+1; DISP(r);
      eltLob(i).ConstructFiniteElement(r);
      size += eltLob(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronLobatto = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronHierarchic> eltHp(10);
  size = 0;
  for (int i = 0; i < eltHp.GetM(); i++)
    {
      int r = rand()%10+1; DISP(r);
      eltHp(i).ConstructFiniteElement(r);
      size += eltHp(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHierarchic = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<HexahedronHcurlLobatto> eltOrtho(10);
  size = 0;
  for (int i = 0; i < eltOrtho.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltOrtho(i).ConstructFiniteElement(r);
      size += eltOrtho(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHcurlLobatto = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronDgLegendre> eltLeg(10);
  size = 0;
  for (int i = 0; i < eltLeg.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltLeg(i).ConstructFiniteElement(r);
      size += eltLeg(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronDgLegendre = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronHcurlFirstFamily> eltHcurl(10);
  size = 0;
  for (int i = 0; i < eltHcurl.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHcurl(i).ConstructFiniteElement(r);
      size += eltHcurl(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHcurlFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronHcurlHpFirstFamily> eltHcurlHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlHp.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHcurlHp(i).ConstructFiniteElement(r);
      size += eltHcurlHp(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHcurlHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronHcurlOptimalFirstFamily> eltHcurlOpt(10);
  size = 0;
  for (int i = 0; i < eltHcurlOpt.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHcurlOpt(i).ConstructFiniteElement(r);
      size += eltHcurlOpt(i).GetMemorySize();
    }


  cout << "Memory taken by HexahedronHcurlOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<HexahedronHcurlOptimalHpFirstFamily> eltHcurlOptHp(10);
  size = 0;
  for (int i = 0; i < eltHcurlOptHp.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHcurlOptHp(i).ConstructFiniteElement(r);
      size += eltHcurlOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHcurlOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
  
  Vector<HexahedronHdivFirstFamily> eltHdiv(10);
  size = 0;
  for (int i = 0; i < eltHdiv.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHdiv(i).ConstructFiniteElement(r);
      size += eltHdiv(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHdivFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronHdivHpFirstFamily> eltHdivHp(10);
  size = 0;
  for (int i = 0; i < eltHdivHp.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHdivHp(i).ConstructFiniteElement(r);
      size += eltHdivHp(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHdivHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Vector<HexahedronHdivOptimalFirstFamily> eltHdivOpt(10);
  size = 0;
  for (int i = 0; i < eltHdivOpt.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHdivOpt(i).ConstructFiniteElement(r);
      size += eltHdivOpt(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHdivOptimalFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;
 
  Vector<HexahedronHdivOptimalHpFirstFamily> eltHdivOptHp(10);
  size = 0;
  for (int i = 0; i < eltHdivOptHp.GetM(); i++)
    {
      int r = rand()%6+1; DISP(r);
      eltHdivOptHp(i).ConstructFiniteElement(r);
      size += eltHdivOptHp(i).GetMemorySize();
    }

  cout << "Memory taken by HexahedronHdivOptimalHpFirstFamily = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */

  /*****************************
   * structures in folder Mesh *
   *****************************/

  /* TriangleLobatto tri;
  QuadrangleLobatto quad;
  TetrahedronClassical tetra;
  tetra.ConstructFiniteElement(3);
  tri.ConstructFiniteElement(6);
  quad.ConstructFiniteElement(6);

  cout << "Waiting... " << endl; cin >> test_input;
  
  int64_t size = 0;
  Mesh<Dimension1> mesh1d;
  mesh1d.CreateRegularMesh(-2.0, 2.0, 10000, 1);
  mesh1d.SetOrder(4);
  mesh1d.NumberMesh();
  size = mesh1d.GetMemorySize();
  
  cout << "Memory taken by Mesh1D = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  Mesh<Dimension2> mesh2d;
  mesh2d.SetGeometryOrder(6);
  mesh2d.Read("cacahuete.mesh");
  
  //VectReal_wp step_subdiv(5); step_subdiv.Fill(); Mlt(0.25, step_subdiv);
  //mesh2d.SubdivideMesh(step_subdiv);
  
  size = mesh2d.GetMemorySize();

  cout << "Memory taken by Mesh2D = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  MeshNumbering<Dimension2> mesh_num2d(mesh2d);
  tri.ConstructNumberMap(mesh_num2d.number_map, 0);
  quad.ConstructNumberMap(mesh_num2d.number_map, 0);
  mesh_num2d.SetOrder(6);
  TinyVector<Vector<int>, 4> order;
  order(0).Reallocate(1); order(0)(0) = 6;
  order(1).Reallocate(1); order(1)(0) = 6;
  mesh_num2d.number_map.ConstructQuadrature2D(order, 0);
  mesh_num2d.NumberMesh();
  size = mesh_num2d.GetMemorySize();

  cout << "Memory taken by MeshNumbering2D = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;


  Mesh<Dimension3> mesh3d;
  mesh3d.SetGeometryOrder(3);
  mesh3d.Read("lung.meshb");
  
  //VectReal_wp step_subdiv(5); step_subdiv.Fill(); Mlt(0.25, step_subdiv);
  //mesh2d.SubdivideMesh(step_subdiv);
  
  size = mesh3d.GetMemorySize();

  cout << "Memory taken by Mesh3D = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  MeshNumbering<Dimension3> mesh_num3d(mesh3d);
  tetra.ConstructNumberMap(mesh_num3d.number_map, 0);
  mesh_num3d.SetOrder(3);
  order(0)(0) = 3; order(1)(0) = 3;
  mesh_num3d.number_map.ConstructQuadrature3D(order, 0, 0);
  mesh_num3d.NumberMesh();
  size = mesh_num3d.GetMemorySize();

  cout << "Memory taken by MeshNumbering3D = " << GetHumanReadableMemory(size) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */
  
  /******************************
   * structures in folder Share *
   ******************************/

  /* int N = 2;
  for (int k = 0; k < 20; k++)
    N *= 2;

  DISP(N);
  Vector<Real_wp> xr(N);
  Vector<Complex_wp> x(N), xc(N/2+1);
  x.Zero(); xr.Zero(); xc.Zero();

  cout << "Waiting... " << endl; cin >> test_input;
  

  MultivariatePolynomial<Real_wp> P;
  P.SetOrder(4, 100);
  P.Zero();
      
  cout << "Memory taken by MultivariatePolynomial = " << GetHumanReadableMemory(P.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;

  
  FftInterface<Complex_wp> fft;
  
  DISP(N);
  fft.Init(N);
  fft.ApplyForward(x);
  
  cout << "Memory taken by FftInterface = " << GetHumanReadableMemory(fft.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input;


  FftRealInterface fft_r;
  
  fft_r.Init(N);
  fft_r.ApplyForward(xr, xc);
  
  cout << "Memory taken by FftRealInterface = " << GetHumanReadableMemory(fft_r.GetMemorySize()) << endl;
  cout << "Waiting... " << endl; cin >> test_input; */


  /***********************************
   * structures in folder Quadrature *
   ***********************************/

  

}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  CheckMemorySize();

  return FinalizeMontjoie();
}
