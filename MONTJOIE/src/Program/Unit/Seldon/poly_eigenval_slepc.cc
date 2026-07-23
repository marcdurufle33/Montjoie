#include "Algebra/MontjoieAlgebra.hxx"

using namespace Montjoie;

// (K + lambda S + lambda^2 M) x = 0
template<class T>
class MyNonLinearEigenProblem : public NonLinearEigenProblem_Base<T>
{
  DistributedMatrix<T, Symmetric, ArrayRowSymSparse>* K, *S, *M;
  SparseDistributedSolver<T> mat_lu;
  
public:
  void SetOperators(DistributedMatrix<T, Symmetric, ArrayRowSymSparse>& K2,
                    DistributedMatrix<T, Symmetric, ArrayRowSymSparse>& S2,
                    DistributedMatrix<T, Symmetric, ArrayRowSymSparse>& M2)
  {
    K = &K2; S = &S2; M = &M2;
    this->Init(K2.GetM());

    this->SetSplitMatrices(3);

    /*
    Vector<T> coef(1); coef(0) = 1.0;
    this->SetNumeratorSplitFct(0, coef);
    
    coef.Reallocate(2); coef.Zero(); coef(0) = 1.0;
    this->SetNumeratorSplitFct(1, coef);
    
    coef.Reallocate(3); coef.Zero(); coef(0) = 1.0;
    this->SetNumeratorSplitFct(2, coef); */

    Vector<FN> fct(3); PetscScalar coef0[1], coef1[2], coef2[3];
    FNCreate(MPI_COMM_WORLD, &fct(0));
    FNSetType(fct(0), FNRATIONAL); coef0[0] = 1.0;
    FNRationalSetNumerator(fct(0), 1, coef0);
    FNRationalSetDenominator(fct(0), 0, NULL);
    
    FNCreate(MPI_COMM_WORLD, &fct(1));
    FNSetType(fct(1), FNSQRT);
    //FNSetType(fct(1), FNRATIONAL); coef1[0] = 1.0; coef1[1] = 0.0;
    //FNRationalSetNumerator(fct(1), 2, coef1);
    //FNRationalSetDenominator(fct(1), 0, NULL);

    FNCreate(MPI_COMM_WORLD, &fct(2));
    //FNSetType(fct(2), FNRATIONAL); coef2[0] = 1.0; coef2[1] = 0.0; coef2[2] = 0.0; 
    FNSetType(fct(2), FNRATIONAL); coef1[0] = 1.0; coef1[1] = 0.0;
    FNRationalSetNumerator(fct(2), 2, coef1);
    FNRationalSetDenominator(fct(2), 0, NULL);
    
    this->SetSlepcFunction_NonRational(fct);
  }

  void ComputeOperator(const T& L)
  {}
  
  void MltOperator(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    T one; SetComplexOne(one);
    Mlt(*K, X, Y);
    MltAdd(L, *S, X, one, Y);
    MltAdd(L*L, *M, X, one, Y);
  }

  void ComputeOperatorExplicit(const T& L, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    Copy(*K, A);
    Add(L, *S, A);
    Add(L*L, *M, A);
  }

  void ComputeJacobian(const T& L)
  {}
 
  void MltJacobian(const T& L, const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    T one; SetComplexOne(one);
    Mlt(*S, X, Y);
    MltAdd(Real_wp(2)*L, *M, X, one, Y);
  }

  void ComputeJacobianExplicit(const T& L, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    Copy(*S, A);
    Add(2.0*L, *M, A);
  }

  void ComputePreconditioning(const T& L)
  {
    DistributedMatrix<T, Symmetric, ArrayRowSymSparse> A;
    A = *K;
    Add(L, *S, A);
    Add(L*L, *M, A);
    
    mat_lu.Factorize(A);
  }
  
  void ComputePreconditioning(const Vector<T>& L, const Vector<T>& coef)
  {
    DistributedMatrix<T, Symmetric, ArrayRowSymSparse> A, B;
    for (int k = 0; k < L.GetM(); k++)
      {
        B = *K;
        Add(L(k), *S, B);
        Add(L(k)*L(k), *M, B);

        if (k == 0)
          {
            A = B; Mlt(coef(0), A);
          }
        else
          Add(coef(k), B, A);
      }
    
    mat_lu.Factorize(A);
  }


  void ComputeSplitPreconditioning(const Vector<int>& num, const Vector<T>& coef)
  {
    DistributedMatrix<T, Symmetric, ArrayRowSymSparse> A;
    A.Reallocate(K->GetM(), K->GetM());
    for (int k = 0; k < num.GetM(); k++)
      {
        int p = num(k);
        if (p == 0)
          Add(coef(k), *K, A);
        else if (p == 1)
          Add(coef(k), *S, A);
        else
          Add(coef(k), *M, A);
      }
    
    mat_lu.Factorize(A);
  }
  
  void ApplyPreconditioning(const SeldonTranspose&, const Vector<T>& X, Vector<T>& Y)
  {
    mat_lu.Solve(Y, X);
  }

  void ComputeOperatorSplitExplicit(int i, DistributedMatrix<T, General, ArrayRowSparse>& A)
  {
    if (i == 0)
      Copy(*K, A);
    else if (i == 1)
      Copy(*S, A);
    else
      Copy(*M, A);
  }
  
  void MltOperatorSplit(int i, const SeldonTranspose& trans, const Vector<T>& X, Vector<T>& Y)
  {
    if (i == 0)
      Mlt(*K, X, Y);
    else if (i == 1)
      Mlt(*S, X, Y);
    else
      Mlt(*M, X, Y);
  }
  
};

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  int ncv = 0;
  if (argc > 1)
    ncv = atoi(argv[1]);
  
  if (false)
    {
      // cas de matrices pleines
      int nb_eigenval = 10;
      
      DISP(sizeof(Petsc_Scalar));
      
      int n = 50;
      Matrix<Petsc_Scalar> K(n, n), S(n, n), M(n, n);
      
      srand(0);
      K.FillRand(); K *= 1e-9;
      S.FillRand(); S *= 1e-9;
      M.FillRand(); M *= 1e-9;
      
      M.WriteText("M.dat");
      K.WriteText("K.dat");
      S.WriteText("S.dat");
      
      Vector<Matrix<Petsc_Scalar>* > list_op(3);
      list_op(0) = &K;
      list_op(1) = &S;
      list_op(2) = &M;
      
      PolynomialDenseEigenProblem<Petsc_Scalar, General, RowMajor> var_eig;
      var_eig.SetStoppingCriterion(1e-12);
      var_eig.SetNbAskedEigenvalues(nb_eigenval);
      var_eig.InitMatrix(list_op);
      //var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0.0);
      var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, 0.6);
      var_eig.SetDiagonalMass(false);
      var_eig.SetSpectralTransformation(true);
      //SlepcParamPep& param = var_eig.GetSlepcParameters();
      //param.SetEigensolverType(param.QARNOLDI);
      
      Vector<Petsc_Scalar> lambda, lambda_imag;
      Matrix<Petsc_Scalar, General, ColMajor> eigen_vec;
      
      FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);
      
      DISP(lambda); DISP(lambda_imag);
    }

  if (false)
    {
      int nb_eigenval = 3;
      // cas de matrices creuses

      // (lambda^2 M + lambda S + K) x = 0
      // M : matrice de masse
      // lambda = sigma + 1/mu => pour calculer les valeurs propres proches de sigma, on calcule les grandes valeurs propres de mu
      DistributedMatrix<Petsc_Scalar, General, ArrayRowSparse> K, S, M;
      K.ReadText("test/matrix/KhPol.dat");
      S.ReadText("test/matrix/ShPol.dat");
      M.ReadText("test/matrix/MhPol.dat");
      
      Vector<DistributedMatrix<Petsc_Scalar, General, ArrayRowSparse>* > list_op(2);
      list_op(0) = &K;
      list_op(1) = &S;
      
      PolynomialSparseEigenProblem<Petsc_Scalar, DistributedMatrix<Petsc_Scalar, General, ArrayRowSparse>, DistributedMatrix<Petsc_Scalar, General, ArrayRowSparse> > var_eig;

      var_eig.SetStoppingCriterion(1e-12);
      var_eig.SetNbAskedEigenvalues(nb_eigenval);
      var_eig.InitMatrix(list_op, M); // K et S sont dans list_op et M dans M
      var_eig.SetPrintLevel(4);
      
      complex<double> center(0.4, -0.15);
      var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, center);
      //FeastParamPep& param = var_eig.GetFeastParameters();
      //param.SetCircleSpectrum(center, 0.08);
      var_eig.SetSpectralTransformation(true);

      Vector<Petsc_Scalar> lambda, lambda_imag;
      Matrix<Petsc_Scalar, General, ColMajor> eigen_vec;
      
      FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);
      //FindEigenvaluesFeast(var_eig, lambda, lambda_imag, eigen_vec);
      
      DISP(lambda); DISP(lambda_imag);
    }

  if (true)
    {
      // cas d'un probleme aux valeurs propres non-lineaire
      int nb_eigenval = 1;
      
      DistributedMatrix<Petsc_Scalar, Symmetric, ArrayRowSymSparse> K, S, M;
      K.ReadText("KhRef.dat");
      S.ReadText("ShRef.dat");
      M.ReadText("MhRef.dat");
      
      MyNonLinearEigenProblem<Petsc_Scalar> var_eig;
      var_eig.SetOperators(K, S, M);
      //var_eig.SetExplicitMatrix();
      
      var_eig.SetStoppingCriterion(1e-12);
      var_eig.SetNbAskedEigenvalues(nb_eigenval);
      if (ncv > 0)
        var_eig.SetNbArnoldiVectors(ncv);
      
      //var_eig.SetTypeSpectrum(var_eig.LARGE_EIGENVALUES, 0.0);
      //var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, Complex_wp(0.4, -0.15));
      var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, Complex_wp(0.1375, -0.12));
      SlepcParamNep& param = var_eig.GetSlepcParameters();
      //param.SetEigensolverType(param.RII);
      //param.SetEigensolverType(param.SLP);
      //param.SetEigensolverType(param.NARNOLDI);
      param.SetEigensolverType(param.NLEIGS);
      //param.SetIntervalRegion(0.3, 0.5, -0.2, -0.1);
      param.SetIntervalRegion(0.1, 0.2, -0.2, -0.1);
      
      Vector<Petsc_Scalar> lambda, lambda_imag;
      Matrix<Petsc_Scalar, General, ColMajor> eigen_vec;
      
      DISP(var_eig.GetM());
      FindEigenvaluesSlepc(var_eig, lambda, lambda_imag, eigen_vec);
      
      DISP(lambda); DISP(lambda_imag);
    }
      
  return FinalizeMontjoie();
}
