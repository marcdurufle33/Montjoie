#ifndef SELDON_FILE_SLEPC_CXX

#include <slepceps.h>
#include <slepcpep.h>
#include <slepcnep.h>
#include <slepc/private/stimpl.h>
#include "petsc/src/mat/impls/aij/seq/aij.h"

typedef struct {
  PetscInt    nmat,maxnmat;
  PetscScalar *coeff;
  Mat         *A;
  Vec         t;
} NEP_NLEIGS_MATSHELL;


namespace Seldon
{

  const char* SlepcParam::GetEigensolverChar() const
  {
    switch (type_solver)
      {
      case POWER : return EPSPOWER;
      case SUBSPACE : return EPSSUBSPACE;
      case ARNOLDI : return EPSARNOLDI;
      case LANCZOS : return EPSLANCZOS;
      case KRYLOVSCHUR : return EPSKRYLOVSCHUR;
      case GD : return EPSGD;
      case JD : return EPSJD;
      case RQCG : return EPSRQCG;
      case LOBPCG : return EPSLOBPCG;
      case CISS : return EPSCISS;
      case LAPACK : return EPSLAPACK;
      case ARPACK : return EPSARPACK;
      case TRLAN : return EPSTRLAN;
      case BLOPEX : return EPSBLOPEX;
      case PRIMME : return EPSPRIMME;
        //case FEAST : return EPSFEAST;
      }
    
    return "";
  }


  void SetParametersSlepc(const SlepcParam& param, EPS& eps)
  {
    int ierr = EPSSetType(eps, param.GetEigensolverChar());
    if (ierr != 0)
      {
	cout << "Chosen type = " << param.GetEigensolverChar() << " Not enabled in Slepc ?" << endl;
	abort();
      }

    int solver = param.GetEigensolverType();
    if (solver == param.BLOPEX)
      {
#ifdef SLEPC_HAVE_BLOPEX
	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSBLOPEXSetBlockSize(eps, bs);
#else
	cout << "Slepc not compiled with BLOPEX" << endl;
	abort();
#endif
      }
    else if (solver == param.CISS)
      {
	int type_extraction = param.GetExtractionType();
	if (type_extraction >= 0)
	  {
	    if (type_extraction == param.EXTRACT_RITZ)
	      EPSCISSSetExtraction(eps, EPS_CISS_EXTRACTION_RITZ);
	    else
	      EPSCISSSetExtraction(eps, EPS_CISS_EXTRACTION_HANKEL);
	  }
	
	int type_quad = param.GetQuadratureRuleType();
	if (type_quad >= 0)
	  {
	    if (type_quad == param.QUADRULE_TRAPEZE)
	      EPSCISSSetQuadRule(eps, EPS_CISS_QUADRULE_TRAPEZOIDAL);
	    else
	      EPSCISSSetQuadRule(eps, EPS_CISS_QUADRULE_CHEBYSHEV);
	  }
	
	PetscScalar a; PetscBool complex_number = PETSC_FALSE;
	if (IsComplexNumber(a))
	  complex_number = PETSC_TRUE;
	
	if (param.GetInnerSteps() > 0)
	  EPSCISSSetRefinement(eps, param.GetInnerSteps(), param.GetOuterSteps());
	
	if (param.GetNumberIntegrationPoints() > 0)
	  EPSCISSSetSizes(eps, param.GetNumberIntegrationPoints(), param.GetBlockSize(),
			  param.GetMomentSize(), param.GetNumberPartitions(),
			  param.GetMaximumBlockSize(), complex_number);
	
	if (param.GetThresholdRank() > 0)
	  EPSCISSSetThreshold(eps, param.GetThresholdRank(), param.GetThresholdSpurious());
      }
    else if (solver == param.FEAST)
      {
#ifdef SLEPC_HAVE_FEAST
	if (param.GetNumberIntegrationPoints() > 0)
	  EPSFEASTSetNumPoints(eps, param.GetNumberIntegrationPoints());
#else
	cout << "Slepc not compiled with FEAST" << endl;
	abort();
#endif
      }    
    else if (solver == param.GD)
      {
	if (param.GetBorthogonalization() >= 0)
	  {
	    PetscBool borth = PETSC_FALSE;
	    if (param.GetBorthogonalization() >= 1)
	      borth = PETSC_TRUE;
	    
	    EPSGDSetBOrth(eps, borth);
	  }

	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSGDSetBlockSize(eps, bs);

	if (param.GetDoubleExpansion() >= 0)
	  {
	    PetscBool exp = PETSC_FALSE;
	    if (param.GetDoubleExpansion() >= 1)
	      exp = PETSC_TRUE;
	    
	    EPSGDSetDoubleExpansion(eps, exp);
	  }
	
	if (param.GetInitialSize() > 0)
	  EPSGDSetInitialSize(eps, param.GetInitialSize());

	if (param.GetKrylovRestart() >= 0)
	  {
	    PetscBool restart = PETSC_FALSE;
	    if (param.GetKrylovRestart() >= 1)
	      restart = PETSC_TRUE;

	    EPSGDSetKrylovStart(eps, restart);
	  }

	if (param.GetRestartNumber() > 0)
	  EPSGDSetRestart(eps, param.GetRestartNumber(), param.GetRestartNumberAdd());
	
        // deprecated function : SetWindowSizes 
	//if (param.GetNumberConvergedVectors() > 0)
        // EPSGDSetWindowSizes(eps, param.GetNumberConvergedVectors(),
	//		      param.GetNumberConvergedVectorsProjected());
      }
    else if (solver == param.JD)
      {
	if (param.GetBorthogonalization() >= 0)
	  {
	    PetscBool borth = PETSC_FALSE;
	    if (param.GetBorthogonalization() >= 1)
	      borth = PETSC_TRUE;

	    EPSJDSetBOrth(eps, borth);
	  }

	PetscInt bs = param.GetBlockSize();
	if (bs > 0)
	  EPSJDSetBlockSize(eps, bs);

	if (param.GetInitialSize() > 0)
	  EPSJDSetInitialSize(eps, param.GetInitialSize());

	if (param.GetKrylovRestart() >= 0)
	  {
	    PetscBool restart = PETSC_FALSE;
	    if (param.GetKrylovRestart() >= 1)
	      restart = PETSC_TRUE;

	    EPSJDSetKrylovStart(eps, restart);
	  }

	if (param.GetRestartNumber() > 0)
	  EPSJDSetRestart(eps, param.GetRestartNumber(), param.GetRestartNumberAdd());

        // deprecated function : SetWindowSizes 	
	//if (param.GetNumberConvergedVectors() > 0)
	//  EPSJDSetWindowSizes(eps, param.GetNumberConvergedVectors(),
	//		      param.GetNumberConvergedVectorsProjected());
      }
    else if (solver == param.KRYLOVSCHUR)
      {
	if (param.UseNonLockingVariant())
	  EPSKrylovSchurSetLocking(eps, PETSC_FALSE);
	else
	  EPSKrylovSchurSetLocking(eps, PETSC_TRUE);

	if (param.GetRestartRatio() > 0)
	  EPSKrylovSchurSetRestart(eps, param.GetRestartRatio());	
      }
    else if (solver == param.LOBPCG)
      {
	if (param.GetBlockSize() > 0)
	  EPSLOBPCGSetBlockSize(eps, param.GetBlockSize());
	
	if (param.UseNonLockingVariant())
	  EPSLOBPCGSetLocking(eps, PETSC_FALSE);
	else
	  EPSLOBPCGSetLocking(eps, PETSC_TRUE);

	if (param.GetRestartRatio() > 0)
	  EPSLOBPCGSetRestart(eps, param.GetRestartRatio());	
      }
    else if (solver == param.PRIMME)
      {
#ifdef SLEPC_HAVE_PRIMME
	if (param.GetBlockSize() > 0)
	  EPSPRIMMESetBlockSize(eps, param.GetBlockSize());

	if (param.GetMethod().size() > 1)
	  {
	    if (param.GetMethod() == "DYNAMIC")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_DYNAMIC);
	    else if (param.GetMethod() == "DEFAULT_MIN_TIME")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_DEFAULT_MIN_TIME);
	    else if (param.GetMethod() == "DEFAULT_MIN_MATVECS")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_DEFAULT_MIN_MATVECS);
	    else if (param.GetMethod() == "ARNOLDI")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_ARNOLDI);
	    else if (param.GetMethod() == "GD")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_GD);
	    else if (param.GetMethod() == "GD_PLUSK")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_GD_PLUSK);
	    else if (param.GetMethod() == "GD_OLSEN_PLUSK")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_GD_OLSEN_PLUSK);
	    else if (param.GetMethod() == "JD_OLSEN_PLUSK")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JD_OLSEN_PLUSK);
	    else if (param.GetMethod() == "RQI")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_RQI);
	    else if (param.GetMethod() == "JDQR")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JDQR);
	    else if (param.GetMethod() == "JDQMR")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JDQMR);
	    else if (param.GetMethod() == "JDQMR_ETOL")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_JDQMR_ETOL);
	    else if (param.GetMethod() == "SUBSPACE_ITERATION")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_SUBSPACE_ITERATION);
	    else if (param.GetMethod() == "LOBPCG_ORTHOBASIS")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_LOBPCG_ORTHOBASIS);
	    else if (param.GetMethod() == "LOBPCG_ORTHOBASISW")
	      EPSPRIMMESetMethod(eps, EPS_PRIMME_LOBPCG_ORTHOBASISW);
	  }
#else
	cout << "Slepc not compiled with PRIMME" << endl;
	abort();
#endif
      }
    else if (solver == param.POWER)
      {
	if (param.GetShiftType() >= 0)
	  {
	    if (param.GetShiftType() == param.SHIFT_CONSTANT)
	      EPSPowerSetShiftType(eps, EPS_POWER_SHIFT_CONSTANT);
	    else if (param.GetShiftType() == param.SHIFT_RAYLEIGH)
	      EPSPowerSetShiftType(eps, EPS_POWER_SHIFT_RAYLEIGH);
	    else if (param.GetShiftType() == param.SHIFT_WILKINSON)
	      EPSPowerSetShiftType(eps, EPS_POWER_SHIFT_WILKINSON);
	  }
      }
    else if (solver == param.RQCG)
      {
	if (param.GetNumberOfSteps() > 0)
	  EPSRQCGSetReset(eps, param.GetNumberOfSteps());	
      }
  }

  
  //! filling the vector y from pointer contained in x (Petsc vector)
  void CopyPointerPetsc(const Vec& x, Vector<PetscScalar>& y)
  {
    // it is assumed that the vector x is stored in a contiguous array
    PetscInt n;
    VecGetLocalSize(x, &n);
    
    PetscScalar* x_array;
    VecGetArrayRead(x, (const PetscScalar**)&x_array);
    
    y.SetData(n, x_array);
  }


  void AllocatePetscVector(const MPI_Comm& comm, Vec& Vr, int n, int nglob,
			   Vector<PetscScalar>& Vr_vec)
  {
    VecCreate(comm, &Vr);
    VecSetSizes(Vr, n, nglob);
    VecSetFromOptions(Vr);

    CopyPointerPetsc(Vr, Vr_vec);
  }

  
  //! matrix-vector product y = mat x (mat : stiffness operator)
  PetscErrorCode MatMult_Matrix(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec0, yvec;
    CopyPointerPetsc(x, xvec0);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    EigenProblem_Base<PetscScalar>& var_eig
      = *reinterpret_cast<EigenProblem_Base<PetscScalar>* >(ctx);

    var_eig.IncrementProdMatVect();
    Vector<PetscScalar> xvec(xvec0);
    if (var_eig.DiagonalMass() || (var_eig.UseCholeskyFactoForMass()))
      {
	// standard eigenvalue problem
	if (var_eig.GetComputationalMode() == var_eig.REGULAR_MODE)
	  {
	    if (var_eig.DiagonalMass())
	      var_eig.MltInvSqrtDiagonalMass(xvec);
	    else
	      var_eig.SolveCholeskyMass(Seldon::SeldonTrans, xvec);
	    
	    var_eig.MltStiffness(xvec, yvec);
            var_eig.IncrementProdMatVect();
	    
	    if (var_eig.DiagonalMass())
	      var_eig.MltInvSqrtDiagonalMass(yvec);
	    else
	      {
                var_eig.SolveCholeskyMass(Seldon::SeldonNoTrans, yvec);
                var_eig.IncrementLinearSolves();
              }
	  }
	else
	  {
	    if (var_eig.DiagonalMass())
	      var_eig.MltSqrtDiagonalMass(xvec);
	    else
	      var_eig.MltCholeskyMass(Seldon::SeldonNoTrans, xvec);
	    
	    var_eig.ComputeSolution(xvec, yvec);
            var_eig.IncrementLinearSolves();
            
	    if (var_eig.DiagonalMass())
	      var_eig.MltSqrtDiagonalMass(yvec);
	    else
	      {
                var_eig.MltCholeskyMass(Seldon::SeldonTrans, yvec);
                var_eig.IncrementLinearSolves();
              }
	  }	
      }
    else
      {
	if (var_eig.GetComputationalMode() == var_eig.INVERT_MODE)
	  {
	    if (var_eig.GetTypeSpectrum() != var_eig.CENTERED_EIGENVALUES)
	      {
		var_eig.MltStiffness(xvec0, xvec);
		var_eig.ComputeSolution(xvec, yvec);
	      }
	    else
	      {
		var_eig.MltMass(xvec0, xvec);
		var_eig.ComputeSolution(xvec, yvec);
	      }

            var_eig.IncrementLinearSolves();
            var_eig.IncrementProdMatVect();
	  }
	else if (var_eig.GetComputationalMode() == var_eig.REGULAR_MODE)
	  {
	    var_eig.MltStiffness(xvec, yvec);
            var_eig.IncrementProdMatVect();
	  }
	else
	  {
	    cout << "Not implemented" << endl;
	    abort();
	  }
      }

    xvec0.Nullify(); yvec.Nullify();
    return 0;
  }


  //! matrix-vector product y = mat x (mat : mass operator)
  PetscErrorCode MatMult_MassMatrix(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    EigenProblem_Base<PetscScalar>& var_eig
      = *reinterpret_cast<EigenProblem_Base<PetscScalar>* >(ctx);

    var_eig.MltMass(xvec, yvec);
    var_eig.IncrementProdMatVect();
 
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }


  //! matrix-vector product with transpose, y = mat^T x
  PetscErrorCode MatMultTranspose_Matrix(Mat A, Vec x, Vec y)
  {
    throw Undefined("Function MatMultTranspose_Matrix not implemented");
    return 0;
  }

  
  //! retrieves diagonal of A, d = diag(A)
  PetscErrorCode MatGetDiagonal_Matrix(Mat A, Vec d)
  {
    throw Undefined("Function MatGetDiagonal_Matrix not implemented");
    return 0;
  }


  template<class T>
  bool PutEigenpairLapackForm(int num, int nev, T& Lr, T& Li, Vector<T>& Vr, Vector<T>& Vi,
			      T& Lr_next, T& Li_next, Vector<T>& Vr_next, Vector<T>& Vi_next,
			      Vector<T>& eigen_values, Vector<T>& lambda_imag,
			      Matrix<T, General, ColMajor>& eigen_vectors)
  {
    bool eigen_pair = false;
    if ((Li != T(0)) && (num < nev-1))
      eigen_pair = true;

    int n = Vr.GetM();
    if (eigen_pair)
      {
	eigen_values(num) = Lr;
	lambda_imag(num) = Li;
	eigen_values(num+1) = Lr_next;
	lambda_imag(num+1) = Li_next;
	for (int j = 0; j < n; j++)	  
	  {
	    eigen_vectors(j, num) = Vr(j);
	    eigen_vectors(j, num+1) = Vi(j);
	  }
      }
    else
      {
	eigen_values(num) = Lr;
	lambda_imag(num) = Li;
	for (int j = 0; j < n; j++)
	  eigen_vectors(j, num) = Vr(j);
      }

    return eigen_pair;
  }


  template<class T>
  bool PutEigenpairLapackForm(int num, int nev, complex<T>& Lr, complex<T>& Li,
			      Vector<complex<T> >& Vr, Vector<complex<T> >& Vi,
			      complex<T>& Lr_next, complex<T>& Li_next,
			      Vector<complex<T> >& Vr_next, Vector<complex<T> >& Vi_next,
			      Vector<complex<T> >& eigen_values, Vector<complex<T> >& lambda_imag,
			      Matrix<complex<T>, General, ColMajor>& eigen_vectors)
  {
    int n = Vr.GetM();
    eigen_values(num) = Lr;
    lambda_imag(num) = Li;
    for (int j = 0; j < n; j++)
      eigen_vectors(j, num) = Vr(j);
    
    return false;
  }


  void FindEigenvaluesSlepc_(EigenProblem_Base<PetscScalar>& var,
			     Vector<PetscScalar>& eigen_values,
			     Vector<PetscScalar>& lambda_imag,
			     Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    // initializing of computation
    PetscScalar shiftr = var.GetShiftValue(), shifti = var.GetImagShiftValue();    
    PetscScalar zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
        
    int print_level = var.GetPrintLevel();
    SlepcParam& param = var.GetSlepcParameters();
    int rank_proc(0);
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(var.GetCommunicator(), &rank_proc);
#endif

    Mat stiff, mass;
    int nev = var.GetNbAskedEigenvalues();
    int ncv = var.GetNbArnoldiVectors();
    int n = var.GetM();
    int nglob = var.GetGlobalM();

    // creation of a matrix-free Petsc structure
#ifdef SELDON_WITH_MPI
    MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
		   reinterpret_cast<void*>(&var), &stiff);

    MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
		   reinterpret_cast<void*>(&var), &mass);
#else
    MatCreateShell(PETSC_COMM_SELF, n, n, nglob, nglob,
		   reinterpret_cast<void*>(&var), &stiff);

    MatCreateShell(PETSC_COMM_SELF, n, n, nglob, nglob,
		   reinterpret_cast<void*>(&var), &mass);
#endif    

    //MatSetFromOptions(stiff);     MatSetFromOptions(mass);
    MatShellSetOperation(stiff, MATOP_MULT, (void(*)())MatMult_Matrix);
    MatShellSetOperation(stiff, MATOP_MULT_TRANSPOSE, (void(*)())MatMultTranspose_Matrix);
    //MatShellSetOperation(stiff, MATOP_GET_DIAGONAL, (void(*)())MatGetDiagonal_Matrix);

    MatShellSetOperation(mass, MATOP_MULT, (void(*)())MatMult_MassMatrix);
    
    // creation of the eigensolver
    EPS eps;
#ifdef SELDON_WITH_MPI
    EPSCreate(var.GetCommunicator(), &eps);
#else
    EPSCreate(PETSC_COMM_SELF, &eps);
#endif
    SetParametersSlepc(param, eps);
    
    // type of eigenproblem (hermitian/generalized)
    bool generalized = true;
    bool isherm = var.IsHermitianProblem();
    if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
      generalized = false;
    else if (var.GetComputationalMode() == var.INVERT_MODE)
      {
	generalized = false;
	isherm = false;
      }
    
    if (generalized)
      EPSSetOperators(eps, stiff, mass);
    else
      EPSSetOperators(eps, stiff, NULL);

    if (isherm)
      {
	if (generalized)
	  EPSSetProblemType(eps, EPS_GHEP);
	else
	  EPSSetProblemType(eps, EPS_HEP);
      }
    else
      {
	if (generalized)
	  EPSSetProblemType(eps, EPS_PGNHEP);
	else
	  EPSSetProblemType(eps, EPS_NHEP);
      }

    EPSWhich which(EPS_LARGEST_MAGNITUDE);
    switch (var.GetTypeSorting())
      {
      case EigenProblem_Base<PetscScalar>::SORTED_REAL : which = EPS_LARGEST_REAL; break;
      case EigenProblem_Base<PetscScalar>::SORTED_IMAG : which = EPS_LARGEST_IMAGINARY; break;
      case EigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = EPS_LARGEST_MAGNITUDE; break;
      case EigenProblem_Base<PetscScalar>::SORTED_USER : which = EPS_WHICH_USER; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case EigenProblem_Base<PetscScalar>::SORTED_REAL : which = EPS_SMALLEST_REAL; break;
          case EigenProblem_Base<PetscScalar>::SORTED_IMAG : which = EPS_SMALLEST_IMAGINARY; break;
          case EigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = EPS_SMALLEST_MAGNITUDE; break;
          }
      }

    if ((var.GetComputationalMode() == var.REGULAR_MODE) &&
	(var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES))
      {
	PetscScalar target = shiftr;
        switch (var.GetTypeSorting())
          {
          case EigenProblem_Base<PetscScalar>::SORTED_REAL : which = EPS_TARGET_REAL; break;
          case EigenProblem_Base<PetscScalar>::SORTED_IMAG : which = EPS_TARGET_IMAGINARY; target = shifti; break;
          case EigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = EPS_TARGET_MAGNITUDE; break;
          }
	
	EPSSetTarget(eps, target);
      }
    
    EPSSetWhichEigenpairs(eps, which);
    if (which == EPS_WHICH_USER)
      EPSSetEigenvalueComparison(eps, &EigenProblem_Base<PetscScalar>::GetComparisonEigenvalueSlepc, &var);
    
    double tol = var.GetStoppingCriterion();
    int nb_max_iter = var.GetNbMaximumIterations();
    EPSSetTolerances(eps, tol, nb_max_iter);
    EPSSetDimensions(eps, nev, ncv, PETSC_DEFAULT);
    EPSSetFromOptions(eps);
    
    // computing needed operators
    if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
      {
	if (var.DiagonalMass())
	  {
            // computation of M
            var.ComputeDiagonalMass();
	    
            // computation of M^{-1/2}
            var.FactorizeDiagonalMass();
          }
        else 
          {
            // computation of M for Cholesky factorisation
            var.ComputeMassForCholesky();
            
            // computation of Cholesky factorisation M = L L^T
            var.FactorizeCholeskyMass();
          }

	if (var.GetComputationalMode() != var.REGULAR_MODE)
	  var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);
      }
    else
      {
	if (var.GetComputationalMode() == var.INVERT_MODE)
          {
	    if (var.GetTypeSpectrum() != var.CENTERED_EIGENVALUES)
              {
                // computation and factorisation of mass matrix
                var.ComputeAndFactorizeStiffnessMatrix(one, zero);
                
                // computation of stiffness matrix
                var.ComputeStiffnessMatrix();
              }
            else
              {
                // computation and factorization of (K - sigma M)
                var.ComputeAndFactorizeStiffnessMatrix(-shiftr, one);

                // computation of M
                var.ComputeMassMatrix();
	      }
	  }
	else if (var.GetComputationalMode() == var.REGULAR_MODE)
          {
	    // factorization of the mass matrix
            var.ComputeAndFactorizeStiffnessMatrix(one, zero);
					 
            // computation of stiffness and mass matrix
            var.ComputeStiffnessMatrix();
            var.ComputeMassMatrix();

	    cout << "Case not implemented" << endl;
	    abort();
	  }
	else
	  {
	    cout << "Case not implemented" << endl;
	    abort();
	  }
      }

    // the eigenvalue problem is solved
    int ierr = EPSSolve(eps);
    if (ierr != 0)
      {
	cout << "Error during solution of eigensystem =  " << ierr << endl;
	abort();
      }

    if (print_level >= 2)
      {
	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
	EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
	PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
      }
    
    EPSConvergedReason reason;
    EPSGetConvergedReason(eps, &reason);
    if (reason < 0)
      {
	if (rank_proc == 0)
	  cout << "The solver did not converge" << endl;
        
        if (reason == EPS_DIVERGED_ITS)
          throw SolverMaximumIterationError("FindEigenvaluesSlepc", "Maximum number of iterations reached");
        else
          throw SolverDivergenceError("FindEigenvaluesSlepc", "The solver diverged");        
      }
    
    // eigenvalues and eigenvectors are extracted
    Vec Vr, Vi, Vr_next, Vi_next;
    Vector<PetscScalar> Vr_vec, Vi_vec, Vr_vec_next, Vi_vec_next;
#ifdef SELDON_WITH_MPI
    AllocatePetscVector(var.GetCommunicator(), Vr, n, nglob, Vr_vec);
    AllocatePetscVector(var.GetCommunicator(), Vi, n, nglob, Vi_vec);
    AllocatePetscVector(var.GetCommunicator(), Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(var.GetCommunicator(), Vi_next, n, nglob, Vi_vec_next);
#else
    AllocatePetscVector(PETSC_COMM_SELF, Vr, n, nglob, Vr_vec);
    AllocatePetscVector(PETSC_COMM_SELF, Vi, n, nglob, Vi_vec);
    AllocatePetscVector(PETSC_COMM_SELF, Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(PETSC_COMM_SELF, Vi_next, n, nglob, Vi_vec_next);
#endif
    
    eigen_values.Reallocate(nev);
    lambda_imag.Reallocate(nev);    
    eigen_vectors.Reallocate(n, nev);
    int num = 0;
    PetscScalar Lr, Li, Lr_next, Li_next;
    bool eigen_pair = true;
    while (num < nev)
      {
	if (eigen_pair)
	  EPSGetEigenpair(eps, num, &Lr, &Li, Vr, Vi);

	if (num < nev-1)
	  EPSGetEigenpair(eps, num+1, &Lr_next, &Li_next, Vr_next, Vi_next);
	
	eigen_pair = PutEigenpairLapackForm(num, nev, Lr, Li, Vr_vec, Vi_vec,
					    Lr_next, Li_next, Vr_vec_next, Vi_vec_next,
					    eigen_values, lambda_imag, eigen_vectors);
	
	if (eigen_pair)
	  num += 2;
	else
	  {
	    Lr = Lr_next; Li = Li_next;
	    Vr_vec = Vr_vec_next; Vi_vec = Vi_vec_next;
	    num++;
	  }
      }

    Vr_vec.Nullify();
    Vi_vec.Nullify();
    Vr_vec_next.Nullify();
    Vi_vec_next.Nullify();

    // modifies eigenvalues and eigenvectors if needed
    ApplyScalingEigenvec(var, eigen_values, lambda_imag, eigen_vectors,
                         shiftr, shifti);
    
    // temporary objects are destroyed
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    EPSDestroy(&eps);
    MatDestroy(&stiff);
    MatDestroy(&mass);

    // clears eigenproblem
    var.Clear();

  }
  
  /*****************************
   * Interface with PEP solver *
   *****************************/


  void ApplySpectralTransform(double shift, Vector<double>& eigen_values, Vector<double>& lambda_imag)
  {
    for (int i = 0; i < eigen_values.GetM(); i++)
      {
        if (lambda_imag(i) == 0.0)
          eigen_values(i) = shift + 1.0/eigen_values(i);
        else
          {
            complex<double> val(eigen_values(i), lambda_imag(i));
            complex<double> z = shift + 1.0/val;
            eigen_values(i) = real(z); lambda_imag(i) = imag(z);
          }
      }
  }

  void ApplySpectralTransform(complex<double> shift, Vector<complex<double> >& eigen_values, Vector<complex<double> >& lambda_imag)
  {
    for (int i = 0; i < eigen_values.GetM(); i++)
      eigen_values(i) = shift + 1.0/eigen_values(i);
  }

  
  void SetParametersSlepc(const SlepcParamPep& param, PEP& pep)
  {
    switch(param.GetEigensolverType())
      {
      case SlepcParamPep::TOAR :  PEPSetType(pep, PEPTOAR); break;
      case SlepcParamPep::STOAR :  PEPSetType(pep, PEPSTOAR); break;
      case SlepcParamPep::QARNOLDI :  PEPSetType(pep, PEPQARNOLDI); break;
      case SlepcParamPep::LINEAR :  PEPSetType(pep, PEPLINEAR); break;
      case SlepcParamPep::JD :  PEPSetType(pep, PEPJD); break;
      }        
  }
  
  struct MySlepcOperator
  {
    PolynomialEigenProblem_Base<PetscScalar>* var;
    int num_op;    
  };

  //! matrix-vector product y = op x
  PetscErrorCode MatMult_PepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    MySlepcOperator& op
      = *reinterpret_cast<MySlepcOperator*>(ctx);
    
    op.var->MltOperator(op.num_op, SeldonNoTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! matrix-vector product y = op^T x
  PetscErrorCode MatMult_PepOperatorTrans(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    MySlepcOperator& op
      = *reinterpret_cast<MySlepcOperator*>(ctx);
    
    op.var->MltOperator(op.num_op, SeldonTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! solving op y = x
  PetscErrorCode MatSolve_PepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    MySlepcOperator& op
      = *reinterpret_cast<MySlepcOperator*>(ctx);
    
    if (op.var->UseSpectralTransformation())
      op.var->SolveOperator(SeldonNoTrans, xvec, yvec);
    else
      op.var->SolveMass(SeldonNoTrans, xvec, yvec);
    
    op.var->IncrementLinearSolves();
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! solving op y = x
  PetscErrorCode MatSolveTrans_PepOperator(Mat mat, Vec x, Vec y)
  {
    cout << "Not implemented" << endl;
    abort();
  }
  
  //! functions to compute eigenvalues with PEP solver of Slepc
  void FindEigenvaluesSlepc_(PolynomialEigenProblem_Base<PetscScalar>& var,
                             Vector<PetscScalar>& eigen_values,
                             Vector<PetscScalar>& lambda_imag,
                             Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    // initializing of computation
    PetscScalar shift = var.GetShiftValue();
    PetscScalar zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);

    // degree of polynom 
    int deg_pol = var.GetPolynomialDegree();
    int nev = var.GetNbAskedEigenvalues();
    int n = var.GetM();
    int nglob = var.GetGlobalM();
    
    // creation of shell matrices
    Vector<Mat> op(deg_pol+1);
    Vector<MySlepcOperator> my_op(deg_pol+1);
    for (int k = 0; k <= deg_pol; k++)
      {
        my_op(k).var = &var;
        my_op(k).num_op = k;
#ifdef SELDON_WITH_MPI
        MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                       reinterpret_cast<void*>(&my_op(k)), &op(k));
#else
        MatCreateShell(PETSC_COMM_SELF, n, n, nglob, nglob,
                       reinterpret_cast<void*>(&my_op(k)), &op(k));
#endif
        
        MatShellSetOperation(op(k), MATOP_MULT, (void(*)())MatMult_PepOperator);
        MatShellSetOperation(op(k), MATOP_MULT_TRANSPOSE, (void(*)())MatMult_PepOperatorTrans);        
      }

    // creation of the eigensolver
    PEP solver; ST st;
#ifdef SELDON_WITH_MPI
    PEPCreate(var.GetCommunicator(), &solver);
#else
    PEPCreate(PETSC_COMM_SELF, &solver);
#endif
    
    PEPSetOperators(solver, deg_pol+1, op.GetData());
    PEPGetST(solver, &st);
    st->matsolve = (PetscErrorCode(*)(Mat, Vec, Vec))MatSolve_PepOperator;
    st->matsolve_trans = (PetscErrorCode(*)(Mat, Vec, Vec))MatSolveTrans_PepOperator;
    //st->D = NULL;
    STSetMatMode(st, ST_MATMODE_SHELL);
    //STSetType(st, STSINVERT);
    
    // sorting and selection of spectrum
    PEPWhich which(PEP_LARGEST_MAGNITUDE);
    switch (var.GetTypeSorting())
      {
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_REAL : which = PEP_LARGEST_REAL; break;
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = PEP_LARGEST_IMAGINARY; break;
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = PEP_LARGEST_MAGNITUDE; break;
      case PolynomialEigenProblem_Base<PetscScalar>::SORTED_USER : which = PEP_WHICH_USER; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_REAL : which = PEP_SMALLEST_REAL; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = PEP_SMALLEST_IMAGINARY; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = PEP_SMALLEST_MAGNITUDE; break;
          }
      }

    if ((!var.UseSpectralTransformation()) &&
        (var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES))
      {
        switch (var.GetTypeSorting())
          {
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_REAL : which = PEP_TARGET_REAL; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = PEP_TARGET_IMAGINARY; break;
          case PolynomialEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = PEP_TARGET_MAGNITUDE; break;
          }
	
	PEPSetTarget(solver, shift);
      }
    
    PEPSetWhichEigenpairs(solver, which);
    if (which == PEP_WHICH_USER)
      PEPSetEigenvalueComparison(solver, &PolynomialEigenProblem_Base<PetscScalar>::GetComparisonEigenvalueSlepc, &var);
    
    // tolerance and number of iterations
    double tol = var.GetStoppingCriterion();
    int nb_max_iter = var.GetNbMaximumIterations();
    PEPSetTolerances(solver, tol, nb_max_iter);    
    PEPSetDimensions(solver, nev, PETSC_DECIDE, PETSC_DECIDE);
    
    if (!var.UseSpectralTransformation())
      var.FactorizeMass();
    else
      {
        // binomial coefficients are computed
        Matrix<int> binom_coef(deg_pol+1, deg_pol+1);
        binom_coef.Zero();
        binom_coef(0, 0) = 1;
        binom_coef(1, 0) = 1; binom_coef(1, 1) = 1;
        for (int n = 2; n <= deg_pol; n++)
          {
            binom_coef(n, 0) = 1; binom_coef(n, n) = 1;
            for (int k = 1; k < n; k++)
              binom_coef(n, k) = binom_coef(n-1, k-1) + binom_coef(n-1, k);
          }
        
        // setting operators to compute
        Vector<PetscScalar> coef(deg_pol+1);
        for (int k = 0; k <= deg_pol; k++)
          {
            coef.Zero(); PetscScalar pow_shift = 1.0;
            for (int j = 0; j <= deg_pol-k; j++)
              {
                coef(j+k) = double(binom_coef(j+k, k))*pow_shift;
                pow_shift *= shift;
              }
            
            var.ComputeOperator(deg_pol-k, coef);
            if (k == 0)
              var.FactorizeOperator(coef);
          }
      }
    
    PEPSetScale(solver, PEP_SCALE_SCALAR, PETSC_DECIDE, PETSC_NULL, PETSC_NULL,
                PETSC_DECIDE, PETSC_DECIDE);
    
    // other parameters of PEP
    SetParametersSlepc(var.GetSlepcParameters(), solver);
    
    // the eigenvalue problem is solved
    int ierr = PEPSolve(solver);
    if (ierr != 0)
      {
	cout << "Error during solution of eigensystem =  " << ierr << endl;
	abort();
      }

    PEPConvergedReason reason;
    PEPGetConvergedReason(solver, &reason);
    if (reason < 0)
      {
        cout << "Failed to converged " << reason << endl;
        abort();
      }

    int print_level = var.GetPrintLevel();
    if (print_level >= 4)
      {
	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
	PEPErrorView(solver, PEP_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
	PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
      }

    int rank_proc = 0;
#ifdef SELDON_WITH_MPI
    MPI_Comm_rank(var.GetCommunicator(), &rank_proc);
#endif
    
    if ((print_level >= 1) && (rank_proc == 0))
      cout << "Number of linear solves = " << var.GetNbLinearSolves() << endl;
    
    // eigenvalues and eigenvectors are extracted
    Vec Vr, Vi, Vr_next, Vi_next;
    Vector<PetscScalar> Vr_vec, Vi_vec, Vr_vec_next, Vi_vec_next;
#ifdef SELDON_WITH_MPI
    AllocatePetscVector(var.GetCommunicator(), Vr, n, nglob, Vr_vec);
    AllocatePetscVector(var.GetCommunicator(), Vi, n, nglob, Vi_vec);
    AllocatePetscVector(var.GetCommunicator(), Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(var.GetCommunicator(), Vi_next, n, nglob, Vi_vec_next);
#else
    AllocatePetscVector(PETSC_COMM_SELF, Vr, n, nglob, Vr_vec);
    AllocatePetscVector(PETSC_COMM_SELF, Vi, n, nglob, Vi_vec);
    AllocatePetscVector(PETSC_COMM_SELF, Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(PETSC_COMM_SELF, Vi_next, n, nglob, Vi_vec_next);
#endif
    
    eigen_values.Reallocate(nev);
    lambda_imag.Reallocate(nev);    
    eigen_vectors.Reallocate(n, nev);
    int num = 0;
    PetscScalar Lr, Li, Lr_next, Li_next;
    bool eigen_pair = true;
    while (num < nev)
      {
	if (eigen_pair)
	  PEPGetEigenpair(solver, num, &Lr, &Li, Vr, Vi);

	if (num < nev-1)
	  PEPGetEigenpair(solver, num+1, &Lr_next, &Li_next, Vr_next, Vi_next);
	
	eigen_pair = PutEigenpairLapackForm(num, nev, Lr, Li, Vr_vec, Vi_vec,
					    Lr_next, Li_next, Vr_vec_next, Vi_vec_next,
					    eigen_values, lambda_imag, eigen_vectors);
	
	if (eigen_pair)
	  num += 2;
	else
	  {
	    Lr = Lr_next; Li = Li_next;
	    Vr_vec = Vr_vec_next; Vi_vec = Vi_vec_next;
	    num++;
	  }
      }
    
    Vr_vec.Nullify();
    Vi_vec.Nullify();
    Vr_vec_next.Nullify();
    Vi_vec_next.Nullify();
    
    if (var.UseSpectralTransformation())
      ApplySpectralTransform(shift, eigen_values, lambda_imag);

    // temporary objects are destroyed
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    PEPDestroy(&solver);
    for (int k = 0; k <= deg_pol; k++)
      MatDestroy(&op(k));

    var.DistributeEigenvectors(eigen_vectors);
  }
  
  
  /*****************************
   * Interface with NEP solver *
   *****************************/

  //! class for matshell operator for NEP solver
  struct NepSlepcOperator
  {
    // if true, the jacobian T'(L) is required
    bool jacobian;
    PetscScalar L; // value of L
    int num_operator; // for split formulation (operator number)
    NonLinearEigenProblem_Base<PetscScalar>* var; // pointer to the non-linear problem

    NepSlepcOperator()
    { jacobian = false; num_operator = -1; var = NULL; SetComplexZero(L); }
    
  };

  //! class for preconditioning for NEP solver 
  struct NepSlepcPreconditioning
  {
    PetscScalar L; // value of L
    NonLinearEigenProblem_Base<PetscScalar>* var; // pointer to the non-linear problem
    int type_solver; // eigensolver to be used
    Vector<PetscScalar> Ltab, coef_tab; // list of values for L and related coefficients
    Vector<int> numL_tab; // operator numbers
    
    NepSlepcPreconditioning()
    { SetComplexZero(L); type_solver = 0; }
  };

  //! preconditioning is updated
  void UpdatePrecond_Nep(PC pc, bool first_construct)
  {
    NepSlepcPreconditioning* shell;
    
    PCShellGetContext(pc, (void**)&shell);
    
    Mat Amat, Pmat;
    PCGetOperators(pc, &Amat, &Pmat);

    MatType type;
    MatGetType(Amat, &type);
    
    if (strcmp(type,MATSHELL) == 0)
      {
        // shell-matrix 
        void* ctxF;
        MatShellGetContext(Amat, &ctxF);

        PetscScalar L;
        if (shell->type_solver == SlepcParamNep::NLEIGS)
          {
            NEP_NLEIGS_MATSHELL* ctx
              = reinterpret_cast<NEP_NLEIGS_MATSHELL*>(ctxF);

            // nleigs solver => we retrieve the values of L and related coefficients
            Vector<PetscScalar> L(ctx->nmat), coef(ctx->nmat);
            Vector<int> numL(ctx->nmat);
            for (int i = 0; i < ctx->nmat; i++)
              {
                MatShellGetContext(ctx->A[i], &ctxF);

                NepSlepcOperator& op
                  = *reinterpret_cast<NepSlepcOperator*>(ctxF);

                L(i) = op.L; coef(i) = ctx->coeff[i];
                numL(i) = op.num_operator;
                shell->var->CheckValueL(L(i));
              }

            // new_coef will be true if the preconditioning changed
            bool new_coef = false;
            if (ctx->nmat != shell->Ltab.GetM())
              new_coef = true;
            else
              {                
                for (int i = 0; i < ctx->nmat; i++)
                  if ((abs(shell->Ltab(i) - L(i)) > 1e-12) || (abs(shell->coef_tab(i) - coef(i)) > 1e-12) || (shell->numL_tab(i) != numL(i)))
                    new_coef = true;
              }

            // preconditioning is recomputed if needed
            if (new_coef || first_construct)
              {
                //cout << " Value of L for preconditioning = " << L << endl;
                shell->Ltab = L; shell->coef_tab = coef; shell->numL_tab = numL;
                if (shell->var->UseSplitMatrices())
                  shell->var->ComputeSplitPreconditioning(numL, coef);
                else
                  shell->var->ComputePreconditioning(L, coef);
              }
          }
        else
          {
            // other solvers => a single value of L
            NepSlepcOperator& op
              = *reinterpret_cast<NepSlepcOperator*>(ctxF);
            
            L = op.L;
            shell->var->CheckValueL(L);
            
            //cout << " Value of L for preconditioning = " << L << endl;
            // preconditioning is recomputed only if coefficients changed
            if ((abs(L - shell->L) >= 1e-12) || first_construct)
              {
                shell->L = L;
                shell->var->ComputePreconditioning(L);
              }
          }
      }
    else
      {
        // case where the matrix is stored
        PetscInt m, n;
        MatGetSize(Amat, &m, &n);
        
        Mat_SeqAIJ* A = (Mat_SeqAIJ*) Amat->data;

        // the matrix is duplicated
        DistributedMatrix<PetscScalar, General, ArrayRowSparse> Ad;
        Ad.Reallocate(m, n);

        for (int i = 0; i < m; i++)
          {
            int nz = A->i[i+1] - A->i[i];
            Ad.ReallocateRow(i, nz);
            for (int j = 0; j < nz; j++)
              {
                Ad.Index(i, j) = A->j[A->i[i]+j];
                Ad.Value(i, j) = A->a[A->i[i]+j];
              }
          }

        // and preconditioning computed
        shell->var->ComputeExplicitPreconditioning(Ad);
      }
  }
  
  // first construction of preconditioning
  PetscErrorCode ConstructPreconditioning_NepOperator(PC pc)
  {
    UpdatePrecond_Nep(pc, true);
    
    return 0;
  }

  // preconditioning may be updated (if L is different)
  PetscErrorCode UpdatePreconditioning_NepOperator(PC pc, KSP ksp, Vec b, Vec x)
  {
    UpdatePrecond_Nep(pc, false);
    return 0;
  }

  // applies the preconditioning
  PetscErrorCode ApplyPreconditioning_NepOperator(PC pc, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);
    
    NepSlepcPreconditioning* shell;

    PCShellGetContext(pc, (void**)&shell);
    shell->var->IncrementLinearSolves();
    shell->var->ApplyPreconditioning(SeldonNoTrans, xvec, yvec);
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  // provides the singular points of T(L)
  PetscErrorCode SetSingularities_NepOperator(NEP nep, PetscInt *maxnp,
                                              PetscScalar *xi, void *ctx)
  {
    NonLinearEigenProblem_Base<PetscScalar>& var
      = *reinterpret_cast<NonLinearEigenProblem_Base<PetscScalar>* >(ctx);

    const Vector<PetscScalar>& s = var.GetSingularities();
    *maxnp = min(*maxnp, PetscInt(s.GetM()));
    for (int k = 0; k < *maxnp; k++)
      xi[k] = s(k);

    return 0;
  }
  
  // sets other parameters of slepc for nep
  void SetParametersSlepc(const SlepcParamNep& param, NEP& nep,
                          NonLinearEigenProblem_Base<PetscScalar>& var,
                          NepSlepcPreconditioning& prec)
  {
    switch(param.GetEigensolverType())
      {
      case SlepcParamNep::RII :  NEPSetType(nep, NEPRII); break;
      case SlepcParamNep::SLP :  NEPSetType(nep, NEPSLP); break;
      case SlepcParamNep::NARNOLDI :  NEPSetType(nep, NEPNARNOLDI); break;
      case SlepcParamNep::CISS :  NEPSetType(nep, NEPCISS); break;
      case SlepcParamNep::INTERPOL :  NEPSetType(nep, NEPINTERPOL); break;
      case SlepcParamNep::NLEIGS :  NEPSetType(nep, NEPNLEIGS); break;
      }        
    
    KSP            ksp; PC pc; bool precond = false;
    if (param.GetEigensolverType() == SlepcParamNep::RII)
      {
        precond = true;
        NEPRIISetLagPreconditioner(nep, 0);
        NEPRIIGetKSP(nep, &ksp);
        if (var.IsHermitianProblem())
          NEPRIISetHermitian(nep, PETSC_TRUE);
      }
    else if (param.GetEigensolverType() == SlepcParamNep::SLP)
      {
        precond = true;
        NEPSLPGetKSP(nep, &ksp);
      }
    else if (param.GetEigensolverType() == SlepcParamNep::NARNOLDI)
      {
        precond = true;
        NEPNArnoldiSetLagPreconditioner(nep, 0);
        NEPNArnoldiGetKSP(nep, &ksp);
        NEPNArnoldiSetLagPreconditioner(nep, 0);
      }
    else if (param.GetEigensolverType() == SlepcParamNep::NLEIGS)
      {
        precond = true;
        if (!var.UseSplitMatrices())
          NEPNLEIGSSetSingularitiesFunction(nep, SetSingularities_NepOperator, reinterpret_cast<void*>(&var));

        if (param.GetRKShifts().GetM() > 0)
          {
            int ns = param.GetRKShifts().GetM();
            PetscScalar shifts[ns];
            for (int k = 0; k < ns; k++)
              shifts[k] = param.GetRKShifts()(k);
            
            NEPNLEIGSSetRKShifts(nep, ns, shifts);
          }
        
        int nsolve; KSP* vec_ksp;
        NEPNLEIGSGetKSPs(nep, &nsolve, &vec_ksp);
        ksp = vec_ksp[0];
        if (nsolve > 1)
          {
            cout << "Too many solves required in nleigs = " << nsolve << endl;
            cout << "Not implemented" << endl;
            abort();
          }
        
        RG rg;
        NEPGetRG(nep, &rg);
        RGSetType(rg, RGINTERVAL);
        RGIntervalSetEndpoints(rg, param.GetLrMin(), param.GetLrMax(),
                               param.GetLiMin(), param.GetLiMax());

        if (!param.InsideRegion(var.GetShiftValue()))
          {
            cout << "Target " << var.GetShiftValue() << " not inside the region"
                 << " [ " << param.GetLrMin() << ", " << param.GetLrMax() << "] x ["
                 << param.GetLiMin() << ", " << param.GetLiMax() << "]" << endl;
            
            abort();
          }
          
        if (param.FullBasis())
          NEPNLEIGSSetFullBasis(nep, PETSC_TRUE);
        else
          NEPNLEIGSSetFullBasis(nep, PETSC_FALSE);
        
        if (param.LockingVariant())
          NEPNLEIGSSetLocking(nep, PETSC_TRUE);
        else
          NEPNLEIGSSetLocking(nep, PETSC_FALSE);

        NEPNLEIGSSetRestart(nep, param.GetRestartNleigs());
        
        if ((param.GetInterpolationDegree() >= 0) && (param.GetInterpolationTolerance() > 0))
          NEPNLEIGSSetInterpolation(nep, param.GetInterpolationTolerance(),
                                    param.GetInterpolationDegree());
      }
    else
      {
        cout << "Solver " << param.GetEigensolverType() << " not interfaced in nep" << endl;
        abort();
      }

    if (var.ExplicitMatrix() && param.UseDefaultPetscSolver())
      precond = false;
    
    if (precond)
      {
        prec.type_solver = param.GetEigensolverType();
        prec.var = &var;
        if (var.ExactPreconditioning())
          KSPSetType(ksp, KSPPREONLY); // (dans le cas ou le preconditionneur est exact)
        else
          KSPSetType(ksp, KSPBCGS);

        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCSHELL);
        PCShellSetContext(pc, reinterpret_cast<void*>(&prec));
        PCShellSetApply(pc, ApplyPreconditioning_NepOperator);
        //PCShellSetSetUp(pc, ConstructPreconditioning_NepOperator);
        PCShellSetPreSolve(pc, UpdatePreconditioning_NepOperator);
      }

    // options from the command line
    if (param.UseCommandLineOptions())
      NEPSetFromOptions(nep);
  }

  // computation of operator T(lambda)
  PetscErrorCode FormFunctionNEP(NEP nep, PetscScalar lambda, Mat A, Mat B, void* ctx)
  {
    void* ctxF;
    MatShellGetContext(A, &ctxF);
    
    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctxF);
    
    op.L = lambda;
    op.var->CheckValueL(lambda);
    op.var->ComputeOperator(op.L);
    
    return 0;
  }

  // computation of operator T'(lambda)
  PetscErrorCode FormJacobianNEP(NEP nep, PetscScalar lambda, Mat A, void* ctx)
  {
    void* ctxF;
    MatShellGetContext(A, &ctxF);
    
    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctxF);
    
    op.L = lambda;
    op.var->CheckValueL(lambda);
    op.var->ComputeJacobian(op.L);
    
    return 0;
  }

  // conversion from a Seldon matrix to a PetscMatrix
  void ConvertToPetsc(DistributedMatrix<PetscScalar, General, ArrayRowSparse>& A, Mat Ap)
  {
    Vector<int> nnz(A.GetM());
    for (int i = 0; i < A.GetM(); i++)
      nnz(i) = A.GetRowSize(i);

    MatSeqAIJSetPreallocation(Ap, PETSC_DEFAULT, nnz.GetData());
    MatSetUp(Ap);

    for (int i = 0; i < A.GetM(); i++)
      MatSetValues(Ap, 1, &i, A.GetRowSize(i), A.GetIndex(i), A.GetData(i), INSERT_VALUES);

    MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY);
  }

  // computes T(L) explicitely
  PetscErrorCode FormFunctionExpNEP(NEP nep, PetscScalar lambda, Mat Ap, Mat Bp, void* ctx)
  {
    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);

    DistributedMatrix<PetscScalar, General, ArrayRowSparse> A;
    op.var->CheckValueL(lambda);
    op.var->ComputeOperatorExplicit(lambda, A);

    ConvertToPetsc(A, Ap);
    
    return 0;
  }

  // computes T'(L) explicitely
  PetscErrorCode FormJacobianExpNEP(NEP nep, PetscScalar lambda, Mat Ap, void* ctx)
  {
    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);
    
    DistributedMatrix<PetscScalar, General, ArrayRowSparse> A;
    op.var->CheckValueL(lambda);
    op.var->ComputeJacobianExplicit(lambda, A);

    ConvertToPetsc(A, Ap);
    
    return 0;
  }

  //! matrix-vector product y = op x
  PetscErrorCode MatMult_NepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);

    if (op.num_operator == -1)
      {
        op.var->CheckValueL(op.L);
        op.var->IncrementProdMatVect();
        if (op.jacobian)
          op.var->MltJacobian(op.L, SeldonNoTrans, xvec, yvec);
        else
          op.var->MltOperator(op.L, SeldonNoTrans, xvec, yvec);
      }
    else
      {
        op.var->IncrementProdMatVect();
        op.var->MltOperatorSplit(op.num_operator, SeldonNoTrans, xvec, yvec);
      }
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  //! matrix-vector product y = op^T x
  PetscErrorCode MatMultTrans_NepOperator(Mat mat, Vec x, Vec y)
  {
    Vector<PetscScalar> xvec, yvec;
    CopyPointerPetsc(x, xvec);
    CopyPointerPetsc(y, yvec);

    void* ctx;
    MatShellGetContext(mat, &ctx);

    NepSlepcOperator& op
      = *reinterpret_cast<NepSlepcOperator*>(ctx);
    
    if (op.num_operator == -1)
      {
        op.var->IncrementProdMatVect();
        op.var->CheckValueL(op.L);
        if (op.jacobian)
          op.var->MltJacobian(op.L, SeldonTrans, xvec, yvec);
        else
          op.var->MltOperator(op.L, SeldonTrans, xvec, yvec);
      }
    else
      {
        op.var->MltOperatorSplit(op.num_operator, SeldonTrans, xvec, yvec);
        op.var->IncrementProdMatVect();
      }
    
    xvec.Nullify(); yvec.Nullify();
    return 0;
  }

  // operator is duplicated
  PetscErrorCode MatDuplicate_NepOperator(Mat A, MatDuplicateOption op, Mat* B)
  {
    NepSlepcOperator         *actx,*bctx;
    PetscInt n, nloc;
#ifdef SELDON_WITH_MPI
    MPI_Comm comm;
#endif

    MatShellGetContext(A, (void**)&actx);
    MatGetLocalSize(A, &nloc, NULL);
    MatGetSize(A, &n, NULL);
    
    bctx = new NepSlepcOperator;
    bctx->jacobian      = actx->jacobian;
    bctx->L  = actx->L;
    bctx->var = actx->var;
    bctx->num_operator = actx->num_operator;
    
#ifdef SELDON_WITH_MPI
    PetscObjectGetComm((PetscObject)A, &comm);
    MatCreateShell(comm, nloc, nloc, n, n,(void*)bctx, B);
#else
    MatCreateShell(PETSC_COMM_SELF, nloc, nloc, n, n,(void*)bctx, B);
#endif        
    MatShellSetOperation(*B, MATOP_MULT, (void(*)(void))MatMult_NepOperator);
    MatShellSetOperation(*B, MATOP_MULT_TRANSPOSE, (void(*)(void))MatMultTrans_NepOperator);
    MatShellSetOperation(*B, MATOP_DUPLICATE, (void(*)(void))MatDuplicate_NepOperator);

    return 0;
  }
  
  //! functions to compute eigenvalues with NEP solver of Slepc
  void FindEigenvaluesSlepc_(NonLinearEigenProblem_Base<PetscScalar>& var,
                             Vector<PetscScalar>& eigen_values,
                             Vector<PetscScalar>& lambda_imag,
                             Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    // initializing of computation
    PetscScalar shift = var.GetShiftValue();
    PetscScalar zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    // dimensions
    int nev = var.GetNbAskedEigenvalues();
    int n = var.GetM();
    int nglob = var.GetGlobalM();
    
    // creation of the eigensolver
    NEP solver;
#ifdef SELDON_WITH_MPI
    NEPCreate(var.GetCommunicator(), &solver);
#else
    NEPCreate(PETSC_COMM_SELF, &solver);
#endif
    
    // creation of shell matrices
    Mat EvalF, EvalJacob;
    NepSlepcOperator operT;
    operT.jacobian = false; operT.var = &var;
    NepSlepcOperator operTp;
    operTp.jacobian = true; operTp.var = &var;

    
    // matrices for the split formulation
    Vector<Mat> EvalF_Split(var.GetNbSplitMatrices());
    Vector<NepSlepcOperator> operSplit(var.GetNbSplitMatrices());
    Vector<FN> FunctionF_Split(var.GetNbSplitMatrices());
    
    if (var.UseSplitMatrices())
      {
        if (var.ExplicitMatrix())
          {
            for (int i = 0; i < var.GetNbSplitMatrices(); i++)
              {
#ifdef SELDON_WITH_MPI
                MatCreate(var.GetCommunicator(), &EvalF_Split(i));
#else
                MatCreate(PETSC_COMM_SELF, &EvalF_Split(i));
#endif
                MatSetSizes(EvalF_Split(i), n, n, nglob, nglob);
                MatSetType(EvalF_Split(i), MATSEQAIJ);

                DistributedMatrix<PetscScalar, General, ArrayRowSparse> A;
                var.ComputeOperatorSplitExplicit(i, A);

                ConvertToPetsc(A, EvalF_Split(i));
              }
          }
        else
          {
            for (int i = 0; i < var.GetNbSplitMatrices(); i++)
              {
                operSplit(i).jacobian = false; operSplit(i).var = &var;
                operSplit(i).num_operator = i;
#ifdef SELDON_WITH_MPI
                MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                               reinterpret_cast<void*>(&operSplit(i)), &EvalF_Split(i));
#else
                MatCreateShell(PETSC_COMM_SELF, n, n, nglob, nglob,
                               reinterpret_cast<void*>(&operSplit(i)), &EvalF_Split(i));
#endif
                
                MatShellSetOperation(EvalF_Split(i), MATOP_MULT, (void(*)(void))MatMult_NepOperator);
                MatShellSetOperation(EvalF_Split(i), MATOP_MULT_TRANSPOSE, (void(*)(void))MatMultTrans_NepOperator);
                MatShellSetOperation(EvalF_Split(i), MATOP_DUPLICATE, (void(*)(void))MatDuplicate_NepOperator);
              }
          }

        if (!var.RationalSplitFunctions())
          {
            FunctionF_Split = var.GetSlepcFunction_NonRational();
            
            NEPSetProblemType(solver, NEP_GENERAL);            
          }
        else
          {
            for (int i = 0; i < var.GetNbSplitMatrices(); i++)
              {
                const Vector<PetscScalar>& P = var.GetNumeratorSplitFct(i);
                const Vector<PetscScalar>& Q = var.GetDenominatorSplitFct(i);
                
#ifdef SELDON_WITH_MPI
                FNCreate(var.GetCommunicator(), &FunctionF_Split(i));
#else
                FNCreate(PETSC_COMM_SELF, &FunctionF_Split(i));
#endif
                FNSetType(FunctionF_Split(i), FNRATIONAL);
                FNRationalSetNumerator(FunctionF_Split(i), P.GetM(), P.GetData());
                FNRationalSetDenominator(FunctionF_Split(i), Q.GetM(), Q.GetData());
              }
            
            NEPSetProblemType(solver, NEP_RATIONAL);
          }

        NEPSetSplitOperator(solver, var.GetNbSplitMatrices(), EvalF_Split.GetData(),
                                FunctionF_Split.GetData(), DIFFERENT_NONZERO_PATTERN);
      }
    else
      {
        if (var.ExplicitMatrix())
          {
#ifdef SELDON_WITH_MPI
            MatCreate(var.GetCommunicator(), &EvalF);
#else
            MatCreate(PETSC_COMM_SELF, &EvalF);
#endif
            MatSetSizes(EvalF, n, n, nglob, nglob);
            MatSetType(EvalF, MATSEQAIJ);
            MatSetUp(EvalF);
          }
        else
          {
#ifdef SELDON_WITH_MPI
            MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                           reinterpret_cast<void*>(&operT), &EvalF);
#else
            MatCreateShell(PETSC_COMM_SELF, n, n, nglob, nglob,
                           reinterpret_cast<void*>(&operT), &EvalF);
#endif
            
            MatShellSetOperation(EvalF, MATOP_MULT, (void(*)(void))MatMult_NepOperator);
            MatShellSetOperation(EvalF, MATOP_MULT_TRANSPOSE, (void(*)(void))MatMultTrans_NepOperator);
            MatShellSetOperation(EvalF, MATOP_DUPLICATE, (void(*)(void))MatDuplicate_NepOperator);
          }
        
        if (var.ExplicitMatrix())
          {
#ifdef SELDON_WITH_MPI
            MatCreate(var.GetCommunicator(), &EvalJacob);
#else
            MatCreate(PETSC_COMM_SELF, &EvalJacob);
#endif
            MatSetSizes(EvalJacob, n, n, nglob, nglob);
            MatSetType(EvalJacob, MATSEQAIJ);
            MatSetUp(EvalJacob);
          }
        else
          {
#ifdef SELDON_WITH_MPI
            MatCreateShell(var.GetCommunicator(), n, n, nglob, nglob,
                           reinterpret_cast<void*>(&operTp), &EvalJacob);
#else
            MatCreateShell(PETSC_COMM_SELF, n, n, nglob, nglob,
                           reinterpret_cast<void*>(&operTp), &EvalJacob);
#endif
            
            MatShellSetOperation(EvalJacob, MATOP_MULT, (void(*)(void))MatMult_NepOperator);
            MatShellSetOperation(EvalJacob, MATOP_MULT_TRANSPOSE, (void(*)(void))MatMultTrans_NepOperator);
            MatShellSetOperation(EvalJacob, MATOP_DUPLICATE, (void(*)(void))MatDuplicate_NepOperator);
          }
        
        if (var.ExplicitMatrix())
          NEPSetFunction(solver, EvalF, EvalF, FormFunctionExpNEP, reinterpret_cast<void*>(&operT));
        else
          NEPSetFunction(solver, EvalF, EvalF, FormFunctionNEP, NULL);
        
        if (var.ExplicitMatrix())
          NEPSetJacobian(solver, EvalJacob, FormJacobianExpNEP, reinterpret_cast<void*>(&operTp));
        else
          NEPSetJacobian(solver, EvalJacob, FormJacobianNEP, NULL);
      }
    
    // sorting and selection of spectrum
    NEPWhich which(NEP_LARGEST_MAGNITUDE);
    switch (var.GetTypeSorting())
      {
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_REAL : which = NEP_LARGEST_REAL; break;
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = NEP_LARGEST_IMAGINARY; break;
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = NEP_LARGEST_MAGNITUDE; break;
      case NonLinearEigenProblem_Base<PetscScalar>::SORTED_USER : which = NEP_WHICH_USER; break;
      }
    
    if (var.GetTypeSpectrum() == var.SMALL_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_REAL : which = NEP_SMALLEST_REAL; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = NEP_SMALLEST_IMAGINARY; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = NEP_SMALLEST_MAGNITUDE; break;
          }
      }

    if (var.GetTypeSpectrum() == var.CENTERED_EIGENVALUES)
      {
        switch (var.GetTypeSorting())
          {
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_REAL : which = NEP_TARGET_REAL; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_IMAG : which = NEP_TARGET_IMAGINARY; break;
          case NonLinearEigenProblem_Base<PetscScalar>::SORTED_MODULUS : which = NEP_TARGET_MAGNITUDE; break;
          }

	NEPSetTarget(solver, shift);
      }

    NEPSetWhichEigenpairs(solver, which);
    if (which == NEP_WHICH_USER)
      NEPSetEigenvalueComparison(solver, &NonLinearEigenProblem_Base<PetscScalar>::GetComparisonEigenvalueSlepc, &var);
    
    // tolerance and number of iterations
    double tol = var.GetStoppingCriterion();
    int nb_max_iter = var.GetNbMaximumIterations();
    NEPSetTolerances(solver, tol, nb_max_iter);
    int ncv = PETSC_DEFAULT;
    if (var.GetNbArnoldiVectors() > 0)
      ncv = var.GetNbArnoldiVectors();

    if (nev > 0)
      NEPSetDimensions(solver, nev, ncv, PETSC_DEFAULT);

    // other parameters of NEP
    NepSlepcPreconditioning prec;
    SetParametersSlepc(var.GetSlepcParameters(), solver, var, prec);
    
    // the eigenvalue problem is solved
    PetscErrorCode ierr = NEPSolve(solver);
    if (ierr != 0)
      {
	cout << "Error during solution of eigensystem =  " << ierr << endl;
	abort();
      }

    int print_level = var.GetPrintLevel();
    if (print_level >= 4)
      {
	PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL);
        NEPConvergedReasonView(solver, PETSC_VIEWER_STDOUT_WORLD);
	NEPErrorView(solver, NEP_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD);
	PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
      }

    // nev is modified if less eigenvalues are converged
    NEPGetConverged(solver, &nev);
    
    // eigenvalues and eigenvectors are extracted
    Vec Vr, Vi, Vr_next, Vi_next;
    Vector<PetscScalar> Vr_vec, Vi_vec, Vr_vec_next, Vi_vec_next;
#ifdef SELDON_WITH_MPI
    AllocatePetscVector(var.GetCommunicator(), Vr, n, nglob, Vr_vec);
    AllocatePetscVector(var.GetCommunicator(), Vi, n, nglob, Vi_vec);
    AllocatePetscVector(var.GetCommunicator(), Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(var.GetCommunicator(), Vi_next, n, nglob, Vi_vec_next);
#else
    AllocatePetscVector(PETSC_COMM_SELF, Vr, n, nglob, Vr_vec);
    AllocatePetscVector(PETSC_COMM_SELF, Vi, n, nglob, Vi_vec);
    AllocatePetscVector(PETSC_COMM_SELF, Vr_next, n, nglob, Vr_vec_next);
    AllocatePetscVector(PETSC_COMM_SELF, Vi_next, n, nglob, Vi_vec_next);
#endif

    eigen_values.Reallocate(nev);
    lambda_imag.Reallocate(nev);    
    eigen_vectors.Reallocate(n, nev);
    int num = 0;
    PetscScalar Lr, Li, Lr_next, Li_next;
    bool eigen_pair = true;
    while (num < nev)
      {
	if (eigen_pair)
	  NEPGetEigenpair(solver, num, &Lr, &Li, Vr, Vi);

	if (num < nev-1)
	  NEPGetEigenpair(solver, num+1, &Lr_next, &Li_next, Vr_next, Vi_next);
	
	eigen_pair = PutEigenpairLapackForm(num, nev, Lr, Li, Vr_vec, Vi_vec,
					    Lr_next, Li_next, Vr_vec_next, Vi_vec_next,
					    eigen_values, lambda_imag, eigen_vectors);
	
	if (eigen_pair)
	  num += 2;
	else
	  {
	    Lr = Lr_next; Li = Li_next;
	    Vr_vec = Vr_vec_next; Vi_vec = Vi_vec_next;
	    num++;
	  }
      }
    
    Vr_vec.Nullify();
    Vi_vec.Nullify();
    Vr_vec_next.Nullify();
    Vi_vec_next.Nullify();
    
    // temporary objects are destroyed
    VecDestroy(&Vr);
    VecDestroy(&Vi);
    NEPDestroy(&solver);
    if (var.UseSplitMatrices())
      {
        for (int i = 0; i < var.GetNbSplitMatrices(); i++)
          {
            MatDestroy(&EvalF_Split(i));
            if (var.RationalSplitFunctions())
              FNDestroy(&FunctionF_Split(i));
          }
      }
    else
      {
        MatDestroy(&EvalF);
        MatDestroy(&EvalJacob);
      }

    var.DistributeEigenvectors(eigen_vectors);
  }
  
#ifdef SELDON_WITH_VIRTUAL
  void FindEigenvaluesSlepc(EigenProblem_Base<PetscScalar>& var,
			    Vector<PetscScalar>& eigen_values,
			    Vector<PetscScalar>& lambda_imag,
			    Matrix<PetscScalar, General, ColMajor>& eigen_vectors)
  {
    FindEigenvaluesSlepc_(var, eigen_values, lambda_imag, eigen_vectors);
  }  

  void FindEigenvaluesSlepc(PolynomialEigenProblem_Base<Petsc_Scalar>& var,
                            Vector<Petsc_Scalar>& eigen_values,
                            Vector<Petsc_Scalar>& lambda_imag,
                            Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors)
  {
    FindEigenvaluesSlepc_(var, eigen_values, lambda_imag, eigen_vectors);
  }

  void FindEigenvaluesSlepc(NonLinearEigenProblem_Base<Petsc_Scalar>& var,
                            Vector<Petsc_Scalar>& eigen_values,
                            Vector<Petsc_Scalar>& lambda_imag,
                            Matrix<Petsc_Scalar, General, ColMajor>& eigen_vectors)
  {
    FindEigenvaluesSlepc_(var, eigen_values, lambda_imag, eigen_vectors);
  }
#else
  template<class EigenProblem, class T, class Allocator1,
           class Allocator2, class Allocator3>
  void FindEigenvaluesSlepc(EigenProblem& var,
			    Vector<T, VectFull, Allocator1>& eigen_values,
			    Vector<T, VectFull, Allocator2>& lambda_imag,
			    Matrix<T, General, ColMajor, Allocator3>& eigen_vectors)
  {
    cout << "Recompile with SELDON_WITH_VIRTUAL" << endl;
    abort();
  }
#endif
  
}

#define SELDON_FILE_SLEPC_HXX
#endif
