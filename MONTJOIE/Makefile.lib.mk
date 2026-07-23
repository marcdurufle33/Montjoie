# Bessel functions
lib/libbessel.a : lib/Bessel.o lib/d1mach.o lib/i1mach.o
	ar rv lib/libbessel.a lib/Bessel.o lib/d1mach.o lib/i1mach.o

lib/libbessel.so: lib/Bessel.o lib/i1mach.o lib/d1mach.o
	$(CC) -shared -Wl,-soname,libbessel.so -rdynamic -o lib/libbessel.so lib/Bessel.o lib/d1mach.o lib/i1mach.o

# Seldon part

# main part of Seldon
OBJ_SELDON := lib/Compil/Seldon/Array3D$(SELD_SUFFIX).o lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o lib/Compil/Seldon/Common$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixPointers$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixPacked$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixHermSymTriang$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixSparse$(SELD_SUFFIX).o lib/Compil/Seldon/Vector$(SELD_SUFFIX).o lib/Compil/Seldon/Vector2$(SELD_SUFFIX).o lib/Compil/Seldon/FunctionsMatrixVector$(SELD_SUFFIX).o lib/Compil/Seldon/FunctionsMatrix$(SELD_SUFFIX).o lib/Compil/Seldon/FunctionsMatrixDense$(SELD_SUFFIX).o lib/Compil/Seldon/FunctionsMatrixArray$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixConversion$(SELD_SUFFIX).o lib/Compil/Seldon/PermutationScalingMatrix$(SELD_SUFFIX).o lib/Compil/Seldon/RelaxationMatrixVector$(SELD_SUFFIX).o lib/Compil/Seldon/IOMatrixMarket$(SELD_SUFFIX).o lib/Compil/Seldon/IlutPreconditioning$(SELD_SUFFIX).o lib/Compil/Seldon/IterativeSolver$(SELD_SUFFIX).o lib/Compil/Seldon/check_dim$(SELD_SUFFIX).o

# complex matrices
OBJ_SELDON := $(OBJ_SELDON) lib/Compil/Seldon/FunctionsMatVectComplex$(SELD_SUFFIX).o lib/Compil/Seldon/FunctionsMatrixComplex$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixComplexConversions$(SELD_SUFFIX).o lib/Compil/Seldon/MatrixComplexSparse$(SELD_SUFFIX).o

# Blas functions
ifeq ($(USE_BLAS),YES)
  OBJ_BLAS := lib/Compil/Seldon/Blas$(SELD_SUFFIX).o lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o
else
  OBJ_BLAS := lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o
endif
OBJ_SELDON := $(OBJ_SELDON) $(OBJ_BLAS)
OBJ_BLAS := lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o lib/Compil/Seldon/Blas$(SELD_SUFFIX).o lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o

# MPI functions
ifeq ($(USE_MPI),YES)
  OBJ_SELDON := $(OBJ_SELDON) lib/Compil/Seldon/DistributedVector$(SELD_SUFFIX).o lib/Compil/Seldon/DistributedMatrix$(SELD_SUFFIX).o
endif

# interface with eigenvalue solvers
# Arpack is always needed in that case
OBJ_EIG := lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o
ifeq ($(USE_ARPACK),YES)
  OBJ_EIG := $(OBJ_EIG) lib/Compil/Seldon/Arpack$(SELD_SUFFIX).o
  ifeq ($(USE_FEAST),YES)
    OBJ_EIG := $(OBJ_EIG) lib/Compil/Seldon/Feast$(SELD_SUFFIX).o
  endif
  ifeq ($(USE_ANASAZI),YES)
    OBJ_EIG := $(OBJ_EIG) lib/Compil/Seldon/Anasazi$(SELD_SUFFIX).o
  endif
  ifeq ($(USE_SLEPC),YES)
    OBJ_EIG := $(OBJ_EIG) lib/Compil/Seldon/Slepc$(SELD_SUFFIX).o
  endif
endif

OBJ_SELDON := $(OBJ_SELDON) $(OBJ_EIG)

# interface with hypre
ifeq ($(USE_HYPRE),YES)
  OBJ_SELDON := $(OBJ_SELDON) lib/Compil/Seldon/Hypre$(SELD_SUFFIX).o
endif

# interfaces with direct solvers
OBJ_SELDON_SOLVE := lib/Compil/Seldon/CholeskySolver$(SELD_SUFFIX).o lib/Compil/Seldon/SparseSeldonSolver$(SELD_SUFFIX).o lib/Compil/Seldon/SparseDirectSolver$(SELD_SUFFIX).o lib/Compil/Seldon/DistributedSolver$(SELD_SUFFIX).o

ifeq ($(USE_MUMPS),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o
endif

ifeq ($(USE_MUMPS),INT64)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o
endif

ifeq ($(USE_PARDISO),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Pardiso$(SELD_SUFFIX).o
endif

ifeq ($(USE_PASTIX),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o
endif

ifeq ($(USE_PASTIX),INT64)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o
endif

ifeq ($(USE_UMFPACK),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/UmfPack$(SELD_SUFFIX).o lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o
endif

ifeq ($(USE_SUPERLU),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/SuperLU$(SELD_SUFFIX).o
endif

ifeq ($(USE_CHOLMOD),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o
endif

ifeq ($(USE_WSMP),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) lib/Compil/Seldon/Wsmp$(SELD_SUFFIX).o
endif

OBJ_SELDON := $(OBJ_SELDON) $(OBJ_SELDON_SOLVE)

# list of dependances for each file .cpp contained in lib/Compil/Seldon
lib/Compil/Seldon/Array3D$(SELD_SUFFIX).o : lib/Compil/Seldon/Array3D.cpp $(CHESELDON)/array/Array3D.cxx $(CHESELDON)/array/Array4D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o : lib/Compil/Seldon/BandMatrix.cpp $(CHESELDON)/matrix_sparse/BandMatrix.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Common$(SELD_SUFFIX).o : lib/Compil/Seldon/Common.cpp $(CHESELDON)/share/Errors.cxx $(CHESELDON)/share/MatrixFlag.cxx $(CHESELDON)/share/Common.cxx $(CHESELDON)/share/Allocator.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/DistributedVector$(SELD_SUFFIX).o : lib/Compil/Seldon/DistributedVector.cpp $(CHESELDON)/share/MpiCommunication.cxx $(CHESELDON)/vector/DistributedVector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/DistributedMatrix$(SELD_SUFFIX).o : lib/Compil/Seldon/DistributedMatrix.cpp $(CHESELDON)/matrix_sparse/DistributedMatrix.cxx $(CHESELDON)/matrix_sparse/DistributedMatrixFunction.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixPointers$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixPointers.cpp $(CHESELDON)/matrix/Matrix_Pointers.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixPacked$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixPacked.cpp $(CHESELDON)/matrix/Matrix_HermPacked.cxx $(CHESELDON)/matrix/Matrix_SymPacked.cxx $(CHESELDON)/matrix/Matrix_TriangPacked.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixHermSymTriang$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixHermSymTriang.cpp $(CHESELDON)/matrix/Matrix_Hermitian.cxx $(CHESELDON)/matrix/Matrix_Symmetric.cxx $(CHESELDON)/matrix/Matrix_Triangular.cxx $(CHESELDON)/matrix/Matrix_Base.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixSparse$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixSparse.cpp $(CHESELDON)/matrix_sparse/Matrix_Sparse.cxx $(CHESELDON)/matrix_sparse/Matrix_SymSparse.cxx $(CHESELDON)/matrix_sparse/Matrix_ArraySparse.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Vector$(SELD_SUFFIX).o : lib/Compil/Seldon/Vector.cpp $(CHESELDON)/vector/Vector.cxx $(CHESELDON)/vector/SparseVector.cxx $(CHESELDON)/vector/Functions_Arrays.cxx $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Vector2$(SELD_SUFFIX).o : lib/Compil/Seldon/Vector2.cpp $(CHESELDON)/vector/Vector2.cxx $(CHESELDON)/vector/Vector3.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/FunctionsMatrixVector$(SELD_SUFFIX).o : lib/Compil/Seldon/FunctionsMatrixVector.cpp $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx $(CHESELDON)/computation/interfaces/Mkl_Sparse.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/FunctionsMatrix$(SELD_SUFFIX).o : lib/Compil/Seldon/FunctionsMatrix.cpp $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx $(CHESELDON)/computation/interfaces/Mkl_Sparse.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/FunctionsMatrixDense$(SELD_SUFFIX).o : lib/Compil/Seldon/FunctionsMatrixDense.cpp $(CHESELDON)/matrix/Functions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/FunctionsMatrixArray$(SELD_SUFFIX).o : lib/Compil/Seldon/FunctionsMatrixArray.cpp $(CHESELDON)/matrix_sparse/Functions_MatrixArray.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/FunctionsMatVectComplex$(SELD_SUFFIX).o : lib/Compil/Seldon/FunctionsMatVectComplex.cpp $(CHESELDON)/matrix_sparse/complex/Functions_MatVectComplex.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/FunctionsMatrixComplex$(SELD_SUFFIX).o : lib/Compil/Seldon/FunctionsMatrixComplex.cpp $(CHESELDON)/matrix_sparse/complex/Functions_MatrixComplex.cxx $(CHESELDON)/matrix_sparse/Functions_MatrixArray.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixComplexConversions$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixComplexConversions.cpp $(CHESELDON)/matrix_sparse/complex/Matrix_ComplexConversions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixComplexSparse$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixComplexSparse.cpp $(CHESELDON)/matrix_sparse/complex/Matrix_ArrayComplexSparse.cxx $(CHESELDON)/matrix_sparse/complex/Matrix_ComplexSparse.cxx $(CHESELDON)/matrix_sparse/complex/Matrix_SymComplexSparse.cxx $(CHESELDON)/matrix_sparse/Matrix_Sparse.cxx $(CHESELDON)/matrix_sparse/IOMatrixMarket.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/MatrixConversion$(SELD_SUFFIX).o : lib/Compil/Seldon/MatrixConversion.cpp $(CHESELDON)/matrix_sparse/Matrix_Conversions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/PermutationScalingMatrix$(SELD_SUFFIX).o : lib/Compil/Seldon/PermutationScalingMatrix.cpp $(CHESELDON)/matrix_sparse/Permutation_ScalingMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/RelaxationMatrixVector$(SELD_SUFFIX).o : lib/Compil/Seldon/RelaxationMatrixVector.cpp $(CHESELDON)/matrix_sparse/Relaxation_MatVect.cxx $(CHESELDON)/matrix_sparse/complex/Functions_MatVectComplex.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Blas$(SELD_SUFFIX).o : lib/Compil/Seldon/Blas.cpp $(CHESELDON)/computation/interfaces/Blas_1.cxx $(CHESELDON)/computation/interfaces/Blas_2.cxx $(CHESELDON)/computation/interfaces/Blas_3.cxx $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o : lib/Compil/Seldon/BlasMpfr.cpp $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx 
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/CholeskySolver$(SELD_SUFFIX).o : lib/Compil/Seldon/CholeskySolver.cpp $(CHESELDON)/computation/solver/SparseCholeskyFactorisation.cxx $(CHESELDON)/computation/solver/DistributedCholeskySolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o : lib/Compil/Seldon/EigenvalueSolver.cpp $(CHESELDON)/computation/interfaces/eigenvalue/VirtualEigenvalueSolver.cxx $(CHESELDON)/computation/interfaces/eigenvalue/PolynomialEigenvalueSolver.cxx $(CHESELDON)/computation/interfaces/eigenvalue/NonLinearEigenvalueSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Arpack$(SELD_SUFFIX).o : lib/Compil/Seldon/Arpack.cpp $(CHESELDON)/computation/interfaces/eigenvalue/ArpackSolver.cxx $(CHESELDON)/computation/interfaces/eigenvalue/Arpack.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Slepc$(SELD_SUFFIX).o : lib/Compil/Seldon/Slepc.cpp $(CHESELDON)/computation/interfaces/eigenvalue/Slepc.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Feast$(SELD_SUFFIX).o : lib/Compil/Seldon/Feast.cpp $(CHESELDON)/computation/interfaces/eigenvalue/Feast.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Anasazi$(SELD_SUFFIX).o : lib/Compil/Seldon/Anasazi.cpp $(CHESELDON)/computation/interfaces/eigenvalue/Anasazi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/IlutPreconditioning$(SELD_SUFFIX).o : lib/Compil/Seldon/IlutPreconditioning.cpp $(CHESELDON)/computation/solver/preconditioner/IlutPreconditioning.cxx $(CHESELDON)/computation/solver/preconditioner/SymmetricIlutPreconditioning.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/IOMatrixMarket$(SELD_SUFFIX).o : lib/Compil/Seldon/IOMatrixMarket.cpp $(CHESELDON)/matrix_sparse/IOMatrixMarket.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/IterativeSolver$(SELD_SUFFIX).o : lib/Compil/Seldon/IterativeSolver.cpp $(CHESELDON)/computation/solver/iterative/Iterative.cxx $(CHESELDON)/computation/solver/preconditioner/Precond_Ssor.cxx $(CHESELDON)/computation/solver/iterative/Cg.cxx $(CHESELDON)/computation/solver/iterative/Cgne.cxx $(CHESELDON)/computation/solver/iterative/Lsqr.cxx $(CHESELDON)/computation/solver/iterative/Cgs.cxx $(CHESELDON)/computation/solver/iterative/BiCg.cxx $(CHESELDON)/computation/solver/iterative/BiCgStab.cxx $(CHESELDON)/computation/solver/iterative/BiCgStabl.cxx $(CHESELDON)/computation/solver/iterative/BiCgcr.cxx $(CHESELDON)/computation/solver/iterative/Gcr.cxx $(CHESELDON)/computation/solver/iterative/CoCg.cxx $(CHESELDON)/computation/solver/iterative/Gmres.cxx $(CHESELDON)/computation/solver/iterative/MinRes.cxx $(CHESELDON)/computation/solver/iterative/Qmr.cxx $(CHESELDON)/computation/solver/iterative/QmrSym.cxx $(CHESELDON)/computation/solver/iterative/QCgs.cxx $(CHESELDON)/computation/solver/iterative/TfQmr.cxx $(CHESELDON)/computation/solver/iterative/Symmlq.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o : lib/Compil/Seldon/Lapack.cpp $(CHESELDON)/computation/interfaces/Lapack_LinearEquations.cxx $(CHESELDON)/computation/interfaces/Lapack_LeastSquares.cxx $(CHESELDON)/computation/interfaces/Lapack_Eigenvalues.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o : lib/Compil/Seldon/LapackMpfr.cpp $(CHESELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx src/Algebra/Eigenvalue.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o : lib/Compil/Seldon/Cholmod.cpp $(CHESELDON)/computation/interfaces/direct/Cholmod.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o : lib/Compil/Seldon/Mumps.cpp $(CHESELDON)/computation/interfaces/direct/Mumps.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o : lib/Compil/Seldon/Pastix.cpp $(CHESELDON)/computation/interfaces/direct/Pastix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Pardiso$(SELD_SUFFIX).o : lib/Compil/Seldon/Pardiso.cpp $(CHESELDON)/computation/interfaces/direct/Pardiso.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/UmfPack$(SELD_SUFFIX).o : lib/Compil/Seldon/UmfPack.cpp $(CHESELDON)/computation/interfaces/direct/UmfPack.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/SuperLU$(SELD_SUFFIX).o : lib/Compil/Seldon/SuperLU.cpp $(CHESELDON)/computation/interfaces/direct/SuperLU.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Wsmp$(SELD_SUFFIX).o : lib/Compil/Seldon/Wsmp.cpp $(CHESELDON)/computation/interfaces/direct/Wsmp.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/SparseSeldonSolver$(SELD_SUFFIX).o : lib/Compil/Seldon/SparseSeldonSolver.cpp $(CHESELDON)/computation/solver/Ordering.cxx $(CHESELDON)/computation/solver/SparseSolver.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/SparseDirectSolver$(SELD_SUFFIX).o : lib/Compil/Seldon/SparseDirectSolver.cpp $(CHESELDON)/computation/interfaces/direct/SparseDirectSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/DistributedSolver$(SELD_SUFFIX).o : lib/Compil/Seldon/DistributedSolver.cpp $(CHESELDON)/computation/solver/DistributedSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/Hypre$(SELD_SUFFIX).o : lib/Compil/Seldon/Hypre.cpp $(CHESELDON)/computation/solver/preconditioner/Hypre.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Seldon/check_dim$(SELD_SUFFIX).o : lib/Compil/Seldon/check_dim.cpp $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)


$(LIB_SELDON_STATIC) :  $(OBJ_SELDON)
	ar rv $(LIB_SELDON_STATIC) $(OBJ_SELDON)

$(LIB_SELDON) :  $(OBJ_SELDON)
	$(CC) -shared -Wl,-soname,libseldon$(SELD_SUFFIX).so -rdynamic -o $(LIB_SELDON) $(OBJ_SELDON)

# files in folder Algebra and Regularization
SRC_ALGEBRA := lib/Compil/Algebra/additional_function.cpp lib/Compil/Algebra/matrix_blockdiagonal.cpp lib/Compil/Algebra/skyline_matrix.cpp lib/Compil/Algebra/tiny_band_matrix.cpp lib/Compil/Algebra/tiny_vector.cpp lib/Compil/Algebra/distributed_block_diag.cpp lib/Compil/Algebra/scalapack.cpp lib/Compil/Algebra/itreg.cpp
LIB_ALGEBRA := $(SRC_ALGEBRA:.cpp=$(MONT_SUFFIX).o)

# dependances for each file in folder lib/Compil/Algebra
lib/Compil/Algebra/additional_function$(MONT_SUFFIX).o : lib/Compil/Algebra/additional_function.cpp src/Algebra/Eigenvalue.cxx src/Algebra/FactorisationLU.cxx src/Algebra/FunctionMatrixExtraction.cxx $(CHESELDON)/matrix_sparse/Functions_MatrixArray.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/tiny_band_matrix$(MONT_SUFFIX).o : lib/Compil/Algebra/tiny_band_matrix.cpp src/Algebra/TinyBandMatrix.cxx src/Algebra/TinyBlockSolver1D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/tiny_vector$(MONT_SUFFIX).o : lib/Compil/Algebra/tiny_vector.cpp $(CHESELDON)/vector/Vector.cxx $(CHESELDON)/matrix/TinyMatrix.cxx $(CHESELDON)/vector/TinyVector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/skyline_matrix$(MONT_SUFFIX).o : lib/Compil/Algebra/skyline_matrix.cpp src/Algebra/GeneralSkyLineMatrix.cxx src/Algebra/SymmetricSkyLineMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/matrix_blockdiagonal$(MONT_SUFFIX).o : lib/Compil/Algebra/matrix_blockdiagonal.cpp src/Algebra/MatrixBlockDiagonal.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/distributed_block_diag$(MONT_SUFFIX).o : lib/Compil/Algebra/distributed_block_diag.cpp src/Algebra/DistributedBlockDiagonalMatrix.cxx $(CHESELDON)/matrix_sparse/DistributedMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/scalapack$(MONT_SUFFIX).o : lib/Compil/Algebra/scalapack.cpp src/Algebra/Scalapack.cxx src/Algebra/DistributedDenseMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/itreg$(MONT_SUFFIX).o : lib/Compil/Algebra/itreg.cpp src/Regularization/ForwardOperator.cxx src/Regularization/RegularizationMethod.cxx src/Regularization/Landweber.cxx src/Regularization/NewtonInner.cxx src/Regularization/NewtonTikh.cxx src/Regularization/IrgnmCg.cxx src/Regularization/MinimizationMethod.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Share
SRC_SHARE := lib/Compil/Share/bessel_function.cpp lib/Compil/Share/common.cpp lib/Compil/Share/fft_interface.cpp lib/Compil/Share/polynomial.cpp lib/Compil/Share/random.cpp lib/Compil/Share/timer.cpp
LIB_SHARE := $(SRC_SHARE:.cpp=$(MONT_SUFFIX).o)

# dependances for each file in folder lib/Compil/Share
lib/Compil/Share/bessel_function$(MONT_SUFFIX).o : lib/Compil/Share/bessel_function.cpp src/Share/BesselFunctionsInterface.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Share/common$(MONT_SUFFIX).o : lib/Compil/Share/common.cpp src/Share/MontjoieTypes.cxx src/Share/CommonMontjoie.cxx src/Share/Precision.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Share/fft_interface$(MONT_SUFFIX).o : lib/Compil/Share/fft_interface.cpp src/Share/FFT.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Share/polynomial$(MONT_SUFFIX).o : lib/Compil/Share/polynomial.cpp src/Share/UnivariatePolynomial.cxx src/Share/MultivariatePolynomial.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Share/random$(MONT_SUFFIX).o : lib/Compil/Share/random.cpp src/Share/RandomGenerator.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Share/timer$(MONT_SUFFIX).o : lib/Compil/Share/timer.cpp src/Share/Timer.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Quadrature
SRC_QUAD := $(wildcard lib/Compil/Quadrature/*.cpp)
LIB_QUAD := $(SRC_QUAD:.cpp=$(MONT_SUFFIX).o)

# dependances for each file in folder lib/Compil/Quadrature
lib/Compil/Quadrature/gauss_jacobi$(MONT_SUFFIX).o : lib/Compil/Quadrature/gauss_jacobi.cpp src/Quadrature/GaussJacobi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/gauss_lobatto_points$(MONT_SUFFIX).o : lib/Compil/Quadrature/gauss_lobatto_points.cpp src/Quadrature/GaussLobattoPoints.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/gauss_lobatto$(MONT_SUFFIX).o : lib/Compil/Quadrature/gauss_lobatto.cpp src/Quadrature/GaussLobatto.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/singular_integration2D$(MONT_SUFFIX).o : lib/Compil/Quadrature/singular_integration2D.cpp src/Quadrature/SingularIntegration2D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/singular_integration3D$(MONT_SUFFIX).o : lib/Compil/Quadrature/singular_integration3D.cpp src/Quadrature/SingularIntegration3D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/triangle_quadrature$(MONT_SUFFIX).o : lib/Compil/Quadrature/triangle_quadrature.cpp src/Quadrature/TriangleQuadrature.cxx src/Quadrature/QuadrangleQuadrature.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/tetrahedron_quadrature$(MONT_SUFFIX).o : lib/Compil/Quadrature/tetrahedron_quadrature.cpp src/Quadrature/TetrahedronQuadrature.cxx src/Quadrature/PyramidQuadrature.cxx src/Quadrature/HexahedronQuadrature.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)


# files in folder Mesh
SRC_MESH := $(wildcard lib/Compil/Mesh/*.cpp)
LIB_MESH := $(SRC_MESH:.cpp=$(MONT_SUFFIX).o)

lib/Compil/Mesh/function_mesh$(MONT_SUFFIX).o : lib/Compil/Mesh/function_mesh.cpp src/Mesh/FunctionsMesh.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/mesh1d$(MONT_SUFFIX).o : lib/Compil/Mesh/mesh1d.cpp src/Mesh/Mesh1D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/mesh2d$(MONT_SUFFIX).o : lib/Compil/Mesh/mesh2d.cpp src/Mesh/Mesh2D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/mesh3d$(MONT_SUFFIX).o : lib/Compil/Mesh/mesh3d.cpp src/Mesh/Mesh3D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/mesh_base$(MONT_SUFFIX).o : lib/Compil/Mesh/mesh_base.cpp src/Mesh/MeshBase.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/mesh_boundaries$(MONT_SUFFIX).o : lib/Compil/Mesh/mesh_boundaries.cpp src/Mesh/MeshBoundaries.cxx src/Mesh/Mesh2DBoundaries.cxx src/Mesh/Mesh3DBoundaries.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/mesh_element$(MONT_SUFFIX).o : lib/Compil/Mesh/mesh_element.cpp src/Mesh/MeshElement.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/number_mesh$(MONT_SUFFIX).o : lib/Compil/Mesh/number_mesh.cpp src/Mesh/NumberMap.cxx src/Mesh/NumberMesh.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Mesh/split_mesh$(MONT_SUFFIX).o : lib/Compil/Mesh/split_mesh.cpp src/Mesh/ParallelMeshFunctions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Output
SRC_OUTPUT := $(wildcard lib/Compil/Output/*.cpp)
LIB_OUTPUT := $(SRC_OUTPUT:.cpp=$(MONT_SUFFIX).o)

LIB_OUTPUT_EXT := lib/Compil/Output/common_output$(MONT_SUFFIX).o lib/Compil/Output/grid_interpolation$(MONT_SUFFIX).o lib/Compil/Output/mesh_interpolation$(MONT_SUFFIX).o lib/Compil/Output/parameter_output$(MONT_SUFFIX).o

LIB_OUTPUT_HARMO := lib/Compil/Output/var_output_problem$(MONT_SUFFIX).o

# dependances for each file in folder lib/Compil/Output
lib/Compil/Output/common_output$(MONT_SUFFIX).o : lib/Compil/Output/common_output.cpp src/Output/CommonInputOutput.cxx src/Output/SplineInterpolation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/grid_interpolation$(MONT_SUFFIX).o : lib/Compil/Output/grid_interpolation.cpp src/Output/GridInterpolation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/mesh_interpolation$(MONT_SUFFIX).o : lib/Compil/Output/mesh_interpolation.cpp src/Output/MeshInterpolation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/mesh_interpolation_fem$(MONT_SUFFIX).o : lib/Compil/Output/mesh_interpolation_fem.cpp src/Output/MeshInterpolationFEM.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/parameter_output$(MONT_SUFFIX).o : lib/Compil/Output/parameter_output.cpp src/Output/ParameterOutput.cxx src/Output/ParameterOutputReprise.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/var_output_problem$(MONT_SUFFIX).o : lib/Compil/Output/var_output_problem.cpp src/Output/OutputHarmonic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder FiniteElement
SRC_FINITE := lib/Compil/FiniteElement/edge_reference.cpp lib/Compil/FiniteElement/face_geom.cpp lib/Compil/FiniteElement/points_reference.cpp lib/Compil/FiniteElement/volume_geom.cpp lib/Compil/FiniteElement/projector.cpp lib/Compil/FiniteElement/fe_h1.cpp lib/Compil/FiniteElement/fe_hcurl2D.cpp lib/Compil/FiniteElement/fe_hcurl3D.cpp lib/Compil/FiniteElement/fe_hdiv.cpp lib/Compil/FiniteElement/element_reference.cpp lib/Compil/FiniteElement/face_reference.cpp lib/Compil/FiniteElement/volume_reference.cpp lib/Compil/FiniteElement/triangle_h1.cpp lib/Compil/FiniteElement/quadrangle_h1.cpp lib/Compil/FiniteElement/triangle_hcurl.cpp lib/Compil/FiniteElement/quadrangle_hcurl.cpp lib/Compil/FiniteElement/triangle_hdiv.cpp lib/Compil/FiniteElement/quadrangle_hdiv.cpp lib/Compil/FiniteElement/tetrahedron_h1.cpp lib/Compil/FiniteElement/pyramid_h1.cpp lib/Compil/FiniteElement/wedge_h1.cpp lib/Compil/FiniteElement/hexahedron_h1.cpp lib/Compil/FiniteElement/tetrahedron_hcurl.cpp lib/Compil/FiniteElement/pyramid_hcurl.cpp lib/Compil/FiniteElement/wedge_hcurl.cpp lib/Compil/FiniteElement/hexahedron_hcurl.cpp lib/Compil/FiniteElement/tetrahedron_hdiv.cpp lib/Compil/FiniteElement/pyramid_hdiv.cpp lib/Compil/FiniteElement/wedge_hdiv.cpp lib/Compil/FiniteElement/hexahedron_hdiv.cpp lib/Compil/FiniteElement/var_finite_element.cpp
LIB_FINITE := $(SRC_FINITE:.cpp=$(MONT_SUFFIX).o)

lib/Compil/FiniteElement/edge_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/edge_reference.cpp src/FiniteElement/Edge/EdgeReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/face_geom$(MONT_SUFFIX).o : lib/Compil/FiniteElement/face_geom.cpp src/FiniteElement/ElementGeomReference.cxx src/FiniteElement/FaceGeomReference.cxx src/FiniteElement/Triangle/TriangleGeomReference.cxx src/FiniteElement/Quadrangle/QuadrangleGeomReference.cxx $(CHESELDON)/vector/Vector.cxx $(CHESELDON)/vector/Functions_Arrays.cxx $(CHESELDON)/matrix/Matrix_Base.cxx $(CHESELDON)/matrix/Matrix_Pointers.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/points_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/points_reference.cpp src/FiniteElement/PointsReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/projector$(MONT_SUFFIX).o : lib/Compil/FiniteElement/projector.cpp src/FiniteElement/ProjectionOperator.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/volume_geom$(MONT_SUFFIX).o : lib/Compil/FiniteElement/volume_geom.cpp src/FiniteElement/VolumeGeomReference.cxx src/FiniteElement/Tetrahedron/TetrahedronGeomReference.cxx src/FiniteElement/Pyramid/PyramidGeomReference.cxx src/FiniteElement/Wedge/WedgeGeomReference.cxx src/FiniteElement/Hexahedron/HexahedronGeomReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/fe_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/fe_h1.cpp src/FiniteElement/FiniteElementH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/fe_hcurl2D$(MONT_SUFFIX).o : lib/Compil/FiniteElement/fe_hcurl2D.cpp src/FiniteElement/FiniteElementHcurl2D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/fe_hcurl3D$(MONT_SUFFIX).o : lib/Compil/FiniteElement/fe_hcurl3D.cpp src/FiniteElement/FiniteElementHcurl3D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/fe_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/fe_hdiv.cpp src/FiniteElement/FiniteElementHdiv.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/element_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/element_reference.cpp src/FiniteElement/ElementReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/face_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/face_reference.cpp src/FiniteElement/FaceReference.cxx src/FiniteElement/Triangle/TriangleReference.cxx src/FiniteElement/Quadrangle/QuadrangleReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/volume_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/volume_reference.cpp src/FiniteElement/VolumeReference.cxx src/FiniteElement/Tetrahedron/TetrahedronReference.cxx src/FiniteElement/Pyramid/PyramidReference.cxx src/FiniteElement/Wedge/WedgeReference.cxx src/FiniteElement/Hexahedron/HexahedronReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/triangle_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/triangle_h1.cpp src/FiniteElement/Triangle/TriangleClassical.cxx src/FiniteElement/Triangle/TriangleQuasiLumped.cxx src/FiniteElement/Triangle/TriangleHierarchic.cxx src/FiniteElement/Triangle/TriangleDgMassLumped.cxx src/FiniteElement/Triangle/TriangleDgOrtho.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/quadrangle_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/quadrangle_h1.cpp src/FiniteElement/Quadrangle/QuadrangleGauss.cxx src/FiniteElement/Quadrangle/QuadrangleLobatto.cxx src/FiniteElement/Quadrangle/QuadrangleRadau.cxx src/FiniteElement/Quadrangle/QuadrangleHierarchic.cxx src/FiniteElement/Quadrangle/QuadrangleDgGauss.cxx src/FiniteElement/Quadrangle/QuadrangleDgOrtho.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/triangle_hcurl$(MONT_SUFFIX).o : lib/Compil/FiniteElement/triangle_hcurl.cpp src/FiniteElement/Triangle/TriangleHcurlFirstFamily.cxx src/FiniteElement/Triangle/TriangleHcurlSecondFamily.cxx src/FiniteElement/Triangle/TriangleHcurlHierarchic.cxx src/FiniteElement/Triangle/TriangleHcurlOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/quadrangle_hcurl$(MONT_SUFFIX).o : lib/Compil/FiniteElement/quadrangle_hcurl.cpp src/FiniteElement/Quadrangle/QuadrangleHcurlFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHcurlGaussFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHcurlOptimalFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHcurlLobatto.cxx src/FiniteElement/Quadrangle/QuadrangleHcurlOptimalHpFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHcurlHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/triangle_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/triangle_hdiv.cpp src/FiniteElement/Triangle/TriangleHdivFirstFamily.cxx src/FiniteElement/Triangle/TriangleHdivOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/quadrangle_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/quadrangle_hdiv.cpp src/FiniteElement/Quadrangle/QuadrangleHdivFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHdivOptimalFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHdivOptimalHpFirstFamily.cxx src/FiniteElement/Quadrangle/QuadrangleHdivHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/tetrahedron_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/tetrahedron_h1.cpp src/FiniteElement/Tetrahedron/TetrahedronClassical.cxx src/FiniteElement/Tetrahedron/TetrahedronHierarchic.cxx src/FiniteElement/Tetrahedron/TetrahedronDgOrtho.cxx src/FiniteElement/Tetrahedron/TetrahedronQuasiLumped.cxx src/FiniteElement/Tetrahedron/TetrahedronMassLumped.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/pyramid_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/pyramid_h1.cpp src/FiniteElement/Pyramid/PyramidClassical.cxx src/FiniteElement/Pyramid/PyramidHierarchic.cxx src/FiniteElement/Pyramid/PyramidDgOrtho.cxx src/FiniteElement/Pyramid/PyramidDgLegendre.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/wedge_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/wedge_h1.cpp src/FiniteElement/Wedge/WedgeClassical.cxx src/FiniteElement/Wedge/WedgeHierarchic.cxx src/FiniteElement/Wedge/WedgeDgClassical.cxx src/FiniteElement/Wedge/WedgeDgOrtho.cxx src/FiniteElement/Wedge/WedgeDgLegendre.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/hexahedron_h1$(MONT_SUFFIX).o : lib/Compil/FiniteElement/hexahedron_h1.cpp src/FiniteElement/Hexahedron/HexahedronGauss.cxx src/FiniteElement/Hexahedron/HexahedronLobatto.cxx src/FiniteElement/Hexahedron/HexahedronHierarchic.cxx src/FiniteElement/Hexahedron/HexahedronDgGauss.cxx src/FiniteElement/Hexahedron/HexahedronDgLegendre.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/tetrahedron_hcurl$(MONT_SUFFIX).o : lib/Compil/FiniteElement/tetrahedron_hcurl.cpp src/FiniteElement/Tetrahedron/TetrahedronHcurlFirstFamily.cxx src/FiniteElement/Tetrahedron/TetrahedronHcurlOptimalHpFirstFamily.cxx src/FiniteElement/Tetrahedron/TetrahedronHcurlLobatto.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/pyramid_hcurl$(MONT_SUFFIX).o : lib/Compil/FiniteElement/pyramid_hcurl.cpp src/FiniteElement/Pyramid/PyramidHcurlFirstFamily.cxx src/FiniteElement/Pyramid/PyramidHcurlOptimalFirstFamily.cxx src/FiniteElement/Pyramid/PyramidHcurlHpFirstFamily.cxx src/FiniteElement/Pyramid/PyramidHcurlOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/wedge_hcurl$(MONT_SUFFIX).o : lib/Compil/FiniteElement/wedge_hcurl.cpp src/FiniteElement/Wedge/WedgeHcurlFirstFamily.cxx src/FiniteElement/Wedge/WedgeHcurlOptimalFirstFamily.cxx src/FiniteElement/Wedge/WedgeHcurlHpFirstFamily.cxx src/FiniteElement/Wedge/WedgeHcurlOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/hexahedron_hcurl$(MONT_SUFFIX).o : lib/Compil/FiniteElement/hexahedron_hcurl.cpp src/FiniteElement/Hexahedron/HexahedronHcurlFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHcurlOptimalFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHcurlHpFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHcurlOptimalHpFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHcurlLobatto.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/tetrahedron_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/tetrahedron_hdiv.cpp src/FiniteElement/Tetrahedron/TetrahedronHdivFirstFamily.cxx src/FiniteElement/Tetrahedron/TetrahedronHdivOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/pyramid_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/pyramid_hdiv.cpp src/FiniteElement/Pyramid/PyramidHdivFirstFamily.cxx src/FiniteElement/Pyramid/PyramidHdivOptimalFirstFamily.cxx src/FiniteElement/Pyramid/PyramidHdivHpFirstFamily.cxx src/FiniteElement/Pyramid/PyramidHdivOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/wedge_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/wedge_hdiv.cpp src/FiniteElement/Wedge/WedgeHdivFirstFamily.cxx src/FiniteElement/Wedge/WedgeHdivOptimalFirstFamily.cxx src/FiniteElement/Wedge/WedgeHdivHpFirstFamily.cxx src/FiniteElement/Wedge/WedgeHdivOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/hexahedron_hdiv$(MONT_SUFFIX).o : lib/Compil/FiniteElement/hexahedron_hdiv.cpp src/FiniteElement/Hexahedron/HexahedronHdivFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHdivOptimalFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHdivHpFirstFamily.cxx src/FiniteElement/Hexahedron/HexahedronHdivOptimalHpFirstFamily.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/var_finite_element$(MONT_SUFFIX).o : lib/Compil/FiniteElement/var_finite_element.cpp src/Harmonic/VarFiniteElement.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Solver
SRC_SOLVER := $(wildcard lib/Compil/Solver/*.cpp)
LIB_SOLVER := $(SRC_SOLVER:.cpp=$(MONT_SUFFIX).o)

LIB_SOLVER_EXT := lib/Compil/Solver/non_linear_solver$(MONT_SUFFIX).o lib/Compil/Solver/linear_solver$(MONT_SUFFIX).o

LIB_SOLVER_HARMO := lib/Compil/Solver/solver_harmonic$(MONT_SUFFIX).o lib/Compil/Solver/eigenvalue_harmonic$(MONT_SUFFIX).o

# dependances for each file in folder lib/Compil/Solver
lib/Compil/Solver/non_linear_solver$(MONT_SUFFIX).o : lib/Compil/Solver/non_linear_solver.cpp src/Solver/NonLinearEquations.cxx src/Solver/NewtonSolver.cxx src/Solver/NonLinearLeastSquares.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Solver/linear_solver$(MONT_SUFFIX).o : lib/Compil/Solver/linear_solver.cpp src/Solver/SolveSystem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Solver/solver_harmonic$(MONT_SUFFIX).o : lib/Compil/Solver/solver_harmonic.cpp src/Solver/SolveHarmonic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Solver/preconditioner$(MONT_SUFFIX).o : lib/Compil/Solver/preconditioner.cpp src/Solver/Preconditioner.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Solver/eigenvalue_harmonic$(MONT_SUFFIX).o : lib/Compil/Solver/eigenvalue_harmonic.cpp src/Solver/EigenvaluesHarmonic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)


# files in folder Harmonic
SRC_HARMONIC := $(wildcard lib/Compil/Harmonic/*.cpp)
LIB_HARMONIC := $(SRC_HARMONIC:.cpp=$(MONT_SUFFIX).o)

lib/Compil/Harmonic/assemble_matrix$(MONT_SUFFIX).o : lib/Compil/Harmonic/assemble_matrix.cpp src/Computation/AssembleMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/fem_matrix_free_class$(MONT_SUFFIX).o : lib/Compil/Harmonic/fem_matrix_free_class.cpp src/Computation/FemMatrixFreeClass.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/generic_equation$(MONT_SUFFIX).o : lib/Compil/Harmonic/generic_equation.cpp src/Harmonic/GenericEquation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_problem1d$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_problem1d.cpp src/Elliptic/VarProblem1D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_problem_base$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_problem_base.cpp src/Harmonic/VarProblemBase.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_axisym_problem$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_axisym_problem.cpp src/Harmonic/VarAxisymProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_geometry_problem$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_geometry_problem.cpp src/Harmonic/VarGeometryProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/distributed_problem$(MONT_SUFFIX).o : lib/Compil/Harmonic/distributed_problem.cpp src/Harmonic/DistributedProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_harmonic_base$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_harmonic_base.cpp src/Harmonic/VarHarmonicBase.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_problem$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_problem.cpp src/Harmonic/VarProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_boundary_condition$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_boundary_condition.cpp src/Harmonic/BoundaryConditionHarmonic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/multi_freq$(MONT_SUFFIX).o : lib/Compil/Harmonic/multi_freq.cpp src/Harmonic/MultiFrequencyProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/transmission_model$(MONT_SUFFIX).o : lib/Compil/Harmonic/transmission_model.cpp src/Harmonic/TransmissionModel.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/gibc_model$(MONT_SUFFIX).o : lib/Compil/Harmonic/gibc_model.cpp src/Harmonic/GeneralizedImpedanceModel.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/transparency_solver$(MONT_SUFFIX).o : lib/Compil/Harmonic/transparency_solver.cpp src/Harmonic/TransparencyCondition.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/source_spatiale$(MONT_SUFFIX).o : lib/Compil/Harmonic/source_spatiale.cpp src/Source/SourceSpatiale.cxx $(CHESELDON)/vector/Vector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/modal_source$(MONT_SUFFIX).o : lib/Compil/Harmonic/modal_source.cpp src/Source/ModalSource.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/var_source_problem$(MONT_SUFFIX).o : lib/Compil/Harmonic/var_source_problem.cpp src/Source/DefineSourceElliptic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/user_source$(MONT_SUFFIX).o : lib/Compil/Harmonic/user_source.cpp src/Source/UserSource.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/physical_constant$(MONT_SUFFIX).o : lib/Compil/Harmonic/physical_constant.cpp src/Elliptic/Maxwell/PhysicalConstant.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Harmonic/physical_index$(MONT_SUFFIX).o : lib/Compil/Harmonic/physical_index.cpp src/Elliptic/PhysicalProperty.cxx src/Elliptic/OneDimPhysicalIndex.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Instationary
SRC_UNSTEADY := $(wildcard lib/Compil/Instationary/*.cpp)
LIB_UNSTEADY := $(SRC_UNSTEADY:.cpp=$(MONT_SUFFIX).o)

LIB_UNSTEADY_EXT := lib/Compil/Instationary/explicit_time_scheme$(MONT_SUFFIX).o lib/Compil/Instationary/implicit_time_scheme$(MONT_SUFFIX).o lib/Compil/Instationary/modified_equation$(MONT_SUFFIX).o lib/Compil/Instationary/local_time_scheme$(MONT_SUFFIX).o

LIB_UNSTEADY_HARMO := lib/Compil/Instationary/var_instationary$(MONT_SUFFIX).o lib/Compil/Instationary/continuous_mass_matrix$(MONT_SUFFIX).o lib/Compil/Instationary/discontinuous_mass_matrix$(MONT_SUFFIX).o lib/Compil/Instationary/time_source$(MONT_SUFFIX).o lib/Compil/Instationary/evaluation_cfl_wave$(MONT_SUFFIX).o

lib/Compil/Instationary/explicit_time_scheme$(MONT_SUFFIX).o : lib/Compil/Instationary/explicit_time_scheme.cpp src/Instationary/ExplicitTimeSchemes.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/implicit_time_scheme$(MONT_SUFFIX).o : lib/Compil/Instationary/implicit_time_scheme.cpp src/Instationary/ImplicitTimeSchemes.cxx src/Instationary/TimeSchemes.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/local_time_scheme$(MONT_SUFFIX).o : lib/Compil/Instationary/local_time_scheme.cpp src/Instationary/LocalTimeSchemes.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/modified_equation$(MONT_SUFFIX).o : lib/Compil/Instationary/modified_equation.cpp src/Instationary/ModifiedEquation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/var_instationary$(MONT_SUFFIX).o : lib/Compil/Instationary/var_instationary.cpp src/Instationary/VarInstationary.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/continuous_mass_matrix$(MONT_SUFFIX).o : lib/Compil/Instationary/continuous_mass_matrix.cpp src/Instationary/ContinuousMassMatrix.cxx src/Instationary/ContinuousUnsteadyMassMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/discontinuous_mass_matrix$(MONT_SUFFIX).o : lib/Compil/Instationary/discontinuous_mass_matrix.cpp src/Instationary/DiscontinuousMassMatrix.cxx src/Instationary/DiscontinuousUnsteadyMassMatrix.cxx $(CHESELDON)/vector/Vector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/time_source$(MONT_SUFFIX).o : lib/Compil/Instationary/time_source.cpp src/Source/TimeSource.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Instationary/evaluation_cfl_wave$(MONT_SUFFIX).o : lib/Compil/Instationary/evaluation_cfl_wave.cpp src/Instationary/EvaluationCflWave.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Elliptic
SRC_ELLIPTIC := $(wildcard lib/Compil/Elliptic/*/*.cpp) $(wildcard lib/Compil/Nathan/*.cpp)

LIB_ELLIPTIC := $(SRC_ELLIPTIC:.cpp=$(MONT_SUFFIX).o)

lib/Compil/Elliptic/Aeroacoustic/aeroacoustic$(MONT_SUFFIX).o : lib/Compil/Elliptic/Aeroacoustic/aeroacoustic.cpp src/Elliptic/Aeroacoustic/AeroAcoustic.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Aeroacoustic/harmonic_galbrun$(MONT_SUFFIX).o : lib/Compil/Elliptic/Aeroacoustic/harmonic_galbrun.cpp src/Elliptic/Aeroacoustic/HarmonicGalbrun.cxx src/Elliptic/Aeroacoustic/HarmonicGalbrunEquation.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Aeroacoustic/symmetric_harmonic_galbrun$(MONT_SUFFIX).o : lib/Compil/Elliptic/Aeroacoustic/symmetric_harmonic_galbrun.cpp src/Elliptic/Aeroacoustic/SymmetricHarmonicGalbrun.cxx src/Elliptic/Aeroacoustic/HarmonicGalbrunEquation.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Aeroacoustic/lee_axi$(MONT_SUFFIX).o : lib/Compil/Elliptic/Aeroacoustic/lee_axi.cpp src/Elliptic/Aeroacoustic/LinearizedEulerEquationAxisym.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Aeroacoustic/axisym_galbrun$(MONT_SUFFIX).o : lib/Compil/Elliptic/Aeroacoustic/axisym_galbrun.cpp src/Elliptic/Aeroacoustic/AxiSymGalbrun.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)


lib/Compil/Nathan/poisson$(MONT_SUFFIX).o: lib/Compil/Nathan/poisson.cpp src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx 
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Elastic/var_elastic$(MONT_SUFFIX).o : lib/Compil/Elliptic/Elastic/var_elastic.cpp src/Elliptic/Elastic/VarElastic.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Elastic/model_elas$(MONT_SUFFIX).o : lib/Compil/Elliptic/Elastic/model_elas.cpp src/Elliptic/Elastic/TransmissionModelElastic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/helmholtz1d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/helmholtz1d.cpp src/Elliptic/Helmholtz/Helmholtz1D.cxx src/Elliptic/VarProblem1D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/helmholtz_radial$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/helmholtz_radial.cpp src/Elliptic/Helmholtz/HelmholtzRadial.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/helmholtz_polar$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/helmholtz_polar.cpp src/Elliptic/Helmholtz/HelmholtzPolar.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/source_helmholtz$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/source_helmholtz.cpp src/Elliptic/Helmholtz/DefineSourceHelmholtz.cxx src/Elliptic/Helmholtz/LaplacianModalSource.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/var_helmholtz$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/var_helmholtz.cpp src/Elliptic/Helmholtz/VarHelmholtz.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/prod_var_helmholtz$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/prod_var_helmholtz.cpp src/Elliptic/Helmholtz/ProdMatVectHelmholtz.cxx src/Computation/ProdMatVectScalarH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/var_laplace$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/var_laplace.cpp src/Elliptic/Helmholtz/VarLaplace.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectScalarH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/model_helmholtz$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/model_helmholtz.cpp src/Elliptic/Helmholtz/ThinSlotHelmholtzModel.cxx src/Elliptic/Helmholtz/TransmissionModelHelmholtz.cxx src/Elliptic/Helmholtz/ImpedanceHelmholtz.cxx $(CHESELDON)/vector/Vector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/axisym_helmholtz$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/axisym_helmholtz.cpp src/Elliptic/Helmholtz/AxiSymHelmholtz.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ProdMatVectScalarH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/cellulaire_periodique$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/cellulaire_periodique.cpp src/Elliptic/Helmholtz/CellulairePeriodique.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/helmholtz_hdiv$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/helmholtz_hdiv.cpp src/Elliptic/Helmholtz/HelmholtzHdiv.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/schrodinger_non_linear1d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/schrodinger_non_linear1d.cpp src/Elliptic/Helmholtz/SchrodingerNonLinear1D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Helmholtz/non_linear_maxwell1d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Helmholtz/non_linear_maxwell1d.cpp src/Elliptic/Helmholtz/NonLinearMaxwell1D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/non_linear_optics$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/non_linear_optics.cpp src/Elliptic/Maxwell/NonLinearOpticsProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/source_maxwell2d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/source_maxwell2d.cpp src/Elliptic/Maxwell/2D/DefineSourceHarmonicMaxwell2D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/harmonic_maxwell2d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/harmonic_maxwell2d.cpp src/Elliptic/Maxwell/2D/HarmonicMaxwell2D.cxx src/Elliptic/Maxwell/2D/StaticMaxwell2D.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ElementaryMatrixHcurl2D.cxx src/Computation/ProdMatVectHcurl2D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/source_maxwell_axi$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/source_maxwell_axi.cpp src/Elliptic/Maxwell/Axi/DefineSourceMaxwellAxi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/harmonic_maxwell_axi$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/harmonic_maxwell_axi.cpp src/Elliptic/Maxwell/Axi/MaxwellAxiSymHarmonic.cxx src/Elliptic/Maxwell/Axi/AxiSymHcurlMaxwell.cxx src/Elliptic/Maxwell/Axi/StaticMaxwellAxi.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/source_maxwell3d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/source_maxwell3d.cpp src/Elliptic/Maxwell/3D/DefineSourceHarmonicMaxwell3D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/harmonic_maxwell3d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/harmonic_maxwell3d.cpp src/Elliptic/Maxwell/3D/HarmonicMaxwell3D.cxx src/Elliptic/Maxwell/3D/StaticMaxwell3D.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectHcurl3D.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/mode_maxwell$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/mode_maxwell.cpp src/Elliptic/Maxwell/3D/ComputationModeMaxwell.cxx src/Harmonic/VarHarmonic.cxx src/Computation/ElementaryMatrixHcurl2D.cxx src/Computation/ProdMatVectHcurl2D.cxx src/Harmonic/GenericEquation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/precond_maxwell3d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/precond_maxwell3d.cpp src/Elliptic/Maxwell/3D/PrecondMaxwell.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Elliptic/Maxwell/model_maxwell3d$(MONT_SUFFIX).o : lib/Compil/Elliptic/Maxwell/model_maxwell3d.cpp src/Elliptic/Maxwell/3D/TransmissionModelMaxwell3D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Hyperbolic
SRC_HYPERBOLIC := $(wildcard lib/Compil/Hyperbolic/*/*.cpp) $(wildcard lib/Compil/Corde/*.cpp)
LIB_HYPERBOLIC := $(SRC_HYPERBOLIC:.cpp=$(MONT_SUFFIX).o)

# dependences for files in folder Hyperbolic
lib/Compil/Hyperbolic/Acoustic/time_acoustic$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Acoustic/time_acoustic.cpp src/Hyperbolic/Acoustic/AxiSymAcoustic.cxx src/Hyperbolic/Acoustic/TimeAcoustic.cxx src/Hyperbolic/Acoustic/ProdMatVectAcoustic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Acoustic/time_reversal$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Acoustic/time_reversal.cpp src/Hyperbolic/Acoustic/TimeReversal.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Acoustic/cross_correlation$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Acoustic/cross_correlation.cpp src/CrossCorrelation/CrossCorrelationFunction.cxx src/CrossCorrelation/CrossCorrelation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Acoustic/advection$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Acoustic/advection.cpp src/Hyperbolic/Acoustic/VarAdvection.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Elastic/time_elastic$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Elastic/time_elastic.cpp src/Hyperbolic/Elastic/TimeElastic.cxx src/Hyperbolic/Elastic/ProdMatVectElastic.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Elastic/reissner_mindlin$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Elastic/reissner_mindlin.cpp src/Hyperbolic/Elastic/ReissnerMindlin.cxx src/Harmonic/VarHarmonic.cxx src/Harmonic/GenericEquation.cxx src/Computation/ElementaryMatrixH1.cxx src/Computation/ProdMatVectH1.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Elastic/fluid_structure$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Elastic/fluid_structure.cpp src/Hyperbolic/Elastic/FluidStructureInteraction.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Elastic/vibro_acoustic$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Elastic/vibro_acoustic.cpp src/Hyperbolic/Elastic/VibroAcousticSession.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Maxwell/time_maxwell2d$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Maxwell/time_maxwell2d.cpp src/Hyperbolic/Maxwell/2D/TimeMaxwell2D.cxx src/Hyperbolic/Maxwell/2D/ProdMatVectMaxwell2D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Maxwell/time_maxwell_axi$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Maxwell/time_maxwell_axi.cpp src/Hyperbolic/Maxwell/Axi/TimeMaxwellAxi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Hyperbolic/Maxwell/time_maxwell3d$(MONT_SUFFIX).o : lib/Compil/Hyperbolic/Maxwell/time_maxwell3d.cpp src/Hyperbolic/Maxwell/3D/TimeMaxwell3D.cxx src/Hyperbolic/Maxwell/3D/ProdMatVectMaxwell3D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Corde/donnees_string$(MONT_SUFFIX).o : lib/Compil/Corde/donnees_string.cpp src/Corde/DonneesString.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Inverse
SRC_INVERSE := $(wildcard lib/Compil/Inverse/*.cpp)
LIB_INVERSE := $(SRC_INVERSE:.cpp=$(MONT_SUFFIX).o)

# dependences for files in folder Inverse
lib/Compil/Inverse/multi_dim_index$(MONT_SUFFIX).o : lib/Compil/Inverse/multi_dim_index.cpp src/Elliptic/MultiDimPhysicalIndex.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Inverse/var_migration$(MONT_SUFFIX).o : lib/Compil/Inverse/var_migration.cpp src/Inverse/VarMigration.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Inverse/boundary_inverse_problem$(MONT_SUFFIX).o : lib/Compil/Inverse/boundary_inverse_problem.cpp src/Inverse/BoundaryInverseProblem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# Montjoie library
OBJ_MONTJOIE :=  $(LIB_SHARE) $(LIB_ALGEBRA) $(LIB_QUAD) $(LIB_OUTPUT) $(LIB_FINITE) $(LIB_MESH) $(LIB_SOLVER) $(LIB_HARMONIC) $(LIB_UNSTEADY)

# object files for only finite element part (not EllipticProblem, VarHarmonic, etc)
LIB_MONTJOIE_FEM := $(LIB_ALGEBRA) $(LIB_SHARE) $(LIB_QUAD) $(LIB_MESH) $(LIB_OUTPUT_EXT) $(LIB_FINITE) $(LIB_SOLVER_EXT) $(LIB_UNSTEADY_EXT)

# object files with dependences on VarHarmonic
LIB_MONTJOIE_HARMO := $(LIB_OUTPUT_HARMO) $(LIB_SOLVER_HARMO) $(LIB_UNSTEADY_HARMO) $(LIB_HARMONIC)

# files related to applications (usually in folders Elliptic, Hyperbolic)
OBJ_MONTJOIE_F := $(LIB_ELLIPTIC) $(LIB_HYPERBOLIC) $(LIB_INVERSE)

$(LIB_MONTJOIE_STATIC) : $(OBJ_MONTJOIE)
	ar rv $(LIB_MONTJOIE_STATIC) $(OBJ_MONTJOIE)

$(LIB_MONTJOIE) : $(OBJ_MONTJOIE)
	$(CC) -shared -Wl,-soname,libmontjoie$(MONT_SUFFIX).so -rdynamic -o $(LIB_MONTJOIE) $(OBJ_MONTJOIE)

$(LIB_MONTJOIE_F_STATIC) : $(OBJ_MONTJOIE_F)
	ar rv $(LIB_MONTJOIE_F_STATIC) $(OBJ_MONTJOIE_F)

$(LIB_MONTJOIE_F) : $(OBJ_MONTJOIE_F)
	$(CC) -shared -Wl,-soname,libmontjoie_full$(MONT_SUFFIX).so -rdynamic -o $(LIB_MONTJOIE_F) $(OBJ_MONTJOIE_F)

# target to remove all the object files
cleanlib :
	rm -f $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_SELDON_STATIC) $(LIB_MONTJOIE_STATIC) $(OBJ_SELDON) $(OBJ_MONTJOIE) $(LIB_MONTJOIE_F) $(OBJ_MONTJOIE_F) $(LIB_MONTJOIE_F_STATIC) src/*.o lib/*.o

# target to remove all the object files, except Seldon
cleanmontjoie :
	rm -f $(LIB_MONTJOIE) $(LIB_MONTJOIE_STATIC) $(OBJ_MONTJOIE) $(LIB_MONTJOIE_F) $(OBJ_MONTJOIE_F) $(LIB_MONTJOIE_F_STATIC)

# target to remove files related to the interface with the direct solvers
cleansolve :
	rm -f $(OBJ_SELDON_SOLVE) lib/Compil/Solver/linear_solver$(MONT_SUFFIX).o lib/Compil/Algebra/scalapack$(MONT_SUFFIX).o

# target to remove files related to the interface with the eigenvalue solvers
cleaneig :
	rm -f $(OBJ_EIG) lib/Compil/Share/common$(MONT_SUFFIX).o lib/Compil/Seldon/Common$(SELD_SUFFIX).o

# target to switch to SELDON_WITH_ABORT
cleanabort :
	rm -f lib/Compil/Seldon/Common$(SELD_SUFFIX).o

# target for flag SELDON_WITH_BLAS/SELDON_WITH_LAPACK
cleanblas :
	rm -f lib/Compil/Algebra/itreg$(MONT_SUFFIX).o lib/Compil/Algebra/additional_function$(MONT_SUFFIX).o lib/Compil/Algebra/Lapack$(MONT_SUFFIX).o $(OBJ_BLAS) lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o

# target for flag MONTJOIE_WITH_GSL
cleangsl :
	rm -f lib/Compil/Algebra/itreg$(MONT_SUFFIX).o lib/Compil/Share/bessel_function$(MONT_SUFFIX).o lib/Compil/Share/random$(MONT_SUFFIX).o lib/Compil/Solver/non_linear_solver$(MONT_SUFFIX).o lib/Compil/Share/fft_interface$(MONT_SUFFIX).o

# target for flag MONTJOIE_WITH_FFTW
cleanfft :
	rm -f lib/Compil/Share/fft_interface$(MONT_SUFFIX).o

# target for finite-element objects
cleanfem :
	rm -f $(LIB_MONTJOIE_FEM) 

# target for harmonic dependances
cleanharmo :
	rm -f $(LIB_MONTJOIE_HARMO) $(LIB_MONTJOIE_F) $(OBJ_MONTJOIE_F)

# elliptic targets
cleanelliptic :
	rm -f $(LIB_ELLIPTIC)

# hyperboic targets
cleanhyperbolic :
	rm -f $(LIB_HYPERBOLIC)

# Uncomment the following line in order to detect the file .o which gives bad result
#LIB_SELDON := $(OBJ_SELDON)

# Uncomment the two following lines in order to detect the file .o which gives bad result
#LIB_MONTJOIE := $(OBJ_MONTJOIE)
#LIB_MONTJOIE_F := $(OBJ_MONTJOIE_F)

ifeq ($(STATIC_COMPILATION),YES)
  LIB_MONTJOIE_F := $(LIB_MONTJOIE_F_STATIC)
  LIB_MONTJOIE := $(LIB_MONTJOIE_STATIC)
  LIB_SELDON := $(LIB_SELDON_STATIC)
endif

ifeq ($(SEPARED_COMPIL),YES)
  LIB := $(LIB_MONTJOIE_F) $(LIB_MONTJOIE) $(LIB_SELDON) $(LIB)
endif
