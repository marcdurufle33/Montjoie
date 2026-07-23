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

lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o : lib/Compil/Seldon/EigenvalueSolver.cpp $(CHESELDON)/computation/interfaces/eigenvalue/VirtualEigenvalueSolver.cxx
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
#lib/Compil/Algebra/tiny_band_matrix.cpp 
SRC_ALGEBRA := lib/Compil/Algebra/additional_function.cpp lib/Compil/Algebra/matrix_blockdiagonal.cpp lib/Compil/Algebra/skyline_matrix.cpp lib/Compil/Algebra/tiny_vector_mesh.cpp lib/Compil/Algebra/distributed_block_diag.cpp lib/Compil/Algebra/tiny_band_mesh.cpp
LIB_ALGEBRA := $(SRC_ALGEBRA:.cpp=$(MONT_SUFFIX).o)

# dependances for each file in folder lib/Compil/Algebra
lib/Compil/Algebra/additional_function$(MONT_SUFFIX).o : lib/Compil/Algebra/additional_function.cpp src/Algebra/Eigenvalue.cxx src/Algebra/FactorisationLU.cxx src/Algebra/FunctionMatrixExtraction.cxx $(CHESELDON)/matrix_sparse/Functions_MatrixArray.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/tiny_band_mesh$(MONT_SUFFIX).o : lib/Compil/Algebra/tiny_band_mesh.cpp src/Algebra/TinyBandMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/tiny_vector_mesh$(MONT_SUFFIX).o : lib/Compil/Algebra/tiny_vector_mesh.cpp $(CHESELDON)/vector/Vector.cxx $(CHESELDON)/matrix/TinyMatrix.cxx $(CHESELDON)/vector/TinyVector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/skyline_matrix$(MONT_SUFFIX).o : lib/Compil/Algebra/skyline_matrix.cpp src/Algebra/GeneralSkyLineMatrix.cxx src/Algebra/SymmetricSkyLineMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/matrix_blockdiagonal$(MONT_SUFFIX).o : lib/Compil/Algebra/matrix_blockdiagonal.cpp src/Algebra/MatrixBlockDiagonal.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Algebra/distributed_block_diag$(MONT_SUFFIX).o : lib/Compil/Algebra/distributed_block_diag.cpp src/Algebra/DistributedBlockDiagonalMatrix.cxx $(CHESELDON)/matrix_sparse/DistributedMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Share
SRC_SHARE := lib/Compil/Share/bessel_function.cpp lib/Compil/Share/common_mesh.cpp lib/Compil/Share/fft_interface.cpp lib/Compil/Share/polynomial.cpp lib/Compil/Share/random.cpp lib/Compil/Share/timer.cpp
LIB_SHARE := $(SRC_SHARE:.cpp=$(MONT_SUFFIX).o)

# dependances for each file in folder lib/Compil/Share
lib/Compil/Share/bessel_function$(MONT_SUFFIX).o : lib/Compil/Share/bessel_function.cpp src/Share/BesselFunctionsInterface.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Share/common_mesh$(MONT_SUFFIX).o : lib/Compil/Share/common_mesh.cpp src/Share/MontjoieTypes.cxx src/Share/CommonMontjoie.cxx src/Share/Precision.cxx
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
SRC_QUAD := lib/Compil/Quadrature/gauss_jacobi.cpp lib/Compil/Quadrature/gauss_lobatto_points.cpp lib/Compil/Quadrature/gauss_lobatto.cpp lib/Compil/Quadrature/triangle_quadrature.cpp lib/Compil/Quadrature/tetrahedron_quadrature.cpp
LIB_QUAD := $(SRC_QUAD:.cpp=$(MONT_SUFFIX).o)

# dependances for each file in folder lib/Compil/Quadrature
lib/Compil/Quadrature/gauss_jacobi$(MONT_SUFFIX).o : lib/Compil/Quadrature/gauss_jacobi.cpp src/Quadrature/GaussJacobi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/gauss_lobatto_points$(MONT_SUFFIX).o : lib/Compil/Quadrature/gauss_lobatto_points.cpp src/Quadrature/GaussLobattoPoints.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Quadrature/gauss_lobatto$(MONT_SUFFIX).o : lib/Compil/Quadrature/gauss_lobatto.cpp src/Quadrature/GaussLobatto.cxx
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

# files in folder FiniteElement
SRC_FINITE := lib/Compil/FiniteElement/edge_reference.cpp lib/Compil/FiniteElement/face_geom_mesh.cpp lib/Compil/FiniteElement/points_reference.cpp lib/Compil/FiniteElement/volume_geom.cpp lib/Compil/FiniteElement/projector_mesh.cpp lib/Compil/FiniteElement/element_reference_mesh.cpp 
LIB_FINITE := $(SRC_FINITE:.cpp=$(MONT_SUFFIX).o)

lib/Compil/FiniteElement/edge_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/edge_reference.cpp src/FiniteElement/Edge/EdgeReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/element_reference_mesh$(MONT_SUFFIX).o : lib/Compil/FiniteElement/element_reference_mesh.cpp src/FiniteElement/ElementReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/face_geom_mesh$(MONT_SUFFIX).o : lib/Compil/FiniteElement/face_geom_mesh.cpp src/FiniteElement/ElementGeomReference.cxx src/FiniteElement/FaceGeomReference.cxx src/FiniteElement/Triangle/TriangleGeomReference.cxx src/FiniteElement/Quadrangle/QuadrangleGeomReference.cxx $(CHESELDON)/vector/Vector.cxx $(CHESELDON)/vector/Functions_Arrays.cxx $(CHESELDON)/matrix/Matrix_Base.cxx $(CHESELDON)/matrix/Matrix_Pointers.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/points_reference$(MONT_SUFFIX).o : lib/Compil/FiniteElement/points_reference.cpp src/FiniteElement/PointsReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/volume_geom$(MONT_SUFFIX).o : lib/Compil/FiniteElement/volume_geom.cpp src/FiniteElement/VolumeGeomReference.cxx src/FiniteElement/Tetrahedron/TetrahedronGeomReference.cxx src/FiniteElement/Pyramid/PyramidGeomReference.cxx src/FiniteElement/Wedge/WedgeGeomReference.cxx src/FiniteElement/Hexahedron/HexahedronGeomReference.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/FiniteElement/projector_mesh$(MONT_SUFFIX).o : lib/Compil/FiniteElement/projector_mesh.cpp src/FiniteElement/ProjectionOperator.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# files in folder Solver
LIB_SOLVER_EXT := lib/Compil/Solver/non_linear_solver$(MONT_SUFFIX).o lib/Compil/Solver/linear_solver$(MONT_SUFFIX).o lib/Compil/Output/common_output$(MONT_SUFFIX).o lib/Compil/Output/grid_interpolation$(MONT_SUFFIX).o lib/Compil/Output/mesh_interpolation$(MONT_SUFFIX).o

# dependances for each file in folder lib/Compil/Solver
lib/Compil/Solver/non_linear_solver$(MONT_SUFFIX).o : lib/Compil/Solver/non_linear_solver.cpp src/Solver/NonLinearEquations.cxx src/Solver/NewtonSolver.cxx src/Solver/NonLinearLeastSquares.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Solver/linear_solver$(MONT_SUFFIX).o : lib/Compil/Solver/linear_solver.cpp src/Solver/SolveSystem.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/common_output$(MONT_SUFFIX).o : lib/Compil/Output/common_output.cpp src/Output/CommonInputOutput.cxx src/Output/SplineInterpolation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/grid_interpolation$(MONT_SUFFIX).o : lib/Compil/Output/grid_interpolation.cpp src/Output/GridInterpolation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

lib/Compil/Output/mesh_interpolation$(MONT_SUFFIX).o : lib/Compil/Output/mesh_interpolation.cpp src/Output/MeshInterpolation.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_LIB) $< -o $@ $(AGRESSIVE_OPTIM)

# Montjoie library
OBJ_MONTJOIE :=  $(LIB_SHARE) $(LIB_ALGEBRA) $(LIB_QUAD) $(LIB_FINITE) $(LIB_MESH) $(LIB_SOLVER_EXT)

# object files for only finite element part (not EllipticProblem, VarHarmonic, etc)
LIB_MONTJOIE_FEM := $(LIB_ALGEBRA) $(LIB_SHARE) $(LIB_QUAD) $(LIB_MESH) $(LIB_FINITE) $(LIB_SOLVER_EXT)

$(LIB_MONTJOIE_STATIC) : $(OBJ_MONTJOIE)
	ar rv $(LIB_MONTJOIE_STATIC) $(OBJ_MONTJOIE)

$(LIB_MONTJOIE) : $(OBJ_MONTJOIE)
	$(CC) -shared -Wl,-soname,libmontjoie$(MONT_SUFFIX).so -rdynamic -o $(LIB_MONTJOIE) $(OBJ_MONTJOIE)

# target to remove all the object files
cleanlib :
	rm -f $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_SELDON_STATIC) $(LIB_MONTJOIE_STATIC) $(OBJ_SELDON) $(OBJ_MONTJOIE) src/*.o lib/*.o

# Uncomment the following line in order to detect the file .o which gives bad result
#LIB_SELDON := $(OBJ_SELDON)

# Uncomment the two following lines in order to detect the file .o which gives bad result
#LIB_MONTJOIE := $(OBJ_MONTJOIE)

ifeq ($(STATIC_COMPILATION),YES)
  LIB_MONTJOIE := $(LIB_MONTJOIE_STATIC)
  LIB_SELDON := $(LIB_SELDON_STATIC)
endif

ifeq ($(SEPARED_COMPIL),YES)
  LIB := $(LIB_MONTJOIE) $(LIB_SELDON) $(LIB)
endif
