# main part of Seldon
OBJ_SELDON := $(DIR_SELDON)/lib/Compil/Seldon/Array3D$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixPointers$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixPacked$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixHermSymTriang$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixSparse$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Vector$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Vector2$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixVector$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrix$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixDense$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixArray$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixConversion$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/PermutationScalingMatrix$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/RelaxationMatrixVector$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/IOMatrixMarket$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/IlutPreconditioning$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/IterativeSolver$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/check_dim$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/TinyVector$(SELD_SUFFIX).o

# complex matrices
OBJ_SELDON := $(OBJ_SELDON) $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatVectComplex$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixComplex$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixComplexConversions$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/MatrixComplexSparse$(SELD_SUFFIX).o

# Blas functions
ifeq ($(USE_BLAS),YES)
  OBJ_BLAS := $(DIR_SELDON)/lib/Compil/Seldon/Blas$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o
endif
OBJ_SELDON := $(OBJ_SELDON) $(OBJ_BLAS)
OBJ_BLAS := $(DIR_SELDON)/lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Blas$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o

# MPI functions
ifeq ($(USE_MPI),YES)
  OBJ_SELDON := $(OBJ_SELDON) $(DIR_SELDON)/lib/Compil/Seldon/DistributedVector$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/DistributedMatrix$(SELD_SUFFIX).o
endif

# interface with eigenvalue solvers
# Arpack is always needed in that case
ifeq ($(USE_BLAS),YES)
OBJ_EIG := $(DIR_SELDON)/lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o
ifeq ($(USE_ARPACK),YES)
  OBJ_EIG := $(OBJ_EIG) $(DIR_SELDON)/lib/Compil/Seldon/Arpack$(SELD_SUFFIX).o
  ifeq ($(USE_FEAST),YES)
    OBJ_EIG := $(OBJ_EIG) $(DIR_SELDON)/lib/Compil/Seldon/Feast$(SELD_SUFFIX).o
  endif
  ifeq ($(USE_ANASAZI),YES)
    OBJ_EIG := $(OBJ_EIG) $(DIR_SELDON)/lib/Compil/Seldon/Anasazi$(SELD_SUFFIX).o
  endif
  ifeq ($(USE_SLEPC),YES)
    OBJ_EIG := $(OBJ_EIG) $(DIR_SELDON)/lib/Compil/Seldon/Slepc$(SELD_SUFFIX).o
  endif
endif
endif

OBJ_SELDON := $(OBJ_SELDON) $(OBJ_EIG)

# interfaces with direct solvers
ifeq ($(USE_BLAS),YES)
OBJ_SELDON_SOLVE := $(DIR_SELDON)/lib/Compil/Seldon/CholeskySolver$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/SparseSeldonSolver$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/SparseDirectSolver$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/DistributedSolver$(SELD_SUFFIX).o
else
OBJ_SELDON_SOLVE :=
endif

ifeq ($(USE_MUMPS),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o
endif

ifeq ($(USE_PARDISO),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/Pardiso$(SELD_SUFFIX).o
endif

ifeq ($(USE_PASTIX),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o
endif

ifeq ($(USE_PASTIX),INT64)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o
endif

ifeq ($(USE_UMFPACK),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/UmfPack$(SELD_SUFFIX).o $(DIR_SELDON)/lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o
endif

ifeq ($(USE_SUPERLU),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/SuperLU$(SELD_SUFFIX).o
endif

ifeq ($(USE_CHOLMOD),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o
endif

ifeq ($(USE_WSMP),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(DIR_SELDON)/lib/Compil/Seldon/Wsmp$(SELD_SUFFIX).o
endif

OBJ_SELDON := $(OBJ_SELDON) $(OBJ_SELDON_SOLVE)

# list of dependances for each file .cpp contained in $(DIR_SELDON)/lib/Compil/Seldon
$(DIR_SELDON)/lib/Compil/Seldon/Array3D$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Array3D.cpp $(DIR_SELDON)/array/Array3D.cxx $(DIR_SELDON)/array/Array4D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/BandMatrix.cpp $(DIR_SELDON)/matrix_sparse/BandMatrix.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx $(DIR_SELDON)/computation/basic_functions/Functions_MatVect.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Common.cpp $(DIR_SELDON)/share/Errors.cxx $(DIR_SELDON)/share/MatrixFlag.cxx $(DIR_SELDON)/share/Common.cxx $(DIR_SELDON)/share/Allocator.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/DistributedVector$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/DistributedVector.cpp $(DIR_SELDON)/share/MpiCommunication.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/DistributedMatrix$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/DistributedMatrix.cpp $(DIR_SELDON)/matrix_sparse/DistributedMatrix.cxx $(DIR_SELDON)/matrix_sparse/DistributedMatrixFunction.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixPointers$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixPointers.cpp $(DIR_SELDON)/matrix/Matrix_Pointers.cxx $(DIR_SELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixPacked$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixPacked.cpp $(DIR_SELDON)/matrix/Matrix_HermPacked.cxx $(DIR_SELDON)/matrix/Matrix_SymPacked.cxx $(DIR_SELDON)/matrix/Matrix_TriangPacked.cxx $(DIR_SELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixHermSymTriang$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixHermSymTriang.cpp $(DIR_SELDON)/matrix/Matrix_Hermitian.cxx $(DIR_SELDON)/matrix/Matrix_Symmetric.cxx $(DIR_SELDON)/matrix/Matrix_Triangular.cxx $(DIR_SELDON)/matrix/Matrix_Base.cxx $(DIR_SELDON)/computation/basic_functions/Functions_MatVect.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixSparse$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixSparse.cpp $(DIR_SELDON)/matrix_sparse/Matrix_Sparse.cxx $(DIR_SELDON)/matrix_sparse/Matrix_SymSparse.cxx $(DIR_SELDON)/matrix_sparse/Matrix_ArraySparse.cxx $(DIR_SELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Vector$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Vector.cpp $(DIR_SELDON)/vector/Vector.cxx $(DIR_SELDON)/vector/SparseVector.cxx $(DIR_SELDON)/vector/Functions_Arrays.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Vector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Vector2$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Vector2.cpp $(DIR_SELDON)/vector/Vector2.cxx $(DIR_SELDON)/vector/Vector3.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixVector$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixVector.cpp $(DIR_SELDON)/computation/basic_functions/Functions_MatVect.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx $(DIR_SELDON)/computation/interfaces/Mkl_Sparse.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrix$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrix.cpp $(DIR_SELDON)/computation/basic_functions/Functions_Matrix.cxx $(DIR_SELDON)/computation/interfaces/Mkl_Sparse.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixDense$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixDense.cpp $(DIR_SELDON)/matrix/Functions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixArray$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixArray.cpp $(DIR_SELDON)/matrix_sparse/Functions_MatrixArray.cxx $(DIR_SELDON)/computation/basic_functions/Functions_MatVect.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatVectComplex$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatVectComplex.cpp $(DIR_SELDON)/matrix_sparse/complex/Functions_MatVectComplex.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixComplex$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/FunctionsMatrixComplex.cpp $(DIR_SELDON)/matrix_sparse/complex/Functions_MatrixComplex.cxx $(DIR_SELDON)/matrix_sparse/Functions_MatrixArray.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixComplexConversions$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixComplexConversions.cpp $(DIR_SELDON)/matrix_sparse/complex/Matrix_ComplexConversions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixComplexSparse$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixComplexSparse.cpp $(DIR_SELDON)/matrix_sparse/complex/Matrix_ArrayComplexSparse.cxx $(DIR_SELDON)/matrix_sparse/complex/Matrix_ComplexSparse.cxx $(DIR_SELDON)/matrix_sparse/complex/Matrix_SymComplexSparse.cxx $(DIR_SELDON)/matrix_sparse/Matrix_Sparse.cxx $(DIR_SELDON)/matrix_sparse/IOMatrixMarket.cxx $(DIR_SELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/MatrixConversion$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/MatrixConversion.cpp $(DIR_SELDON)/matrix_sparse/Matrix_Conversions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/PermutationScalingMatrix$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/PermutationScalingMatrix.cpp $(DIR_SELDON)/matrix_sparse/Permutation_ScalingMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/RelaxationMatrixVector$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/RelaxationMatrixVector.cpp $(DIR_SELDON)/matrix_sparse/Relaxation_MatVect.cxx $(DIR_SELDON)/matrix_sparse/complex/Functions_MatVectComplex.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Blas$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Blas.cpp $(DIR_SELDON)/computation/interfaces/Blas_1.cxx $(DIR_SELDON)/computation/interfaces/Blas_2.cxx $(DIR_SELDON)/computation/interfaces/Blas_3.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Vector.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/BlasMpfr.cpp $(DIR_SELDON)/computation/basic_functions/Functions_Vector.cxx $(DIR_SELDON)/computation/basic_functions/Functions_MatVect.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Matrix.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx 
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/CholeskySolver$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/CholeskySolver.cpp $(DIR_SELDON)/computation/solver/SparseCholeskyFactorisation.cxx $(DIR_SELDON)/computation/solver/DistributedCholeskySolver.cxx 
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/EigenvalueSolver.cpp $(DIR_SELDON)/computation/interfaces/eigenvalue/VirtualEigenvalueSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Arpack$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Arpack.cpp $(DIR_SELDON)/computation/interfaces/eigenvalue/ArpackSolver.cxx $(DIR_SELDON)/computation/interfaces/eigenvalue/Arpack.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Feast$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Feast.cpp $(DIR_SELDON)/computation/interfaces/eigenvalue/Feast.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Anasazi$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Anasazi.cpp $(DIR_SELDON)/computation/interfaces/eigenvalue/Anasazi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Slepc$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Slepc.cpp $(DIR_SELDON)/computation/interfaces/eigenvalue/Slepc.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/IlutPreconditioning$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/IlutPreconditioning.cpp $(DIR_SELDON)/computation/solver/preconditioner/IlutPreconditioning.cxx $(DIR_SELDON)/computation/solver/preconditioner/SymmetricIlutPreconditioning.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/IOMatrixMarket$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/IOMatrixMarket.cpp $(DIR_SELDON)/matrix_sparse/IOMatrixMarket.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/IterativeSolver$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/IterativeSolver.cpp $(DIR_SELDON)/computation/solver/iterative/Iterative.cxx $(DIR_SELDON)/computation/solver/preconditioner/Precond_Ssor.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Lapack.cpp $(DIR_SELDON)/computation/interfaces/Lapack_LinearEquations.cxx $(DIR_SELDON)/computation/interfaces/Lapack_LeastSquares.cxx $(DIR_SELDON)/computation/interfaces/Lapack_Eigenvalues.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/LapackMpfr.cpp $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx src/Algebra/Eigenvalue.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Cholmod.cpp $(DIR_SELDON)/computation/interfaces/direct/Cholmod.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Mumps.cpp $(DIR_SELDON)/computation/interfaces/direct/Mumps.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Pastix.cpp $(DIR_SELDON)/computation/interfaces/direct/Pastix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Pardiso$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Pardiso.cpp $(DIR_SELDON)/computation/interfaces/direct/Pardiso.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/UmfPack$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/UmfPack.cpp $(DIR_SELDON)/computation/interfaces/direct/UmfPack.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/SuperLU$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/SuperLU.cpp $(DIR_SELDON)/computation/interfaces/direct/SuperLU.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/Wsmp$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/Wsmp.cpp $(DIR_SELDON)/computation/interfaces/direct/Wsmp.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/SparseSeldonSolver$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/SparseSeldonSolver.cpp $(DIR_SELDON)/computation/solver/Ordering.cxx $(DIR_SELDON)/computation/solver/SparseSolver.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/SparseDirectSolver$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/SparseDirectSolver.cpp $(DIR_SELDON)/computation/interfaces/direct/SparseDirectSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/DistributedSolver$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/DistributedSolver.cpp $(DIR_SELDON)/computation/solver/DistributedSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/check_dim$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/check_dim.cpp $(DIR_SELDON)/computation/basic_functions/Functions_Vector.cxx $(DIR_SELDON)/computation/basic_functions/Functions_MatVect.cxx $(DIR_SELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(DIR_SELDON)/lib/Compil/Seldon/TinyVector$(SELD_SUFFIX).o : $(DIR_SELDON)/lib/Compil/Seldon/TinyVector.cpp $(DIR_SELDON)/vector/TinyVector.cxx $(DIR_SELDON)/matrix/TinyMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)


$(LIB_SELDON_STATIC) :  $(OBJ_SELDON)
	ar rv $(LIB_SELDON_STATIC) $(OBJ_SELDON)

$(LIB_SELDON) :  $(OBJ_SELDON)
	$(CC) -shared -Wl,-soname,libseldon$(SELD_SUFFIX).so -rdynamic -o $(LIB_SELDON) $(OBJ_SELDON)

# target to remove all the object files
cleanlib :
	rm -f $(LIB_SELDON) $(LIB_SELDON_STATIC) $(OBJ_SELDON)

# target to remove files related to the interface with the direct solvers
cleansolve :
	rm -f $(OBJ_SELDON_SOLVE)

# target to remove files related to the interface with the eigenvalue solvers
cleaneig :
	rm -f $(OBJ_EIG) $(DIR_SELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o

# target to switch to SELDON_WITH_ABORT
cleanabort :
	rm -f $(DIR_SELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o

# Uncomment the following line in order to detect the file .o which gives bad result
#LIB_SELDON := $(OBJ_SELDON)

ifeq ($(STATIC_COMPILATION),YES)
  LIB_SELDON := $(LIB_SELDON_STATIC)
endif

ifeq ($(SEPARED_COMPIL),YES)
  LIB := $(LIB_SELDON) $(LIB)
endif
