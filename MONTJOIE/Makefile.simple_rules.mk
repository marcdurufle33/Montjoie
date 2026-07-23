# Makefile rules
#
# This file contains rules for compiling various Montjoie programs
#
# To be able to use this Makefile rules with your own Makefile, you have to :
# a) set up your compilation variables
# -> set your c++ compilator into the CC variable
# -> set your fortran 90 compilator into the FCC variable
# -> set your total include path into the INCLUDE variable
# -> set your total library path into the LIB variable
# -> set your compilator flags into the FLAGS variable
# -> set your compilator optimisation into the AGRESSIVE_OPTIM variable
# b) set up your defaut compiled program (optional)
# -> set your default binary filename into the BIN variable
# -> set yout default source filename into the LSTBIN variable
#
# Author : SEMIN Adrien
# Author of the initial Makefile (under AFS filesystem) : DURUFLE Marc

ifeq ($(USE_BLAS),YES)
  BLAS_MESSAGE="Utilisation de Blas"
else
  BLAS_MESSAGE="Blas non utilise"
endif

ifeq ($(USE_ARPACK),YES)
  ARPACK_MESSAGE="Utilisation de Arpack"
else
  ARPACK_MESSAGE="Arpack non utilise"
endif

ifeq ($(USE_MUMPS),YES)
  MUMPS_MESSAGE="Utilisation de Mumps"
else
  MUMPS_MESSAGE="Mumps non utilise"
endif

ifeq ($(USE_PASTIX),YES)
  PASTIX_MESSAGE="Utilisation de Pastix"
else
  PASTIX_MESSAGE="Pastix non utilise"
endif

ifeq ($(USE_MPFR),YES)
  MPFR_MESSAGE="Utilisation de Mpfr"
else
  MPFR_MESSAGE="Mpfr non utilise"
endif

ifeq ($(OPTIMIZATION),YES)
  OPTIM_MESSAGE="Compilation en mode optimise"
else
  OPTIM_MESSAGE="Compilation en mode debogage"
endif

ifeq ($(USE_MPI),YES)
  MPI_MESSAGE="Compilation en mode parallele"
else
  MPI_MESSAGE="Compilation en mode sequentiel"
endif

ifeq ($(USE_VERBOSE),YES)
  VERBOSE=
else
  VERBOSE=@
endif

# Default executable
regular: lib $(BIN)

ifeq ($(STATIC_COMPILATION),NO)
  LIB_BESSEL := lib/libbessel.so
else
  LIB_BESSEL := lib/libbessel.a
endif

# Semi-external libraries
lib: messages $(LIB_BESSEL) lib/etime.o $(MPFR_LIB)

messages :
	@echo $(BLAS_MESSAGE)
	@echo $(ARPACK_MESSAGE)
	@echo $(MUMPS_MESSAGE)
	@echo $(PASTIX_MESSAGE)
	@echo $(MPFR_MESSAGE)
	@echo
	@echo $(MPI_MESSAGE)
	@echo $(OPTIM_MESSAGE)
	@echo

lib/Bessel.o : lib/Bessel.f
	$(VERBOSE)$(FCC) -c $< -o $@ -O3 -fPIC

lib/d1mach.o : lib/d1mach.f
	$(VERBOSE)$(FCC) -c $< -o $@ -O3 -fPIC

lib/i1mach.o : lib/i1mach.f
	$(VERBOSE)$(FCC) -c $< -o $@ -O3 -fPIC

lib/etime.o : lib/etime.f
	$(VERBOSE)$(FCC) -c $< -o $@ -O3 -fPIC

lib/mpreal.o : src/Share/mpreal.cpp
	$(VERBOSE)$(CC) -c -Ilib $< -o $@ -O3 -fPIC -std=c++11

include Makefile.lib_simple.mk

# Default binary
$(BIN) : lib $(LSTBIN) $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	@echo "a partir de $(LSTBIN)"
	$(VERBOSE)$(CC) $(FLAGS) $(LSTBIN) -c -o $(BIN).o $(INCLUDE) \
	  $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) $(BIN).o -o $(BIN) $(AGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv $(BIN) .
	@echo "Regle $@ compilee avec succes."

# clean .o files in src directory
clean : 
	@rm -f src/*.o

cleanall : 
	@rm -f src/*.o lib/*.o lib/*.so lib/*.a lib/*/*/*.o lib/*/*/*/*.o

.PHONY : all lib clean

	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

.PHONY : all lib clean
