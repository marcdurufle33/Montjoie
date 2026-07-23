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

ifeq ($(OPTIMIZATION),FASTER)
  OPTIM_MESSAGE="Compilation en mode optimise"
else
  OPTIM_MESSAGE="Compilation en mode debogage"
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

include Makefile.libmesh.mk

# Default binary
$(BIN) : lib $(LSTBIN) $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	@echo "a partir de $(LSTBIN)"
	$(VERBOSE)$(CC) $(FLAGS) $(LSTBIN) -c -o $(BIN).o $(INCLUDE) \
	  $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) $(BIN).o -o $(BIN) $(AGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv $(BIN) .
	@echo "Regle $@ compilee avec succes."

mesh : convert manipule manipule2D sym_mesh per_mesh tetmesh

# List of all binaries

######################
# Targets for Meshes #
######################

# conversion of meshes (2-D and 3-D)
convert : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/convert2mesh.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/convert_mesh.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/convert_mesh.x .
	@echo "Regle $@ compilee avec succes."

# manipulation of meshes (3-D only)
manipule : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/manipule.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/manipule.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/manipule.x .
	@echo "Regle $@ compilee avec succes."

# manipulation of meshes (2-D only)
manipule2D : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/manipule2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/manipule2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/manipule2D.x .
	@echo "Regle $@ compilee avec succes."

# creation of films (2-D only)
film2D : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/film2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/film2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/film2D.x .
	@echo "Regle $@ compilee avec succes."

# symetrization of meshes (2-D and 3-D)
sym_mesh : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# periodization of meshes (2-D and 3-D)
per_mesh : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/periodize_mesh.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# adds circular layers
add_circle : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/add_circular_layer.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# adds spherical layers
add_sphere : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/add_spherical_layer.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# tetmesh tool for using ghs3D
tetmesh : lib $(LIB_SELDON) $(LIB_MONTJOIE)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# clean .o files in src directory
clean : 
	@rm -f src/*.o

cleanall : 
	@rm -f src/*.o lib/*.o lib/*.so lib/*.a lib/*/*/*.o lib/*/*/*/*.o

# clean Montjoie input / output files and .o files in src directory
cleanresult : 
	@rm -f *.dat *.txt *.ini
	@rm -f *.bb
	@rm -f *.mesh
	@rm -f mu[0-9]* rho[0-9]*

# clean Montjoie input / output files and binary files
mrproper : clean
	@rm -f *.x
	@rm -f *.ini
	@rm -f convert_mesh
	@rm -f tetmesh
	@rm -f lib/*.o
	@rm -f src/*.o
	@rm -f *~

.PHONY : all lib clean mrproper

	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

.PHONY : all lib clean mrproper
