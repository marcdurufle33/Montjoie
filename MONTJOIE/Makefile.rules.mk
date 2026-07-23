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

ifeq ($(PLAFRIM1),YES)
  PLAF_MESSAGE="Nous sommes sur Plafrim 1"
else
  PLAF_MESSAGE="Nous sommes sur Plafrim 2"
endif

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

ifeq ($(OPTIMIZATION),FASTER)
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
	@echo $(PLAF_MESSAGE)
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

include Makefile.lib.mk

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

# targets for non-regression tests
all : vibro fluid_structure laplace helm2D helm_div helm3D helm_radial maxwell_axi maxwell2D static_maxwell maxwell3D mode_maxwell time_maxwell time3D acous2D acous3D harmonic_elas2D harmonic_elas3D elas2D elas3D reissner static_elastic soundboard schrodinger1D aero2D aero3D galbrun harmonic_aero2D acous_axi aero_axi galbrun_axi solution_helm solution solution_disc CellulairePeriodique CellulairePeriodique3D kdv camassa helm_axi write_nodal piano

piano : multistringStiffNL PianoStiffNL StringImpedanceStiffLIN StringImpedanceStiffNL StringSoundboardStiffNL 

mesh : convert manipule manipule2D sym_mesh per_mesh tetmesh

# List of all binaries

####################
# Various Problems #
####################

# Lung problem (Adrien)
lung : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/lung.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# vibro acoustic (Juliette)
vibro : lib  $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/vibro.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Viscous thermal equation (Marie)
visco : lib  $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/visco.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Viscous cell problem for homogenized viscothermal acoustics (Alexis)
viscous_cell : lib  $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/viscous_cell.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# modeEz
modeEz : lib  $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/mode_Ez.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# modeEz
modeEs : lib  $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Augustin/mode_Es.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Cross Correlation post-computing (Adrien Semin)
crosscorrelation : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/CrossCorrelation.cc \
	  -c -o src/$@.o -fopenmp $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) -fopenmp 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Cross Correlation post-computing - compute directly picture (Adrien Semin)
crosspicture : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/CrossPicture.cc \
	  -c -o src/$@.o -fopenmp $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) -fopenmp 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Cross Correlation post-computing - compute directly picture (Adrien Semin)
crosspicturesecond: lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/CrossPictureSecond.cc \
	  -c -o src/$@.o -fopenmp $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) -fopenmp 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Map of speed post-computing (Adrien Semin)
mapspeed : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/MapSpeed.cc \
	  -c -o src/$@.o -fopenmp $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) -fopenmp 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Append files
appendfiles : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/AppendFiles.cc \
	  -c -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Differentiate files
differentiatefiles : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/DifferentiateFiles.cc \
	  -c -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."


# coupling between acoustics and elastics (fluid-structure interaction)
fluid_structure : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/fluid_structure.cc -c -o src/$@.o \
	  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# writing nodal points of a data file
write_nodal : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/write_nodal_points.cc \
	  -c -o src/$@.o -fopenmp $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/write_nodal.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) -fopenmp 
	$(VERBOSE)mv src/write_nodal.x .
	@echo "Regle $@ compilee avec succes."


#################
# Piano/Strings #
#################


# MultiString (version objet)
multistringStiffLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING src/Program/Corde/test_multistring.cc \
	  -c -o src/multistringStiffLIN.o  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/multistringStiffLIN.o -o src/multistringStiffLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/multistringStiffLIN.x .
	@echo "Regle $@ compilee avec succes."


# MultiString (version objet)
multistringLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING src/Program/Corde/test_multistring.cc \
	  -c -o src/multistringLIN.o  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/multistringLIN.o -o src/multistringLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/multistringLIN.x .
	@echo "Regle $@ compilee avec succes."


# MultiString (version objet)
multistringLINTL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING_TL src/Program/Corde/test_multistring.cc \
	  -c -o src/$@   $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) src/$@ -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# MultiString (version objet)
multistringNL2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING_2T src/Program/Corde/test_multistring.cc \
	  -c -o src/$@   $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@ -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# MultiString (version objet)
multistringStiffNL2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING_2T src/Program/Corde/test_multistring.cc \
	  -c -o src/$@   $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@ -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."	
# MultiString (version objet)
multistringStiffLIN2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING_2T src/Program/Corde/test_multistring.cc \
	  -c -o src/$@   $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@ -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."	
# MultiString (version objet)
multistringLIN2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING_2T src/Program/Corde/test_multistring.cc \
	  -c -o src/$@   $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@ -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."		

# MultiString (version objet)
multistringStiff : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_STIFF_STRING src/Program/Corde/test_multistring.cc \
	  -c -o src/multistringStiff.o  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/multistringStiff.o -o src/multistringStiff.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/multistringStiff.x .
	@echo "Regle $@ compilee avec succes."


# MultiString (version objet)
multistringStiffNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING src/Program/Corde/test_multistring.cc \
	  -c -o src/multistringStiffNL.o  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/multistringStiffNL.o -o src/multistringStiffNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/multistringStiffNL.x .
	@echo "Regle $@ compilee avec succes."

# MultiString (version objet)
multistringNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING src/Program/Corde/test_multistring.cc \
	  -c -o src/multistringNL.o  $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/multistringNL.o -o src/multistringNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/multistringNL.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plate-fluid
PianoFluide : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Corde/PianoFluide.cc \
	  -c -o src/PianoFluide.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/PianoFluide.o -o src/PianoFluide.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoFluide.x .
	@echo "Regle $@ compilee avec succes."

# Couplage piano total 
PianoStiffNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING src/Program/Corde/CouplagePiano.cc \
	  -c -o src/CouplagePianoStiffNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CouplagePianoStiffNL.o -o src/PianoStiffNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoStiffNL.x PianoStiffNL.x
	@echo "Regle $@ compilee avec succes."

# Couplage piano total 
PianoStiffNL2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING_2T src/Program/Corde/CouplagePiano.cc \
	  -c -o src/CouplagePianoStiffNL2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CouplagePianoStiffNL2T.o -o src/PianoStiffNL2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoStiffNL2T.x PianoStiffNL2T.x
	@echo "Regle $@ compilee avec succes."

# Couplage piano total
PianoStiffLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING src/Program/Corde/CouplagePiano.cc \
	  -c -o src/CouplagePianoStiffLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CouplagePianoStiffLIN.o -o src/PianoStiffLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoStiffLIN.x .
	@echo "Regle $@ compilee avec succes."


# Couplage piano total
PianoStiff : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_STIFF_STRING src/Program/Corde/CouplagePiano.cc \
	  -c -o src/CouplagePiano.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CouplagePiano.o -o src/PianoStiff.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoStiff.x .
	@echo "Regle $@ compilee avec succes."

# Couplage piano total
PianoNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING src/Program/Corde/CouplagePiano.cc \
	  -c -o src/CouplagePianoNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CouplagePianoNL.o -o src/PianoNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoNL.x .
	@echo "Regle $@ compilee avec succes."

# Couplage piano total
PianoLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING src/Program/Corde/CouplagePiano.cc \
	  -c -o src/CouplagePianoLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CouplagePianoLIN.o -o src/PianoLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/PianoLIN.x .
	@echo "Regle $@ compilee avec succes."



# Test couplage plaque corde avec impedance
StringImpedanceStiffNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING src/Program/Corde/CordeImpedance.cc \
	  -c -o src/CordeImpedanceStiffNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeImpedanceStiffNL.o -o src/StringImpedanceStiffNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/StringImpedanceStiffNL.x .
	@echo "Regle $@ compilee avec succes."

# Test couplage plaque corde avec impedance Linearise
StringImpedanceStiffLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING src/Program/Corde/CordeImpedance.cc \
	  -c -o src/CordeImpedanceStiffLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeImpedanceStiffLIN.o -o src/StringImpedanceStiffLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/StringImpedanceStiffLIN.x .
	@echo "Regle $@ compilee avec succes."

# Test couplage plaque corde avec impedance Linearise
StringImpedanceLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING src/Program/Corde/CordeImpedance.cc \
	  -c -o src/CordeImpedanceLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeImpedanceLIN.o -o src/StringImpedanceLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/StringImpedanceLIN.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringImpedanceLIN2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING_2T src/Program/Corde/CordeImpedance.cc \
	  -c -o src/CordeImpedanceLIN2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeImpedanceLIN2T.o -o src/StringImpedanceLIN2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringImpedanceLIN2T.x .
	@echo "Regle $@ compilee avec succes."	

# Test couplage plaque corde avec impedance NL
StringImpedanceNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING src/Program/Corde/CordeImpedance.cc \
	  -c -o src/CordeImpedanceNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeImpedanceNL.o -o src/StringImpedanceNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)

	$(VERBOSE)mv src/StringImpedanceNL.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringSoundboardStiffLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardStiffLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardStiffLIN.o -o src/StringSoundboardStiffLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardStiffLIN.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringSoundboardLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardLIN.o -o src/StringSoundboardLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardLIN.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringSoundboardStiff : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_STIFF_STRING src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardStiff.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardStiff.o -o src/StringSoundboardStiff.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardStiff.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringSoundboardStiffNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardStiffNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardStiffNL.o -o src/StringSoundboardStiffNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardStiffNL.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringSoundboardStiffNL2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING_2T src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardStiffNL2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardStiffNL2T.o -o src/StringSoundboardStiffNL2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardStiffNL2T.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringSoundboardStiffLIN2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING_2T src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardStiffLIN2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardStiffLIN2T.o -o src/StringSoundboardStiffLIN2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardStiffLIN2T.x .
	@echo "Regle $@ compilee avec succes."

StringSoundboardNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING src/Program/Corde/CordeSoundboard.cc \
	  -c -o src/CordeSoundboardNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeSoundboardNL.o -o src/StringSoundboardNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringSoundboardNL.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeStiffLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STIFF_STRING src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeStiffLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeStiffLIN.o -o src/StringBridgeStiffLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeStiffLIN.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeLIN : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeLIN.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeLIN.o -o src/StringBridgeLIN.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeLIN.x .
	@echo "Regle $@ compilee avec succes."


# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeLINTL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING_TL src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeLINTL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeLINTL.o -o src/StringBridgeLINTL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeLINTL.x .
	@echo "Regle $@ compilee avec succes."


# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeLIN2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_LINEAR_STRING_2T src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeLIN2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeLIN2T.o -o src/StringBridgeLIN2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeLIN2T.x .
	@echo "Regle $@ compilee avec succes."	

# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeNL.o -o src/StringBridgeNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeNL.x .
	@echo "Regle $@ compilee avec succes."	

# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeNL2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STRING_2T src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeNL2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeNL2T.o -o src/StringBridgeNL2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeNL2T.x .
	@echo "Regle $@ compilee avec succes."	

# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeStiffNL : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeStiffNL.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeStiffNL.o -o src/StringBridgeStiffNL.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeStiffNL.x .
	@echo "Regle $@ compilee avec succes."

# Couplage plaque corde avec vraie plaque diagonalisee
StringBridgeStiffNL2T : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) -DPIANO_NONLINEAR_STIFF_STRING_2T src/Program/Corde/CordeBridge.cc \
	  -c -o src/CordeBridgeStiffNL2T.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/CordeBridgeStiffNL2T.o -o src/StringBridgeStiffNL2T.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/StringBridgeStiffNL2T.x .
	@echo "Regle $@ compilee avec succes."	


# Test de la resolution exacte de la plaque
Test_Plaque : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Corde/test_soundboard.cc \
	  -c -o src/Test_Soundboard.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Test_Soundboard.o \
	  -o src/Test_Soundboard.x $(FLAGS) $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/Test_Soundboard.x .
	@echo "Regle $@ compilee avec succes."


#########################
# Targets for Beginners #
#########################

# Problems in Beginner directory with only one or two instanciations of EllipticProblem/HyperbolicProblem
# Usually the instanciation may vary in these files since they are often used, so check in these files
# to see which finite element is used

# Laplace equation (Poisson)
claplace : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/laplace.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/laplace.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/laplace.x .
	@echo "Regle $@ compilee avec succes."

# 1-D Helmholtz equation
helm1D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/helmholtz1D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz1D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz1D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D Helmholtz equation
chelm2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/helmholtz2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D Helmholtz equation
chelm3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/helmholtz3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D static Maxwell equation
cstatic_maxwell : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/static_maxwell.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/static_maxwell.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/static_maxwell.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-harmonic Maxwell equation
cmaxwell3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/maxwell3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/maxwell3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/maxwell3D.x .
	@echo "Regle $@ compilee avec succes."

# 1-D non-linear Maxwell's equations
maxwell1D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/maxwell1D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/maxwell1D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/maxwell1D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-domain Maxwell equation
ctime3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/maxwell3D_time.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/time3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/time3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D acoustic equation
cacous2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/acous2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/acous2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/acous2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D acoustic equation
cacous3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/acous3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/acous3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/acous3D.x .
	@echo "Regle $@ compilee avec succes."

# advection equation
advection : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/advection.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/advec.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/advec.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-harmonic elastic equation
charmonic_elas2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/harmonic_elastic2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/harmonic_elas2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/harmonic_elas2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-harmonic elastic equation
charmonic_elas3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/harmonic_elastic3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/harmonic_elas3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/harmonic_elas3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-domain elastic equation
celas2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/elastic2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/elas2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/elas2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-domain elastic equation
celas3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/elastic3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/elas3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/elas3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-domain aero-acoustic equation
caero2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/aero2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/aero2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/aero2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-domain aero-acoustic equation
caero3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/aero3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/aero3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/aero3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-harmonic aero-acoustic equation
charmonic_aero2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/harmonic_aero2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/harmonic_aero2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/harmonic_aero2D.x .
	@echo "Regle $@ compilee avec succes."

# Difference between two matrices
diff_matrix : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/diff_matrix.cc \
	  -c -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# eigenvalues of finite element matrix
eigenval_matrix : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/eigenval_aero.cc \
	  -c -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB) 
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."


##############################
# Targets for Advanced users #
##############################

# Problems in directory Advanced with several instanciations of EllipticProblem/HyperbolicProblem

# Laplace equation (Poisson)
laplace : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/laplace.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/laplace.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/laplace.x .
	@echo "Regle $@ compilee avec succes."

# 2-D Helmholtz equation
helm2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/helmholtz2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz2D.x .
	@echo "Regle $@ compilee avec succes."

# Helmholtz equation solved with H(div) formulation
helm_div : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/helmholtz_sipg.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helm_hdiv.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helm_hdiv.x .
	@echo "Regle $@ compilee avec succes."

# 3-D Helmholtz equation
helm3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/helmholtz3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz3D.x .
	@echo "Regle $@ compilee avec succes."

# 2.5-D Helmholtz equation (axisymmetric domains)
helm_axi : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/helmholtz_axi.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz_axi.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz_axi.x .
	@echo "Regle $@ compilee avec succes."

# 1.5-D Helmholtz equation (radial domains)
helm_radial : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Beginner/helm_radial.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/helmholtz_radial.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/helmholtz_radial.x .
	@echo "Regle $@ compilee avec succes."

# 2.5-D time-harmonic Maxwell equation (axisymmetric domains)
maxwell_axi : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/maxwell_axi.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/maxwell_axi.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/maxwell_axi.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-harmonic Maxwell equation
maxwell2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/maxwell2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/maxwell2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/maxwell2D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D and 3-D static Maxwell equation
static_maxwell : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/static_maxwell.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/static_maxwell.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/static_maxwell.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-harmonic Maxwell equation
maxwell3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/maxwell3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/maxwell3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/maxwell3D.x .
	@echo "Regle $@ compilee avec succes."

# computation of modes propagating along guides for 3-D Maxwell equations
mode_maxwell : lib  $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/mode_maxwell.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/mode_maxwell.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/mode_maxwell.x .
	@echo "Regle $@ compilee avec succes."

# 2.5-D time-domain Maxwell equation (axisymmetric domains)
time_maxwell : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/time_maxwell.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/time_maxwell.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/time_maxwell.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-domain Maxwell equation
time3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/maxwell3D_time.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/time3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/time3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D acoustic equation
acous2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/acous2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/acous2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/acous2D.x .
	@echo "Regle $@ compilee avec succes."

# 2.5-D acoustic equation (axisymmetric domains)
acous_axi : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/acous_axi.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/acous_axi.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/acous_axi.x .
	@echo "Regle $@ compilee avec succes."

# 3-D acoustic equation
acous3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/acous3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/acous3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/acous3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-harmonic elastic equation
harmonic_elas2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/harmonic_elastic2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/harmonic_elas2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/harmonic_elas2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-harmonic elastic equation
harmonic_elas3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/harmonic_elastic3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/harmonic_elas3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/harmonic_elas3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-domain elastic equation
elas2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/elastic2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/elas2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/elas2D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-domain elastic equation
elas3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/elastic3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/elas3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/elas3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-domain Reissner-Mindlin equation
reissner : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/reissner_mindlin.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/reissner.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/reissner.x .
	@echo "Regle $@ compilee avec succes."

# Computation of eigenmodes for 3-D elastic equation
static_elastic : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/static_elastic3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/static_elastic.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/static_elastic.x .
	@echo "Regle $@ compilee avec succes."

# Computation of eigenmodes for Reissner-Mindlin equation 
soundboard : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/soundboard.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/soundboard.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/soundboard.x .
	@echo "Regle $@ compilee avec succes."

# 1-D non-linear schrodinger equation
schrodinger1D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/schrodinger1D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/schrodinger1D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/schrodinger1D.x .
	@echo "Regle $@ compilee avec succes."

# 3-D non-linear schrodinger equation
schrodinger3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/schrodinger3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/schrodinger3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/schrodinger3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-domain aero-acoustic equation
aero2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/aero2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/aero2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/aero2D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-domain and time-harmonic Galbrun equation
galbrun : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/galbrun.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/galbrun.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/galbrun.x .
	@echo "Regle $@ compilee avec succes."

# 3-D time-domain aero-acoustic equation
aero3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/aero3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/aero3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/aero3D.x .
	@echo "Regle $@ compilee avec succes."

# 2-D time-harmonic aero-acoustic equation
harmonic_aero2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/harmonic_aero2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/harmonic_aero2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/harmonic_aero2D.x .
	@echo "Regle $@ compilee avec succes."

# axisymmmetric time-harmonic aero-acoustic equation
aero_axi : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/aero_axi.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/aero_axi.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/aero_axi.x .
	@echo "Regle $@ compilee avec succes."

# axisymmmetric time-harmonic Galbrun equation
galbrun_axi : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/galbrun_axi.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/galbrun_axi.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/galbrun_axi.x .
	@echo "Regle $@ compilee avec succes."

# Vlasov-Maxwell equations
vlasov : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/vlasov.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/vlasov.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/vlasov.x .
	@echo "Regle $@ compilee avec succes."

# computation of analytical solutions for sperical layers and Helmholtz equation
solution_helm : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/SolutionHelmholtz.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)  
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/solution_helm.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/solution_helm.x .
	@echo "Regle $@ compilee avec succes."

# computation of analytical solutions for sperical layers and time-harmonic Maxwell equations
solution : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/SolutionSphere.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/solution.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/solution.x .
	@echo "Regle $@ compilee avec succes."

# analytical solution for circular layers and Helmholtz equation
solution_disc : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/SolutionDisc.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/solution_disc.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/solution_disc.x .
	@echo "Regle $@ compilee avec succes."


# migration
migration : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/migration.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/migration.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/migration.x .
	@echo "Regle $@ compilee avec succes."


############################
# Targets for Expert users #
############################

# Programs in directory Expert, they are not documented and not very easy to understand

# dispersion analysis in 2-D or 3-D
dispersion : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# computation of eigenmodes with PML layers (for laser cavities)
mode_pml : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/mode_pml.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/mode_pml.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/mode_pml.x .
	@echo "Regle $@ compilee avec succes."

# computation of eigenmodes with PML layers (for laser cavities)
mode1D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/mode1D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/mode1D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/mode1D.x .
	@echo "Regle $@ compilee avec succes."

# computation of eigenmodes with PML layers (for laser cavities)
mode2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/mode2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/mode2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/mode2D.x .
	@echo "Regle $@ compilee avec succes."

# piezoelectricity
piezo : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/piezo.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/piezo.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/piezo.x .
	@echo "Regle $@ compilee avec succes."

# boundary inverse problems
inv_harmonic : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/inverse_harmonic.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/inv_harmonic.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/inv_harmonic.x .
	@echo "Regle $@ compilee avec succes."

# boundary inverse problems (3-D)
inv_harmonic3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/inverse_harmonic3D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/inv_harmonic3D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/inv_harmonic3D.x .
	@echo "Regle $@ compilee avec succes."

# cellulairePeriodic (Berangere)
CellulairePeriodique : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# cellulairePeriodic3D (Berangere)
CellulairePeriodique3D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

######################
# Targets for Test   #
######################


# program for checking 3-D Hcurl finite element
checkHcurl : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/check_finite_elementHcurl.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/check.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/check.x .
	@echo "Regle $@ compilee avec succes."

# program for checking 3-D Hcurl finite element
checkHdiv : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/check_finite_elementHdiv.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/check.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/check.x .
	@echo "Regle $@ compilee avec succes."

# program for checking 2-D Hcurl finite element
checkHcurl2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/check_finite_elementEdge.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/check.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/check.x .
	@echo "Regle $@ compilee avec succes."

# program for checking H1 finite element
checkH1 : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/check_finite_element.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/check.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/check.x .
	@echo "Regle $@ compilee avec succes."

# KdV equation
kdv : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# Camassa-Holm equation
camassa : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Advanced/camassa_holm.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# subdomain
subdomain : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Test/iterative_subdomain.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# audio
audio : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/audiofils.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# multi_freq
multi_freq : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Expert/multi_freq.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."


######################
# Targets for Meshes #
######################

# conversion of meshes (2-D and 3-D)
convert : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/convert2mesh.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/convert_mesh.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/convert_mesh.x .
	@echo "Regle $@ compilee avec succes."

# manipulation of meshes (3-D only)
manipule : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/manipule.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/manipule.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/manipule.x .
	@echo "Regle $@ compilee avec succes."

# manipulation of meshes (2-D only)
manipule2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/manipule2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/manipule2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/manipule2D.x .
	@echo "Regle $@ compilee avec succes."

# creation of films (2-D only)
film2D : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/film2D.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/film2D.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/film2D.x .
	@echo "Regle $@ compilee avec succes."

# symetrization of meshes (2-D and 3-D)
sym_mesh : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# periodization of meshes (2-D and 3-D)
per_mesh : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/periodize_mesh.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# adds circular layers
add_circle : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/add_circular_layer.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# adds spherical layers
add_sphere : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/add_spherical_layer.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

# tetmesh tool for using ghs3D
tetmesh : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Mesh/$@.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/$@.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/$@.x .
	@echo "Regle $@ compilee avec succes."

#AJOUT NATHAN
poisson : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Nathan/poisson.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/poisson.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/poisson.x .
	@echo "Regle $@ compilee avec succes."


#AJOUT NATHAN
# Time-harmonic Goldstein equation
goldstein : lib $(LIB_SELDON) $(LIB_MONTJOIE) $(LIB_MONTJOIE_F)
	@echo "Compilation de la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/Program/Nathan/goldstein.cc -c \
	  -o src/$@.o $(INCLUDE) $(AGRESSIVE_OPTIM)
	@echo "Edition de liens pour la regle $@"
	$(VERBOSE)$(CC) $(FLAGS) src/$@.o -o src/goldstein.x $(FLAGS) \
	  $(AGGRESSIVE_OPTIM) $(MPFR_LIB) $(LIB)
	$(VERBOSE)mv src/goldstein.x .
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
