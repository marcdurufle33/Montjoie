#!/usr/bin/env python

import sys, pwd, os

import string, subprocess

liste_fichier = ["src/Program/Unit/Seldon/array3d_test.cc", \
                 "src/Program/Unit/Seldon/band_test.cc", \
                 "src/Program/Unit/Seldon/cholesky_test.cc", \
                 "src/Program/Unit/Seldon/conversion_matrix.cc", \
                 "src/Program/Unit/Seldon/direct_test.cc", \
                 "src/Program/Unit/Seldon/precond_test.cc", \
                 "src/Program/Unit/Seldon/eigenvalue_solver_test.cc", \
                 "src/Program/Unit/Seldon/eigenvalue_test.cc", \
                 "src/Program/Unit/Seldon/eigenvalue_sparse.cc", \
                 "src/Program/Unit/Seldon/least_squares_test.cc", \
                 "src/Program/Unit/Seldon/lu_test.cc", \
                 "src/Program/Unit/Seldon/function_matrix_complex_test.cc", \
                 "src/Program/Unit/Seldon/function_matrix_test.cc", \
                 "src/Program/Unit/Seldon/function_matvect_test.cc", \
                 "src/Program/Unit/Seldon/function_vector.cc", \
                 "src/Program/Unit/Seldon/matrix_complex_sparse_test.cc", \
                 "src/Program/Unit/Seldon/matrix_sparse_test.cc", \
                 "src/Program/Unit/Seldon/matrix_test.cc", \
                 "src/Program/Unit/Seldon/matrix_product_test.cc", \
                 "src/Program/Unit/Seldon/relaxation_test.cc", \
                 "src/Program/Unit/Seldon/vector_test.cc", \
                 "src/Program/Unit/Seldon/vector2_test.cc", \
                 "src/Program/Unit/Algebra/diagonal_matrix.cc", \
                 "src/Program/Unit/Algebra/function_matrix_extraction_test.cc", \
                 "src/Program/Unit/Algebra/MatrixBlockDiagTest.cc", \
                 "src/Program/Unit/Algebra/skyline_matrix_test.cc", \
                 "src/Program/Unit/Algebra/tiny_vector_test.cc", \
                 "src/Program/Unit/Algebra/tiny_matrix_test.cc", \
                 "src/Program/Unit/Algebra/tiny_array3d_test.cc", \
                 "src/Program/Unit/Algebra/tiny_band_test.cc", \
                 "src/Program/Unit/Algebra/tiny_block_solver.cc", \
                 "src/Program/Unit/Share/fft_test.cc", \
                 "src/Program/Unit/Hyperbolic/test_finite_difference_1d.cc", \
                 "src/Program/Unit/Instationary/explicit_scheme.cc", \
                 "src/Program/Unit/Instationary/implicit_scheme.cc", \
                 "src/Program/Unit/Instationary/mass_matrix_test.cc", \
                 "src/Program/Unit/physical_indices_test.cc", \
                 "src/Program/Unit/Share/random_test.cc" ]

liste_fichier_par = ["src/Program/Unit/Algebra/distributed_matrix_test.cc", \
                     "src/Program/Unit/Algebra/distributed_block_diagonal.cc", \
                     "src/Program/Unit/Algebra/distributed_cholesky.cc"]

# repertoire courrant
current_dir = os.getcwd()
    
def SetSimulation(index):
    fid = open("../Makefile", "r")
    lignes = fid.readlines()
    fid.close()

    fid = open("../UnitMakefile"+str(index), "w")
    Nseq = len(liste_fichier)
    for ligne in lignes:
        if (len(lignes) > 1):
            if (ligne[0] != '#'):
                if (ligne.startswith("LSTBIN ")):
                    if (index >= Nseq):
                        fid.write("LSTBIN = "+liste_fichier_par[index-Nseq]+ '\n')
                    else:
                        fid.write("LSTBIN = "+liste_fichier[index]+ '\n')
                elif (ligne.startswith("BIN ")):
                    fid.write("BIN = src/test"+str(index)+".x\n")
                else:
                    fid.write(ligne)
                    #fid.write('\n')
        else:
            fid.write('\n')
    
    fid.close()
    #print "Compiling ", nom, " ..."
    command = "cd " + current_dir + "/.. && "
    command += "make -f UnitMakefile" + str(index) + " >& sortMake" + str(index)
    if (index >= Nseq):
            command += " && mpirun -np 4 ./test" + str(index) + ".x >& sortie" + str(index) \
                + " && rm UnitMakefile"+str(index)+ " test" + str(index) + ".x"
    else:
        command += " && ./test" + str(index) + ".x >& sortie" + str(index) \
            + " && rm UnitMakefile"+str(index)+ " test" + str(index) + ".x"
    return command


# on ignore les premieres simulations si demande
debut = 0
if (len(sys.argv) >= 2):
    debut = int(sys.argv[1])

success_test = True
os.environ['SHELL'] = 'bash'
fid = open("/home/durufle/log_compil_test.txt", "w")
for run_index in range(debut, len(liste_fichier)):
    print "Running test #" + str(debut + run_index)
    
    # Getting parameters of the simulation
    command_line = SetSimulation(run_index)
    
    # Launches the simulation.
    #sortie = os.system(command_line)    
    p = subprocess.Popen(command_line, shell=True, executable="/bin/bash")
    sortie = p.wait()
    if (sortie != 0):
        print "Test ", run_index, " : ", liste_fichier[run_index], " failed "

if (success_test):
    print "All tests passed successfully"
    fid.write("All tests passed successfully\n")
else:
    sys.exit(-1)
