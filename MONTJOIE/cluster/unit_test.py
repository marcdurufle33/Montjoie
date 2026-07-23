#!/usr/bin/env python

import sys, pwd, os

import string, network
import subprocess

# Network.
net = network.Network("all")

# Network status.
network_busy = False
quota_reached = False
all_launched = False
all_done = False

# Waiting time.
waiting_command = "sleep 5"


liste_fichier = ["src/Program/Unit/Seldon/array3d_test.cc", \
                 "src/Program/Unit/Seldon/band_test.cc", \
                 "src/Program/Unit/Seldon/cholesky_test.cc", \
                 "src/Program/Unit/Seldon/conversion_matrix.cc", \
                 "src/Program/Unit/Seldon/direct_test.cc", \
                 "src/Program/Unit/Seldon/precond_test.cc", \
                 "src/Program/Unit/Seldon/eigenvalue_solver_test.cc", \
                 "src/Program/Unit/Seldon/eigenvalue_test.cc", \
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
                 "src/Program/Unit/Algebra/function_matrix_extraction_test.cc", \
                 "src/Program/Unit/Algebra/MatrixBlockDiagTest.cc", \
                 "src/Program/Unit/Algebra/skyline_matrix_test.cc", \
                 "src/Program/Unit/Algebra/tiny_vector_test.cc", \
                 "src/Program/Unit/Algebra/tiny_matrix_test.cc", \
                 "src/Program/Unit/Algebra/tiny_array3d_test.cc", \
                 "src/Program/Unit/Algebra/tiny_band_test.cc", \
                 "src/Program/Unit/Share/fft_test.cc", \
                 "src/Program/Unit/Hyperbolic/test_finite_difference_1d.cc", \
                 "src/Program/Unit/Instationary/explicit_scheme.cc", \
                 "src/Program/Unit/Instationary/implicit_scheme.cc", \
                 "src/Program/Unit/Instationary/mass_matrix_test.cc", \
                 "src/Program/Unit/physical_indices_test.cc", \
                 "src/Program/Unit/random_test.cc" ]

liste_fichier_par = ["src/Program/Unit/Algebra/distributed_matrix_test.cc", \
                     "src/Program/Unit/Algebra/distributed_block_diagonal.cc" ]

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
    command = "cd " + current_dir + "/..; source ~/.bash_profile;"
    command += "make -f UnitMakefile" + str(index) + " >& sortMake" + str(index)
    if (index >= Nseq):
            command += " && mpirun -np 4 ./test" + str(index) + ".x >& sortie" + str(index) \
                + " && rm UnitMakefile"+str(index)+ " test" + str(index) + ".x"
    else:
        command += " && ./test" + str(index) + ".x >& sortie" + str(index) \
            + " && rm UnitMakefile"+str(index)+ " test" + str(index) + ".x"
    return command


# liste de machines
machine_list = net.hosts
nb_machines = len(network.comp_host_tuples)*8

# on ignore les premieres simulations si demande
debut = 0
if (len(sys.argv) >= 2):
    debut = int(sys.argv[1])

nb_simul_ef = len(liste_fichier) + len(liste_fichier_par) - debut
# nombre de simulations terminee
simul_done = [0] * nb_simul_ef     #simulation state
# pour chaque simulation, numero de machine
simul_machines_index = [0] * nb_simul_ef
# processus de chaque  simulation 
simul_process = [0] * nb_simul_ef
# nombre de simulation sur chaque machine
simul_machines = [0] * len(machine_list)

# retourne une machine libre
def GetMachine():
    i = 0
    for machine in machine_list:
        load_machine = machine.LoadAverage(1.0)
        if  load_machine[0] < machine.cpu - 0.7 and simul_machines[i] < machine.cpu:
            break
        
        i = i + 1
        
    network_busy = (i == len(machine_list))
    return network_busy, machine, i


success_test = True
fid = open("/home/durufle/log_compil_test.txt", "w")
while min(simul_done) != 2:
    nb_running_simulations = simul_done.count(1)
    try:
        run_index = simul_done.index(0)
    except:
        run_index = -1
    
    quota_reached = nb_running_simulations >= nb_machines
    all_launched = run_index == -1
    if not all_launched and not quota_reached:
        # Searches for a machine.
        [network_busy, machine, i] = GetMachine()
        if not network_busy:
            print "Running test #" + str(debut + run_index) \
                  + " on " + machine.name + "."
            simul_machines[i] += 8
            simul_machines_index[run_index] = i

            # Getting parameters of the simulation
            command_line = SetSimulation(debut+run_index)
            
            # Launches the simulation.
            simul_process[run_index] = net.LaunchBG(command_line, machine)
            simul_done[run_index] = 1
    
    if all_launched or quota_reached or network_busy:

        # Selects "running" simulations that have just ended.
        for i in range(nb_simul_ef):
            if (simul_done[i] == 1 \
                and simul_process[i].poll() != None):
                if (simul_process[i].poll() != 0):
                    if (i+debut >= len(liste_fichier)):
                        j = i + debut - len(liste_fichier)
                        print "Test parallele ", j, " : ", liste_fichier_par[j], " failed."
                        fid.write("Test parallel " + str(j) + " : " +  liste_fichier_par[j] + " failed.\n")
                    else:
                        print "Test ", i, " : ", liste_fichier[i+debut], " failed."
                        fid.write("Test " +  str(i) + " : " + liste_fichier[i+debut] + " failed.\n")
                    success_test = False
                
                simul_done[i] = 2
                simul_machines[simul_machines_index[i]] -= 8
        
        os.system(waiting_command)
        
    all_done = min(simul_done) >= 2


if (success_test):
    print "All tests passed successfully"
    fid.write("All tests passed successfully\n")
else:
    sys.exit(-1)
