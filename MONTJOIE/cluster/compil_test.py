#!/usr/bin/python

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

# liste des regles a compiler
target_compil = ["lung", "vibro", "fluid_structure", "multistringStiffNL", "PianoStiffNL", \
         "laplace", "helm2D", "helm3D", "maxwell2D","maxwell_axi", "static_maxwell", "maxwell3D", \
         "mode_maxwell", "time_maxwell", "time3D", "acous2D", "acous3D", "harmonic_elas2D", \
         "harmonic_elas3D", "elas2D", "elas3D", "reissner", "static_elastic", \
         "soundboard", "aero2D", "galbrun", "aero3D", "harmonic_aero2D", "StringImpedanceStiffLIN", \
         "solution_helm", "solution", "dispersion", "CellulairePeriodique", \
         "CellulairePeriodique3D", "kdv", "camassa", "convert", "manipule", "manipule2D", \
         "sym_mesh", "per_mesh", "tetmesh", "schrodinger1D", "helm_radial", "helm_axi", \
         "StringImpedanceStiffNL", "StringSoundboardStiffNL", "subdomain", "write_nodal", \
         "galbrun_axi", "aero_axi", "helm_div", "schrodinger3D"]


# repertoire courrant
current_dir = os.getcwd()

# liste de machines
machine_list = net.hosts
nb_machines = len(network.comp_host_tuples)*8

nb_simul_ef = len(target_compil)

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

# pour tous les index on donne la commande a realiser
def SetSimulation(index):
    return "cd "+current_dir+"/..; source ~/.bash_profile; export OMP_NUM_THREADS=1; make " + target_compil[index] + " >& cluster/errMake" + str(index)

# tant que toutes les compilations ne sont pas terminees
success_compil = True
nom_fichier = "/home/durufle/log_compil_test.txt"
fid = open(nom_fichier, "w")
fid.close()
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
            print "Running compilation #" + str(run_index) \
                  + " on " + machine.name + "."
            simul_machines[i] += 8
            simul_machines_index[run_index] = i

            # Getting parameters of the simulation
            command_line = SetSimulation(run_index)
            
            # Launches the simulation.
            simul_process[run_index] = net.LaunchBG(command_line, machine)
            simul_done[run_index] = 1
    
    if all_launched or quota_reached or network_busy:

        # Selects "running" simulations that have just ended.
        for i in range(nb_simul_ef):
            if (simul_done[i] == 1 \
                and simul_process[i].poll() != None):
                if (simul_process[i].poll() != 0):
                    fid = open(nom_fichier, "a")
                    print "Compilation of target ", i, " : ",  target_compil[i], " failed"
                    fid.write("Compilation of target " +  target_compil[i] + " failed.\n")
                    fid.close()
                    success_compil = False
                
                simul_done[i] = 2
                simul_machines[simul_machines_index[i]] -= 8
        
        os.system(waiting_command)
        
    all_done = min(simul_done) >= 2

if (success_compil):
    fid = open(nom_fichier, "a")
    print "All compilations successfully passed\n"
    fid.write("All compilations successfully passed\n")
    fid.close()
else:
    sys.exit(1)


