#!/usr/bin/env python

import sys, pwd, os, math

import string, network
import subprocess

# Network.
net = network.Network("all")

# Network status.
network_busy = False
all_launched = False
all_done = False

# Waiting time.
waiting_command = "sleep 5"

# repertoire courrant
current_dir = os.getcwd()

from test_regression import *

# liste de machines
machine_list = net.hosts

# on ignore les premieres simulations si demande
debut = 0
if (len(sys.argv) >= 2):
    debut = int(sys.argv[1])

nb_simul_ef = len(liste_data_files) - debut
# nombre de simulations terminee
simul_done = [0] * nb_simul_ef     #simulation state
# pour chaque simulation, numero de machine
simul_machines_index = [0] * nb_simul_ef
# processus de chaque  simulation 
simul_process = [0] * nb_simul_ef
# nombre de simulation sur chaque machine
simul_machines = [0] * (len(machine_list)+1)

# retourne une machine libre
def GetMachine(n):
    network_busy = True
    if (n > 8):
        # at most 10 jobs in the cluster
        if (simul_machines[len(machine_list)] < 10):
            network_busy = False
        
        return network_busy, "cluster", len(machine_list)
    
    i = 0
    for machine in machine_list:
        load_machine = machine.LoadAverage(1.0)
        # if  ( (load_machine[0]+n < (machine.cpu + 0.3)) and (simul_machines[i]+n <= machine.cpu) ):
        if (simul_machines[i] < machine.cpu):
            network_busy = False
            break
        
        i = i + 1
        
    return network_busy, machine, i

def LaunchParallelJob(command_line, index_simu, nb_proc, pbs_job):
    if (pbs_job):
        nb_nodes = nb_proc/8
        fid = open("job_temp", "w")
        fid.write("#PBS -N time3D\n")
        fid.write("#PBS -l walltime=1:59:00\n")
        fid.write("#PBS -l mem="+str(nb_nodes*20)+"gb\n")
        fid.write("#PBS -l nodes="+str(nb_nodes)+":ppn=8\n");
        fid.write("export OMP_NUM_THREADS=1\n");
        #    fid.write("source /home/durufle/.bash_profile\n")
        fid.write(command_line)
        fid.write('\n')
        fid.close()
        os.system("qsub job_temp")
    else:
        if (nb_proc%24 == 0):
            nb_cores = 24
            nb_nodes = nb_proc/24
        elif (nb_proc%16 == 0):
            nb_cores = 16
            nb_nodes = nb_proc/16
        else:
            nb_cores = 8
            nb_nodes = nb_proc/8
        
        fid = open("job_testN"+str(index_simu), "w")
        fid.write("#!/bin/bash\n")
        fid.write("#SBATCH --job-name=testN"+str(index_simu)+"\n")
        fid.write("#SBATCH --time=01:59:59\n")
        fid.write("#SBATCH -N " + str(nb_nodes) + "\n")
        #fid.write("#SBATCH --cpus-per-task "+ str(nb_cores) + "\n")
        fid.write("#SBATCH -n "+ str(nb_nodes*nb_cores) + "\n")
        fid.write("#SBATCH -p court\n")
        fid.write("cd /home/durufle/montjoie/MONTJOIE\n")
        fid.write("export OMP_NUM_THREADS=1\n")
        # pour que ca tourne sur des coeurs separes
        fid.write("unset MV2_ENABLE_AFFINITY\n")
        fid.write(command_line)
        fid.write('\n')
        fid.close()
        os.system("sbatch job_testN"+str(index_simu))

def FindJobName(nom):
    err_squeue = os.system("squeue >& sortQ")
    if (err_squeue != 0):
        return True
    else:
        fid = open("sortQ", "r")
        lignes = fid.readlines()
        fid.close()
        for ligne in lignes:
            mots = ligne.split()
            if (len(mots) > 3):
                if (mots[2] == nom):
                    return True
        
        return False

success_test = True
fid = open("/home/durufle/log_regression_test.txt", "w")
while min(simul_done) != 2:
    nb_running_simulations = simul_done.count(1)
    try:
        run_index = simul_done.index(0)
    except:
        run_index = -1
    
    all_launched = run_index == -1
    if not all_launched:
        # Searches for a machine.
        nb_proc = nb_proc_for_simulation[debut+run_index]
        [network_busy, machine, i] = GetMachine(nb_proc)
        if not network_busy:
            # Getting parameters of the simulation
            command_line = SetSimulation(debut+run_index)
            
            # Launches the simulation.
            simul_done[run_index] = 1            
            if (machine == "cluster"):
                print("Running test #" + str(debut + run_index) + " on the cluster")
                simul_machines[i] += 1
                simul_machines_index[run_index] = i
                LaunchParallelJob(command_line, debut+run_index, nb_proc, False)
            else:
                print("Running test #" + str(debut + run_index) \
                    + " on " + machine.name + ".")
                simul_machines[i] += nb_proc
                simul_machines_index[run_index] = i
                command_line = "source ~/.bash_profile; export OMP_NUM_THREADS=1; " + command_line
                simul_process[run_index] = net.LaunchBG(command_line, machine)
            
    if all_launched or network_busy:

        # Selects "running" simulations that have just ended.
        for i in range(nb_simul_ef):
            nb_proc = nb_proc_for_simulation[debut+i]
            if (simul_done[i] == 1):
                if (simul_machines_index[i] == len(machine_list)):
                    # parallel job
                    job_found = FindJobName("testN"+str(debut+i))
                    if (not job_found):
                        simul_done[i] = 2
                        simul_machines[len(machine_list)] -= 1
                else:
                    if (simul_process[i].poll() != None):
                        if (simul_process[i].poll() != 0):
                            fid.write("Test " +  liste_data_files[i+debut] + " failed.\n")
                            success_test = False
                            
                        print("Job ", i, " completed")
                        simul_done[i] = 2
                        simul_machines[simul_machines_index[i]] -= nb_proc
        
        os.system(waiting_command)
        
    all_done = min(simul_done) >= 2


if (success_test):
    fid.write("All tests passed successfully")
else:
    sys.exit(-1)
