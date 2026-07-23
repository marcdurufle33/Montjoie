#!/usr/bin/env python
import os, string, sys, math, array, numpy


threshold = 1e-6
threshold_L1 = 1e-3

output_keywords = ["FileOutputPlane","FileOutputGrille","OutputFile","FileOutputGrille3D","FileOutputLine","FileOutputCircle", \
                   "FileOutputPoint","FileOutputPlaneAxi","FileOutputCircleAxi","FileOutputLineAxi","FileOutputPointAxi",\
                   "FileOutputMeshVolumetric","FileOutputMeshSurfacic","FileRCS","FluidFileOutputGrille","SolidFileOutputGrille", \
                   "AcousticFileOutputGrille", "FileOutputGrid", "FileOutputDisplayZ", "FileOutputString", "FileOutputSismoString", \
                   "FileOutputShank", "FileOutputEnergyPiano", "AcousticFileOutputPoint", "AcousticFileOutputGrille", "MecanicFileOutputPlane", "DisplayEigenvectors"]

from module_regression import *
from generation_example import *

def CheckSimulation(index):
    num = liste_code[index]
    data_file = liste_data_files[index]
    nb_snap = liste_nb_snapshots[index]
    name_elt = liste_element[num]
    commande = liste_command[num]
    dim_N = dimension_element[num]
    nb_u = liste_nb_components_unknown[num]
    nb_du = liste_nb_components_gradient[num]
    complex_eq = complexe_solution[num]
    # composantes a lire ? (component = -1 si toutes les composantes)
    component = 0; param = []
    [keyword_found,parameters] = ReadDataFile(data_file, "ElectricOrMagnetic");
    #print "param ",parameters
    if (keyword_found):
        component = int(parameters[0]);
    
    error_regression = False
    comparison_one_file = False
    directory_ref = "/home/durufle/NonRegressionTest/" + GetBaseString(data_file)[11:]+"/"
    name_output = []
    # on liste toutes les sorties
    for j in range(len(output_keywords)):
        #print commande
        if (commande.startswith("./maxwell_axi.x")):
            dim_N = 3
        
        [keyword_found, parameters] = ReadDataFile(data_file, output_keywords[j]);
        if (keyword_found):
            n0 = 0
            if (output_keywords[j] == "DisplayEigenvectors"):
                n0 = 1
            
            old_param = parameters;
            parameters = [0]*(len(old_param)-n0)
            for n in range(n0, len(old_param)):
                parameters[n-n0] = name_elt+'_'+old_param[n]
                            
            extensions = GetOutputExtensions(output_keywords[j], complex_eq, component, nb_u, nb_du, data_file)
            #print("extensions", extensions, nb_snap)
            name_output = GetOutputFiles(output_keywords[j], parameters, extensions, nb_snap)
            #print("parameters", parameters,"files ",name_output)
            norme_max = 0.0
            for k in range(len(name_output)):
               # on lit le fichier de resultat et le fichier de reference
               name_file_ref = directory_ref + name_output[k]
               
               [success1, sol1] = ReadOutputFile("../"+name_output[k], output_keywords[j], dim_N)
               [success2, sol2] = ReadOutputFile(name_file_ref, output_keywords[j], dim_N)
               #print("file ","../"+name_output[k],"ref",name_file_ref)
               #print output_keywords[j], dim_N
               #print("success1",success1,"success2",success2)
               #print("Fichier teste ", name_output[k])
               # si les deux fichiers existent, et que la taille est la meme
               if (success1) and (success2) and (len(sol1)==len(sol2)):
                  # comparaison
                  diff = numpy.zeros(len(sol1))
                  norme_sol = numpy.zeros(len(sol1))
                  nb_pts = 0
                  for n in range(len(sol1)):
                      #if (sol1[n] != 0) and (sol2[n] != 0):
                      diff[nb_pts] = abs(sol1[n]-sol2[n])
                      norme_sol[nb_pts] = abs(sol1[n])
                      nb_pts = nb_pts + 1
                  
                  norme = 0
                  normeL1 = 0
                  # calcul de l'erreur mediane
                  if (nb_pts == 0):
                      print("Solution nulle ",name_output[k])
                  else:
                      # erreur mediane
                      diff = numpy.sort(diff[0:nb_pts])
                      norme_sol_L1 = sum(norme_sol)
                      norme_sol_max = norme_sol.max()
                      p = int(0.95*len(diff))
                      if (norme_sol_max <= threshold*norme_max):
                          print("Solution nulle ",name_output[k])
                      else:
                          norme = diff[p] / norme_sol_max
                          normeL1 = sum(diff) / norme_sol_L1
                          norme_max = max(norme_max, norme_sol_max)
                          
                  #print "Norme = ", norme_sol, norme_max
                  print("Erreur = ", norme, normeL1)
                  comparison_one_file = True
                  if ((norme > threshold) or (normeL1 > threshold_L1)):
                     print("ERROR - Regression in the code ")
                     print("The data file is ",data_file)
                     print("The output file is ",name_output[k])
                     print("The reference file is ",name_file_ref)
                     print("Relative error is equal to ",norme)
                     error_regression = True
                
               else:
                  if (success1 != success2):
                     if (success1 == False):
                         print("No data file found ",name_output[k])
                     if (success2 == False):
                         print("No reference file found ",name_file_ref)
    
    
    if not comparison_one_file:
        print("No comparison tests performed, problem of path ?")
        print(directory_ref, name_output)
        error_regression = True
    
    if not error_regression:
       print("Data file ",data_file," successfully tested")
       

def CopySimulation(index):
    num = liste_code[index]
    data_file = liste_data_files[index]
    nb_snap = liste_nb_snapshots[index]
    dim_N = dimension_element[num]
    name_elt = liste_element[num]
    nb_u = liste_nb_components_unknown[num]
    nb_du = liste_nb_components_gradient[num]
    complex_eq = complexe_solution[num]
    # composantes a lire ? (component = -1 si toutes les composantes)
    component = 0; param = []
    [keyword_found,parameters] = ReadDataFile(data_file, "ElectricOrMagnetic");
    #print "param ",parameters
    if (keyword_found):
        component = int(parameters[0]);
    
    error_regression = False
    # on liste toutes les sorties
    for j in range(len(output_keywords)):
        #print "keyword => ", output_keywords[j]
        [keyword_found,parameters] = ReadDataFile(data_file, output_keywords[j]);
        if (keyword_found):
            n0 = 0
            if (output_keywords[j] == "DisplayEigenvectors"):
                n0 = 1

            old_param = parameters;
            parameters = [0]*(len(old_param)-n0)
            for n in range(n0, len(old_param)):
                parameters[n-n0] = name_elt+'_'+old_param[n]
            
            extensions = GetOutputExtensions(output_keywords[j], complex_eq, component, nb_u, nb_du, data_file)
            name_output = GetOutputFiles(output_keywords[j], parameters, extensions, nb_snap)
            #print "files ",name_output
            directory_ref = "/home/durufle/NonRegressionTest/" + GetBaseString(data_file)[11:]+"/"
            os.system("mkdir -p " + directory_ref)
            for k in range(len(name_output)):
               # on lit le fichier de resultat et le fichier de reference
               name_file_ref = directory_ref + name_output[k]
               
               [success1,sol1] = ReadOutputFile("../"+name_output[k], output_keywords[j], dim_N)
               print(name_file_ref, name_output[k])
               # si le fichier existe, on le copie
               if (success1):
                  os.system("cp "+"../"+name_output[k]+' '+name_file_ref)


def SetSimulation(index):
    num = liste_code[index]
    data_file = liste_data_files[index]
    exec_command = liste_command[num]
    name_elt = liste_element[num]
    nb_proc = nb_proc_for_simulation[num]
    print("fichier de donnees : ", data_file)
    command_line = TestRegression(exec_command, nb_proc, data_file, name_elt, index)
    #print "ligne commande", command_line
    return command_line

def TestRegression(exec_command, nb_proc, data_file, name_elt, index):
    
    # on modifie le fichier de donnee pour mettre le type d'element fini
    parameters = [name_elt]
    data_file_out = "temp" + str(index) + ".ini";
    sortie_file = "sortRegression" + str(index)
    ModifyDataFile(data_file, data_file_out, "TypeElement", parameters)
    
    # on change les noms des fichiers de sortie pour mettre
    # le type d'element fini devant
    for j in range(len(output_keywords)):
        [keyword_found,parameters] = ReadDataFile(data_file, output_keywords[j]);
        if (keyword_found):
            n0 = 0
            if (output_keywords[j] == "DisplayEigenvectors"):
                n0 = 1
                
            #print(output_keywords[j], n0)
            for n in range(n0, len(parameters)):
                parameters[n] = name_elt+'_'+parameters[n]
                
            ModifyDataFile(data_file_out, data_file_out, output_keywords[j], parameters);


    # commande a lancer pour Montjoie
    command_line = "cd ~/montjoie/MONTJOIE; mpirun -np " + str(nb_proc) + " " + exec_command + " cluster/" + \
        data_file_out + " >& " + sortie_file;
    
    return command_line

