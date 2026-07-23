#!/usr/bin/env python

# fichier de generation des tableaux python pour les tests de non-regression
# pour rajouter un test de non-regression, placer le dans un sous-repertoire de example
# (par exemple example/laplace/MonFichier.ini), et rajouter un champ :
# DataNonRegressionTest = nom_executable nb_snapshots nb_proc
# (par exemple DataNonRegressionTest = acous2D.x 4 8)
# le nombre d'instantanes sera egal a 0 si la simulation demande est une
# simulation harmonique

from module_regression import *

import os

liste_data_files = []
liste_nb_snapshots = []
liste_element = []
liste_command = []
dimension_element = []
liste_nb_components_unknown = []
liste_nb_components_gradient = []
complexe_solution = []
nb_proc_for_simulation = []

# on boucle sur les sous-repertoires de example
for rep in sorted(os.listdir("../example")):
    if (os.path.isdir("../example/"+rep)):
        for fichier in sorted(os.listdir("../example/"+rep)):
            nom_complet = "../example/" + rep + "/" + fichier
            # on ne considere que les fichiers ini 
            if (fichier.endswith(".ini")):
                keyword_found, param = ReadDataFile(nom_complet, "DataNonRegressionTest")
                if (keyword_found):
                    #print fichier, param
                    # on rajoute le fichier de donnees a la liste
                    liste_data_files.append(nom_complet)
                    if (len(param) < 3):
                        print("You need to enter at least three parameters in file", fichier)
                        print("Actual parameters", param)
                        
                    executable = param[0]
                    if (not executable.startswith("./")):
                        executable = "./"+executable
                    
                    liste_command.append(executable)                    
                    
                    # nombre de snapshots
                    try:
                        liste_nb_snapshots.append(int(param[1]))
                    except:
                        print("Probleme de cast en entier sur le fichier ", fichier)
                        print(param[1], "n'est pas un entier")
                    
                    # nombre de processeurs
                    try:
                        nb_proc_for_simulation.append(int(param[2]))
                    except:
                        print("Probleme de cast en entier sur le fichier ", fichier)
                        print(param[2], "n'est pas un entier")
                    
                    # element et dimension
                    keyword_found, element = ReadDataFile(nom_complet, "TypeElement")
                    if (len(element) > 0):
                        liste_element.append(element[0])
                    else:
                        liste_element.append('')
                        
                    dim = GetDimensionElement(element)
                    dimension_element.append(dim)
                    
                    # nombre de composantes et reel ou complexe
                    keyword_found, equation = ReadDataFile(nom_complet, "TypeEquation")
                    nb_comp, nb_grad = GetNbComponentsEquation(dim, element, equation)
                    liste_nb_components_unknown.append(nb_comp)
                    liste_nb_components_gradient.append(nb_grad)
                    complex_sol = IsEquationSolvedWithComplex(equation)
                    complexe_solution.append(complex_sol)

#print len(liste_data_files)
#print liste_data_files
#print liste_command
#print liste_nb_snapshots
#print liste_element
#print complexe_solution
#print liste_nb_components_unknown
#print dimension_element
#print nb_proc_for_simulation

liste_code = range(len(liste_data_files))
