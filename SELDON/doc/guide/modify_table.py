import os, sys
import string

# pour alterner les lignes 1/2 dans les category-table
def alterne_numeros(fichier_dox):
    f = open(fichier_dox, "r")
    lignes_dox = f.readlines()
    f.close()

    # on alterne les numeros
    level_table = -1
    longueur = len('category-table-tr-')
    for i in range(len(lignes_dox)):
        if (lignes_dox[i].find('<table class="category-table">') >= 0):
            level_table = 0
            num = 1
        else:
            if (level_table >= 0):
                if (lignes_dox[i].find('<table') >= 0):
                    level_table += 1
                elif (lignes_dox[i].find('</table') >= 0):
                    level_table -= 1;
                elif (lignes_dox[i].find('category-table-tr-') >= 0):
                    pos = lignes_dox[i].find('category-table-tr-')
                    lignes_dox[i] = lignes_dox[i][:pos+longueur] + str(num) + lignes_dox[i][pos+longueur+1:]
                    num = 3-num

    # et on ecrase le fichier
    f = open(fichier_dox, "w")
    for i in range(len(lignes_dox)):
        f.write(lignes_dox[i])
    
    f.close()

liste_fichiers = sorted(os.listdir("."))
for f in liste_fichiers:
    if (f.endswith('.dox')):
        alterne_numeros(f)



