import os, string, sys, math, array, numpy

# fonction pour lire les parametres d'un fichier de donnees
# entree : file_name, nom du fichier de donnees
#          keyword, mot-cle a rechercher
# sortie : keyword_found : True si le mot cle a ete trouve
#          param : liste des parametres associes au mot-cle
def ReadDataFile(file_name, keyword, num = 0):
    # on lit le fichier de donnees
    file_in = open(file_name, mode='r');
    ligne_datafile = file_in.readlines();
    file_in.close();

    # on le cherche dans chaque ligne
    keyword_found = False
    param = []
    ind = 0
    for i in range(len(ligne_datafile)):
        ligne_val = ligne_datafile[i]
        # on zappe le caractere entree '\n'
        ligne_val = ligne_val[0:(len(ligne_val)-1)]
        liste_mots = ligne_val.split(" ")
        header = liste_mots[0]        
        if (header==keyword):
            param = liste_mots[2:len(liste_mots)]
            keyword_found = True
            if (ind == num):
                return keyword_found, param
            
            ind = ind + 1
    
    if (num == -1):
        return ind
    
    return keyword_found, param

# fonction pour modifier un fichier de donnees
# entree : file_name : nom du fichier d'entree
#          keyword : mot cle a modifier
#          param : parametres associe
#          file_name_out : nom du fichier de sortie
def ModifyDataFile(file_name, file_name_out, keyword, param, num = 0):
    # on lit le fichier d'entree
    file_in = open(file_name, mode='r');
    ligne_datafile = file_in.readlines();
    file_in.close();

    # on ouvre le fichier de sortie 
    file_out = open(file_name_out, mode='w');
    
    keyword_found = False
    # boucle sur les lignes
    ind = 0
    for i in range(len(ligne_datafile)):
        ligne_val = ligne_datafile[i]
        # on zappe le caractere '\n'
        ligne_val = ligne_val[0:(len(ligne_val)-1)]
        liste_mots = ligne_val.split(" ")
        header = liste_mots[0]
        if (header==keyword):
            if (ind == num):
                file_out.write(keyword+" =")
                for j in range(len(param)):
                    file_out.write(' '+param[j])
                    
                file_out.write('\n')
            else:
                file_out.write(ligne_datafile[i])
                
            ind = ind + 1
        else:
            file_out.write(ligne_datafile[i])
    
    
    
    file_out.close()

# retourne l'extension d'un nom de fichier
def GetExtension(nom):
    parts = []
    if (nom.rfind('.') != -1):
      parts = nom[nom.rfind('.'):len(nom)];
    return parts

# retourne la racine d'un nom de fichier
def GetBaseString(nom):
    parts = nom
    if (nom.rfind('.') != -1):
      parts = nom[0:nom.rfind('.')]
    
    return parts

# convertit un entier en string de 4 lettres : 0 -> '0000', 1 -> '0001', 23 -> '0023'
def NumberToString(n):
    n0 = n
    chaine = ""
    for i in range(4):
        r = n0%10; n0 = n0//10
        chaine = str(r)+chaine

    return chaine

# retourne la liste des fichiers de sortie
# exemple : si name_output = 'toto.dat', que extensions = ['_U0.dat', '_U1.dat']
# et que nb_snap = 2, ca va donner ['toto0001_U0.dat', 'toto0001_U1.dat', 'toto0002_U0.dat', 'toto0002_U1.dat']
def GetOutputFiles(keyword, name_output, extensions, nb_snapshot):
    base_name = name_output
    for i in range(len(name_output)):
        base_name[i] = GetBaseString(name_output[i])
            
    nb_snap = nb_snapshot
    if (keyword == "FileOutputSismoString"):
        nb_snap = 0
    elif (keyword == "FileOutputString"):
        nb_snap = 0
    elif (keyword == "FileOutputShank"):
        nb_snap = 0
    elif (keyword == "FileOutputEnergyPiano"):
        nb_snap = 0
    
    if (nb_snap == 0):
        nb = 0
        nb_files = len(base_name)
        if (keyword == "FileOutputSismoString"):
            nb_files = 2
        
        name_output = list(range(nb_files*len(extensions)))
        for i in range(nb_files):
            for j in range(len(extensions)):
                name_output[nb] = base_name[i]+extensions[j]
                nb = nb + 1;
    else:
        nb = 0
        name_output = list(range(len(base_name)*len(extensions)*nb_snap))
        for i in range(len(base_name)):
            for j in range(len(extensions)):
                for k in range(nb_snap):
                    if (GetExtension(extensions[j]) == ".bb"):
                        name_output[nb] = base_name[i]+GetBaseString(extensions[j])+ \
                        '.'+str(k+1)+GetExtension(extensions[j])
                    else:
                        name_output[nb] = base_name[i]+NumberToString(k+1)+extensions[j]
                    
                    nb = nb + 1;
    return name_output

def GetOutputExtensions(keyword, complex_eq, component, nb_u, nb_du, data_file):
    extensions = []
    if (keyword.startswith("FileOutputMesh")):
        if (complex_eq):
            if ((component == -1) or (component == -2)):
                for i in range(nb_u):
                    extensions.append("_U"+str(i)+"_real.bb")
                    extensions.append("_U"+str(i)+"_imag.bb")
                    extensions.append("_U"+str(i)+"_abs.bb")
                    
                if (component == -1):
                    for i in range(nb_du):
                        extensions.append("_dU"+str(i)+"_real.bb")
                        extensions.append("_dU"+str(i)+"_imag.bb")
                        extensions.append("_dU"+str(i)+"_abs.bb")
            else:
                extensions.append("_real.bb")
                extensions.append("_imag.bb")
                extensions.append("_abs.bb")
               
        else: 
            if ((component == -1) or (component == -2)):
                for i in range(nb_u):
                    extensions.append("_U"+str(i)+".bb")
                 
                if (component == -1):
                    for i in range(nb_du):
                        extensions.append("_dU"+str(i)+".bb")
                  
            else:
                extensions.append(".bb")
               
    elif (keyword == "FileRCS") or (keyword == "FileOutputPoint"):
        extensions.append(".dat")
    elif (keyword == "FileOutputDisplayZ"):
        [fourier_found, param_fourier] = ReadDataFile(data_file, "FourierTransform");
        [polar_found, param_polar] = ReadDataFile(data_file, "Polarization");
        elliptic = False
        if (polar_found):
            if (param_polar[0] == "Elliptic"):
                elliptic = True
        
        if (fourier_found and (param_fourier[0] == "ODD_FREQUENCIES")):
            n = int(param_fourier[1])
            if (elliptic):
                for k in range(n):
                    extensions.append("_X_Freq"+str(k)+".dat")
                    extensions.append("_Y_Freq"+str(k)+".dat")
            else:
                for k in range(n):
                    extensions.append("_Freq"+str(k)+".dat")
        else:
            if (elliptic):
                extensions.append("_X.dat")
                extensions.append("_Y.dat")
                extensions.append("_timeX.dat")
                extensions.append("_timeY.dat")
            else:
                extensions.append(".dat")
                extensions.append("_time.dat")
    elif (keyword == "FileOutputSismoString"):
        nb_string = ReadDataFile(data_file, "String", -1)
        for i in range(nb_string):
            extensions.append(str(i)+".txt")
    elif (keyword == "FileOutputString"):
        extensions.append(".txt")
    elif (keyword == "FileOutputShank"):
        extensions.append(".txt")
    elif (keyword == "FileOutputEnergyPiano"):
        extensions.append(".txt")
    elif (keyword == "DisplayEigenvectors"):
        if ((component == -1) or (component == -2)):
            for i in range(nb_u):
                extensions.append("G0_U"+str(i)+".dat")
                  
            if (component == -1):
                for i in range(nb_du):
                    extensions.append("G0_dU"+str(i)+".dat")
        else:
            extensions.append("G0.dat");
    else:
        if ((component == -1) or (component == -2)):
            for i in range(nb_u):
                extensions.append("_U"+str(i)+".dat")
                  
            if (component == -1):
                for i in range(nb_du):
                    extensions.append("_dU"+str(i)+".dat")
                  
        else:
            extensions.append(".dat")
    
    return extensions

# lecture d'un fichier de sortie Montjoie
# en sortie le vecteur solution U
# en entree : name_file nom du fichier de sortie Montjoie
#            keyword : mot-cle associe au fichier de sortie (FileOutputPlane, etc)
#            dim_N : cas 2-D ou 3-D
def ReadOutputFile(name_file, keyword, dim_N2):
   xsol = []
   success = True  
   dim_N = dim_N2
   if (keyword == "FileOutputPlane"):
       dim_N = max(dim_N, 2)

   if (keyword == "MecanicFileOutputPlane"):
       dim_N = 2

   if (keyword == "FileOutputPlaneAxi"):
       dim_N = 3

   if (keyword == "FileOutputLineAxi"):
       dim_N = 3

   if (keyword == "FileOutputCircleAxi"):
       dim_N = 3

   if (keyword == "FileOutputPointsFileAxi"):
       dim_N = 3

   if (keyword == "FileOutputGrille"):
       dim_N = 3

   if (keyword == "AcousticFileOutputGrille"):
       dim_N = 3
       
   if (dim_N == 1):
       if (keyword == "FileOutputDisplayZ"):
           # fichier binaire lisible par matlab
           try:
               fileobj = open(name_file, mode='rb')
               ntmp = array.array('i')
               ntmp.fromfile(fileobj, 4)
               prec = ntmp[2]
               if (prec == 3):
                   prec = 2
               
               xtmp = array.array('d')
               xtmp.fromfile(fileobj, 2)
               itmp = array.array('i')
               itmp.fromfile(fileobj, 1)
               taille = itmp[0]
               itmp = array.array('i')
               itmp.fromfile(fileobj, 1)
                              
               xtmp = array.array('d')
               try:
                   xtmp.fromfile(fileobj, prec*taille)
               except EOFError:
                   ntmp.append(0)
               except MemoryError:
                   print("name_file",name_file,"taille",taille)
                   
               xsol = numpy.array(xtmp.tolist())
               fileobj.close()
           except IOError as xxx_todo_changeme:
               (errno, strerror) = xxx_todo_changeme.args
               success = False                          
       else:
           # using loadtxt
           try:
               #print "nom fichier", name_file
               V = numpy.loadtxt(name_file)
               xsol = numpy.reshape(V, numpy.size(V))
           except IOError as xxx_todo_changeme1:
               (errno, strerror) = xxx_todo_changeme1.args
               success = False           
   elif (dim_N == 2):
       
       if (keyword == "FileOutputMeshVolumetric"):
           # fichier bb
           try:
               file_in = open(name_file,"r")
               lignes = file_in.readlines()
               header = lignes[0].split()
               n = int(header[2])
               xsol = list(range(n))
               for i in range(n):
                   xsol[i] = float(lignes[i+1])
                   
               file_in.close()
           except IOError as xxx_todo_changeme2:
               (errno, strerror) = xxx_todo_changeme2.args
               success = False
       elif ((keyword == "FileRCS") or (keyword == "FileOutputPoint")):
           try:
               file_in = open(name_file,"r")
               lignes = file_in.readlines()
               nb = 0
               rcs_file = False
               if (keyword == "FileRCS"):
                   rcs_file = True
               xsol = list(range(len(lignes)*len(lignes[0].split())))
               for i in range(len(lignes)):
                   reels = lignes[i].split()
                   for j in range(len(reels)):
                       xsol[nb] = float(reels[j])
                       if (rcs_file and (j > 0) and ((xsol[nb] < -100) or (math.isinf(xsol[nb])))):
                           xsol[nb] = 0
                       nb = nb+1
         
               file_in.close()
           except IOError as xxx_todo_changeme3:
               (errno, strerror) = xxx_todo_changeme3.args
               success = False
       else:
           #fichier binaire lisible par matlab
           try:
               fileobj = open(name_file, mode='rb')
               ntmp = array.array('i')
               ntmp.fromfile(fileobj, 4)
               nb_grids = ntmp[0]
               prec = ntmp[2]
               if (prec >= 2):
                   prec = 2
               else:
                   prec = 1
                   
               prec_type = 'd'
               if ((ntmp[2] == 0) or (ntmp[2] == 2)):
                   prec_type = 'f';
               
               taille = 0
               if (ntmp[3] == 0):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj, 2)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj, 2)
                       ntmp.fromfile(fileobj, 1)
               elif (ntmp[3] == 1):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,4)
                       ntmp.fromfile(fileobj, 1)
               elif (ntmp[3] == 2):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,2)
               elif (ntmp[3] == 3):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,3)
                       ntmp.fromfile(fileobj, 1)
               elif (ntmp[3] == 6):
                   for k in range(nb_grids):
                       ntmp.fromfile(fileobj, 1)
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,2*ntmp[4])
               
               itmp = array.array('i')
               itmp.fromfile(fileobj, 1)
               taille = itmp[0]
               xtmp = array.array(prec_type)
               try:
                   xtmp.fromfile(fileobj,prec*taille)
               except EOFError:
                   ntmp.append(0)
               except MemoryError:
                   print("name_file",name_file,"taille",taille)
            
               xsol = xtmp.tolist()
               fileobj.close()
           except IOError as xxx_todo_changeme4:
               (errno, strerror) = xxx_todo_changeme4.args
               success = False
   else:
       if (keyword.startswith("FileOutputMesh")):
           # fichier bb
           try:
               file_in = open(name_file,"r")
               lignes = file_in.readlines()
               header = lignes[0].split()
               n = int(header[2])
               xsol = list(range(n))
               for i in range(n):
                   xsol[i] = float(lignes[i+1])
                   
                   file_in.close()
           except IOError as xxx_todo_changeme5:
               (errno, strerror) = xxx_todo_changeme5.args
               success = False
       elif ((keyword == "FileRCS") or (keyword == "FileOutputPoint")):
           try:
               file_in = open(name_file,"r")
               lignes = file_in.readlines()
               nb = 0
               rcs_file = False
               if (keyword == "FileRCS"):
                   rcs_file = True

               xsol = list(range(len(lignes)*len(lignes[0].split())))
               for i in range(len(lignes)):
                   reels = lignes[i].split()
                   for j in range(len(reels)):
                       xsol[nb] = float(reels[j])
                       if (rcs_file and (j > 0) and ((xsol[nb] < -100) or (math.isinf(xsol[nb])))):
                           xsol[nb] = 0

                       nb = nb+1
                       
               file_in.close()
           except IOError as xxx_todo_changeme6:
               (errno, strerror) = xxx_todo_changeme6.args
               success = False
       else:
           #fichier binaire lisible par matlab
           try :
               fileobj = open(name_file, mode='rb')
               ntmp = array.array('i')
               ntmp.fromfile(fileobj, 4)
               nb_grids = ntmp[0]
               prec = ntmp[2]
               if (prec >= 2):
                   prec = 2
               else:
                   prec = 1
                   
               prec_type = 'd'
               if ((ntmp[2] == 0) or (ntmp[2] == 2)):
                   prec_type = 'f';
                   
               taille = 0
               if (ntmp[3] == 0):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,2)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj,2)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj,2)
                       ntmp.fromfile(fileobj, 1)                   
               elif (ntmp[3] == 1):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,2)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj,2)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj,2)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj,3)
               elif (ntmp[3] == 2):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj, 9)
                       ntmp.fromfile(fileobj, 2)
               elif (ntmp[3] == 3):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,6)
                       ntmp.fromfile(fileobj, 1)
               elif (ntmp[3] == 4):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj,3)
               elif (ntmp[3] == 5):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       ntmp.fromfile(fileobj, 1)
                       xtmp.fromfile(fileobj,3*ntmp[4])
               elif (ntmp[3] == 6):
                   for k in range(nb_grids):
                       xtmp = array.array(prec_type)
                       xtmp.fromfile(fileobj, 8)
                       ntmp.fromfile(fileobj, 1)
               
               itmp = array.array('i')
               itmp.fromfile(fileobj, 1)
               taille = itmp[0]
               xtmp = array.array(prec_type)
               try:
                   xtmp.fromfile(fileobj,prec*taille)
               except EOFError:
                   ntmp.append(0)
               except MemoryError:
                   print("name_file",name_file,"taille",taille)
                   
               xsol = xtmp.tolist()
               fileobj.close()
           except IOError as xxx_todo_changeme7:
               (errno, strerror) = xxx_todo_changeme7.args
               success = False
   
   return success,xsol

#return dimension of finite element
def GetDimensionElement(param):
    dim = 3
    if (len(param) > 0):
        if (param[0].startswith("QUADRANGLE") or param[0].startswith("TRIANGLE")):
            dim = 2
        elif (param[0].startswith("EDGE")):
            return 1
    else:
        return 1
    
    return dim

def GetNbComponentsEquation(dim, element, param):
    nb_comp = 0
    nb_grad = 0
    if (len(param) > 0):
        if ( (param[0] == "HELMHOLTZ") or (param[0] == "ACOUSTIC") or (param[0] == "LAPLACE") or (param[0] == "HARMONIC_ADVECTION_DG")):
            nb_comp = 1
            nb_grad = dim
        elif ( (param[0] == "HELMHOLTZ_SIPG") or (param[0] == "ACOUSTIC_SIPG") or (param[0] == "LAPLACE_SIPG") ):
            nb_comp = 1
            nb_grad = dim
        elif ( (param[0] == "HELMHOLTZ_DIV") or (param[0] == "LAPLACE_DIV") ):
            nb_comp = dim
            nb_grad = 1
        elif ((param[0] == "HELMHOLTZ_DG") or (param[0] == "HELMHOLTZ_HDG") or (param[0] == "LAPLACE_DG") \
              or (param[0] == "ACOUSTIC_DG") or (param[0] == "ACOUSTIC_HDG")):
            nb_comp = 1+dim
            nb_grad = 0
        elif ( (param[0] == "TIME_ELASTIC") or (param[0] == "TIME_GALBRUN") or \
                   (param[0] == "HARMONIC_ELASTIC") or (param[0] == "STATIC_ELASTIC")):
            nb_comp = dim
            nb_grad = dim*dim
        elif ( (param[0] == "TIME_VIBRO") or (param[0] == "TIME_FLUID_STRUCTURE")):
            nb_comp = 1
            nb_grad = dim
        elif (param[0] == "HARMONIC_MAXWELL2D"):
            nb_comp = 2
            nb_grad = 1
        elif ( (param[0] == "HARMONIC_MAXWELL") or (param[0] == "TIME_MAXWELL3D") or (param[0] == "TIME_MAXWELL_HDG") or (param[0] == "STATIC_MAXWELL") ):
            nb_comp = 3
            nb_grad = 3
        elif ( (param[0] == "TIME_REISSNER_MINDLIN") or (param[0] == "HARMONIC_REISSNER_MINDLIN")):
            nb_comp = 3
            nb_grad = 6
        elif ( (param[0] == "MODE_MAXWELL")):
            nb_comp = 4
            nb_grad = 2
        elif ( (param[0] == "AERO_ACOUSTIC") or (param[0] == "HARMONIC_AEROACOUSTIC")):
            nb_comp = dim + 1
            nb_grad = 0
        elif (param[0] == "HARMONIC_LEE"):
            nb_comp = dim+2
            nb_grad = 0
        elif ( (param[0] == "HARMONIC_GALBRUN") or (param[0] == "HARMONIC_GALBRUN_MODEL_DIV")):
            nb_comp = 2*dim + 1
            nb_grad = nb_comp*dim
        elif (param[0] == "HARMONIC_ELASTIC_DG"):
            nb_comp = dim*(dim+3)//2
            nb_grad = 0
        elif (param[0] == "HARMONIC_GALBRUN_SIPG"):
            nb_comp = dim
            nb_grad = nb_comp*dim
        elif (param[0] == "HARMONIC_GALBRUN_H1"):
            nb_comp = dim+1
            nb_grad = nb_comp*(dim+1)
        elif (param[0] == 'CELLULAIRE_PERIODIQUE'):
            nb_comp = 3
            nb_grad = 0
    
    return nb_comp, nb_grad

def IsEquationSolvedWithComplex(param):
    if (len(param) > 0):
        if ( (param[0] == "ACOUSTIC") or (param[0] == "LAPLACE") or \
                 (param[0] == "AERO_ACOUSTIC")):
            return False

        if ( (param[0] == "CELLULAIRE_PERIODIQUE") or (param[0] == "STATIC_MAXWELL") \
                 or (param[0] == "TIME_MAXWELL3D") or (param[0] == "TIME_MAXWELL_HDG") or (param[0] == "TIME_REISSNER_MINDLIN") \
                 or (param[0] == "TIME_VIBRO") or (param[0] == "TIME_FLUID_STRUCTURE") \
                 or (param[0] == "MODE_MAXWELL") ):
            return False

        if ( (param[0] == "TIME_ELASTIC") or (param[0] == "TIME_GALBRUN") ):
            return False

    return True
