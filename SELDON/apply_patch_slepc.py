#!/usr/bin/env python

import os

# - we add matsolve in include/slepc/private/stimpl.h in the structure _p_ST :
 #   PetscErrorCode (*matsolve)(Mat A, Vec b, Vec x);
#    PetscErrorCode (*matsolve_trans)(Mat A, Vec b, Vec x);
    
try:
    fid = open("include/slepc/private/stimpl.h")
    lignes = fid.readlines(); fid.close();
    test_struct = False;
    for i in range(len(lignes)):
        if (lignes[i].startswith('struct _p_ST')):
            test_struct = True;

        if (test_struct):
            # pour detecter si le patch a deja ete applique
            if (lignes[i].find("matsolve") >= 0):
                test_struct = False;
            
            if (lignes[i].startswith('};')):
                lignes.insert(i, "  PetscErrorCode (*matsolve_trans)(Mat A, Vec b, Vec x);\n");
                lignes.insert(i, "  PetscErrorCode (*matsolve)(Mat A, Vec b, Vec x);\n");
                break;
    
    fid = open("include/slepc/private/stimpl.h", 'w')
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify include/slepc/private/stimpl.h")


# - we modify STMatSolve in src/sys/classes/st/interface/stsles.c to put
#    if (st->matsolve != NULL)
#    return (*st->matsolve)(st->A[0], b, x);

#  et in STMatSolveTranspose
#    if (st->matsolve_trans != NULL)  
#    return (*st->matsolve_trans)(st->A[0], b, x);
try:
    fid = open("src/sys/classes/st/interface/stsles.c")
    lignes = fid.readlines(); fid.close();
    test_func = False; i = 0;
    while (i < len(lignes)):
        if (lignes[i].startswith("PetscErrorCode STMatSolve(ST st")):
            if (lignes[i+2].find("matsolve") == -1):
                lignes.insert(i+2, "    return (*st->matsolve)(st->A[0], b, x);\n");
                lignes.insert(i+2, "  if (st->matsolve != NULL)\n");

        if (lignes[i].startswith("PetscErrorCode STMatSolveTranspose(ST st")):
            if (lignes[i+2].find("matsolve") == -1):
                lignes.insert(i+2, "    return (*st->matsolve_trans)(st->A[0], b, x);\n");
                lignes.insert(i+2, "  if (st->matsolve_trans != NULL)\n");

        i += 1
    
    fid = open("src/sys/classes/st/interface/stsles.c", 'w')
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify src/sys/classes/st/interface/stsles.c")
                
# - we initialize matsolve in STCreate (stfunc.c)
#    st->matsolve = NULL;
#    st->matsolve_trans = NULL;
try:
    fid = open("src/sys/classes/st/interface/stfunc.c")
    lignes = fid.readlines(); fid.close();
    i = 0;
    while (i < len(lignes)):
        if (lignes[i].startswith("PetscErrorCode STCreate")):
            j = i; test_solve = False; pos_solve = i+4;
            while (j < len(lignes)):
                if (lignes[j].startswith('}')):
                    break;

                if (lignes[j].find("st->data") >= 0):
                    pos_solve = j+1;
                
                if (lignes[j].find("matsolve") >= 0):
                    test_solve = True

                j += 1
            
            if (not test_solve):
                lignes.insert(pos_solve, "  st->matsolve_trans = NULL;\n");
                lignes.insert(pos_solve, "  st->matsolve = NULL;\n");

            break;

        i += 1
                
    fid = open("src/sys/classes/st/interface/stfunc.c", 'w')
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify src/sys/classes/st/interface/stfunc.c")

# we modify nleigs.c (corrections for matrix-shell implementation)
try:
    fid = open("src/nep/impls/nleigs/nleigs.c")
    lignes = fid.readlines(); fid.close()
    i = 0;
    while (i < len(lignes)):
        if (lignes[i].startswith("static PetscErrorCode NEPNLEIGSDividedDifferences_callback")):
            j = i; pos_bv_create = -1
            pres_ligne_rand_ctx = False
            while (j < len(lignes)):
                if (lignes[j].startswith('}')):
                    break;
                
                if (lignes[j].find("BVGetRandomContext") >= 0):
                    pres_ligne_rand_ctx = True
                    ligne_rand_ctx = lignes[j]
                    del lignes[j]
                    j -= 1
            
                if (lignes[j].find("MatCreateVecs(D[0]") >= 0):
                    pos_bv_create = j;
                    break;
            
                j += 1
            
            if (pos_bv_create >= 0):
                if (lignes[pos_bv_create].find("BVCreate") < 0):
                    lignes.insert(pos_bv_create, "    ierr = BVCreate(PetscObjectComm((PetscObject)nep), &nep->V);CHKERRQ(ierr);\n");
                    if (pres_ligne_rand_ctx):
                        lignes.insert(pos_bv_create+1, ligne_rand_ctx)
            
            break;
                
        i += 1;

    i = 0
    while (i < len(lignes)):
        if (lignes[i].find("MatCreateShell") >= 0):
            mots = lignes[i].split(',')
            if (mots[1] == 'n'):
                chaine = "PetscObject)";
                pos = lignes[i].find(chaine)
                name_matrix = lignes[i][pos+len(chaine)];
                lignes.insert(i, "PetscInt nloc; ierr = MatGetLocalSize("+name_matrix+",&nloc,NULL);CHKERRQ(ierr);\n");
                lignes[i+1] = mots[0] +", nloc, nloc";
                for j in range(3, len(mots)):
                    lignes[i+1] += "," + mots[j]

                i += 1
        i += 1
    
    fid = open("src/nep/impls/nleigs/nleigs.c", "w")
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify src/nep/impls/nleigs/nleigs.c")

# we add a symbolic link to petsc sources
chemin = os.environ.get("PETSC_DIR")
presence_src = False
for rep in os.listdir(chemin + "/include/petsc"):
    if (rep == "src"):
        presence_src = True

if (not presence_src):
    commande_ln = "ln -s " + chemin +"/src " + chemin + "/include/petsc/src"
    os.system(commande_ln)
