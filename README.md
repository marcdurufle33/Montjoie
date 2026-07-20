/************
* PREREQUIS *
*************/

Pour installer Montjoie, il faut installer des packages de base de
linux (pour avoir un compilateur C++/Fortran et gerer les Makefile).
Il est conseille d'avoir des notions de base du fonctionnement d'un Makefile.
Pour installer les packages de base, vous pouvez taper :

sudo apt install make cmake
sudo apt install gfortran g++

Pour une compilation/execution en parallele, installer une librairie MPI:

sudo apt install libopenmpi-dev

Il est conseille de mettre les librairies exterieures dans un repertoire independant
(par exemple un repertoire qu'on appellera Solve):

cd
mkdir Solve
cd Solve

/***************
 * BLAS/LAPACK *
 ***************/
 
Pour Blas/Lapack, si c'est possible, il est plus efficace d'utiliser
une version optimisee comme la MKL. Sur les clusters (comme Plafrim),
ce type de librairie est installee de base. On indique ci-dessous les commandes
pour installer Blas/Cblas/Lapack avec la version de base (non-optimisee) :

wget http://www.netlib.org/blas/blas.tgz
tar zxvf blas.tgz
cd BLAS-3.8.0
make
mv blas_LINUX.a libblas.a
cd ..

wget http://www.netlib.org/blas/blast-forum/cblas.tgz
tar zxvf cblas.tgz
cd CBLAS
# Modifier le chemin ou est Blas (variable BLLIB) dans Makefile.in
# par exemple, on ecrira BLLIB = /home/durufle/Solve/BLAS-3.8.0/libblas.a
make
mv lib/cblas_LINUX.a libcblas.a
ar rv libcblas_mkl.a src/zdotcsub.o src/zdotusub.o src/cdotcsub.o src/cdotusub.o
cd ..

wget http://www.netlib.org/lapack/lapack-3.5.0.tgz
tar zxvf lapack-3.5.0.tgz
cd lapack-3.5.0
cmake .
make
mv lib/liblapack.a liblapack.a


/*********
 * MUMPS *
 *********/

Pour une execution rapide de Montjoie, il est fortement recommande d'avoir une librairie
de solveur direct (Mumps, Pastix, Pardiso, etc). On indique ici comment installer Mumps :

Metis (utilise par Mumps)

installer metis
- soit en telechargeant le clone git : https://github.com/KarypisLab/METIS
- soit en installant le paquet linux : sudo apt install libmetis-dev

Mumps

wget graal.ens-lyon.fr/MUMPS/MUMPS_5.6.0.tar.gz
tar zxvf MUMPS_5.6.0.tar.gz
cd MUMPS_5.6.0
cp Make.inc/Makefile.debian.SEQ Makefile.inc (ou .PAR en parallele)

Modifier Makefile.inc (avec un editeur de texte comme emacs) pour specifier les renumerotation
On vous propose ici deux options (compilation en sequentiel ou parallele)

# En sequentiel (on choisit Metis et Pord):
# renseigner ISCOTCH a rien et IMETIS le repertoire ou vous avez installe Metis :
ISCOTCH =
IMETIS = -I/usr/include  (ou le chemin ou vous avez installe metis)
ORDERINGSF = -Dmetis -Dpord 

# En parallele (on prend tous les ordonnancements)
# Pour Scotch et ParMetis, voir MONTJOIE/INSTALL.txt :
ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
# En parallele, ne pas oublier de modifier IMETIS et ISCOTCH par exemple:
LSCOTCHDIR = /home/durufle/Solve/scotch_6.0.4/
ISCOTCH   = -I$(LSCOTCHDIR)/include # only needed for ptscotch
LMETISDIR = /home/durufle/Solve/parmetis-4.0.3/
IMETIS    = -I$(LMETISDIR)/include -I$(LMETISDIR)/metis/include

# en parallele ou sequentiel
make alllib
cd ..

/*************
* MONTJOIE *
************/

Pour installer d'autres librairies (non-necessaires, mais parfois utiles), une liste exhaustive est
presente dans le fichier MONTJOIE/INSTALL.txt.  Le repertoire SELDON contient les sources de Seldon,
une librairie C++ d'algebre lineaire pour gerer les matrices et vecteurs
(documentation sur https://www.math.u-bordeaux.fr/~durufle/seldon/ pour Seldon).
La documention pour Montjoie est sur https://www.math.u-bordeaux.fr/~durufle/montjoie/
Le repertoire MONTJOIE contient les sources de Montjoie, un code elements finis.
Une fois Blas/Lapack et Mumps installes, vous pouvez compiler Montjoie en tapant :

cd MONTJOIE
cp Makefile.LINUX Makefile

Il faut ensuite modifier le fichier Makefile pour specifier ou sont les
repertoires d'installation de Blas/Lapack et Mumps.
Dans ce fichier, il faut chercher la ligne commencant par EXTERNAL := ...
Vous mettez le repertoire ou vous avez installe ces librairies. Vous pouvez verifier que CHEMUMPS_SEQ
est correct aussi (dans le cas d'une compilation sequentielle).
Une fois les chemins precises, on peut compiler en tapant :

make -j4 helm2D

Il est ici conseille de compiler sur plusieurs coeurs en meme temps (option -j4 pour compiler sur 4 coeurs),
vous pouvez mettre -j8 si votre machine a 8 coeurs par exemple. Ca permet d'avoir une compilation
plus rapide, mais ca peut ralentir votre machine.
Une fois le code compile, vous pouvez verifier qu'il fonctionne en tapant

./helmholtz2D.x example/helm2D/carre.ini

Les fichiers de resultat sont lisible avec python
(sous Ubuntu, taper : sudo apt install python3-numpy python3-scipy python3-matplotlib).
On tapera par exemple :

ipython3 --pylab
from MATLAB.visuND import *
X, Y, Z, coor, V = loadND('totalCarre_U0.dat')
plot2dinst(X, Y, real(V), -0.1, 0.1)

Il y a aussi des scripts MATLAB (dans le repertoire MATLAB) si vous etes plus a l'aise avec ce logiciel.

La 