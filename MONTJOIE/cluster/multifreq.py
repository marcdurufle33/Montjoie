from module_regression import *
import os
from numpy import *
from visuND import *

f0=0.01;
fmax=4;
N=30;
freq=linspace(f0,fmax,N);

Err = zeros(len(freq))

for i in range(len(freq)):
#for i in range(2):
    f=freq[i]
    ModifyDataFile("../sol_num.ini","../temp.ini","Frequency",[str(freq[i]),"0.0"])
    os.system("../helmholtz2D.x ../temp.ini > sortHelm")
    [X,Y,Z,coor,V]=loadND('totalCarre_U0.dat')
    ModifyDataFile("../sol_exacte.ini","../tempex.ini","Frequency",[str(freq[i])])
    os.system("../solution_disc.x ../tempex.ini > sortEx")
    [X1,Y1,Z1,coor1,V1]=loadND('UnExacte.dat')
    Err[i]  = erreurL2(V,V1.T)
    print "Erreur a freq", freq[i], " = ", Err[i]

print Err

savetxt("ErrPMLNonDispersive001.dat",Err)
savetxt("Freq.dat",freq)
