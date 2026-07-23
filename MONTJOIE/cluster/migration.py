from module_regression import *
from visuND import *

import os, string, sys, numpy

Freq = numpy.array([2.0])
Sources = numpy.array([-2.0, 0.0, 3.0])
#Sources = numpy.array([0.0])

image = numpy.zeros([200, 200])
file_exp = 'totalTest'
file_simu = 'totalNum'
ini_exp = 'EXECUTION/disque_exp.ini'
ini_simu = 'EXECUTION/disque_num.ini'
os.chdir('../')
for f in Freq:
    for x in Sources:
        ModifyDataFile(ini_exp, 'temp_exp.ini', 'Frequency', [str(f), '0.0'])
        ModifyDataFile('temp_exp.ini', 'temp_exp.ini', 'TypeSource', ['SRC_SURFACE', '2', 'GAUSSIAN', str(x), '-5.0', '1.0', '1.5', '1.0'])
        #ModifyDataFile('temp_exp.ini', 'temp_exp.ini', 'TypeSource', ['SRC_SURFACE', '2', 'UNIFORM', '1.0'])
        os.system('./helmholtz2D.x temp_exp.ini')
        [X, Y, Z, coor, Vexp] = loadND(file_exp+"_U0.dat")

        ModifyDataFile(ini_simu, 'temp_simu.ini', 'Frequency', [str(f), '0.0'])
        ModifyDataFile('temp_simu.ini', 'temp_simu.ini', 'TypeSource', ['SRC_SURFACE', '2', 'GAUSSIAN', str(x), '-5.0', '1.0', '1.5', '1.0'])
        #ModifyDataFile('temp_simu.ini', 'temp_simu.ini', 'TypeSource', ['SRC_SURFACE', '2', 'UNIFORM', '1.0'])
        os.system('./helmholtz2D.x temp_simu.ini')
        [X, Y, Z, coor, Vsimu] = loadND(file_simu+"_U0.dat")
        
        x, y = ReadPointsMesh(file_exp+"_U0_real.mesh")
        PT = numpy.zeros([len(x), 2]); PT[:,0] = x; PT[:,1] = y;
        numpy.savetxt('PointsExp.txt', PT)

        x, y = ReadPointsMesh(file_simu+"_U0_real.mesh")
        PT = numpy.zeros([len(x), 2]); PT[:,0] = x; PT[:,1] = y;
        numpy.savetxt('PointsSimu.txt', PT)

        v = ReadComplexBb(file_exp+"_U0")
        A = numpy.zeros([len(v),1]) + 1j*numpy.zeros([len(v), 1]); A[:,0] = numpy.conj(v)
        write_full('ValExp.dat', A)

        v = ReadComplexBb(file_simu+"_U0")
        A = numpy.zeros([len(v),1]) + 1j*numpy.zeros([len(v), 1]); A[:,0] = numpy.conj(v)
        write_full('ValSimu.dat', A)

        ModifyDataFile('temp_simu.ini', 'temp_simu.ini', 'TypeSource', ['SRC_SURFACE', '2', 'VARIABLE', 'PointsExp.txt', 'ValExp.dat'])
        os.system('./helmholtz2D.x temp_simu.ini')
        [X, Y, Z, coor, Vexp_retro] = loadND(file_simu+"_U0.dat")

        ModifyDataFile('temp_simu.ini', 'temp_simu.ini', 'TypeSource', ['SRC_SURFACE', '2', 'VARIABLE', 'PointsSimu.txt', 'ValSimu.dat'])
        os.system('./helmholtz2D.x temp_simu.ini')
        [X, Y, Z, coor, Vsimu_retro] = loadND(file_simu+"_U0.dat")
        
        image = image + (Vexp_retro-Vsimu_retro)*Vsimu
        #image = image + Vexp_retro*Vsimu
        
os.chdir('cluster/')
plot2dinst(X, Y, abs(image))
