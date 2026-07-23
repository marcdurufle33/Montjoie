#!/usr/bin/env python
import os, string, sys, math, array, test_regression

for num_simul in range(0,120):
    command_line = test_regression.SetSimulation(num_simul)
    os.system(command_line)
    test_regression.CheckSimulation(num_simul)
    
