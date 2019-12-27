#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '2/15/2019'
__version__ = '1.0'




R = 8.315e-3   # ideal gas constant, kJ/mol/K 
T = 298.15   # absolute room temperature, K, or 25 C

deftMW = 40   # molecular weight, kDa
deftKcat = 200   # default catalytic rate constant, 1/s
deftKm = 0.2   # default Michaelis constant, mM

deftKcatRelBnds = [0.1, 10]   # relative range of catalytic rate constant
deftKmRelBnds = [0.1, 10]   # relative range of Michaelis constant
deftKeqRelBnds = [1, 10]   # relative range of equilibrium constant

nsteps = 100   # # of integration step
eigThreshold = 1  # threshold of eigenvalues (-1e-6 recommended, if too much system failure in ensemble models, increase gradually to 1 or larger for real values)



