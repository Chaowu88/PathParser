#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '2/15/2019'
__version__ = '1.0'


'''
This script performs thermodynamics analysis: 1 maximizing the minimal driving force; 2 minimizing the totol enzyme protein cost of a given pathway
'''


import argparse
import os
import re
from parse_network import parse_network, get_full_stoichiometric_matrix, get_steady_state_net_fluxes




if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'This script does thermodynamic analysis: 1 maximizing the minimal driving force; 2 minimizing the totol enzyme protein cost of a given pathway')
	parser.add_argument('-o', '--outDir', type = str, required = True, help = 'output directory')
	parser.add_argument('-r', '--reactionFile', type = str, required = True, help = 'reaction list file, required fields: Enzyme ID, Reversibility, ΔrGm, Substrates, Products, and Enzyme MW')
	parser.add_argument('-i', '--iniMetabs', type = str, required = False, help = 'metabolites as initial substrates, sep by ",". By default, they will be detected automatically, sometimes they should be set explicitly, e.g. for cylic pathways')
	parser.add_argument('-f', '--finMetabs', type = str, required = False, help = 'metabolites as end products, sep by ",". By default, they will be detected automatically, sometimes they should be set explicitly, e.g. for cylic pathways')
	parser.add_argument('-eb', '--exBalMetabs', type = str, required = False, help = 'metabolites excluded from mass balance, sep by ","')
	parser.add_argument('-eo', '--exOptMetabs', type = str, required = False, help = 'metabolites excluded from optimization, sep by ","')
	parser.add_argument('-b', '--concBnds', type = str, required = True, help='concentration lower and upper bound (mM) for all metabolites, sep by ","')
	parser.add_argument('-a', '--assignFlux', type = str, required = False, help='assign flux to some enzyme in the format "enzyme ID:value", then flux distribution will be calculated. If not assigned, influx to pathway will be set to 1, flux distribution can also be calculated. NOTE the calculated flux distribution is equivalent to occurance when computing protein cost')
	parser.add_argument('-w', '--runWhich', type = str, required = True, help = "which analysis to run, '1' for maximizing the minimal driving force, '2' for minimizing the totol enzyme protein cost, '12' for both")
	args = parser.parse_args()
	
	outDir = args.outDir
	reactionFile = args.reactionFile
	iniMetabs = args.iniMetabs
	finMetabs = args.finMetabs
	exBalMetabs = args.exBalMetabs
	exOptMetabs = args.exOptMetabs
	concBnds = args.concBnds
	assignFlux = args.assignFlux
	runWhich = args.runWhich
	
	os.makedirs(outDir, exist_ok = True)

	
	## get stoichiometric matrix ---------------------------------------------------------------------------
	print('\n\nParsing network')
	print('.' * 50)	
	
	# get the stoichiometric matrix of a pathway
	iniMetabs = iniMetabs.split(',') if iniMetabs else []
	finMetabs = finMetabs.split(',') if finMetabs else []
	exBalMetabs = exBalMetabs.split(',') if exBalMetabs else []
	exOptMetabs = exOptMetabs.split(',') if exOptMetabs else []
	
	S4Bal, S4Opt, enzymeInfo, metabInfo = parse_network(reactionFile, iniMetabs, finMetabs, exBalMetabs, exOptMetabs)
	
	S4BalFull = get_full_stoichiometric_matrix(S4Bal, metabInfo)   # S4BalFull also includes input and output reactions of the pathway
	
	# get flux distribution in steady state
	if assignFlux and not re.search(r'1', runWhich):
		speEnz, speFlux = assignFlux.split(':')
		speFlux = float(speFlux)
		
		Vss = get_steady_state_net_fluxes(S4BalFull, enzymeInfo, metabInfo, speEnz, speFlux)
		
	else:
		Vss = get_steady_state_net_fluxes(S4BalFull, enzymeInfo, metabInfo)
	
	print('\nDone.')
	
	
	## thermodynamic analysis ------------------------------------------------------------------------------
	# maximize minimal driving force
	if re.search(r'1', runWhich):
	
		from thermodynamics import optimize_minimal_driving_force
		from output import print_driving_force_optimization_results, plot_cumulative_deltaGs, save_driving_force_optimization_results	
			
		print('\n\nMaximize minimal driving force')
		print('.' * 50)

		# maximize minimal driving force
		concLB, concUB = map(float, concBnds.split(','))
			
		optConcs, optDeltaGs, refDeltaGs = optimize_minimal_driving_force(S4Opt, enzymeInfo, concLB, concUB)
			
		# output results
		print_driving_force_optimization_results(optConcs, optDeltaGs, refDeltaGs)
			
		plot_cumulative_deltaGs(optDeltaGs, refDeltaGs, outDir)
		
		save_driving_force_optimization_results(optConcs, optDeltaGs, refDeltaGs, outDir)
	
		print('\nDone.')
	
	
	# minimize enzyme cost	
	if re.search(r'2', runWhich):
		
		from parse_network import get_full_stoichiometric_matrix, get_steady_state_net_fluxes
		from thermodynamics import optimize_enzyme_cost
		from output import print_enzyme_cost_optimization_results, plot_enzyme_costs, save_enzyme_cost_optimization_results
	
		print('\n\nMinimizing enzyme cost')
		print('.' * 50)
		
		# minimize enzyme cost
		concLB, concUB = map(float, concBnds.split(','))
			
		optConcs, optEnzyCosts, optEnzyCostTotal = optimize_enzyme_cost(S4Opt, Vss, enzymeInfo, concLB, concUB)	
			
		# output results
		print_enzyme_cost_optimization_results(optConcs, optEnzyCosts, optEnzyCostTotal)
			
		plot_enzyme_costs(optEnzyCosts, outDir)
		
		save_enzyme_cost_optimization_results(optConcs, optEnzyCosts, optEnzyCostTotal, outDir)

		print('\nDone.')
	


		
		
		
		
