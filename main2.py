#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '2/19/2019'
__version__ = '1.0'


'''
This script performs robustness analysis: 1 calculating robustness index for enzymes; 2 calculating probability of system failure under enzyme perturbation; 3 calculating flux fold change under enzyme perturbation
'''


import argparse
import os
import re
import numpy as np
import pandas as pd
from constants import nsteps
from parse_network import parse_network, get_full_stoichiometric_matrix, get_steady_state_net_fluxes
from ensemble_models import generate_ensemble_models, simulate_perturbation




if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'This script does thermodynamic analysis: 1 maximizing the minimal driving force; 2 minimizing the totol enzyme protein cost of a given pathway')
	parser.add_argument('-o', '--outDir', type = str, required = True, help = 'output directory')
	parser.add_argument('-r', '--reactionFile', type = str, required = True, help = 'reaction list file, required fields: Enzyme ID, Reversibility, ΔrGm, Substrates and Products')
	parser.add_argument('-i', '--iniMetabs', type = str, required = False, help = 'metabolites as initial substrates, sep by ",". By deafult, they will be detected automatically, sometimes they should be set explicitly, e.g. for cylic pathways')
	parser.add_argument('-f', '--finMetabs', type = str, required = False, help = 'metabolites as end products, sep by ",". By deafult, they will be detected automatically, sometimes they should be set explicitly, e.g. for cylic pathways')
	parser.add_argument('-eb', '--exBalMetabs', type = str, required = False, help = 'metabolites excluded from mass balance, sep by ","')
	parser.add_argument('-eo', '--exOptMetabs', type = str, required = False, help = 'metabolites excluded from optimization, sep by ","')
	parser.add_argument('-n', '--nmodels', type = int, required = True, help = 'number of models in an ensemble')
	parser.add_argument('-b', '--enzymeBnds', type = str, required = True, help = 'lower and upper bound of relative enzyme level, sep by ","')
	parser.add_argument('-d', '--ifDump', type = str, required = True, choices = ['yes', 'no'], help = "whether to dump generated models, 'yes' or 'no'")
	parser.add_argument('-w', '--runWhich', type = str, required = True, help = "which analysis to run, '1' for robustmess index, '2' for probability of system failure, '3' for flux fold change, '12', '23', ... for combinations")
	parser.add_argument('-p', '--nprocess', type = int, required = True, help = "number of processes to run simultaneously")
	parser.add_argument('-t', '--ifReal', action = 'store_true', required = True, help = "whether to use the real value of concentrations, Kms and Keqs, 'yes' or 'no'")
	subparsers = parser.add_subparsers(dest = 'ifReal')
	parser_yes = subparsers.add_parser('yes')
	parser_yes.add_argument('-a', '--assignFlux', type = str, required = True, help='assign flux (mmol/gCDW/h) to some enzyme in the format "enzyme ID:value", then flux distribution of reference state will be calculated')	
	parser_yes.add_argument('-mc', '--metabConcFile', type = str, required = True, help = 'file of metabolite concentrations (mM) in reference state')
	parser_yes.add_argument('-ec', '--enzConcFile', type = str, required = True, help = 'file of enzyme concentrations (mmol/gCDW) in reference state')
	parser_no = subparsers.add_parser('no')
	args = parser.parse_args()
	
	outDir = args.outDir
	reactionFile = args.reactionFile
	iniMetabs = args.iniMetabs
	finMetabs = args.finMetabs
	exBalMetabs = args.exBalMetabs
	exOptMetabs = args.exOptMetabs
	nmodels = args.nmodels
	enzymeBnds = args.enzymeBnds
	ifDump = args.ifDump
	runWhich = args.runWhich
	nprocess = args.nprocess
	ifReal = args.ifReal
	if ifReal == 'yes':
		assignFlux = args.assignFlux
		metabConcFile = args.metabConcFile
		enzConcFile = args.enzConcFile
		
	
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
	
	S4BalFull = get_full_stoichiometric_matrix(S4Bal, metabInfo)   

	# get flux distribution in steady state
	if ifReal == 'yes':
		speEnz, speFlux = assignFlux.split(':')
		speFlux = float(speFlux)
		
		Vss = get_steady_state_net_fluxes(S4BalFull, enzymeInfo, metabInfo, speEnz, speFlux)
		
	else:
		Vss = get_steady_state_net_fluxes(S4BalFull, enzymeInfo, metabInfo)
	
	print('\nDone.')
	
	
	## generate ensemble models and simulate perturbation ---------------------------------------------------
	print('\n\nGenerating ensemble models')
	print('.' * 50)
	
	# generate ensemble models

	S4OptFull = get_full_stoichiometric_matrix(S4Opt, metabInfo)   
	
	if ifReal == 'yes':
		from parse_network import read_concentrations
		
		Css = read_concentrations(metabConcFile)
		Ess = read_concentrations(enzConcFile)
		
		EssMean = Ess.mean()
		for enzyme in S4OptFull.columns: 
			Ess.loc[enzyme] = Ess.get(enzyme, EssMean)   
		
		ensembleModels = generate_ensemble_models(S4OptFull, enzymeInfo, Vss, nmodels, Ess, Css)

	else:	
		ensembleModels = generate_ensemble_models(S4OptFull, enzymeInfo, Vss, nmodels)
		
		
	# simulate perturbation (estimate metabolite concentrations at different enzyme levels)
	metabs = S4OptFull.index
	enzymes = S4OptFull.columns
	
	Smetab2rnx = S4OptFull.T / S4OptFull.T.abs()
	Smetab2rnx = Smetab2rnx.replace(np.nan, 0)
	
	enzymeLB, enzymeUB = map(float, enzymeBnds.split(','))
	
	if ifReal == 'yes':
		enzymeLBs = Ess * enzymeLB
		enzymeUBs = Ess * enzymeUB
	
		pertResults = simulate_perturbation(ensembleModels, S4OptFull, Smetab2rnx, enzymes, metabs, nsteps, enzymeLBs, enzymeUBs, nmodels, nprocess, Ess, Css)
	
	else:
		enzymeLBs = pd.Series(np.full(len(enzymes), enzymeLB), index = enzymes)
		enzymeUBs = pd.Series(np.full(len(enzymes), enzymeUB), index = enzymes)
		
		pertResults = simulate_perturbation(ensembleModels, S4OptFull, Smetab2rnx, enzymes, metabs, nsteps, enzymeLBs, enzymeUBs, nmodels, nprocess)	
	
			
	if ifDump == 'yes':
		from output import dump_ensemble_models
			
		dump_ensemble_models(pertResults, outDir)	
	
	print('\nDone.')
	
	
	## estimate robustmess ----------------------------------------------------------------------------------
	innerEnzymes = [enz for enz in enzymes if not re.match(r'.+_(in|out)', enz)]
		
	# calculate robustness index
	if re.search(r'1', runWhich):
		
		from robustness import calculate_robustness_index
		from output import plot_robustness_index, save_robustness_index
		
		print('\n\nCalculating robustness index')
		print('.' * 50)
		
		# calculate robustness index
		robustIdx = calculate_robustness_index(pertResults, innerEnzymes, nsteps)
		
		# output results	
		plot_robustness_index(robustIdx, outDir)
		save_robustness_index(robustIdx, outDir)	
		
		print('\nDone.')	
	
	# calculate probability of system failure under enzyme perturbation
	if re.search(r'2', runWhich):	
		
		from robustness import calculate_system_failure_probability
		from output import plot_system_failure_probability, save_system_failure_probability
		
		print('\n\nCalculating probability of system failure')
		print('.' * 50)
		
		# calculate probability of system failure
		failurePro = calculate_system_failure_probability(pertResults, innerEnzymes, nsteps, nmodels, enzymeLB, enzymeUB)
		
		# output results	
		plot_system_failure_probability(failurePro, outDir)
		save_system_failure_probability(failurePro, outDir)	
			
		print('\nDone.')
	
	# calculate flux fold change under enzyme perturbation
	if re.search(r'3', runWhich):
		
		from robustness import calculate_flux_fold_change, calculate_flux_control_index
		from output import plot_flux_fold_change, plot_flux_control_index, save_flux_fold_change
		
		print('\n\nCalculating flux fold change')
		print('.' * 50)
			
		fluxChangeBnds = (0.2, 5)   # may need to set for plot
		
		# calculate flux fold change and flux control index
		fluxChange = calculate_flux_fold_change(ifReal, Smetab2rnx, ensembleModels, Vss, pertResults, enzymes, innerEnzymes, nsteps, enzymeLB, enzymeUB, nprocess, fluxBnds = fluxChangeBnds)
		
		fluxConIdx = calculate_flux_control_index(fluxChange, innerEnzymes, fluxBnds = fluxChangeBnds)	

		# output results	
		plot_flux_fold_change(innerEnzymes, fluxChange, outDir, fluxBndsShow = fluxChangeBnds)
		plot_flux_control_index(fluxConIdx, outDir)
		save_flux_fold_change(innerEnzymes, fluxChange, outDir)	
		
		print('\nDone.')
	
	
	
	





