#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '10/19/2018'
__version__ = '1.0'




def generate_ensemble_models(S, enzymeInfo, Vss, nmodels, Ess = [], Css = []):
	'''
	Parameters
	S: stoichiometric matrix, metabolite in rows, reaction in columns (including input and output reactions)
	enzymeInfo: df, reaction in rows (same order with S)
	Vss: ser, fluxes in steady state
	nmodels: int, number of ensemble models
	Css: ser, metabolite concentration in steady state
	Ess: ser, enzyme concentration in steady state
	
	Returns
	ensembleModels: lst
	NOTE kcats, ... are in order of enzymes
	NOTE product concs, coes, Kms = [] and Keqs = 0 for irreversible and output reactions
	'''
		
	import numpy as np
	import re
	from numpy.random import rand
	from constants import deftKm, deftKmRelBnds, deftKeqRelBnds
	from common_rate_laws import v_expression
		
	ensembleModels = []   
	for i in range(nmodels):
	
		reverses, kcats, subConcss, subCoess, subKmss, proConcss, proCoess, proKmss, Keqs = [], [], [], [], [], [], [], [], []
		for enzyme in S.columns:
			
			subs = S.index[S.loc[:, enzyme] < 0]
			pros = S.index[S.loc[:, enzyme] > 0]
			reverse = 0 if enzyme not in enzymeInfo.index else enzymeInfo.loc[enzyme, 'rev']   
			
			# generate random Kms and Keq, metabolite concentrations set 1, then calculate kcat
			if subs.size == 0: # no substrate indicates a input reaction (in form of X_in -> X) 
				if len(Css) > 0:    	
					KmLBs = deftKm * deftKmRelBnds[0]
					KmUBs = deftKm * deftKmRelBnds[1]
				
				else:   
					KmLBs = deftKmRelBnds[0]
					KmUBs = deftKmRelBnds[1]
					
				subConcs = np.ones(1)
				subCoes = np.ones(1)
				subKms = np.power(10, np.log10(KmLBs) + (np.log10(KmUBs) - np.log10(KmLBs)) * rand(1))
				
			elif pros.size == 0:   # no product indicates a output reaction (in form of X -> X_out) 	
				if len(Css) > 0:   
					subConcs = Css.loc[subs].values
				
					KmLBs = deftKm * deftKmRelBnds[0]
					KmUBs = deftKm * deftKmRelBnds[1] 

				else:   
					subConcs = np.ones(1)
					
					KmLBs = deftKmRelBnds[0]
					KmUBs = deftKmRelBnds[1]
					
				subCoes = np.ones(1)
				subKms = np.power(10, np.log10(KmLBs) + (np.log10(KmUBs) - np.log10(KmLBs)) * rand(1))
				
			else:
				if len(Css) > 0:   
					subConcs = Css.loc[subs].values
					
					KmLBs = np.array([item[1] for item in enzymeInfo.loc[enzyme, 'subsKm'].loc[subs]])
					KmUBs = np.array([item[2] for item in enzymeInfo.loc[enzyme, 'subsKm'].loc[subs]])
				
				else:   	
					subConcs = np.ones(len(subs))
					
					KmLBs = deftKmRelBnds[0]
					KmUBs = deftKmRelBnds[1]
				
				subKms = np.power(10, np.log10(KmLBs) + (np.log10(KmUBs) - np.log10(KmLBs)) * rand(len(subs)))	
				subCoes = S.loc[subs, enzyme].abs().values
			
			if reverse:
				if len(Css) > 0:   	
					proConcs = Css.loc[pros].values
					
					KmLBs = np.array([item[1] for item in enzymeInfo.loc[enzyme, 'prosKm'].loc[pros]])
					KmUBs = np.array([item[2] for item in enzymeInfo.loc[enzyme, 'prosKm'].loc[pros]])
					
					KeqLBs = enzymeInfo.loc[enzyme, 'Keq'][1]
					KeqUBs = enzymeInfo.loc[enzyme, 'Keq'][2]
				
				else:   
					proConcs = np.ones(len(pros))
					
					KmLBs = deftKmRelBnds[0]
					KmUBs = deftKmRelBnds[1]
					
					KeqLBs = deftKeqRelBnds[0] 
					KeqUBs = deftKeqRelBnds[1]	
					
				proKms = np.power(10, np.log10(KmLBs) + (np.log10(KmUBs) - np.log10(KmLBs)) * rand(len(pros)))	
				proCoes = S.loc[pros, enzyme].abs().values   
				Keq = np.power(10, np.log10(KeqLBs) + (np.log10(KeqUBs) - np.log10(KeqLBs)) * rand())
				
			else:	
				proConcs = []
				proCoes = []
				proKms = []
				
				Keq = 0
			
			if len(Css) > 0:   
				kcat = Vss[enzyme] / Ess[enzyme] / v_expression(reverse, subConcs, subCoes, subKms, proConcs, proCoes, proKms, Keq) / 3600   # NOTE V in mmol/gCDW/h, E in mmol/gCDW, kcat should be in 1/s

			else:   	
				kcat = Vss[enzyme] / v_expression(reverse, subConcs, subCoes, subKms, proConcs, proCoes, proKms, Keq)   # E = 1 in reference state for relative values	
			
			reverses.append(reverse)
			kcats.append(kcat)
			subConcss.append(subConcs)
			subCoess.append(subCoes)
			subKmss.append(subKms)
			proConcss.append(proConcs)
			proCoess.append(proCoes)
			proKmss.append(proKms)
			Keqs.append(Keq)
		
		ensembleModels.append([reverses, kcats, subConcss, subCoess, subKmss, proConcss, proCoess, proKmss, Keqs])	
	
	return ensembleModels


def simulation_worker(i, ensembleModel, S, Smetab2rnx, E, Eini, X, Xini, enzymes, nsteps, enzymeLBs, enzymeUBs):
	'''
	Parameters
	i: int, model #
	ensembleModel: lst
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	E: sym array, enzyme concentrations, in order of enzymes
	Eini: array, initial enzyme concentrations, in order of enzymes
	X: sym array, metabolites concentrations, in order of metabs
	Xini: array, initial metabolites concentrations, in order of metabs
	enzymes: lst, enzyme IDs
	nsteps: int, # of integration steps
	enzymeLBs: ser, lower bounds of enzyme level
	enzymeUBs: ser, upper bounds of enzyme level
	
	Returns
	resultPerModel: dict
	'''
	
	import numpy as np
	import pandas as pd
	from scipy.linalg import eigvals
	from constants import eigThreshold
	from utilities import get_Jacobian, get_dVdE, solve_dXdE, get_lambdify_function
	from common_rate_laws import v_expression
	import platform
	if platform.system() == 'Linux':
		import os
		os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
		
	print('\nprocessing model %s ...' % (i + 1))
	
	reverses, kcats, subConcss, subCoess, subKmss, proConcss, proCoess, proKmss, Keqs = ensembleModel
	
	# calculate the Jacobian matrix of reference state and keep those model with all Jacobian eigenvalues real parts < 0
	J = get_Jacobian(S, Smetab2rnx, v_expression, Eini, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)	
	Jlam = get_lambdify_function(X, J)	
	
	Jss = np.matrix(Jlam(*Xini)).astype(np.float)
	
	if np.any(eigvals(Jss).real >= eigThreshold):
		print('Jacobian matrix singular, model abandoned')
		return
		
	# solve ODE to get relation of X ~ E
	J = get_Jacobian(S, Smetab2rnx, v_expression, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
	dVdE = get_dVdE(Smetab2rnx, v_expression, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
	
	XE = np.concatenate((X, E))
	Jlam = get_lambdify_function(XE, J)
	dVdElam = get_lambdify_function(XE, dVdE)
	
	resultPerModel = {}
	for enzyme in enzymes:
		
		# enzyme concentration increase
		Espan1 = pd.DataFrame(np.array([Eini, Eini]).T, index = enzymes)
		Espan1.loc[enzyme, 1] = enzymeUBs.loc[enzyme]
		
		Eout1, Xout1 = solve_dXdE(Espan1, nsteps, Xini, Jlam, dVdElam, S)
	
		Eout1 = Eout1.dropna(axis = 1)
		Xout1 = Xout1.dropna(axis = 1)

		# enzyme concentration decrease
		Espan2 = pd.DataFrame(np.array([Eini, Eini]).T, index = enzymes)
		Espan2.loc[enzyme, 1] = enzymeLBs.loc[enzyme]
		
		Eout2, Xout2 = solve_dXdE(Espan2, nsteps, Xini, Jlam, dVdElam, S)

		Eout2 = Eout2.dropna(axis = 1)
		Xout2 = Xout2.dropna(axis = 1)
	
		resultPerModel[enzyme] = [Eout2, Eout1, Xout2, Xout1]
	
	return resultPerModel
	
	
def simulate_perturbation(ensembleModels, S, Smetab2rnx, enzymes, metabs, nsteps, enzymeLBs, enzymeUBs, nmodels, nprocess, Eini = [], Xini = []):
	'''
	Parameters
	ensembleModels: lst
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	enzymes: lst, enzyme IDs
	metabs: lst, metabolite IDs
	nsteps: int, # of integration steps
	enzymeLBs: ser, lower bounds of enzyme level
	enzymeUBs: ser, upper bounds of enzyme level
	nmodels: int, number of ensemble models
	nprocess: int, number of processes
	Eini: ser, initial enzyme concentrations, if real values used
	Xini: ser, initial enzyme concentrations if real values used
	
	Returns
	results: dict
	'''
	
	import numpy as np	
	from sympy import symbols	
	from multiprocessing import Pool
	
	X = np.array(symbols(' '.join(metabs)))
	E = np.array(symbols(' '.join(enzymes)))   
		
	if len(Eini) > 0:   
		Xini = Xini.loc[metabs]   
		Eini = Eini.loc[enzymes]   
		
		ifReal = 'yes'
	
	else:
		Xini = np.ones(len(metabs))	
		Eini = np.ones(len(enzymes))
		
		ifReal = 'no'
	
	# multiprocessing
	pool = Pool(processes = nprocess)	
		
	tmp = []   
	for i in range(nmodels):

		res = pool.apply_async(func = simulation_worker, args = (i, ensembleModels[i], S, Smetab2rnx, E, Eini, X, Xini, enzymes, nsteps, enzymeLBs, enzymeUBs))
		
		tmp.append(res)
		
	pool.close()
	pool.join()
	
	for i, res in enumerate(tmp): tmp[i] = res.get()
	
	# get results
	results = {enzyme: [] for enzyme in enzymes}
	for i in range(nmodels): 
		if tmp[i]:
			for enzyme in enzymes:
				results[enzyme].append(tmp[i][enzyme])
	
	return results





