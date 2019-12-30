#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '10/21/2018'
__version__ = '1.0'


import numpy as np
import pandas as pd




def calculate_robustness_index(results, enzymesInner, nsteps):
	'''
	Parameters
	results: dict
	enzymesInner: lst, enzyme IDs with initial and final reaction
	nsteps: int, # of integration steps
	enzymeLBs: ser, lower bounds of enzyme level
	enzymeUBs: ser, upper bounds of enzyme level
		
	Returns
	robustIdx: ser, median of robustness index Si for each enzyme
	'''
	
	from scipy.stats import lognorm
	
	robustIdx = pd.Series(index = enzymesInner)
	for enzyme in robustIdx.index:
	
		Ss = []
		for i in range(len(results[enzyme])):
			
			resulti = results[enzyme][i]
				
			# feasible LB and UB of enzyme level
			Eref = resulti[0].loc[enzyme, resulti[0].columns[0]]
			
			LB = resulti[0].loc[enzyme, resulti[0].columns[-1]]
			UB = resulti[1].loc[enzyme, resulti[1].columns[-1]]
			
			# calculate the probability of maintaining stability
			p = lognorm.cdf(UB, s = 0.5, scale = Eref) - lognorm.cdf(LB, s = 0.5, scale = Eref)   # ln(E) ~ N(ln(Eref), 0.5)
			
			# calculate the robustness index
			if p <= 0: p = 0.0001   
				
			S = -p * np.log(p)	
			#S = p
				
			Ss.append(S)	
				
		robustIdx.loc[enzyme] = np.mean(Ss)
		
	return robustIdx

	
def calculate_system_failure_probability(results, enzymesInner, nsteps, nmodels, enzymeLB, enzymeUB):
	'''
	Parameters
	results: dict
	enzymesInner: lst, enzyme IDs with initial and final reaction
	nsteps: int, # of integration steps
	nmodels: int, # of ensemble models
	enzymeLB: float, lower bound of relative enzyme level
	enzymeUB: float, upper bound of relative enzyme level
	
	Returns
	failurePro: df, probability of system failure, enzyme in rows, enzyme level in columns
	'''

	ERange = np.concatenate((np.linspace(enzymeLB, 1, nsteps + 1), np.linspace(1, enzymeUB, nsteps + 1)[1:]))
	
	failurePro = pd.DataFrame(index = enzymesInner, columns = ERange)	
	for enzyme in failurePro.index:
	
		nmodels = len(results[enzyme])
		
		# count for decreased enzyme level
		for Elevel in failurePro.columns[:nsteps + 1]:
			
			count = 0
			for	i in range(nmodels):

				Elength = results[enzyme][i][0].shape[1]
				feasibleLB = 1 - (Elength - 1) * (1 - enzymeLB) / nsteps
				
				if Elevel >= feasibleLB: count += 1
		
			failurePro.loc[enzyme, Elevel] = 1 - count / nmodels
			
		# count for increased enzyme level
		for Elevel in failurePro.columns[nsteps + 1:]:
			
			count = 0
			for i in range(nmodels):
				
				Elength = results[enzyme][i][1].shape[1]
				feasibleUB = 1 + (Elength - 1) * (enzymeUB - 1) / nsteps
				
				if Elevel <= feasibleUB: count += 1
				
			failurePro.loc[enzyme, Elevel] = 1 - count / nmodels
			
	return failurePro	
	
	
def flux_change_calculation_enzymeDOWN_worker(ifReal, enzyme, enzymes, Smetab2rnx, ensembleModels, Vss, results, fluxRange, ERangeDown, nsteps, enzymeLB, nwindows):
	'''
	Parameters	
	ifReal: str, whether using real values, 'yes' or 'no'	
	enzyme: str, enzyme ID
	enzymes: lst, enzyme IDs
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	ensembleModels: lst
	Vss: ser, fluxes in steady state
	results: dict
	fluxRange: array, range of flux change
	ERangeDown: array, range of enzyme level change (down regulated)
	nsteps: int, # of integration steps
	enzymeLB: float, lower bound of relative enzyme level
	nwindows: int, # of window to get the histogram of flux change. better set a odd number
	
	Returns
	fluxChangeEdown: dict
	'''
	
	from utilities import get_V
	from common_rate_laws import v_expression
	
	fluxChangeEdown = pd.DataFrame(index = range(nwindows), columns = ERangeDown)
			
	nmodels = len(results[enzyme])	

	# stats for decreased enzyme level
	for Elevel in fluxChangeEdown.columns:
		
		colID = int(round((1 - Elevel) * nsteps / (1 - enzymeLB), 0))
		
		# calculate all flux changes	
		fluxChangeThisEnzyme = []
		for i in range(nmodels):
		
			if results[enzyme][i][0].shape[1] < colID + 1: continue    
	
			reverses, kcats, subConcss, subCoess, subKmss, proConcss, proCoess, proKmss, Keqs = ensembleModels[i]	
			
			Enew = results[enzyme][i][0].loc[:, colID]
			
			Xnew = results[enzyme][i][2].loc[:, colID]
			
			Vnew = get_V(Smetab2rnx, v_expression, Enew, Xnew, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
			if ifReal == 'yes': Vnew = Vnew * 3600   # V in mmol/gCDW/h for real values
			Vnew = pd.Series(np.array(Vnew).reshape(-1).astype(np.float), index = enzymes)	
			
			fluxChangeThisEnzyme.append(Vnew[enzyme] / Vss[enzyme])
		
		# get the histogram of flux changes
		fluxChangeEdown.loc[:, Elevel] = np.histogram(fluxChangeThisEnzyme, bins = fluxRange)[0][::-1]   
	
	return fluxChangeEdown
	
	
def flux_change_calculation_enzymeUP_worker(ifReal, enzyme, enzymes, Smetab2rnx, ensembleModels, Vss, results, fluxRange, ERangeUp, nsteps, enzymeUB, nwindows):
	'''
	Parameters	
	ifReal: str, whether using real values, 'yes' or 'no'	
	enzyme: str, enzyme ID
	enzymes: lst, enzyme IDs
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	ensembleModels: lst
	Vss: ser, fluxes in steady state
	results: dict
	fluxRange: array, range of flux change
	ERangeUp: array, range of enzyme level change (up regulated)
	nsteps: int, # of integration steps
	enzymeUB: float, upper bound of relative enzyme level
	nwindows: int, # of window to get the histogram of flux change. better set a odd number
	
	Returns
	fluxChangeEup: dict
	'''
	
	from utilities import get_V
	from common_rate_laws import v_expression
	
	fluxChangeEup = pd.DataFrame(index = range(nwindows), columns = ERangeUp)
			
	nmodels = len(results[enzyme])	
	
	for Elevel in fluxChangeEup.columns:
				
			colID = int(round((Elevel - 1) * nsteps / (enzymeUB - 1), 0))
				
			# calculate all flux changes	
			fluxChangeThisEnzyme = []
			for i in range(nmodels):
			
				if results[enzyme][i][1].shape[1] < colID + 1: continue    
				
				reverses, kcats, subConcss, subCoess, subKmss, proConcss, proCoess, proKmss, Keqs = ensembleModels[i]
				
				Enew = results[enzyme][i][1].loc[:, colID]
				
				Xnew = results[enzyme][i][3].loc[:, colID]
				
				Vnew = get_V(Smetab2rnx, v_expression, Enew, Xnew, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
				if ifReal == 'yes': Vnew = Vnew * 3600   # V in mmol/gCDW/h for real values
				Vnew = pd.Series(np.array(Vnew).reshape(-1).astype(np.float), index = enzymes)	
				
				fluxChangeThisEnzyme.append(Vnew[enzyme] / Vss[enzyme])
				
			# get the histogram of flux changes
			fluxChangeEup.loc[:, Elevel] = np.histogram(fluxChangeThisEnzyme, bins = fluxRange)[0]   
			
	return fluxChangeEup
	
	
def calculate_flux_fold_change(ifReal, Smetab2rnx, ensembleModels, Vss, results, enzymes, enzymesInner, nsteps, enzymeLB, enzymeUB, nprocess, fluxBnds = (0.1, 10), nwindows = 49):
	'''
	Parameters
	ifReal: str, whether using real values, 'yes' or 'no'
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	ensembleModels: lst
	Vss: ser, fluxes in steady state
	results: dict
	enzymes: lst, enzyme IDs
	enzymesInner: lst, enzyme IDs with initial and final reaction
	nsteps: int, # of integration steps
	enzymeLB: float, lower bound of relative enzyme level
	enzymeUB: float, upper bound of relative enzyme level
	nprocess: int, number of processes to run simutaneously
	fluxBnds: 2-tuple, relative bounds of flux change
	nwindows: int, # of window to get the histogram of flux change. better set a odd number, the higher value of nwindows, the higher resolution of figure
	
	Returns
	fluxChange: dict 
	'''
	
	from multiprocessing import Pool
	
	ERangeDown = np.linspace(enzymeLB, 1, nsteps + 1)
	ERangeUp = np.linspace(1, enzymeUB, nsteps + 1)
	fluxRange = np.logspace(np.log10(fluxBnds[0]), np.log10(fluxBnds[1]), nwindows + 1)
	
	# decreased enzyme level
	pool1 = Pool(processes = nprocess)

	fluxChangeEdown = {}
	for enzyme in enzymesInner:
	
		res = pool1.apply_async(func = flux_change_calculation_enzymeDOWN_worker, args = (ifReal, enzyme, enzymes, Smetab2rnx, ensembleModels, Vss, results, fluxRange, ERangeDown, nsteps, enzymeLB, nwindows))
	
		fluxChangeEdown[enzyme] = res
	
	pool1.close()	
	pool1.join()
	
	for enzyme, res in fluxChangeEdown.items(): fluxChangeEdown[enzyme] = res.get()
	
	# increased enzyme level
	pool2 = Pool(processes = nprocess)
	
	fluxChangeEup = {}
	for enzyme in enzymesInner:
	
		res = pool2.apply_async(func = flux_change_calculation_enzymeUP_worker, args = (ifReal, enzyme, enzymes, Smetab2rnx, ensembleModels, Vss, results, fluxRange, ERangeUp, nsteps, enzymeUB, nwindows))
	
		fluxChangeEup[enzyme] = res
	
	pool2.close()	
	pool2.join()
	
	for enzyme, res in fluxChangeEup.items(): fluxChangeEup[enzyme] = res.get()
	
	# combine data
	fluxChange = {}
	for enzyme in enzymesInner:
	
		fluxChange[enzyme] = pd.concat([fluxChangeEdown[enzyme], fluxChangeEup[enzyme].iloc[:, 1:]], axis = 1)

	return fluxChange	
	
	
def calculate_flux_control_index(fluxChange, enzymes, fluxBnds = (0.1, 10)):
	'''
	Parameters
	fluxChange: dict
	enzymes: lst, enzyme IDs
	fluxBnds: 2-tuple, relative bounds of flux change
	
	Returns
	ConIdx: df, enyzmes in rows
	'''
	
	ConIdx = pd.DataFrame(index = enzymes, columns = ['Down regulation', 'Up regulation'])   

	for enzyme in enzymes:
		
		fluxChangeThisEnzyme = fluxChange[enzyme].copy()
		
		fluxChangeThisEnzyme.columns = fluxChangeThisEnzyme.columns.astype('float')
		fluxChangeThisEnzyme.index = np.logspace(np.log10(fluxBnds[1]), np.log10(fluxBnds[0]), fluxChangeThisEnzyme.index.size)

		steps = fluxChangeThisEnzyme.columns.size

		# decreased enzyme level
		fluxChangeEdown = fluxChangeThisEnzyme.iloc[:, :(steps-1)//2]
		
		ConIdxEdown = []
		for v in fluxChangeEdown.index:
			for e in fluxChangeEdown.columns:
				ConIdxEdown.extend([np.log10(v) / np.log10(e)] * fluxChangeEdown.loc[v, e])
		
		ConIdx.loc[enzyme:enzyme, 'Down regulation'] = [ConIdxEdown]
		
		# increased enzyme level
		fluxChangeEup = fluxChangeThisEnzyme.iloc[:, (steps+1)//2:]
		
		ConIdxEup = []
		for v in fluxChangeEup.index:
			for e in fluxChangeEup.columns:
				ConIdxEup.extend([np.log10(v) / np.log10(e)] * fluxChangeEup.loc[v, e])

		ConIdx.loc[enzyme:enzyme, 'Up regulation'] = [ConIdxEup]

	return ConIdx
	
	
	
	
	
	
	
	
	
	
	
