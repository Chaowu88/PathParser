#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '2/14/2019'
__version__ = '1.0'


import numpy as np
import pandas as pd
	
	


def parse_network(reactionFile, iniMetabs = [], finMetabs = [], exBalMetabs = [], exOptMetabs = []):
	'''
	Parameters
	reactionFile: str, reaction list file
	iniMetabs: lst, metabolites as initial substrates
	finMetabs: lst, metabolites as end products
	exBalMetabs: lst, metabolites excluded from mass balance
	exOptMetabs: lst, metabolites excluded from optimization
	
	Returns
	S4Bal: df, stoichiometric matrix for mass balance, metabolite in rows, reaction in columns (same order with S). negative for substrates, positive for products
	S4Opt: df, stoichiometric matrix for optimization, metabolite in rows, reaction in columns (same order with S). negative for substrates, positive for products
	enzymeInfo: df, reaction in rows (same order with S)
	metabInfo: ser, metabolite in index (same order with S), values are 0 (inner metabolite) or -1 (initial substrate) or 1 (final product)
	'''
		
	import re	
	from constants import R, T, deftMW, deftKcat, deftKm, deftKcatRelBnds, deftKmRelBnds, deftKeqRelBnds
	
	
	def parse_reactantStr(reaStr):
		'''
		Parameters
		reaStr: str

		Returns
		coes: lst, reactant coefficients
		reas: lst, reactants
		'''
		
		coe_reas = re.findall(r'(\d+\.?\d*|)(\w+)', reaStr)
	
		coes = [1.0 if not i[0] else float(i[0]) for i in coe_reas]
		reas = [i[1] for i in coe_reas]
	
		return coes, reas
		
	def parse_kineticStr(kinStr, deftVal, deftRelBnds):
		'''
		Parameters
		kinStr: str
		deftVal: float, default values
		deftRelBnds: lst, default lower and upper bound

		Returns
		infos: lst of tuple (value, lb, ub)
		'''
			
		infosRaw = [re.match(r'(\d+\.?\d*)?\(?(\d+\.?\d*)?,?(\d+\.?\d*)?\)?', item).groups() for item in kinStr.split(';')]
		infos = []
		for val, lb, ub in infosRaw:
	
			val = float(val) if val else deftVal
			lb = float(lb) if lb else val * deftRelBnds[0]
			ub = float(ub) if ub else val * deftRelBnds[1]
			
			infos.append([val, lb, ub])
		
		return infos	
		
		
	inputs = pd.read_csv(reactionFile, sep = '\t', header = None, index_col = 0, names = ['id', 'rev', 'deltaGm', 'subs', 'pros', 'subsKm', 'prosKm', 'kcat', 'MW'], comment = '#', na_filter = False, dtype = str)
	
	
	# get S and enzymeInfo
	S = pd.DataFrame(columns = inputs.index)
	
	for enzyme in inputs.index:	
		
		subsStr, prosStr = inputs.loc[enzyme, 'subs':'pros']
		
		coes, subs = parse_reactantStr(subsStr)
		for sub, coe in zip(subs, coes): S.loc[sub, enzyme] = -coe
		
		coes, pros = parse_reactantStr(prosStr)	
		for pro, coe in zip(pros, coes): S.loc[pro, enzyme] = coe
		
	S = S.replace(np.nan, 0)
	
	
	# get S4Bal and S4Opt
	balMetabs = [metab for metab in S.index if metab not in exBalMetabs]
	S4Bal = S.loc[balMetabs, :]
	
	optMetabs = [metab for metab in S.index if metab not in exOptMetabs]
	S4Opt = S.loc[optMetabs, :]
	
	
	# get enzymeInfo
	enzymeInfo = pd.DataFrame(index = inputs.index, columns = ['rev', 'Keq', 'subsKm', 'prosKm', 'kcat', 'MW'], dtype = object)
	
	for enzyme in enzymeInfo.index:
	
		rev, deltaGm, subsStr, prosStr, subsKmStr, prosKmStr, kcatStr, MW = inputs.loc[enzyme, ['rev', 'deltaGm', 'subs', 'pros', 'subsKm', 'prosKm', 'kcat', 'MW']]
		
		enzymeInfo.loc[enzyme, 'rev'] = float(rev)
			
		Keq = np.exp(-float(deltaGm) / R / T)
		enzymeInfo.loc[enzyme:enzyme, 'Keq'] = [[Keq, Keq * deftKeqRelBnds[0], Keq * deftKeqRelBnds[1]]]
		
		subs = parse_reactantStr(subsStr)[1]
		subKmInfos = parse_kineticStr(subsKmStr, deftKm, deftKmRelBnds)
		subKms = pd.Series(index = S.index, dtype = object)
		subKms.loc[subs] = subKmInfos   
		enzymeInfo.loc[enzyme:enzyme, 'subsKm'] = [subKms]   
		
		pros = parse_reactantStr(prosStr)[1]
		proKmInfos = parse_kineticStr(prosKmStr, deftKm, deftKmRelBnds)
		proKms = pd.Series(index = S.index, dtype = object)
		proKms.loc[pros] = proKmInfos   
		enzymeInfo.loc[enzyme:enzyme, 'prosKm'] = [proKms]   
		
		kcatInfos = parse_kineticStr(kcatStr, deftKcat, deftKcatRelBnds)
		enzymeInfo.loc[enzyme:enzyme, 'kcat'] = kcatInfos   
		
		enzymeInfo.loc[enzyme, 'MW'] = float(MW) if MW else deftMW
	
	
	# get metabInfo
	iniSubs = S.loc[np.all(S <= 0, axis = 1), :].index   
	finPros = S.loc[np.all(S >= 0, axis = 1), :].index   
	
	innerMetabs = set(S.index) - set(iniSubs) - set(finPros)
	
	metabInfo = pd.Series(dict(list(dict.fromkeys(iniSubs, -1).items()) +   
	                           list(dict.fromkeys(finPros, 1).items()) + 
					           list(dict.fromkeys(innerMetabs, 0).items()) +
					           list(dict.fromkeys(iniMetabs, -1).items()) + 
							   list(dict.fromkeys(finMetabs, 1).items())))
							   
	
	return S4Bal, S4Opt, enzymeInfo, metabInfo
			
		
def get_full_stoichiometric_matrix(S, metabInfo):
	'''
	Parameters
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns. negative for substrates, positive for products
	metabInfo: ser, metabolite in index, values are 0 (inner metabolite) or -1 (initial substrate) or 1 (final product)
	
	Returns
	SFull: same with S, including input and output reactions of the pathway
	NOTE input and output reactions all in the form X -> X
	'''
	
	iniSubs = [metab for metab in S.index if metabInfo.loc[metab] == -1]
	finPros = [metab for metab in S.index if metabInfo.loc[metab] == 1]
	
	# make new S to include input and output reactions of the pathway
	SFull = S.copy(deep = 'all')
	
	for iniSub in iniSubs:
		SFull.loc[iniSub, iniSub+'_in'] = 1.0
		
	for finPro in finPros:	
		SFull.loc[finPro, finPro+'_out'] = -1.0
		
	SFull = SFull.replace(np.nan, 0)	
	
	return SFull


def get_steady_state_net_fluxes(S, enzymeInfo, metabInfo, speEnz = None, speFlux = None):		
	'''
	Parameters
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns. negative for substrates, positive for products 
	enzymeInfo: df, reaction in rows
	metabInfo: ser, metabolite in index, values are 0 (inner metabolite) or -1 (initial substrate) or 1 (final product)
	speEnz: str, enzyme specified as initial enzyme
	speFlux: float, flux specified as initial flux
	
	Returns
	Vss: ser, net fluxes in steady state (including in and out fluxes)	
	'''
	
	import re
	from scipy.optimize import lsq_linear
	
	# make new S to solve all net fluxes
	Snew = S.copy(deep = 'all')
	
	if speEnz:
		Snew.loc[speEnz, speEnz] = 1.0
		
	else:	
		iniSubs = [metab for metab in S.index if metabInfo.loc[metab] == -1]	
			
		Snew.loc[iniSubs[0] + '_in', iniSubs[0] + '_in'] = 1.0	

	Snew = Snew.replace(np.nan, 0)	
	
	# solve all net fluxes
	beq = np.zeros(Snew.shape[0])
		
	if speEnz:
		beq[-1] = speFlux
	
	else:		
		beq[-1] = 1.0   
	
	lbs = [0 if (r not in enzymeInfo.index or enzymeInfo.loc[r, 'rev'] == 0) else -np.inf for r in Snew.columns]   
	
	bnds = (lbs, np.inf)
	
	fluxes = lsq_linear(Snew.values, beq, bounds = bnds).x   
	
	Vss = pd.Series(fluxes.reshape(-1), index = Snew.columns)
	
	return Vss
	
	
def read_concentrations(concFile):
	'''
	Parameters	
	concFile: str, metabolite concentration file
	
	Returns
	Css: ser, metabolite concentration in steady state
	'''
	
	Concs = pd.read_csv(concFile, sep = '\t', squeeze = True, header = None, index_col = 0, comment = '#')
	
	return Concs
	
	
	
	
	
	
	
	
	
	
	
	
	
