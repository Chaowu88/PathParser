#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '10/18/2018'
__version__ = '1.0'


import numpy as np
import pandas as pd




def optimize_minimal_driving_force(S, Vss, enzymeInfo, concLB, concUB):
	'''
	Parameters
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns. negative for substrates, positive for products
	Vss: ser, net fluxes in steady state (including in and out fluxes)
	enzymeInfo: df, reaction in rows
	concLB: float, concentration lower bound (mM) for all metabolites
	concUB: float, concentration upper bound (mM) for all metabolites
	
	Returns
	optConcs: ser, optimal log(concentrations)
	optDeltaGs: ser, optimal minimal driving forces
	refDeltaGs: float, reference minimal driving forces (all concentrations at 1 mM)
	'''

	from openopt import LP
	from constants import R, T

	f = np.zeros(S.shape[0] + 1)
	f[0] = -1
	
	S = S * Vss[S.columns]
	A = np.concatenate((np.ones((S.shape[1], 1)), R * T * S.T), axis = 1)
	
	b = -np.array([-R * T * np.log(item[0]) for item in enzymeInfo.loc[:, 'Keq']])
	b = b * Vss[S.columns].values
	
	lb = [-np.inf] + [np.log(concLB)] * S.shape[0]
	ub = [np.inf] + [np.log(concUB)] * S.shape[0]
	
	meth = 'cvxopt_lp'
	
	
	p = LP(f = f, A = A, b = b, lb = lb, ub = ub, iprint = -1, name = 'Maximize minimal driving force')
	r = p.solve(meth, plot = 0)
	
	
	optLogConcs = r.xf[1:]
	optConcs = pd.Series(np.exp(optLogConcs), index = S.index)   
	
	optDeltaGs = pd.Series(-b + R * T * np.dot(S.T, optLogConcs), index = S.columns)
	
	refDeltaGs = pd.Series(-b, index = S.columns)
	
	return optConcs, optDeltaGs, refDeltaGs


def optimize_enzyme_cost(S, Vss, enzymeInfo, concLB, concUB):
	'''
	Parameters
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns. negative for substrates, positive for products	
	Vss: ser, net fluxes in steady state (including in and out fluxes)
	enzymeInfo: df, reaction in rows
	concLB: float, concentration lower bound (mM) for all metabolites
	concUB: float, concentration upper bound (mM) for all metabolites
	
	Returns
	optConcs: ser, optimal log(concentrations)
	optEnzyCosts: ser, optimal enzyme costs
	optEnzyCostTotal: float, optimal total enzyme cost
	'''
	
	from sympy import symbols, Matrix, lambdify
	from openopt import LP, NLP
	from constants import R, T
	from common_rate_laws import E_expression
	
	# get initial concentration guess
	def get_initial_guess(S, enzymeInfo, concLB, concUB):
		
		nMetabs = len(S.index)
		
		ini = []
		for i, metab in enumerate(S.index):
		
			f = np.zeros(nMetabs)
			f[i] = 1
			
			A = R * T * S.T
			b = -np.array([-R * T * np.log(item[0]) for item in enzymeInfo.loc[:, 'Keq']])
			
			lb = [np.log(concLB)] * nMetabs
			ub = [np.log(concUB)] * nMetabs
			
			meth = 'pclp'
		
			p = LP(f = f, A = A, b = b, lb = lb, ub = ub, iprint = -1)
			r = p.solve(meth, plot = 0)
			
			ini.append(r.xf[i])
			
		ini = np.array(ini)
		
		return ini
	
	# get f and dfdx
	def get_enzyme_cost(logConcs, S, Vss, enzymeInfo):
			
		enzyCosts = []
		for enzyme in S.columns:
		
			ifRev = enzymeInfo.loc[enzyme, 'rev']
		
			subMasks = S.loc[:, enzyme] < 0
			subLogConcs = logConcs[subMasks]
			subCoes = S.loc[subMasks, enzyme].abs().values   
			subKms = np.array([item[0] for item in enzymeInfo.loc[enzyme, 'subsKm'].loc[subMasks]])
			
			proMasks = S.loc[:, enzyme] > 0
			proLogConcs = logConcs[proMasks]
			proCoes = S.loc[proMasks, enzyme].abs().values   
			proKms = np.array([item[0] for item in enzymeInfo.loc[enzyme, 'prosKm'].loc[proMasks]])
			
			MW = enzymeInfo.loc[enzyme, 'MW']
			v = Vss[enzyme]
			kcat = enzymeInfo.loc[enzyme, 'kcat'][0]
			deltaGm = -R * T * np.log(enzymeInfo.loc[enzyme, 'Keq'][0])
			
			enzyCost = MW * E_expression(ifRev, subLogConcs, subCoes, subKms, proLogConcs, proCoes, proKms, v, kcat, deltaGm)
			
			enzyCosts.append(enzyCost)
		
		return enzyCosts	
		
	def symbolize_f(logConcs, S, Vss, enzymeInfo):
		
		obj = np.sum(get_enzyme_cost(logConcs, S, Vss, enzymeInfo))
		
		return obj
		
	def symbolize_dfdx(logConcs, S, Vss, enzymeInfo):
		
		jac = Matrix([symbolize_f(logConcs, S, Vss, enzymeInfo)]).jacobian(logConcs)
		
		return jac
	
		
	logConcs = np.array(symbols(' '.join(S.index)))
	
	pref = lambdify(logConcs, symbolize_f(logConcs, S, Vss, enzymeInfo), modules = 'numpy')	
	predfdx = lambdify(logConcs, symbolize_dfdx(logConcs, S, Vss, enzymeInfo), modules = 'numpy')	
	
	f = lambda x, pref, predfdx: pref(*x)
	dfdx = lambda x, pref, predfdx: predfdx(*x)
	
	# optimize enzyme cost
	A = R * T * S.T
	b = -np.array([-R * T * np.log(item[0]) for item in enzymeInfo.loc[:, 'Keq']])
	
	lb = np.full(len(logConcs), np.log(concLB))
	ub = np.full(len(logConcs), np.log(concUB))
		
	meth = 'ralg'	
	
	ini = np.log(concLB) + np.random.rand(len(logConcs)) * (np.log(concUB) - np.log(concLB))
	#ini = get_initial_guess(S, enzymeInfo, concLB, concUB)
	
	
	p = NLP(f = f, x0 = ini, df = dfdx, args = (pref, predfdx), A = A, b = b, lb = lb, ub =ub, iprint = -1, name = 'Minimize enzyme cost')
	r = p.solve(meth, plot = 0)
	
	
	optLogConcs = r.xf
	optEnzyCostTotal = r.ff
		
	optConcs = pd.Series(np.exp(optLogConcs), index = S.index)   

	optEnzyCosts = pd.Series(get_enzyme_cost(optLogConcs, S, Vss, enzymeInfo), index = S.columns)
	
	
	return optConcs, optEnzyCosts, optEnzyCostTotal	
	
	
	
	
	
	




