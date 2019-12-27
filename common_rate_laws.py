#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '2/15/2019'
__version__ = '1.0'




def E_expression(ifRev, sLogConcs, sCoes, sKms, pLogConcs, pCoes, pKms, v, kcat, deltaGm):
	'''
	Parameters
	ifRev: int, i reversible, 0 irreversible
	sLogConcs: array, log(substrate concentrations)
	sCoes: array, substrate coefficients
	sKms: array, substrate Km values
	pLogConcs: array, log(product concentrations)
	pCoes: array, product coefficients
	pKms: array, product Km values
	v: scalar, flux through current enzyme
	kcat: scalar, enzyme kcat value
	deltaGm: scalar, enzyme ΔG'm
	
	Returns
	E: scalar, enzyme concentration, in exponential form
	'''

	from constants import R, T
	from numpy import sum, log
	from sympy import exp
	
	if ifRev:
		E = v / kcat * (exp(sum(pCoes * (pLogConcs - log(pKms))) - sum(sCoes * (sLogConcs - log(sKms)))) + exp(sum(sCoes * (log(sKms) - sLogConcs))) + 1) / (1 - exp(deltaGm / R / T + sum(pCoes * pLogConcs) - sum(sCoes * sLogConcs)))
		
	else:
		E = v / kcat * (1 + exp(sum(sCoes * (log(sKms) - sLogConcs))))	
			
	return E

	
def v_expression(ifRev, sConcs, sCoes, sKms, pConcs, pCoes, pKms, Keq):
	'''
	Parameters
	ifRev: int, i reversible, 0 irreversible
	sConcs: array, substrate concentrations
	sCoes: array, substrate coefficients
	sKms: array, substrate Km values
	pConcs: array, product concentrations
	pCoes: array, product coefficients
	pKms: array, product Km values
	Keq: scalar, equilibrium constant
	
	Returns
	kin: scalar, kinetic part of generalized rate law
	'''
	
	from numpy import product
	
	if ifRev:
		kin = product((1 / sKms)**sCoes) * (product(sConcs**sCoes) - product(pConcs**pCoes) / Keq) / (product((1 + sConcs / sKms)**sCoes) + product((1 + pConcs / pKms)**pCoes) - 1)
		
	else:
		kin = product((sConcs / sKms)**sCoes) / product((1 + sConcs / sKms)**sCoes)	
			
	return kin
	







