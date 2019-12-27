#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '10/20/2018'
__version__ = '1.0'




def assign_metabolites_to_reactions(Smetab2rnx, X):
	'''
	Parameters
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	X: array, metabolite concentrations
	
	Returns
	X4rnxs: array, each entry is tuple of substrate and product concatenations for some reaction (in order of enzymes) 
	'''
	
	import numpy as np
	
	def transformer(row, X):
		
		if X[row == -1].size > 0:
			return np.array(X[row == -1]), np.array(X[row == 1])
		else:
			return np.ones(1), np.array(X[row == 1])   

	X4rnxs = Smetab2rnx.apply(transformer, axis = 1, raw = True, args = (X,))
    
	return np.array(X4rnxs)			
	

def get_V(Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs):
	'''
	Parameters
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	model: func, rate law model
	E: array, enzyme concentrations, in order of enzymes
	X: array, metabolites concentrations, in order of metabs
	reverses: array, whether reversible, in order of enzymes
	kcats: array, kcat, in order of enzymes
	subCoess: array of array, substrate coefficients in order of enzymes
	subKmss: array of array, substrate Kms, in order of enzymes
	proCoess: array of array, product coefficients, in order of enzymes
	proKmss: array of array, product Kms, in order of enzymes
	Keqs: array, Keq, in order of enzymes
	
	Returns
	V: sym mat, fluxes
	'''
	
	from sympy import Matrix
	
	X4rnxs = assign_metabolites_to_reactions(Smetab2rnx, X)

	V = [kcats[i] * E[i] * model(reverses[i], X4rnxs[i][0], subCoess[i], subKmss[i], X4rnxs[i][1], proCoess[i], proKmss[i], Keqs[i]) for i in range(Smetab2rnx.shape[0])]

	return Matrix(V)     
	
	
def get_dVdX(Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs):
	'''
	Parameters
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	model: func, rate law model
	E: array, enzyme concentrations, in order of enzymes
	X: array, metabolites concentrations, in order of metabs
	reverses: array, whether reversible, in order of enzymes
	kcats: array, kcat, in order of enzymes
	subCoess: array of array, substrate coefficients in order of enzymes
	subKmss: array of array, substrate Kms, in order of enzymes
	proCoess: array of array, product coefficients, in order of enzymes
	proKmss: array of array, product Kms, in order of enzymes
	Keqs: array, Keq, in order of enzymes
	
	Returns
	J: sym mat, Jacobian matrix 
	'''

	V = get_V(Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
	
	dVdX = V.jacobian(X)
	
	return dVdX
	
	
def get_dVdE(Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs):
	'''
	Parameters
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	model: func, rate law model
	E: array, enzyme concentrations, in order of enzymes
	X: array, metabolites concentrations, in order of metabs
	reverses: array, whether reversible, in order of enzymes
	kcats: array, kcat, in order of enzymes
	subCoess: array of array, substrate coefficients in order of enzymes
	subKmss: array of array, substrate Kms, in order of enzymes
	proCoess: array of array, product coefficients, in order of enzymes
	proKmss: array of array, product Kms, in order of enzymes
	Keqs: array, Keq, in order of enzymes
	
	Returns
	dVdE: sym mat, dVdE
	'''
	
	V = get_V(Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
		
	dVdE = V.jacobian(E)
	
	return dVdE	
	
	
def get_Jacobian(S, Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs):
	'''
	Parameters
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns
	Smetab2rnx: df, transforme X to metabolites needed in each reaction
	model: func, rate law model
	E: array, enzyme concentrations, in order of enzymes
	X: array, metabolites concentrations, in order of metabs
	reverses: array, whether reversible, in order of enzymes
	kcats: array, kcat, in order of enzymes
	subCoess: array of array, substrate coefficients in order of enzymes
	subKmss: array of array, substrate Kms, in order of enzymes
	proCoess: array of array, product coefficients, in order of enzymes
	proKmss: array of array, product Kms, in order of enzymes
	Keqs: array, Keq, in order of enzymes
	
	Returns
	J: sym mat, Jacobian matrix 	
	'''
	
	from sympy import Matrix
	
	dVdX = get_dVdX(Smetab2rnx, model, E, X, reverses, kcats, subCoess, subKmss, proCoess, proKmss, Keqs)
	
	J = Matrix(S) * dVdX   
	
	return J
	
	
def get_lambdify_function(args, func):
	'''
	Parameters	
	args: list of sym variable
	func: sym function
	
	Returns
	funLam: lambdified function
	'''
	
	from sympy import lambdify
	
	funLam = lambdify(args, func, modules = 'numpy')	
	
	return funLam
	
	
def solve_dXdE(Espan, nsteps, Xini, Jlam, dVdElam, S):
	'''
	Parameters
	Espan: df, 1st and 2nd columns are integration interval, enzyme in rows
	nsteps: int, # of integration steps
	Xini: array, ini values of X
	Jlam: lambdified function, Jacobian matrix
	dVdElam: lambdified function, dVdE
	S: df, stoichiometric matrix, metabolite in rows, reaction in columns
		
	Returns
	Eout: df, enzyme expression range, enzyme in rows, columns are the same with Xout 
	Xout: df, metabolite concentration range, metabolite in rows, columns are the same with Eout (initial input metabolite not included)
	'''

	import numpy as np
	import pandas as pd
	from scipy.linalg import eigvals, pinv2
	from constants import eigThreshold
	
	# prepare initial X, E
	Espan = np.matrix(Espan)
	
	dE = (Espan[:, 1] - Espan[:, 0]) / nsteps
	
	X = np.matrix(Xini[:, np.newaxis])
	E = Espan[:, 0]
	
	# prepare initial Xout, Eout
	Xout = pd.DataFrame(index = S.index, columns = range(nsteps + 1))
	Eout = pd.DataFrame(index = S.columns, columns = range(nsteps + 1))
	
	Xout.iloc[:, 0] = X
	Eout.iloc[:, 0] = E
	
	for i in range(1, nsteps + 1):
	
		XE = np.array(np.concatenate((X, E)))
		
		# update Jacobian matrix and screen
		J = np.matrix(Jlam(*XE)).astype(np.float)

		if np.any(eigvals(J).real >= eigThreshold): break   

		# update X, E and screen
		dVdE = np.matrix(dVdElam(*XE)).astype(np.float)

		dX = -pinv2(J) * np.matrix(S) * dVdE * np.matrix(dE)
		
		X = X + dX
		E = E + dE
		
		if X.min() <= 0: break   
		
		# update Xout, Eout
		Xout.iloc[:, i] = X
		Eout.iloc[:, i] = E
		
	return Eout, Xout	
	
