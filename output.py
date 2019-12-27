#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '10/18/2018'
__version__ = '1.0'


import numpy as np
import pandas as pd




def print_driving_force_optimization_results(optConcs, optDeltaGs, refDeltaGs):
	'''
	Parameters	
	optConcs: ser, optimal concentrations
	optDeltaGs: ser, optimal minimal driving forces
	refDeltaGs: ser, reference minimal driving forces (all concentrations at 1 mM)
	'''
	
	print("\nReaction\tΔG' at 1mM\toptimized ΔG'")
	
	for rnx in optDeltaGs.index: 
		print('%s\t%.1f\t%.1f' % (rnx, refDeltaGs.loc[rnx], optDeltaGs.loc[rnx]))
	
	print('\nMetabolite\toptimized conc. (mM)')	
	
	for metab in optConcs.index:
		print('%s\t%6.3f' % (metab, optConcs.loc[metab]))

			
def save_driving_force_optimization_results(optConcs, optDeltaGs, refDeltaGs, outDir):
	'''
	Parameters	
	optConcs: ser, optimal concentrations
	optDeltaGs: ser, optimal minimal driving forces
	refDeltaGs: ser, reference minimal driving forces (all concentrations at 1 mM)
	outDir: str, output directory
	'''
	
	allDeltaGs = pd.concat([optDeltaGs, refDeltaGs], axis = 1)
	allDeltaGs.to_csv('%s/minimal_driving_forces.tsv' % outDir, sep = '\t', header = ["ΔG' at 1mM", "optimized ΔG'"], index_label = '#Reaction')
	
	optConcs.to_csv('%s/metabConc_MDF.tsv' % outDir, sep = '\t', header = ['Optimized Concentration (mM)'], index_label = '#Metabolite')
	

def plot_cumulative_deltaGs(optDeltaGs, refDeltaGs, outDir):
	'''
	Parameters	
	optDeltaGs: ser, optimal minimal driving forces
	refDeltaGs: float, reference minimal driving forces (all concentrations at 1 mM)
	outDir: str, output directory
	'''
	
	import re
	import platform
	system = platform.system()
	if re.search(r'linux', system, flags = re.I):
		import matplotlib
		matplotlib.use('agg')    
	import matplotlib.pyplot as plt
	
	enzymes = optDeltaGs.index
	optDeltaGs = optDeltaGs.values
	refDeltaGs = refDeltaGs.values

	x = np.arange(1, len(enzymes) + 2)
	
	xticks = list(map(np.mean, zip(x[1:], x[:-1])))   

	cumOptDeltaG = np.cumsum(optDeltaGs)
	cumRefDeltaG = np.cumsum(refDeltaGs)
	
	cumOptDeltaG = np.insert(cumOptDeltaG, 0, 0)
	cumRefDeltaG = np.insert(cumRefDeltaG, 0, 0)
	
	
	#plt.style.use('ggplot')
	
	plt.figure(figsize = (enzymes.size * 1, 6))
	
	plt.plot(x, cumOptDeltaG, marker = 'o', label = 'optimal conc.', color = '#E24A33', linewidth = 2)
	plt.plot(x, cumRefDeltaG, marker = '^', label = 'conc. at 1 mM', color = '#1f77b4', linewidth = 2)
	
	plt.ylabel("Cumulative $\Delta$G' (kJ mol$^{-1}$)", fontsize = 20)
	
	plt.xticks(xticks, enzymes, fontsize = 15)
		
	plt.legend(fontsize = 15)

	plt.savefig('%s/minimal_driving_forces.jpg' % outDir, dpi = 300, bbox_inches = 'tight')
	
	
def print_enzyme_cost_optimization_results(optConcs, optEnzyCosts, optEnzyCostTotal):
	'''
	Parameters
	optConcs: ser, optimal concentrations
	optEnzyCosts: ser, optimal enzyme costs
	optEnzyCostTotal: float, optimal total enzyme cost
	'''
	
	print('\nOptimized total enzyme cost: %.1e g/mol/s' % (optEnzyCostTotal * 1000))   # MW in kDa
	
	print('\nReaction\toptimized enzyme cost (g/mol/s)')
	for rnx in optEnzyCosts.index: 
		print('%s\t%.1e' % (rnx, optEnzyCosts.loc[rnx] * 1000))   # MW in kDa
	
	print('\nMetabolite\toptimized conc. (mM)')
	for metab in optConcs.index:
		print('%s\t%6.3f' % (metab, optConcs.loc[metab]))

			
def save_enzyme_cost_optimization_results(optConcs, optEnzyCosts, optEnzyCostTotal, outDir):
	'''
	Parameters	
	optConcs: ser, optimal concentrations
	optEnzyCosts: ser, optimal enzyme costs
	optEnzyCostTotal: float, optimal total enzyme cost
	outDir: str, output directory
	'''
	
	optEnzyCosts.loc['Total'] = optEnzyCostTotal
	optEnzyCosts = optEnzyCosts * 1000   # MW in kDa
	optEnzyCosts.to_csv('%s/enzyme_protein_costs.tsv' % outDir, sep = '\t', header = ['optimized enzyme cost (g/mol/s)'], index_label = '#Reaction')
	
	optConcs = optConcs
	optConcs.to_csv('%s/metabConc_EPC.tsv' % outDir, sep = '\t', header = ['Optimized Concentration (mM)'], index_label = '#Metabolite')
	
			
def plot_enzyme_costs(optEnzyCosts, outDir):
	'''
	Parameters	
	optEnzyCosts: ser, optimal enzyme costs
	outDir: str, output directory
	'''
	
	import re
	import platform
	system = platform.system()
	if re.search(r'linux', system, flags = re.I):
		import matplotlib
		matplotlib.use('agg')
	import matplotlib.pyplot as plt
	
	enzymes = optEnzyCosts.index
	costs = optEnzyCosts.values
	
	#plt.style.use('ggplot')
	
	plt.figure(figsize = (enzymes.size * 1, 6))
	
	plt.bar(np.arange(len(enzymes)), costs * 1000, color = '#1f77b4', tick_label = enzymes)   # MW in kDa
	
	plt.xticks(fontsize = 15)
	plt.ylabel("Enzyme protein cost (g/(mol s$^{-1}$)", fontsize = 20)
	
	plt.savefig('%s/enzyme_protein_costs.jpg' % outDir, dpi = 300, bbox_inches = 'tight')	
	
	
def dump_ensemble_models(pertResults, outDir):
	'''
	Parameters
	pertResults: dict, simulation results from ensemble models
	outDir: str, output directory
	'''
	
	import pickle
	
	with open('%s/simulation_results.bin' % outDir, 'wb') as f:
			
		pickle.dump(pertResults, f)
	
	
def save_robustness_index(robustIdx, outDir):
	'''
	Parameters
	robustIdx: ser, median of robustness index Si for each enzyme
	outDir: str, output directory
	'''
	
	robustIdx.to_csv('%s/robustness_index.tsv' % outDir, sep = '\t', header = ['S index'], index_label = '#Reaction')
	
			
def plot_robustness_index(robustIdx, outDir):
	'''
	Parameters
	robustIdx: ser, median of robustness index Si for each enzyme
	outDir: str, output directory
	'''
	
	import re
	import platform
	system = platform.system()
	if re.search(r'linux', system, flags = re.I):
		import matplotlib
		matplotlib.use('agg')
	import matplotlib.pyplot as plt	
		
	enzymes = robustIdx.index
	Sindex = robustIdx.values
	
	#plt.style.use('ggplot')
	
	plt.figure(figsize = (enzymes.size * 1, 6))
	
	plt.bar(range(enzymes.size), Sindex, tick_label = enzymes, color = '#E24A33')
	
	plt.xticks(fontsize = 15)
	plt.ylabel('Robustness index (totally %.3f)' % np.sum(Sindex), fontsize = 20)
	
	plt.savefig('%s/robustness_index.jpg' % outDir, dpi = 300, bbox_inches = 'tight')
	

def save_system_failure_probability(failurePro, outDir):
	'''
	Parameters
	failurePro: df, probability of system failure, enzyme in rows, enzyme level in columns
	outDir: str, output directory	
	'''
	
	failurePro.to_csv('%s/system_failure.tsv' % outDir, sep = '\t', index_label = '#Reaction')
	
	
def plot_system_failure_probability(failurePro, outDir):
	'''
	Parameters
	failurePro: df, probability of system failure, enzyme in rows, enzyme level in columns
	outDir: str, output directory
	'''
	
	import re
	import platform
	system = platform.system()
	if re.search(r'linux', system, flags = re.I):
		import matplotlib
		matplotlib.use('agg')
	import matplotlib.pyplot as plt
	
	nEnzyme = failurePro.shape[0]
	nCol = 3
	nRow = np.ceil(nEnzyme / nCol)
	
	#plt.style.use('ggplot')
	
	plt.figure(figsize = (8, nRow * 7/3))

	for i, enzyme in enumerate(failurePro.index):
		
		ind = i + 1
		ax = plt.subplot(nRow, nCol, ind)
		
		x = failurePro.columns
		y = failurePro.loc[enzyme, :]
		
		ax.semilogx(x, y, color = '#E24A33')
		
		ax.set_xlabel('%s fold change' % enzyme, fontsize = 12)
		
		ax.set_xticks([x.min(), 0.5, 1, 2, x.max()])
		ax.set_xticklabels(['%.1f' % x.min(), '0.5', '1', '2', '%.0f' % x.max()])
		
		ax.set_ylim([-0.05, 1.05])
		ax.set_yticks([0, 0.5, 1])
		ax.set_yticklabels(['0', '50%', '100%'])
		
	plt.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.95, top = 0.9, wspace = 0.4, hspace = 0.5)   
	
	plt.suptitle('Probability of system failure', fontsize = 20) 
	
	plt.savefig('%s/system_failure.jpg' % outDir, dpi = 300)
	
	
def save_flux_fold_change(enzymes, fluxChange, outDir):
	'''
	Parameters	
	enzymes: lst, enzyme IDs
	fluxChange: dict
	outDir: str, output directory	
	'''
	
	for enzyme in enzymes:
		fluxChange[enzyme].to_csv('%s/flux_change_%s.tsv' % (outDir, enzyme), sep = '\t', index_label = '#Flux change no.')
	
	
def plot_flux_fold_change(enzymes, fluxChange, outDir, fluxBndsShow = (0.1, 10)):
	'''
	Parameters	
	enzymes: lst, enzyme IDs
	fluxChange: dict
	outDir: str, output directory
	fluxBndsShow: 2-tuple, flux change for show
	'''
	
	import re
	import platform
	system = platform.system()
	if re.search(r'linux', system, flags = re.I):
		import matplotlib
		matplotlib.use('agg')
	import matplotlib.pyplot as plt
	from matplotlib.colors import LinearSegmentedColormap
	import seaborn as sns
	
	nEnzyme = len(enzymes)
	nCol = 3
	nRow = np.ceil(nEnzyme / nCol)
	
	cmap = LinearSegmentedColormap.from_list(name = 'mycolor', colors = [(1,1,1), (31/256,119/256,180/256)], N=10)
	
	plt.figure(figsize = (8, nRow * 7/3))
	
	for i, enzyme in enumerate(enzymes):
	
		data = fluxChange[enzyme]
	
		ind = i + 1
		ax = plt.subplot(nRow, nCol, ind)
		
		
		sns.heatmap(data, cmap = cmap, ax = ax, cbar = True)
		
		ax.set_xlabel('%s fold change' % enzyme, fontsize = 12)
		
		xticks = np.arange(data.shape[1]) + 0.5   
		ax.set_xticks(np.percentile(xticks, [0, (np.log10(0.5)+1)/2 * 100, 50, (np.log10(2)+1)/2 * 100, 100]))
		ax.set_xticklabels(['%.1f' % data.columns.min(), '0.5', '1', '2', '%.0f' % data.columns.max()], rotation = 0)
		
		yticks = np.arange(data.shape[0]) + 0.5   
		ax.set_yticks(np.percentile(yticks, [0, 50, 100]))
		ax.set_yticklabels(['%s' % fluxBndsShow[1], '1', '%s' % fluxBndsShow[0]], rotation = 0)
		
		for side in ['bottom', 'left',  'top', 'right']:
			ax.spines[side].set_visible(True)

	plt.subplots_adjust(left = 0.1, bottom = 0.1, right = 0.95, top = 0.9, wspace = 0.4, hspace = 0.5)   
	
	plt.suptitle('Flux fold change', fontsize = 20) 
	
	plt.savefig('%s/flux_change.jpg' % outDir, dpi = 300)	
	
	
def plot_flux_control_index(ConIdx, outDir):
	'''
	Parameters
	ConIdx: df, enyzmes in rows, each cell is the list of flux_control_index
	outDir: str, output directory
	'''
	
	import re
	import platform
	system = platform.system()
	if re.search(r'linux', system, flags = re.I):
		import matplotlib
		matplotlib.use('agg')
	import matplotlib.pyplot as plt
	
	totalWidth = 0.8
	singleWidth = totalWidth / 2

	nenzymes = ConIdx.shape[0]
	x = np.arange(nenzymes)
	x0 = x - 0.5 * singleWidth

	ConIdxMed = ConIdx.applymap(np.mean)
	xticks = ConIdx.index + '\n\n' + ConIdxMed['Up regulation'].round(2).apply(str) + '\n' + ConIdxMed['Down regulation'].round(2).apply(str)


	#plt.style.use('ggplot')

	plt.figure(figsize = (nenzymes * 1, 6))

	colors = ['lightblue', 'pink']

	for i in range(ConIdx.shape[1]):
		
		plt.boxplot(ConIdx.iloc[:, i], positions = x0 + i * singleWidth, widths = 0.7*singleWidth, notch = True, patch_artist = True, showmeans = True, meanline = True, meanprops = {'color':'k'}, boxprops = {'facecolor':colors[i]}, showfliers = False)

	plt.xlim((x[0]-0.6, x[-1]+0.6))
	plt.xticks(x, xticks, fontsize = 15)
	
	plt.ylabel('Flux control index', fontsize = 20)

	for i in range(ConIdx.shape[1]): 
		plt.scatter([], [], marker = 's', color = colors[::-1][i], label = ConIdx.columns[::-1][i])   
	plt.legend(loc = 'center', bbox_to_anchor = (1.2, 0.5), fontsize = 15)
	
	plt.savefig('%s/flux_control_index.jpg' % outDir, dpi = 300, bbox_inches = 'tight')
	
	
	
	
	
	
