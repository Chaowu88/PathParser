# PathParser
PathParser is a python-based computational tool for the thermodynamics and robustness analysis of native and designed metabolic pathways. The following functionalities are provided:   
   
__1 Max-min driving force (MDF) optimization.__ The Gibbs free energy change of the least favorable reaction will be maximized to evaluate thermodynamic feasibility of the entire pathway.  
   
__2 Protein cost estimation.__ The smallest enzyme investment will be assessed by calculating the minimal total enzyme cost supporting a given pathway flux.   
   
__3 Robustness analysis.__ Specifically, an ensemble of models is generated to simulate the system response to enzyme perturbation based on bifurcation theory and a continuation method. Probability of system failure and flux fold change against enzyme perturbation are estimated as well as a flux control index.   
## Dependencies 
PathParser was developed and tested using Python 3.6+ with the following packages:   
   
numpy1.16.1, pandas0.23.4, scipy1.2.1, sympy1.1.1, matplotlib3.0.3, seaborn0.8.1, openopt0.5625
## Usage
__main1.py__ performs MDF optimization and protein cost estimation with the following arguments:
   
>-o, --outDir: output directory   
-r, --reactionFile: reaction file, required fields: Enzyme ID, Substrates, Products, Reversibility, ΔrGm and Enzyme MW. See below as an example   
-b, --concBnds: concentration lower and upper bound (mM) for all metabolites, sep by ","   
-w, --runWhich: which analysis to run, '1' for maximizing the minimal driving force, '2' for minimizing the totol enzyme protein cost, '12' for both   
-i, --iniMetabs: optional, metabolites as initial substrates, sep by ",". By default, they will be detected automatically, sometimes they should be set explicitly, e.g. for cylic pathways   
-f, --finMetabs: optional, metabolites as end products, sep by ",". By default, they will be detected automatically, sometimes they should be set explicitly, e.g. for cylic pathways  
-eb, --exBalMetabs: optional, metabolites excluded from mass balance, sep by ","  
-eo, --exOptMetabs: optional, metabolites excluded from optimization, sep by ","  
-a, --assignFlux: optional, assign flux to some enzyme in the format "enzyme ID:value", then flux distribution will be calculated. By default, influx to pathway will be set to 1  
-h, --help: show help message and exit  
   
example:   
```
python way\to\PathParser\main1.py -o way\to\PathParser\example\CBB -r way\to\PathParser\example\CBB.tsv -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP -b 0.001,10 -w 12
```
__main2.py__ performs robustness analysis with the following arguments:
    
>-o, --outDir: see above  
-r, --reactionFile: reaction file, required fields: Enzyme ID, Substrates, Products, Reversibility, ΔrGm and Enzyme MW. See above as an example  
-n, --nmodels: number of models in an ensemble  
-b, --enzymeBnds: lower and upper bound of relative enzyme level, sep by ","  
-d, --ifDump: whether to dump generated models, "yes" or "no"  
-p, --nprocess: number of processes to run simultaneously  
-w, --runWhich: which analysis to run, '1' for robustmess index, '2' for probability of system failure, '3' for flux fold change, and any other combination of the numbers     
-t, --ifReal: whether to use the real value of concentrations, Kms and Keqs, "yes" or "no"  
-a, --assignFlux: assign flux (mmol/gCDW/h) to some enzyme in the format "enzyme ID:value", then flux distribution of reference state will be calculated, required if --ifReal is "yes"  
-mc, --metabConcFile: file of metabolite concentrations (mM) in reference state, required if --ifReal is "yes"  
-ec, --enzConcFile: file of enzyme concentrations (mmol/gCDW) in reference state, required if --ifReal is "yes"  
-i, --iniMetabs: optional, see above  
-f, --finMetabs: optional, see above  
-eb, --exBalMetabs: optional, see above  
-eo, --exOptMetabs: optional, see above  
 
__NOTE.__   
  
1 It is highly recommended to run this script in a high-performance computer.  
2 Robustness against enzyme perturbation can be evaluated in both relative and absolute manners, real flux values as well as metabolite concentrations and enzyme concentrations should be provided if the latter.  
    
example:
```
python way\to \PathParser\main2.py -o way\to\example\CBB -r example\example\CBB.tsv -f GAP -eb ATP,ADP,Pi,NADH,NAD,NADPH,NADP -eo ATP,ADP,Pi,NADH,NAD,NADPH,NADP -n 1000 -b 0.1,10 -d no -w 123 -p 30 -t no
```
## License
PathParser is released under a GNU General Public [License](https://github.com/Chaowu88/PathParser/blob/master/LICENSE).
## Citation
Chao Wu, Huaiguang Jiang, Isha Kalra, Xin Wang, Melissa Cano, PinChing Maness, Jianping Yu, Wei Xiong. A generalized computational framework to streamline thermodynamics and kinetics analysis of metabolic pathways. Metabolic Engineering. 2019 Aug 8. pii: S1096-7176(19)30258-7.   
   
https://www.sciencedirect.com/science/article/pii/S1096717619302587

    
    
    

