# PathParser
PathParser is a python-based computational tool for the thermodynamics and robustness analysis of native and designed metabolic pathways. The following functionalities are provided:   
   
__1 Max-min driving force (MDF) optimization.__ The Gibbs free energy change of the least favorable reaction will be maximized to evaluate thermodynamic feasibility of the entire pathway.  
   
__2 Protein cost estimation.__ The smallest enzyme investment will be assessed by calculating the minimal total enzyme cost supporting a given pathway flux.   
   
__3 Robustness analysis.__ Specifically, an ensemble of models is generated to simulate the system response to enzyme perturbation based on bifurcation theory and a continuation method. Probability of system failure and flux fold change against enzyme perturbation are estimated as well as a flux control index.   
## Dependencies 
PathParser was developed and tested using Python 3.6+ with the following packages:   
   
numpy1.16.1, pandas0.23.4, scipy1.2.1, sympy1.1.1, matplotlib3.0.3, seaborn0.8.1, openopt0.5625
## Usage
main1.py performs MDF optimization and protein cost estimation with the following arguments:
   
-o, --outDir: output directory
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



    
    
    

