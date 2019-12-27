# PathParser
PathParser is a python-based computational tool for the thermodynamics and robustness analysis of native and designed metabolic pathways. The following functionalities are provided:   
   
1 Max-min driving force (MDF) optimization. The Gibbs free energy change of the least favorable reaction will be maximized to evaluate thermodynamic feasibility of the entire pathway.  
   
2 Protein cost estimation. The smallest enzyme investment will be assessed by calculating the minimal total enzyme cost supporting a given pathway flux.   
   
3 Robustness analysis. Specifically, an ensemble of models is generated to simulate the system response to enzyme perturbation based on bifurcation theory and a continuation method. Probability of system failure and flux fold change against enzyme perturbation are estimated as well as a flux control index.   
## Dependencies 
