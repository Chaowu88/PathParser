# PathParser
PathParser is a python-based computational tool for the thermodynamics and robustness analysis of native and designed metabolic pathways. The following functionalities are provided:   
   
__1 Max-min driving force (MDF) optimization.__ The Gibbs free energy change of the least favorable reaction will be maximized to evaluate thermodynamic feasibility of the entire pathway.  
   
__2 Protein cost estimation.__ The smallest enzyme investment will be assessed by calculating the minimal total enzyme cost supporting a given pathway flux.   
   
__3 Robustness analysis.__ Specifically, an ensemble of models is generated to simulate the system response to enzyme perturbation based on bifurcation theory and a continuation method. Probability of system failure and flux fold change against enzyme perturbation are estimated as well as a flux control index.   
## Dependencies 
PathParser was developed and tested using Python 3.6+ with the following packages:   
   
numpy1.16.1, pandas0.23.4, scipy1.2.1, sympy1.1.1, matplotlib3.0.3, seaborn0.8.1, openopt0.5625
## Usage
