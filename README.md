
##----- Maia and Guimaraes ------

###--- The hierarchical coevolutionary units of ecological networks ---

Repository of data and scripts of manuscript entitled "The hierarchical coevolutionary units of ecological networks" of Kate Pereira Maia and Paulo Roberto Guimaraes Jr. 

REPO DOI: [![DOI](https://zenodo.org/badge/738971666.svg)](https://zenodo.org/doi/10.5281/zenodo.12628129)

The repository is divided in the following four folders:

#### Data:
The dataset file has information on each ecological networks used: network ID (ID), interaction sign (IntSign, either A for antagonism or M for mutualism), interaction type (IntType) and additional natural history information	(IntNH).
The Matrices folder has the 376 ecological networks used in the study as incidence matrices.

#### Scripts:
This folder contains the 14 scripts and additional functions used to manipulate data, produce files in outputs and plot results. Description of code and outputs can be found in each script, but in summary:
Scripts 1-5: compute groups across network scales, and tests for a hierarchical organisation. 
Script 6: plots groups and hierarchical structure results. 
Script 7: produces results on groups and hierarchical structure found in the manuscript. 
Script 8-10: compute matrices of coevolutionary effects and tests if these effects respect groups and hierarchical organisation. 
Script 11: plots coevolutionary effects results. 
Script 12: produces results on coevolutionary effects found in the manuscript and tables S2-S3.
Scripts 13-4: performs sensitivity analyses.

#### Output:
To facilitate (some scripts are time consuming), all outputs produced by scripts 1-5 and 8-10 were already generated. Those include matrices of coevolutionary effects (TE-matrices, Q-matrices), profile of subgraphs found in each network (m6-Profiles), dataframes on node and net-level structure used for plots and analyses, and dataframes on coevolutionary effects used for plots and analyses. Detailed information on each output can be found in their corresponding scripts.

#### Modularity:
The Matrices folder has the matrices (_GC1) for the largest component of each network (generated in Script 01) and used to run the modularity analyses with the MODULAR programme (similarly for the Null_Matrices folder). The folder also contains the modularity results (produced by MODULAR) for two algorithms: Simulated Annealing and Fast Grid. Because MODULAR does not accept row/col names, Script 01 also saves species names in the Sp_Names folder.   
