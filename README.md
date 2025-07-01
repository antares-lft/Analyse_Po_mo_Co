# Analyse_Po_mo_Co

1. Tree modelisation tree :
   - Visualisation of the tree. The colors represent the logarithme of the population size.
2. Simulation script:
   - Modify the necessary scripts to run the scripts based on the parameters (recombination, population size).
   - Contains the Julia scripts for the different population sizes to obtain the equilibrium fitness landscape as well as the bash files to run on the cluster.
4. Data Extraction : the goal is to extract the dN/dS, the pN/pS, the branch length and the population size for each simulation
   - script : the code that allows the extraction
   - table : all the csv files
5. Script R analyses :
   - Data Analysis : script to analyse the different table (from the previous extraction) to make the graph
6. Graphiques and figures : graphs of pN/pS and dN/dS link to the population size
   - dN_dS_constant_popsize_vs._analytical_dN_dS : Boxplot to select the best branches for the next graphs and compare to the analytical value obtained by the Julia script
   - Recombination : Each folder contains a population estimate (using pS). dN/dS vs pop size (implemented or experimental) and pN/pS vs pop size (implemented or experimental)
