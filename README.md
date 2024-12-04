# Dypsidinae_MSc
Here all scripts used during my master thesis are available

There was carried out two analyses the first making a consensus tree of 1000 Astral trees made on random genetrees from the posterior distribution for the data used in Eiserhardt et al. 2022 and the second dataset with more genes and with exclusion of paralogs.

# First analysis
Workflow_with_heyduk.py was used to run the analysis on Genome.au.dk. For this i used the following scripts:
1. MrBayes_script.py is a script used to make a MrBayes block using the resutls from IQtree ModelFinder. The script was changed accroding to the desired number of generations and when genes needed to run longer append=yes was added, so MrBayes continued from the previous analysis.
2. nexus to newick.py is a script used for converting the trees from nexus to newick format and making one file with the 1000 genetrees from the two MrBayes runs.
3. make_random_treesets.py is a script that takes the one random genetree from the posterior distribution from each gene that can be used by IQtree
4. astral to mrBayes
5. make nexus
6. rename
7. run MrBayes
8. rename with bash

# Second analysis
Workflow_phylopalm.py was used to run the analysis on Genome.au.dk. For this I additionally used the following scripts:
1. MrBayes_script.py is a script used to make a MrBayes block using the resutls from IQtree ModelFinder. The script was changed accroding to the desired number of generations and when genes needed to run longer append=yes was added, so MrBayes continued from the previous analysis.
2. read MrBayes results.py is used for cheking if the MrBayes run have converged. If they had not converged MrBayes was run for more generations
3. nexus_to_newick_phylopalm.py is a script used for converting the trees from nexus to newick format and making one file with the 1000 last genetrees from the two MrBayes runs.
4. make_random_treesets_phylopalm.py is a script that takes the one random genetree from the posterior distribution from each gene that can be used by IQtree
5. astral to mrBayes
6. make nexus
7. rename
8. run MrBayes
9. rename with bash

