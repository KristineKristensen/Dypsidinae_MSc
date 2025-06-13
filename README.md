# Dypsidinae_MSc
Here scripts used during my master thesis are available.

There was first made a consensus tree of 1,000 ASTRAL-IV trees made on randomly sampled gene trees from the posterior distribution for the data used in Eiserhardt et al. 2022 and 4 trees using wASTRAL(1. Original). Following this the same method was used on a dataset with more genes and with separation of multicopy genes(2. PhyloPalm). Finally the same analysis were made on a revised version of the the dataset used in 2.PhyloPalm (3_new_dataset).

Most of the computing for this project was performed on the GenomeDK cluster. I would like to thank GenomeDK and Aarhus University for providing computational resources and support that contributed to these research results.

# 1_Original
Workflow_with_heyduk.py was used to run the analysis on the GenomeDK cluster. For this i used the following scripts:

1. MrBayes.py \is a script used to make a MrBayes block using the results from ModelFinder 
2. MrBayes_2.py is like MrBayes.py but used when genes needed to run longer e.g. a difference is append=yes, so MrBayes continues from the previous analysis.
3. nexus_to_newick.py is used for converting the trees from NEXUS to Newick format and gather 1000 gene trees from the two MrBayes runs into one file.
4. random_tree_sets.py takes the one random gene tree from the posterior distribution of each gene to make a file with one random tree from each gene that can be used by ASTRAL-IV.
5. astral_2_mrbayes.py is used for removing LPP in the tree, so only the branch lengths is kept.
6. make_nexus.R is used to convert the ASTRAL-IV species trees without LPP into a nexus format that MrBayes can read.
7. format_nexus.py is a script for changing the format of the NEXUS file, to resample the nexus file that MrBayes creates.
8. mb_writing_sumt.py is used for writing a MrBayes block that is used for employing the sumt command on the last 1502 runs of MrBayes to find one consensus gene tree.
9. one_treefile.py takes the consensus gene trees from MrBayes and collects them into one file and convert the file to Newick format so the file can be given to wASTRAL.
10. genes_converged_IQtree.py makes one file with the gene trees from IQ-TREE that corresponds to genes that converged in MrBayes to make a file that can be given to wASTRAL.


# 2_PhyloPalm
Workflow_phylopalm.py was used to run the analysis on the GenomeDK cluster. For this I used the following scripts:

1. IQtree_2_MrBayes.py is like MrBayes.py  a script used to make a MrBayes block using the results from ModelFinder. 
2. reading_MB_logfiles.py reads the log file of MrBayes and makes as CSV file with the information on ASDF, PSRF, ESS and the plot to evaluate convergence. 
3. IQtree_2_MrBayes_2.py is like MrBayes_2.py script for genes that needed to run longer because they had not converged.
4. nexus_to_newick_phylopalm.py is  like nexus_to_newick.py used for converting the trees from NEXUS to Newick format and making one file with the 1,000 last gene trees from the two MrBayes runs.
5. random_treesets_phylopalm.py is a script that like random_tree_sets.py takes the one random gene tree from the posterior distribution of each gene to make a file with one random tree from each gene that can be used by ASTRAL-IV.
6. astral_for_mrbayes_phylopalm.py is like astral_2_mrbayes.py used for removing LPP in the tree, so only the branch lengths is kept.
7. make_nexus.R in 1_Heyduk were used without changes and can be found there.
8. format_nexus.py  in 1_Heyduk were used without changes and can be found there.
9. mb_writing_sumt.py is used for writing a MrBayes block that is used for employing the sumt command on the last 1502 runs of MrBayes to find one consensus gene tree.
10. Make_one_treefile.py takes the consensus gene trees from MrBayes and collects them into one file and convert the file to Newick format so the file can be given to wASTRAL.
11. genes_converged_IQtree.py makes one file with the gene trees from IQ-TREE that corresponds to genes that converged in MrBayes to make a file that can be given to wASTRAL.

# 3_new_dataset
Workflow.py was used to run the analysis on the GenomeDK cluster. For this I used the following scripts:

1. splitting_species_fasta.py a script used to split the species fasta into fasta files for each genomic region
2. rm_species_MSA.py a script that removes species from the MSA
3. amas_raw.sh and amas_gt1.sh scripts for calculatin summary statistic with AMAS (not written by me)
4. optrimal.R an adapted version of the optrimAl.R script (Shee Z.Q., Frodin D.G., Cámara-Leret R. and Pokorny L. 2020. Reconstructing the Complex Evolutionary History of the Papuasian Schefflera Radiation Through Herbariomics. Front. Plant Sci., 11. https://doi.org/10.3389/fpls.2020.00258)
5. IQtree_2_MrBayes.py is like MrBayes.py  a script used to make a MrBayes block using the results from ModelFinder. 
6. reading_MB_logfiles.py reads the log file of MrBayes and makes as CSV file with the information on ASDF, PSRF, ESS and the plot to evaluate convergence. 
7. IQtree_2_MrBayes_2.py is like MrBayes_2.py script for genes that needed to run longer because they had not converged.
8. reading_MB_logfiles_2.py like MrBayes_2.py script for genes that needed to run longer because they had not converged.
9. IQtree_2_MrBayes_n.py and reading_MB_logfiles_n.py was then used to make genes run longer and check convergence for n=3 and n=4 for 5 and 7.5 milion generations.
10. nexus_to_newick.py is used for converting the trees from NEXUS to Newick format and making one file with the 1,000 last gene trees from the two MrBayes runs.
11. random_tree_sets.py is a script that like random_tree_sets.py takes the one random gene tree from the posterior distribution of each gene to make a file with one random tree from each gene that can be used by ASTRAL-IV.
12. mb_writing_sumt.py is used for writing a MrBayes block that is used for employing the sumt command on the last 1502 runs of MrBayes to find one consensus gene tree.
13. Make_one_treefile.py takes the consensus gene trees from MrBayes and collects them into one file and convert the file to Newick format so the file can be given to wASTRAL.
14. genes_converged_IQtree.py makes one file with the gene trees from IQ-TREE that corresponds to genes that converged in MrBayes to make a file that can be given to wASTRAL.
15. model_selected.py Script that writes a csv with the subsitution model selected for each genomic regions


# Make_figures
Cotains the scripts used to make figures.


