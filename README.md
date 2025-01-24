# Dypsidinae_MSc
Here scripts used during my master thesis are available.

There was first made a consensus tree of 1,000 ASTRAL -IV trees made on random gene trees from the posterior distribution for the data used in Eiserhardt et al. 2022 (1. Heyduk). Following this the same method was used on a dataset with more genes and with separation of multicopy genes and on this dataset species trees using gene trees from IQ-TREE and MrBayes were also inferred using wASTRAL (2. PhyloPalm).
Some/all of the computing for this project was performed on the GenomeDK cluster. We would like to thank GenomeDK and Aarhus University for providing computational resources and support that contributed to these research results.

# 1_Heyduk
Workflow_with_heyduk.py was used to run the analysis on the GenomeDK cluster. For this i used the following scripts:

1. MrBayes.py is a script used to make a MrBayes block using the results from ModelFinder 
2. MrBayes_2.py is like MrBayes.py but used when genes needed to run longer e.g. a difference is append=yes, so MrBayes continues from the previous analysis.
3. nexus_to_newick.py is used for converting the trees from NEXUS to Newick format and gather 1000 gene trees from the two MrBayes runs into one file.
4. random_tree_sets.py takes the one random gene tree from the posterior distribution of each gene to make a file with one random tree from each gene that can be used by ASTRAL-IV.
5. astral_2_mrbayes.py is used for removing LPP in the tree, so only the branch lengths is kept.
6. make_nexus.R is used to convert the ASTRAL-IV species trees without LPP into a nexus format that MrBayes can read.
7. format_nexus.py is a script for changing the format of the NEXUS file, to resample the nexus file that MrBayes creates.

# 2_PhyloPalm
Workflow_phylopalm.py was used to run the analysis on the GenomeDK cluster. For this I used the following scripts:

1. IQtree_2_MrBayes.py is like MrBayes.py  a script used to make a MrBayes block using the results from ModelFinder. 
2. reading_MB_logfiles.py reads the log file of MrBayes and makes as CSV file with the information on ASDF, PSRF, ESS and the plot to evaluate convergence. 
3.IQtree_2_MrBayes_2.py is like MrBayes_2.py script for genes that needed to run longer because they had not converged.
4. nexus_to_newick_phylopalm.py is  like nexus_to_newick.py used for converting the trees from NEXUS to Newick format and making one file with the 1,000 last gene trees from the two MrBayes runs.
5. random_treesets_phylopalm.py is a script that like random_tree_sets.py takes the one random gene tree from the posterior distribution of each gene to make a file with one random tree from each gene that can be used by ASTRAL-IV.
6. astral_for_mrbayes_phylopalm.py is like astral_2_mrbayes.py used for removing LPP in the tree, so only the branch lengths is kept.
7. make_nexus.R in 1_Heyduk were used without changes and can be found there.
8. format_nexus.py  in 1_Heyduk were used without changes and can be found there.
9. mb_writing_sumt.py is used for writing a MrBayes block that is used for employing the sumt command on the last 1502 runs of MrBayes to find one consensus gene tree.
10. Make_one_treefile.py takes the consensus gene trees from MrBayes and collects them into one file and convert the file to Newick format so the file can be given to wASTRAL.
11.genes_converged_IQtree.py makes one file with the gene trees from IQ-TREE that corresponds to genes that converged in MrBayes to make a file that can be given to wASTRAL.

Programs used:
* Borowiec, M.L.,2016.AMAS: a fast tool for alignment manipulation and computing of summary statistics. PeerJ. 4:e1660. https://doi.org/10.7717/peerj.1660
* Cock, P.A., Antao, T., Chang, J.T., Chapman, B.A., Cox, C.J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B. and de Hoon, M.J.L., 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 25, 1422-1423. https://doi.org/10.1093/bioinformatics/btp163
* Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F, von Haeseler, A., Jermiin, L.S., 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. Nat Methods. 14, 587–589. https://doi.org/10.1038/nmeth.4285 
* Minh, B.Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M.D., von Haeseler, A., Lanfear, R., 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. biol. evol. 37, 5, 1530–1534. https://doi.org/10.1093/molbev/msaa015
* Ronquist, F., Teslenko M., van der Mark, P., Ayres, D.L., Darling, A., Höhna, S., Larget, B., Liu, L., Suchard, M.A., Huelsenbeck, J.P., 2012. MrBayes 3.2: Efficient Bayesian Phylogenetic Inference and Model Choice Across a Large Model Space. Syst. Biol. 61, 3, 539–542. https://doi.org/10.1093/sysbio/sys029
* Tabatabaee, Y., Zhang, C., Warnow, T., Mirarab, S., 2023. Phylogenomic branch length estimation using quartets. Bioinformatics. 39, Supplement_1, i185–i193, https://doi.org/10.1093/bioinformatics/btad221
* Zhang, C., Mirarab, S., 2022. Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees. Mol. biol. evol. 39, 12. https://doi.org/10.1093/molbev/msac215



