
from os import O_SYNC, name
from gwf import Workflow
import os.path
import csv

gwf = Workflow()

   
# ########################################################################################################################
# ##############################################---- IQTREE ----##########################################################
# ########################################################################################################################

def iqtree(path_in, genes):
    """Using IQTREE to find the best model with ModelFinder"""
    inputs = [path_in+genes+".fasta", path_in+genes+".part"]
    outputs = ["/home/kris/dypsidinae/B_phylopalm/1.IQtree/"+genes+".part.treefile"]
    options = {'cores': 2, 'memory': "40g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """
     
    cd /home/kris/dypsidinae/data/1.Dypsidiinae/
        
    #Activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate iqtree
        
    iqtree2 -s {genes}.fasta -T AUTO -m TESTONLY -p {genes}.part -mset mrbayes  
   
    mv *.part.ckp.gz /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *aligned_noempty.fasta-out_clean.part.best_model.nex /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *.part.iqtree /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *.part.model.gz /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *.part.treefile /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *.fasta-out_clean.part.best_scheme /home/kris/dypsidinae/B_phylopalm/1.IQtree/ 
    mv *.fasta-out_clean.part.best_scheme.nex /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *.fasta-out_clean.part.log /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    mv *.fasta-out_clean.part.model.gz /home/kris/dypsidinae/B_phylopalm/1.IQtree/
    """.format(path_in = path_in, genes = genes)

    return (inputs, outputs, options, spec) 
 
# ########################################################################################################################
# ##############################################---- AMAS_convert ----##########################################################
# ########################################################################################################################

def amas_c(path_in, genes):
    """Using AMAS to convert files"""
    inputs = [path_in+genes+".fasta"]
    outputs = ["/home/kris/dypsidinae/B_phylopalm/2.Nex_files/"+genes+".fasta-out.nex"]
    options = {'cores': 2, 'memory': "40g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """
     
    cd /home/kris/dypsidinae/B_phylopalm/2.Nex_files/
        
    #Activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate amas
    
    AMAS.py convert -d dna -f fasta -i {genes}".fasta" -u nexus
       
    """.format(path_in = path_in, genes = genes)

    return (inputs, outputs, options, spec) 


# ########################################################################################################################
# ##############################################---- MrBayes ----#########################################################
# ########################################################################################################################

def MrBayes(path_in, genes):
    """Using MrBayes to find gene trees"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/B_phylopalm/3.MrBayes/"+genes+".log"]
    options = {'cores':"2", 'memory': "120g", 'walltime': "36:00:00", 'account':"dypsidinae"}

    spec = """
    
    cd /home/kris/dypsidinae/scripts/

    #activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes.py {genes}
    
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/
        
    #Activate MrBayes
    ./mb  /home/kris/dypsidinae/3.MrBayes/{genes}_MrBayes_block.nex 

    mv *.log /home/kris/dypsidinae/3.MrBayes/
    
    
    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec) 

# ########################################################################################################################
# ##############################################---- MrBayes_2 ----#########################################################
# ########################################################################################################################


def MrBayes_2(path_in, genes):
    """Using MrBayes on genes that did not converge the first time"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/3.MrBayes/"+genes+"_2.log"]
    options = {'cores': "8", 'memory': "100g", 'walltime': "30:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts/

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes_2.py {genes}

    #Activate MrBayes
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/
    ./mb  /home/kris/dypsidinae/3.MrBayes/{genes}_MrBayes_block_2.nex
   
    mv *.log /home/kris/dypsidinae/3.MrBayes/
    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec)



# ########################################################################################################################
# ##############################################---- MrBayes_3 ----#########################################################
# ########################################################################################################################

def MrBayes_3(path_in, genes):
    """Using MrBayes to find the best model for genes that did not converge in the second round"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/"+genes+"_3.log"]
    options = {'cores': "16", 'memory': "100g", 'walltime': "12:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts/
    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes_3.py {genes}

    #Activate MrBayes 
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/  
    ./mb  /home/kris/dypsidinae/3.MrBayes/{genes}_MrBayes_block_3.nex

    mv *.log /home/kris/dypsidinae/3.MrBayes/

    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec)


# ################################################################################################################
# ##############################################---- MrBayes_4 ----#########################################################
# ########################################################################################################################


def MrBayes_4(path_in, genes):
    """Using MrBayes on genes that did not converge the third time"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/"+genes+"_4.log"]
    options = {'cores': "16", 'memory': "100g", 'walltime': "8:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts/

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes_4.py {genes}

    #Activate MrBayes
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/    
    ./mb  /home/kris/dypsidinae/3.MrBayes/{genes}_MrBayes_block_4.nex


    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec)

# ########################################################################################################################
# #####################################---- converting & discarding burnin ----############################################
# ########################################################################################################################

def treefile(path_in, output, genes):
    """Converting from nexus to newick format and making one file with genetrees from both MrBayes runs, where burnin is excluded"""
      
    inputs = [path_in + genes + ".fasta-out.nex.run1.t" , path_in + genes + ".fasta-out.nex.run2.t"]
    outputs = [path_in + output]
    options = {'cores': 2, 'memory': "100g", 'walltime': "6:00:00", 'account':"dypsidinae"}
    spec = """
    
    cd /home/kris/dypsidinae/scripts/

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate base

    python nexus_to_newick.py {genes}

        """.format(path_in = path_in, output=output, genes = genes)
    
    return (inputs, outputs, options, spec)

# ########################################################################################################################
# #####################################---- random treefiles ----#########################################################
# ########################################################################################################################

def random_tree_sets(path_in, number, path_out, output):
    """making random sets of gene trees"""
    inputs = [path_in]
    outputs = [path_out + output]
    options = {'cores': "2", 'memory': "40g", 'walltime': "8:00:00", 'account':"dypsidinae"}

    spec = """
    
    cd /home/kris/dypsidinae/scripts/

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env

    python3 random_tree_sets.py {number}

    cd {path_in}
    mv {number}random_trees.tre {path_out}
  
    
    """.format(path_in = path_in, number=number, path_out=path_out, output=output )
    
    return (inputs, outputs, options, spec)


# ########################################################################################################################
# ################---- Astral4 Tree Search on gene trees from the posteriror distribution ----############################
# ########################################################################################################################

#Showing Posterior Probabilities
def astral4(path_in, gene_tree_file,number, output):
    """Using Astral 4 to construct species trees from the 1000 gene tree files"""
    inputs = [path_in + gene_tree_file]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "100g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """
    cd {path_in}
    /home/kris/ASTER/bin/astral4 -t 16 -o {number}_astral_tree_probabilities.tre {gene_tree_file} 2>{number}_log_astral.out
    

    """.format(path_in = path_in, gene_tree_file = gene_tree_file, number=number, output=output)
    
    return (inputs, outputs, options, spec)



# #######################################################################################################################
# #####################################---- Astral into MrBayes ----#####################################################
# #########################################################################################################################

#Showing number of species trees supporting clades
def astral_2_MrBayes(path_in,input, output):
    """converting the 1000 Astral trees into one file for MrBayes"""
    inputs = [path_in+input]
    outputs = [path_in + output]
    options = {'cores': '1', 'memory': "40g", 'walltime': "4:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env

    python3 astral_for_mrbayes.py

    """.format(path_in = path_in, input=input, output=output)

    return (inputs, outputs, options, spec)

# ########################################################################################################################
# ##############################################---- MrBayes_consensus tree ----##########################################
# ########################################################################################################################

def consensus_tree(path_in, input, output):
    """Using MrBayes to finde consensus tree"""
    inputs = [path_in+input]
    outputs = [output]
    options = {'cores': "16", 'memory': "100g", 'walltime': "8:00:00", 'account':"dypsidinae"}

    spec = """

    #Activate MrBayes
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/    
    ./mb  concensus.nex


    """.format(path_in = path_in, input=input, output=output)

    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################


genes = []

                                                    
#Running IQTREE for trimmed and aligned files                                             
for i in range(0, len(genes)):
   gwf.target_from_template('Iqtree_'+str(i), iqtree(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/data/Dryad/B_main_analysis/A_alignments/"))  

# running AMAS to convert .fasta files to .nex files
for i in range(0, len(genes)):
   gwf.target_from_template('amas_'+str(i), amas_c(genes = genes[i], path_in="/home/kris/dypsidinae/data/nex_files/"))


#running MrBayes for 2.5 milion generationer
for i in range(0, len(genes)):
    gwf.target_from_template('MrBayes_test2'+str(i), MrBayes(genes = genes[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running MrBayes for 3.5 milion generationer
#for i in range(0, len(genes_v2)):
 #   gwf.target_from_template('MrBayes_2.1_'+str(i), MrBayes_2(genes = genes_v2[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running MrBayes for 5 milion generationer
for i in range(0, len(genes_v3)):
    gwf.target_from_template('MrBayes_3__'+str(i), MrBayes_3(genes = genes_v3[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running MrBayes for 7.5 milion generationer
for i in range(0, len(genes_v4)):
    gwf.target_from_template('MrBayes_4__'+str(i), MrBayes_4(genes = genes_v4[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running conversion to newick
for i in range(0, len(genes_converged)):
    gwf.target_from_template('treefile_'+str(i), treefile(genes = genes_converged[i],
                                                    path_in = "/home/kris/dypsidinae/data/nex_files/",
                                                    output= str(genes_converged[i])+"_genetree.tre"))
    
# making 1000 files containing one random gene tree from the posterior distribution of each gene
for i in range(1000):
    gwf.target_from_template('random_tree_sets_'+str(i), random_tree_sets(path_in = "/home/kris/dypsidinae/data/nex_files/",
                                                    number=i,
                                                    path_out="/home/kris/dypsidinae/Astral/",
                                                    output=str(i)+"random_trees.tre"))

# Running ASTRAL on the 1000 files with gene trees
for i in range(1000):
    gwf.target_from_template('astral4'+str(i), astral4(path_in = "/home/kris/dypsidinae/Astral/", number=i, gene_tree_file=str(i)+"random_trees.tre", output=str(i)+"_astral_tree_probabilities.tre"))

# Running ASTRAL4 on the 1000 trees made with ASTRAL
gwf.target_from_template('astral', astral(path_in = "/home/kris/dypsidinae/Astral/",
                                           output="consensus_astral_tree.tre"))


# Running wASTRAL
gwf.target_from_template('wastral', wastral(path_in = "/home/kris/dypsidinae/Astral/",
                                           output="wastral_tree.tre"))


# Running ASTRAL
gwf.target_from_template('astral_gene_supporting', astral_gene_supporting(path_in = "/home/kris/dypsidinae/Astral/",
                                                    input="wastral_tree.tre",
                                                    output="wastral_scored.tre"))
#making file for MrBayes
gwf.target_from_template('astral_2_MrBayes, astral_2_MrBayes(path_in = "/home/kris/dypsidinae/Astral/",
                                                    input= *"_astral_tree_probabilities.tre",
                                                    output="astral_for_mrbayes.run1.t"))

#making consensus tree with MrBayes
gwf.target_from_template('consensus_tree', consensus_tree(path_in= "/home/kris/dypsidinae/Astral/", input= "astral_for_mrbayes.run1.t",
                                                    output="astral_for_mrbayes.run1.t.con.tre"))
