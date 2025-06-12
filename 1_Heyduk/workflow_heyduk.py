
from os import O_SYNC, name
from gwf import Workflow
import os.path
import csv

gwf = Workflow()

   
# ########################################################################################################################
# ##############################################---- IQTREE ----##########################################################
# ########################################################################################################################

def iqtree(path_in, genes):
    """Using IQTREE to find the best model"""
    inputs = [path_in+genes+".fasta", path_in+genes+".part"]
    outputs = ["/home/kris/dypsidinae/1.IQtree_2/"+genes+".part.treefile"]
    options = {'cores': 2, 'memory': "40g", 'walltime': "12:00:00", 'account':"dypsidinae"}

    spec = """
     
    cd /home/kris/dypsidinae/data/Dryad/B_main_analysis/A_alignments

        
    #Activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate iqtree
        
    iqtree2 -s {genes}.fasta -T AUTO -m TESTONLY -p {genes}.part -mset mrbayes  
   
    mv *.part.ckp.gz /home/kris/dypsidinae/1.IQtree_2
    mv *aligned_noempty.fasta-out_clean.part.best_model.nex /home/kris/dypsidinae/1.IQtree_2
    mv *.part.iqtree /home/kris/dypsidinae/1.IQtree_2
    mv *.part.model.gz /home/kris/dypsidinae/1.IQtree_2
    mv *.part.treefile /home/kris/dypsidinae/1.IQtree_2
    mv *.fasta-out_clean.part.best_scheme /home/kris/dypsidinae/1.IQtree_2 
    mv *.fasta-out_clean.part.best_scheme.nex /home/kris/dypsidinae/1.IQtree_2
    mv *.fasta-out_clean.part.log /home/kris/dypsidinae/1.IQtree_2
    mv *.fasta-out_clean.part.model.gz /home/kris/dypsidinae/1.IQtree_2
    """.format(path_in = path_in, genes = genes)

    return (inputs, outputs, options, spec) 
 
# ########################################################################################################################
# ##############################################---- AMAS_convert ----##########################################################
# ########################################################################################################################

def amas_c(path_in, genes):
    """Using AMAS to convert files"""
    inputs = [path_in+genes+".fasta"]
    outputs = ["/home/kris/dypsidinae/data/nex_files/"+genes+".fasta-out.nex"]
    options = {'cores': 2, 'memory': "40g", 'walltime': "12:00:00", 'account':"dypsidinae"}

    spec = """
     
    cd /home/kris/dypsidinae/data/nex_files
        
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
    outputs = ["/home/kris/dypsidinae/2.MrBayes/"+genes+".log"]
    options = {'cores':"2", 'memory': "120g", 'walltime': "36:00:00", 'account':"dypsidinae"}

    spec = """
    
    cd /home/kris/dypsidinae/scripts/

    #activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes.py {genes}
    
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/
        
    #Activate MrBayes
    ./mb  /home/kris/dypsidinae/2.MrBayes/{genes}_MrBayes_block.nex 

    mv *.log /home/kris/dypsidinae/2.MrBayes/
    
    
    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec) 

#Stop and check if genes have converged
# ########################################################################################################################
# ##############################################---- MrBayes_2 ----#########################################################
# ########################################################################################################################


def MrBayes_2(path_in, genes):
    """Using MrBayes on genes that did not converge the first time"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/2.MrBayes/"+genes+"_2.log"]
    options = {'cores': "8", 'memory': "100g", 'walltime': "30:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts/

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes_2.py {genes}

    #Activate MrBayes
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/
    ./mb  /home/kris/dypsidinae/2.MrBayes/{genes}_MrBayes_block_2.nex
   
    mv *.log /home/kris/dypsidinae/2.MrBayes/
    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec)


#Stop and check if genes have converged
# ########################################################################################################################
# ##############################################---- MrBayes_3 ----#########################################################
# ########################################################################################################################

def MrBayes_3(path_in, genes):
    """Using MrBayes to find the best model for genes that did not converge the second time"""
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
    ./mb  /home/kris/dypsidinae/2.MrBayes/{genes}_MrBayes_block_3.nex

    mv *.log /home/kris/dypsidinae/2.MrBayes/

    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec)

#Stop and check if genes have converged
# ################################################################################################################
# ##############################################---- MrBayes_4 ----#########################################################
# ########################################################################################################################


def MrBayes_4(path_in, genes):
    """Using MrBayes on genes that did not converge the third time"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/"+genes+"_4.log"]
    options = {'cores': "16", 'memory': "100g", 'walltime': "8:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/sccripts/

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes_4.py {genes}

    #Activate MrBayes
    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/    
    ./mb  /home/kris/dypsidinae/2.MrBayes/{genes}_MrBayes_block_4.nex


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
    """making random set of genetrees"""
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
# #####################################---- Astral4 Tree Search ----#####################################################
# ########################################################################################################################

#Showing Posterior Probabilities
def astral4(path_in, gene_tree_file,number, output):
    """Using Astral to construct a species tree from the genetrees"""
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

def astral_2_MrBayes(path_in,output):
    """converting the 1000 Astral trees into one file for MrBayes"""
    inputs = [path_in]
    outputs = [path_in + output]
    options = {'cores': '1', 'memory': "40g", 'walltime': "4:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts

    #activating the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate base

    python astral_2_mrbayes.py

    """.format(path_in = path_in, output=output)

    return (inputs, outputs, options, spec)


# ########################################################################################################################
# ##############################################---- IQTREE ----##########################################################
# ########################################################################################################################

def iqtree(path_in, genes):
    """Using IQTREE to find the best model"""
    inputs = [path_in+genes+".fasta", path_in+genes+".part"]
    outputs = ["/home/kris/dypsidinae/4.IQtree/"+genes+".part.treefile"]
    options = {'cores': 2, 'memory': "40g", 'walltime': "12:00:00", 'account':"dypsidinae"}

    spec = """
     
    cd /home/kris/dypsidinae/data/Dryad/B_main_analysis/A_alignments

        
    #Activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate iqtree
        
    iqtree2 -s {genes}.fasta -T AUTO -m TEST  -p {genes}.part -mset mrbayes  -B 1000
   
    mv *.part.ckp.gz /home/kris/dypsidinae/4.IQtree
    mv *aligned_noempty.fasta-out_clean.part.best_model.nex /home/kris/dypsidinae/4.IQtree # slet linje
    mv *.part.iqtree /home/kris/dypsidinae/4.IQtree
    mv *.part.model.gz /home/kris/dypsidinae/4.IQtree
    mv *.part.treefile /home/kris/dypsidinae/4.IQtree
    mv *.fasta-out_clean.part.best_scheme /home/kris/dypsidinae/4.IQtree 
    mv *.fasta-out_clean.part.best_scheme.nex /home/kris/dypsidinae/4.IQtree
    mv *.fasta-out_clean.part.log /home/kris/dypsidinae/4.IQtree
    mv *.fasta-out_clean.part.model.gz /home/kris/dypsidinae/4.IQtree
    """.format(path_in = path_in, genes = genes)

    return (inputs, outputs, options, spec) 

# ########################################################################################################################
# ###########################---- wAstral Tree Search on gene trees from IQtree ----#####################################
# ########################################################################################################################

def ML_wastral(path_in, gene_tree_file,output):
    """Using wAstral to construct ML species trees from the gene tree files"""
    inputs = [path_in + gene_tree_file]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "100g", 'walltime': "48:00:00", 'account':"dypsidinae"}

    spec = """
    cd {path_in}
    /home/kris/ASTER/bin/wastral -t 16 -o astral_ML.tre {gene_tree_file} 2> astral_ML_log_astral.out

    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output)
    
    return (inputs, outputs, options, spec)

 # #######################################################################################################################
# #####################################---- consensus genetrees ----#####################################################
# #######################################################################################################################
#Finding the consensus gene trees
def MrB_sumt(path_in, genes, path_out):
    """Using MrBayes to find gene trees"""
    inputs = [path_in+genes]
    outputs = [path_out+genes]
    options = {'cores':"2", 'memory': "120g", 'walltime': "166:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts/


    #activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 mb_writing_sumt.py {genes}

    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/

    #Activate MrBayes
    ./mb  {path_out}{genes}_mb_sumt.nex
    
    mv {genes}.con.tre {path_out}

    """.format(path_in = path_in, genes= genes, path_out= path_out)

    return (inputs, outputs, options, spec)


# ########################################################################################################################
# ###########################---- making one gene tree file from MrBayes ----#############################################
# ########################################################################################################################
#Making one file in newick format with one gene trees per gene from MrBayes
def MrBayes_one_file(path_in, output):
    """Using wAstral to construct species trees from the gene tree files from MrBayes"""
    inputs = [path_in]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "100g", 'walltime': "48:00:00", 'account':"dypsidinae"}

    spec = """
   
   cd /home/kris/dypsidinae/scripts/

    #activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate base

    python Make_one_treefile.py 

    """.format(path_in = path_in, output=output)
    
    return (inputs, outputs, options, spec)


# ########################################################################################################################
# ###########################---- wAstral Tree Search on gene trees from MrBayes ----#####################################
# ########################################################################################################################

#Making a species tree from the gene trees
def MrBayes_wastral(path_in, gene_tree_file, output):
    """Using wAstral to construct species trees from the gene tree files from MrBayes"""
    inputs = [path_in + gene_tree_file]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "100g", 'walltime': "48:00:00", 'account':"dypsidinae"}

    spec = """
   

    cd {path_in}
    /home/kris/ASTER/bin/wastral -t 16 -o heyduk_astral_Bayes.tre {gene_tree_file} 2> heyduk_astral_Bayes_log.out

    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output)
    
    return (inputs, outputs, options, spec)

# ########################################################################################################################
# ###########################----gene trees from converged genes for IQtree ----##########################################
# ########################################################################################################################
#Making one file  with all the gene trees from IQtree with the genes that converged using MrBayes
def IQtree_converged_genes(path_in, output):
    """selecting gene trees that converged in MrBayes"""
    inputs = [path_in]
    outputs = [path_in + output]
    options = {'cores': '4', 'memory': "100g", 'walltime': "48:00:00", 'account':"dypsidinae"}

    spec = """

   cd /home/kris/dypsidinae/scripts/

    #activate the enviroment
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate base

    python genes_converged_IQtree.py

    """.format(path_in = path_in, output=output)

    return (inputs, outputs, options, spec)




# ########################################################################################################################
# ###########################---- wAstral Tree Search on gene trees from IQtree that converged ----#######################
# ########################################################################################################################

def ML_converged_wastral(path_in, gene_tree_file,output):
    """Using wAstral to construct ML species trees from the gene tree files"""
    inputs = [path_in + gene_tree_file]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "100g", 'walltime': "48:00:00", 'account':"dypsidinae"}

    spec = """
    cd {path_in}
    /home/kris/ASTER/bin/wastral -t 16 -o astral_ML_converged.tre {gene_tree_file} 2> astral_ML_converged_log.out

    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output)

    return (inputs, outputs, options, spec)


########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################


sp = ['95', '2015', '2', '78', '3', '4', '202', '2012', '6', '8', '10', '11', '12', '14', '149', '150', '82', '21', '204', '2033', '28', '91', '2051', '2034', '40', '41', '195', '44', '45', '102', '103', '105', '109', '51', '111', '52', '171', '112', '53', '2052', '174', '113', '57', '59', '115', '61', '63', '181', '120', '67', '123', '2014', '92', '183', '71', '72', '100', '1', '5', '144', '146', '7', '79', '80', '13', '15', '17', '18', '19', '20', '81', '148', '83', '151', '152', '22', '23', '24', '153', '203', '25', '27', '154', '29', '30', '84', '32', '33', '85', '86', '87', '34', '89', '90', '155', '156', '201', '35', '36', '158', '37', '38', '160', '161', '162', '39', '163', '193', '194', '97', '42', '98', '43', '99', '166', '46', '101', '104', '47', '106', '107', '108', '48', '169', '49', '50', '172', '54', '173', '175', '114', '58', '176', '60', '178', '179', '62', '2053', '117', '118', '182', '65', '121', '66', '122', '68', '69', '124', '125', '70', '126', '184', '73', '185', '191', '1011', '1012', '205', '198', '199', '196', '187', '77', '9', '26', '31', '88', '170', '110', '55', '116', '128', '74']

genes = ['reduced_2339_aligned_noempty.fasta-out_clean', 'reduced_514_aligned_noempty.fasta-out_clean' , 'reduced_1007_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_51_aligned_noempty.fasta-out_clean', 'reduced_1013_aligned_noempty.fasta-out_clean', 'reduced_2370_aligned_noempty.fasta-out_clean', 'reduced_52_aligned_noempty.fasta-out_clean', 'reduced_1017_aligned_noempty.fasta-out_clean', 'reduced_2377_aligned_noempty.fasta-out_clean', 'reduced_556_aligned_noempty.fasta-out_clean', 'reduced_1020_aligned_noempty.fasta-out_clean', 'reduced_237_aligned_noempty.fasta-out_clean', 'reduced_563_aligned_noempty.fasta-out_clean', 'reduced_1025_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean', 'reduced_576_aligned_noempty.fasta-out_clean', 'reduced_1035_aligned_noempty.fasta-out_clean', 'reduced_240_aligned_noempty.fasta-out_clean', 'reduced_587_aligned_noempty.fasta-out_clean', 'reduced_1050_aligned_noempty.fasta-out_clean', 'reduced_2459_aligned_noempty.fasta-out_clean', 'reduced_604_aligned_noempty.fasta-out_clean', 'reduced_1052_aligned_noempty.fasta-out_clean', 'reduced_245_aligned_noempty.fasta-out_clean', 'reduced_609_aligned_noempty.fasta-out_clean', 'reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_61_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_1168_aligned_noempty.fasta-out_clean', 'reduced_252e_aligned_noempty.fasta-out_clean', 'reduced_630_aligned_noempty.fasta-out_clean', 'reduced_1171_aligned_noempty.fasta-out_clean', 'reduced_252p_aligned_noempty.fasta-out_clean', 'reduced_637_aligned_noempty.fasta-out_clean', 'reduced_1197_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_673_aligned_noempty.fasta-out_clean', 'reduced_1201_aligned_noempty.fasta-out_clean', 'reduced_2550_aligned_noempty.fasta-out_clean', 'reduced_680_aligned_noempty.fasta-out_clean', 'reduced_120_aligned_noempty.fasta-out_clean', 'reduced_2561_aligned_noempty.fasta-out_clean', 'reduced_717_aligned_noempty.fasta-out_clean', 'reduced_122_aligned_noempty.fasta-out_clean', 'reduced_257_aligned_noempty.fasta-out_clean', 'reduced_727_aligned_noempty.fasta-out_clean', 'reduced_125_aligned_noempty.fasta-out_clean', 'reduced_267_aligned_noempty.fasta-out_clean', 'reduced_732_aligned_noempty.fasta-out_clean', 'reduced_12_aligned_noempty.fasta-out_clean', 'reduced_269_aligned_noempty.fasta-out_clean', 'reduced_736_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_277_aligned_noempty.fasta-out_clean', 'reduced_740_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_280_aligned_noempty.fasta-out_clean', 'reduced_743_aligned_noempty.fasta-out_clean', 'reduced_1484_aligned_noempty.fasta-out_clean', 'reduced_281_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_148_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_758_aligned_noempty.fasta-out_clean', 'reduced_1494_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_785_aligned_noempty.fasta-out_clean', 'reduced_14_aligned_noempty.fasta-out_clean', 'reduced_293_aligned_noempty.fasta-out_clean', 'reduced_790_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean', 'reduced_793_aligned_noempty.fasta-out_clean', 'reduced_1615_aligned_noempty.fasta-out_clean', 'reduced_305_aligned_noempty.fasta-out_clean', 'reduced_7_aligned_noempty.fasta-out_clean', 'reduced_168_aligned_noempty.fasta-out_clean', 'reduced_308_aligned_noempty.fasta-out_clean', 'reduced_807_aligned_noempty.fasta-out_clean', 'reduced_17_aligned_noempty.fasta-out_clean', 'reduced_310_aligned_noempty.fasta-out_clean', 'reduced_808_aligned_noempty.fasta-out_clean', 'reduced_1801_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_1815_aligned_noempty.fasta-out_clean', 'reduced_326_aligned_noempty.fasta-out_clean', 'reduced_825_aligned_noempty.fasta-out_clean', 'reduced_182_aligned_noempty.fasta-out_clean', 'reduced_32e_aligned.fasta-out_clean', 'reduced_82_aligned_noempty.fasta-out_clean', 'reduced_1842_aligned_noempty.fasta-out_clean', 'reduced_32s_aligned_noempty.fasta-out_clean', 'reduced_83_aligned_noempty.fasta-out_clean', 'reduced_1854_aligned_noempty.fasta-out_clean', 'reduced_332_aligned_noempty.fasta-out_clean', 'reduced_84_aligned_noempty.fasta-out_clean', 'reduced_1877_aligned_noempty.fasta-out_clean', 'reduced_357_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_1901_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_191_aligned_noempty.fasta-out_clean', 'reduced_362_aligned_noempty.fasta-out_clean', 'reduced_872_aligned_noempty.fasta-out_clean', 'reduced_194_aligned_noempty.fasta-out_clean', 'reduced_363_aligned_noempty.fasta-out_clean', 'reduced_874_aligned_noempty.fasta-out_clean', 'reduced_197_aligned_noempty.fasta-out_clean', 'reduced_369_aligned_noempty.fasta-out_clean', 'reduced_883e_aligned_noempty.fasta-out_clean', 'reduced_1986_aligned_noempty.fasta-out_clean', 'reduced_378e_aligned_noempty.fasta-out_clean', 'reduced_883n_aligned_noempty.fasta-out_clean', 'reduced_201_aligned_noempty.fasta-out_clean', 'reduced_378s_aligned_noempty.fasta-out_clean', 'reduced_886_aligned_noempty.fasta-out_clean', 'reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_38_aligned_noempty.fasta-out_clean', 'reduced_88_aligned_noempty.fasta-out_clean', 'reduced_204s_aligned_noempty.fasta-out_clean', 'reduced_391_aligned_noempty.fasta-out_clean', 'reduced_897_aligned_noempty.fasta-out_clean', 'reduced_2056_aligned_noempty.fasta-out_clean', 'reduced_392_aligned_noempty.fasta-out_clean', 'reduced_89_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_415_aligned_noempty.fasta-out_clean', 'reduced_938_aligned_noempty.fasta-out_clean', 'reduced_215_aligned_noempty.fasta-out_clean', 'reduced_417_aligned_noempty.fasta-out_clean', 'reduced_948_aligned_noempty.fasta-out_clean', 'reduced_2164_aligned_noempty.fasta-out_clean', 'reduced_421_aligned_noempty.fasta-out_clean', 'reduced_94_aligned_noempty.fasta-out_clean', 'reduced_218_aligned_noempty.fasta-out_clean', 'reduced_449_aligned_noempty.fasta-out_clean', 'reduced_950_aligned_noempty.fasta-out_clean', 'reduced_21_aligned_noempty.fasta-out_clean', 'reduced_464_aligned_noempty.fasta-out_clean', 'reduced_958_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_484_aligned_noempty.fasta-out_clean', 'reduced_964_aligned_noempty.fasta-out_clean', 'reduced_225_aligned_noempty.fasta-out_clean', 'reduced_490_aligned_noempty.fasta-out_clean', 'reduced_977_aligned_noempty.fasta-out_clean', 'reduced_226_aligned_noempty.fasta-out_clean', 'reduced_497_aligned_noempty.fasta-out_clean', 'reduced_982_aligned_noempty.fasta-out_clean', 'reduced_2291_aligned_noempty.fasta-out_clean', 'reduced_4_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_231_aligned_noempty.fasta-out_clean', 'reduced_508_aligned_noempty.fasta-out_clean', 'reduced_989_aligned_noempty.fasta-out_clean']

genes_v2 = ['reduced_191_aligned_noempty.fasta-out_clean', 'reduced_252e_aligned_noempty.fasta-out_clean', 'reduced_21_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_362_aligned_noempty.fasta-out_clean', 'reduced_883e_aligned_noempty.fasta-out_clean', 'reduced_226_aligned_noempty.fasta-out_clean', 'reduced_563_aligned_noempty.fasta-out_clean', 'reduced_89_aligned_noempty.fasta-out_clean']

genes_v3 =['reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean','reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_1020_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_12_aligned_noempty.fasta-out_clean', 'reduced_1007_aligned_noempty.fasta-out_clean', 'reduced_83_aligned_noempty.fasta-out_clean', 'reduced_1168_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_2370_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean','reduced_191_aligned_noempty.fasta-out_clean', 'reduced_391_aligned_noempty.fasta-out_clean', 'reduced_673_aligned_noempty.fasta-out_clean', 'reduced_369_aligned_noempty.fasta-out_clean', 'reduced_1171_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_1201_aligned_noempty.fasta-out_clean', 'reduced_897_aligned_noempty.fasta-out_clean', 'reduced_758_aligned_noempty.fasta-out_clean'] 

genes_v4= ['reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_191_aligned_noempty.fasta-out_clean', 'reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean']

genes_converged=['reduced_2339_aligned_noempty.fasta-out_clean', 'reduced_514_aligned_noempty.fasta-out_clean', 'reduced_1007_aligned_noempty.fasta-out_clean', 'reduced_51_aligned_noempty.fasta-out_clean', 'reduced_1013_aligned_noempty.fasta-out_clean', 'reduced_2370_aligned_noempty.fasta-out_clean', 'reduced_52_aligned_noempty.fasta-out_clean', 'reduced_1017_aligned_noempty.fasta-out_clean', 'reduced_2377_aligned_noempty.fasta-out_clean', 'reduced_556_aligned_noempty.fasta-out_clean', 'reduced_1020_aligned_noempty.fasta-out_clean', 'reduced_237_aligned_noempty.fasta-out_clean', 'reduced_563_aligned_noempty.fasta-out_clean', 'reduced_1025_aligned_noempty.fasta-out_clean', 'reduced_576_aligned_noempty.fasta-out_clean', 'reduced_1035_aligned_noempty.fasta-out_clean', 'reduced_240_aligned_noempty.fasta-out_clean', 'reduced_587_aligned_noempty.fasta-out_clean', 'reduced_1050_aligned_noempty.fasta-out_clean', 'reduced_2459_aligned_noempty.fasta-out_clean', 'reduced_604_aligned_noempty.fasta-out_clean', 'reduced_1052_aligned_noempty.fasta-out_clean', 'reduced_245_aligned_noempty.fasta-out_clean', 'reduced_609_aligned_noempty.fasta-out_clean', 'reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_61_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_1168_aligned_noempty.fasta-out_clean', 'reduced_252e_aligned_noempty.fasta-out_clean', 'reduced_630_aligned_noempty.fasta-out_clean', 'reduced_1171_aligned_noempty.fasta-out_clean', 'reduced_252p_aligned_noempty.fasta-out_clean', 'reduced_637_aligned_noempty.fasta-out_clean', 'reduced_1197_aligned_noempty.fasta-out_clean', 'reduced_673_aligned_noempty.fasta-out_clean', 'reduced_1201_aligned_noempty.fasta-out_clean', 'reduced_2550_aligned_noempty.fasta-out_clean', 'reduced_680_aligned_noempty.fasta-out_clean', 'reduced_120_aligned_noempty.fasta-out_clean', 'reduced_2561_aligned_noempty.fasta-out_clean', 'reduced_717_aligned_noempty.fasta-out_clean', 'reduced_122_aligned_noempty.fasta-out_clean', 'reduced_257_aligned_noempty.fasta-out_clean', 'reduced_727_aligned_noempty.fasta-out_clean', 'reduced_125_aligned_noempty.fasta-out_clean', 'reduced_267_aligned_noempty.fasta-out_clean', 'reduced_732_aligned_noempty.fasta-out_clean', 'reduced_12_aligned_noempty.fasta-out_clean', 'reduced_269_aligned_noempty.fasta-out_clean', 'reduced_736_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_277_aligned_noempty.fasta-out_clean', 'reduced_740_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_280_aligned_noempty.fasta-out_clean', 'reduced_743_aligned_noempty.fasta-out_clean', 'reduced_1484_aligned_noempty.fasta-out_clean', 'reduced_281_aligned_noempty.fasta-out_clean', 'reduced_148_aligned_noempty.fasta-out_clean', 'reduced_758_aligned_noempty.fasta-out_clean', 'reduced_1494_aligned_noempty.fasta-out_clean', 'reduced_785_aligned_noempty.fasta-out_clean', 'reduced_14_aligned_noempty.fasta-out_clean', 'reduced_293_aligned_noempty.fasta-out_clean', 'reduced_790_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_793_aligned_noempty.fasta-out_clean', 'reduced_1615_aligned_noempty.fasta-out_clean', 'reduced_305_aligned_noempty.fasta-out_clean', 'reduced_7_aligned_noempty.fasta-out_clean', 'reduced_168_aligned_noempty.fasta-out_clean', 'reduced_308_aligned_noempty.fasta-out_clean', 'reduced_807_aligned_noempty.fasta-out_clean', 'reduced_17_aligned_noempty.fasta-out_clean', 'reduced_310_aligned_noempty.fasta-out_clean', 'reduced_808_aligned_noempty.fasta-out_clean', 'reduced_1801_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_1815_aligned_noempty.fasta-out_clean', 'reduced_326_aligned_noempty.fasta-out_clean', 'reduced_825_aligned_noempty.fasta-out_clean', 'reduced_182_aligned_noempty.fasta-out_clean', 'reduced_32e_aligned.fasta-out_clean', 'reduced_82_aligned_noempty.fasta-out_clean', 'reduced_1842_aligned_noempty.fasta-out_clean', 'reduced_32s_aligned_noempty.fasta-out_clean', 'reduced_83_aligned_noempty.fasta-out_clean', 'reduced_1854_aligned_noempty.fasta-out_clean', 'reduced_332_aligned_noempty.fasta-out_clean', 'reduced_84_aligned_noempty.fasta-out_clean', 'reduced_1877_aligned_noempty.fasta-out_clean', 'reduced_357_aligned_noempty.fasta-out_clean', 'reduced_1901_aligned_noempty.fasta-out_clean', 'reduced_362_aligned_noempty.fasta-out_clean', 'reduced_872_aligned_noempty.fasta-out_clean', 'reduced_194_aligned_noempty.fasta-out_clean', 'reduced_363_aligned_noempty.fasta-out_clean', 'reduced_874_aligned_noempty.fasta-out_clean', 'reduced_197_aligned_noempty.fasta-out_clean', 'reduced_369_aligned_noempty.fasta-out_clean', 'reduced_883e_aligned_noempty.fasta-out_clean', 'reduced_1986_aligned_noempty.fasta-out_clean', 'reduced_378e_aligned_noempty.fasta-out_clean', 'reduced_883n_aligned_noempty.fasta-out_clean', 'reduced_201_aligned_noempty.fasta-out_clean', 'reduced_378s_aligned_noempty.fasta-out_clean', 'reduced_886_aligned_noempty.fasta-out_clean', 'reduced_38_aligned_noempty.fasta-out_clean', 'reduced_88_aligned_noempty.fasta-out_clean', 'reduced_204s_aligned_noempty.fasta-out_clean', 'reduced_391_aligned_noempty.fasta-out_clean', 'reduced_897_aligned_noempty.fasta-out_clean', 'reduced_2056_aligned_noempty.fasta-out_clean', 'reduced_392_aligned_noempty.fasta-out_clean', 'reduced_89_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_415_aligned_noempty.fasta-out_clean', 'reduced_938_aligned_noempty.fasta-out_clean', 'reduced_215_aligned_noempty.fasta-out_clean', 'reduced_417_aligned_noempty.fasta-out_clean', 'reduced_948_aligned_noempty.fasta-out_clean', 'reduced_2164_aligned_noempty.fasta-out_clean', 'reduced_421_aligned_noempty.fasta-out_clean', 'reduced_94_aligned_noempty.fasta-out_clean', 'reduced_218_aligned_noempty.fasta-out_clean', 'reduced_449_aligned_noempty.fasta-out_clean', 'reduced_950_aligned_noempty.fasta-out_clean', 'reduced_21_aligned_noempty.fasta-out_clean', 'reduced_464_aligned_noempty.fasta-out_clean', 'reduced_958_aligned_noempty.fasta-out_clean', 'reduced_484_aligned_noempty.fasta-out_clean', 'reduced_964_aligned_noempty.fasta-out_clean', 'reduced_225_aligned_noempty.fasta-out_clean', 'reduced_490_aligned_noempty.fasta-out_clean', 'reduced_977_aligned_noempty.fasta-out_clean', 'reduced_226_aligned_noempty.fasta-out_clean', 'reduced_497_aligned_noempty.fasta-out_clean', 'reduced_982_aligned_noempty.fasta-out_clean', 'reduced_2291_aligned_noempty.fasta-out_clean', 'reduced_4_aligned_noempty.fasta-out_clean', 'reduced_231_aligned_noempty.fasta-out_clean', 'reduced_508_aligned_noempty.fasta-out_clean', 'reduced_989_aligned_noempty.fasta-out_clean']



                                                    
#Running IQTREE for trimmed and aligned files                                             
for i in range(0, len(genes)):
   gwf.target_from_template('Iqtree_'+str(i), iqtree(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/data/Dryad/B_main_analysis/A_alignments/"))  

# running AMAS to convert .fasta files to .nex files
for i in range(0, len(genes)):
   gwf.target_from_template('amas_'+str(i), amas_c(genes = genes[i], path_in="/home/kris/dypsidinae/data/nex_files/"))

"""
#running MrBayes
for i in range(0, len(genes)):
    gwf.target_from_template('MrBayes_test2'+str(i), MrBayes(genes = genes[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running MrBayes
#for i in range(0, len(genes_v2)):
 #   gwf.target_from_template('MrBayes_2.1_'+str(i), MrBayes_2(genes = genes_v2[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running MrBayes
for i in range(0, len(genes_v3)):
    gwf.target_from_template('MrBayes_3__'+str(i), MrBayes_3(genes = genes_v3[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

#running MrBayes
for i in range(0, len(genes_v4)):
    gwf.target_from_template('MrBayes_4__'+str(i), MrBayes_4(genes = genes_v4[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

# make treefile
for i in range(len(genes)):
    gwf.target_from_template('rm_burnin_'+str(i),rm_burnin(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/data/nex_files/",
						    output_1="_1_tree_file.nex",
						    output_2= "_2_tree_file.nex"))

#running conversion to newick
for i in range(0, len(genes_converged)):
    gwf.target_from_template('treefile_'+str(i), treefile(genes = genes_converged[i],
                                                    path_in = "/home/kris/dypsidinae/data/nex_files/",
                                                    output= str(genes_converged[i])+"_genetree.tre"))

for i in range(1000):
    gwf.target_from_template('random_tree_sets_'+str(i), random_tree_sets(path_in = "/home/kris/dypsidinae/data/nex_files/",
                                                    number=i,
                                                    path_out="/home/kris/dypsidinae/3.Astral/",
                                                    output=str(i)+"random_trees.tre"))
# Running ASTRAL 
for i in range(1000):
    gwf.target_from_template('astral4'+str(i), astral4(path_in = "/home/kris/dypsidinae/3.Astral/", number=i, gene_tree_file=str(i)+"random_trees.tre", output=str(i)+"_astral_tree_probabilities.tre"))


#making file for MrBayes
for i in range(1000):
    gwf.target_from_template('astral_2_MrBayes_'+str(i), astral_2_MrBayes(path_in = "/home/kris/dypsidinae/Astral/",input=str(i)+"_astral_tree_probabilities.tre", number=i,
                                                    output=str(i)+"astral_for_mrbayes.run1.t"))


#making file for MrBayes
gwf.target_from_template('astral_2_MrBayes', astral_2_MrBayes(path_in = "/home/kris/dypsidinae/3.Astral/",
                                                    output="astral_for_mrbayes.tre"))
"""
#Running IQTREE for trimmed and aligned files to find ML genetrees
for i in range(0, len(genes)):
   gwf.target_from_template('ML_tree_'+str(i), ML_tree(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/data/Dryad/B_main_analysis/A_alignments/",
                                                    path_out="/home/kris/dypsidinae/4.IQtree/")) 

#Running wASTRAL
gwf.target_from_template('ML_wastral_'+str(i), ML_wastral(path_in = "/home/kris/dypsidinae/4.IQtree/", gene_tree_file="ML_gene_trees.nex", output= "astral_ML.tre"))




#Make one file with a gene tree for each gene
gwf.target_from_template('MrBayes_one_file', MrBayes_one_file(path_in = "/home/kris/dypsidinae/2.MrBayes/", output="MB_gene_trees.tre"))

#Running wASTRAL
#gwf.target_from_template('MB_wastral_'+str(i), MrBayes_wastral(path_in = "/home/kris/dypsidinae/B_phylopalm/3.MrBayes/", gene_tree_file="MB_gene_trees.tre", output= "astral_Bayes.tre"))


#Take genes that converged in MrBayes
gwf.target_from_template('IQtree_converged_genes', IQtree_converged_genes(path_in = "/home/kris/dypsidinae/4.IQtree/", output="only_converged_gene_trees.nex"))

#running wASTRAL
gwf.target_from_template('ML_converged_wastral', ML_converged_wastral(path_in = "/home/kris/dypsidinae/4.IQtree/", gene_tree_file="only_converged_gene_trees.nex", output= "astral_ML_converged.tre"))


# #########################################################################################################################
# #####################################---- Scoring ML converged Astral Tree  ----#########################################
# #########################################################################################################################

#Scoring existing tree with ASTRAL
def astral_scoring_ML_converged(path_in,input, output):
    """Using Astral to score existing tree"""
    inputs = [path_in+input]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "40g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/Astral/

    java -D"java.library.path=lib/" -jar astral.5.7.8.jar -q {path_in}{input} -t 2 -o {path_in}scored_wastral_ML_converged.tre -i {path_in}only_converged_gene_trees.nex 2>scored_ML_converged_log.out

    mv scored_wastral_ML_converged.tre /home/kris/dypsidinae/4.IQtree/
    mv scored_ML_converged_log.out /home/kris/dypsidinae/4.IQtree/


    """.format(path_in = path_in, input=input, output=output)

    return (inputs, outputs, options, spec)

# Running ASTRAL
gwf.target_from_template('astral_scoring_ML_converged', astral_scoring_ML_converged(path_in = "/home/kris/dypsidinae/4.IQtree/",
                                                    input="astral_ML_converged.tre",
                                                    output="scored_wastral_ML_converged.tre"))


# #########################################################################################################################
# #####################################---- Scoring ML Astral Tree  ----###################################################
# #########################################################################################################################

#Scoring existing tree with ASTRAL
def astral_scoring_ML(path_in,input, output):
    """Using Astral to score existing tree"""
    inputs = [path_in+input]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "40g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/Astral/

    java -D"java.library.path=lib/" -jar astral.5.7.8.jar -q {path_in}{input} -t 2 -o {path_in}scored_wastral_ML.tre -i {path_in}gene_trees.nex 2>scored_ML_log.out

    mv scored_wastral_ML.tre /home/kris/dypsidinae/4.IQtree/
    mv scored_ML_log.out /home/kris/dypsidinae/4.IQtree/


    """.format(path_in = path_in, input=input, output=output)

    return (inputs, outputs, options, spec)


# Running ASTRAL
gwf.target_from_template('astral_scoring_ML', astral_scoring_ML(path_in = "/home/kris/dypsidinae/4.IQtree/",
                                                    input="astral_ML.tre",
                                                    output="scored_wastral_ML.tre"))


# #########################################################################################################################
# #####################################---- Scoring ML freerate Astral Tree  ----##########################################
# #########################################################################################################################

#Scoring existing tree with ASTRAL
def astral_scoring_ML_freerate(path_in,input, output):
    """Using Astral to score existing tree"""
    inputs = [path_in+input]
    outputs = [path_in + output]
    options = {'cores': '16', 'memory': "40g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/Astral/

    java -D"java.library.path=lib/" -jar astral.5.7.8.jar -q {path_in}{input} -t 2 -o {path_in}scored_wastral_ML_freerate.tre -i {path_in}gene_trees_freerate.nex. 2>scored_ML_freerate_log.out

    mv scored_wastral_ML_freerate.tre /home/kris/dypsidinae/4.IQtree_freerate/
    mv scored_ML_freerate_log.out /home/kris/dypsidinae/4.IQtree_freerate/


    """.format(path_in = path_in, input=input, output=output)

    return (inputs, outputs, options, spec)


# Running ASTRAL
gwf.target_from_template('astral_scoring_ML_freerate', astral_scoring_ML_freerate(path_in = "/home/kris/dypsidinae/4.IQtree_freerate/",
                                                    input="astral_ML_freerate.tre",
                                                    output="scored_wastral_ML_freerate.tre"))

astral_ML_freerate.tre                                 
gene_trees_ML_freerate.nex 
