'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used run Modelfinder with IQtree followed by MrBayes and Astral for pylogenetic inference
------------------------------------------------------------------------------------------------------------------------
This code is inspired by the species_workflow.py by Oscar Wrisberg
Edited by Kristine N Kristensen
------------------------------------------------------------------------------------------------------------------------
'''

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
    """Using MrBayes to find the best model"""
    inputs = [path_in+genes+".fasta-out.nex"]
    outputs = ["/home/kris/dypsidinae/2.MrBayes/"+genes+".log"]
    options = {'cores': 8, 'memory': "100g", 'walltime': "12:00:00", 'account':"dypsidinae"}

    spec = """

    cd /home/kris/dypsidinae/scripts/
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env
    python3 MrBayes_2.py {genes}

    cd /home/kris/dypsidinae/programs/MrBayes-3.2.7a/bin/

    #Activate MrBayes
    ./mb  /home/kris/dypsidinae/2.MrBayes/{genes}_MrBayes_block.nex

    mv *.log /home/kris/dypsidinae/2.MrBayes/
   


    """.format(path_in = path_in, genes=genes)

    return (inputs, outputs, options, spec)




# ########################################################################################################################
# #####################################---- Astral Tree Search ----#####################################################
# ########################################################################################################################

#Here you should stop and go to folder /home/kris/dypsidinae/MrBayes  and do cat *treefile???hedder noget andet > gene_trees.nex

#Showing Posterior Probabilities
def astral_2(path_in, gene_tree_file, output, genes):
    """Using Astral to construct a species tree from the genetrees"""
    inputs = [path_in]
    outputs = [path_in + output]
    options = {'cores': 10, 'memory': "100g", 'walltime': "12:00:00", 'account':"dypsidinae"}

    spec = """
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    cd /home/kris/Astral/
  
    java -D"java.library.path=lib/" -jar astral.5.15.5.jar -i /home/kris/dypsidinae/A/3.Astral/gene_trees.nex -o astral_tree_probabilities.tre 2> log_posterior_probability.out
    
    mv astral_tree_probabilities.tre /home/kris/dypsidinae/3.Astral/
    mv log_posterior_probability.out /home/kris/dypsidinae/3.Astral/
    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output, genes = genes)
    
    return (inputs, outputs, options, spec)
        
#Showing number of genes supporting clades   
def astral_gene_supporting(path_in, gene_tree_file, output, genes):
    """Using Astral to construct a species tree based on the genetrees"""
    inputs = [path_in+"gene_trees.nex"]
    outputs = [path_in + output]
    options = {'cores': 10, 'memory': "40g", 'walltime': "24:00:00", 'account':"dypsidinae"}

    spec = """
    source /home/kris/miniconda3/etc/profile.d/conda.sh
    cd /home/kris/Astral/
  
    java -D"java.library.path=lib/" -jar astral.5.15.5.jar -i /home/kris/dypsidinae/A/3.Astral/gene_trees.nex -o astral_gene_probability.tre -t 2 2> log_gene_effective.out
      
    mv astral_gene_probability.tre /home/kris/Dypsidinae/A/3.Astral/
    mv log_gene_effective.out /home/kris/dypsidinae/A/3.Astral/
    
    """.format(path_in = path_in, gene_tree_file = gene_tree_file, output=output, genes = genes)

    return (inputs, outputs, options, spec)


########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################


sp = ['95', '2015', '2', '78', '3', '4', '202', '2012', '6', '8', '10', '11', '12', '14', '149', '150', '82', '21', '204', '2033', '28', '91', '2051', '2034', '40', '41', '195', '44', '45', '102', '103', '105', '109', '51', '111', '52', '171', '112', '53', '2052', '174', '113', '57', '59', '115', '61', '63', '181', '120', '67', '123', '2014', '92', '183', '71', '72', '100', '1', '5', '144', '146', '7', '79', '80', '13', '15', '17', '18', '19', '20', '81', '148', '83', '151', '152', '22', '23', '24', '153', '203', '25', '27', '154', '29', '30', '84', '32', '33', '85', '86', '87', '34', '89', '90', '155', '156', '201', '35', '36', '158', '37', '38', '160', '161', '162', '39', '163', '193', '194', '97', '42', '98', '43', '99', '166', '46', '101', '104', '47', '106', '107', '108', '48', '169', '49', '50', '172', '54', '173', '175', '114', '58', '176', '60', '178', '179', '62', '2053', '117', '118', '182', '65', '121', '66', '122', '68', '69', '124', '125', '70', '126', '184', '73', '185', '191', '1011', '1012', '205', '198', '199', '196', '187', '77', '9', '26', '31', '88', '170', '110', '55', '116', '128', '74']

genes = ['reduced_2339_aligned_noempty.fasta-out_clean', 'reduced_514_aligned_noempty.fasta-out_clean' , 'reduced_1007_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_51_aligned_noempty.fasta-out_clean', 'reduced_1013_aligned_noempty.fasta-out_clean', 'reduced_2370_aligned_noempty.fasta-out_clean', 'reduced_52_aligned_noempty.fasta-out_clean', 'reduced_1017_aligned_noempty.fasta-out_clean', 'reduced_2377_aligned_noempty.fasta-out_clean', 'reduced_556_aligned_noempty.fasta-out_clean', 'reduced_1020_aligned_noempty.fasta-out_clean', 'reduced_237_aligned_noempty.fasta-out_clean', 'reduced_563_aligned_noempty.fasta-out_clean', 'reduced_1025_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean', 'reduced_576_aligned_noempty.fasta-out_clean', 'reduced_1035_aligned_noempty.fasta-out_clean', 'reduced_240_aligned_noempty.fasta-out_clean', 'reduced_587_aligned_noempty.fasta-out_clean', 'reduced_1050_aligned_noempty.fasta-out_clean', 'reduced_2459_aligned_noempty.fasta-out_clean', 'reduced_604_aligned_noempty.fasta-out_clean', 'reduced_1052_aligned_noempty.fasta-out_clean', 'reduced_245_aligned_noempty.fasta-out_clean', 'reduced_609_aligned_noempty.fasta-out_clean', 'reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_61_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_1168_aligned_noempty.fasta-out_clean', 'reduced_252e_aligned_noempty.fasta-out_clean', 'reduced_630_aligned_noempty.fasta-out_clean', 'reduced_1171_aligned_noempty.fasta-out_clean', 'reduced_252p_aligned_noempty.fasta-out_clean', 'reduced_637_aligned_noempty.fasta-out_clean', 'reduced_1197_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_673_aligned_noempty.fasta-out_clean', 'reduced_1201_aligned_noempty.fasta-out_clean', 'reduced_2550_aligned_noempty.fasta-out_clean', 'reduced_680_aligned_noempty.fasta-out_clean', 'reduced_120_aligned_noempty.fasta-out_clean', 'reduced_2561_aligned_noempty.fasta-out_clean', 'reduced_717_aligned_noempty.fasta-out_clean', 'reduced_122_aligned_noempty.fasta-out_clean', 'reduced_257_aligned_noempty.fasta-out_clean', 'reduced_727_aligned_noempty.fasta-out_clean', 'reduced_125_aligned_noempty.fasta-out_clean', 'reduced_267_aligned_noempty.fasta-out_clean', 'reduced_732_aligned_noempty.fasta-out_clean', 'reduced_12_aligned_noempty.fasta-out_clean', 'reduced_269_aligned_noempty.fasta-out_clean', 'reduced_736_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_277_aligned_noempty.fasta-out_clean', 'reduced_740_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_280_aligned_noempty.fasta-out_clean', 'reduced_743_aligned_noempty.fasta-out_clean', 'reduced_1484_aligned_noempty.fasta-out_clean', 'reduced_281_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_148_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_758_aligned_noempty.fasta-out_clean', 'reduced_1494_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_785_aligned_noempty.fasta-out_clean', 'reduced_14_aligned_noempty.fasta-out_clean', 'reduced_293_aligned_noempty.fasta-out_clean', 'reduced_790_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean', 'reduced_793_aligned_noempty.fasta-out_clean', 'reduced_1615_aligned_noempty.fasta-out_clean', 'reduced_305_aligned_noempty.fasta-out_clean', 'reduced_7_aligned_noempty.fasta-out_clean', 'reduced_168_aligned_noempty.fasta-out_clean', 'reduced_308_aligned_noempty.fasta-out_clean', 'reduced_807_aligned_noempty.fasta-out_clean', 'reduced_17_aligned_noempty.fasta-out_clean', 'reduced_310_aligned_noempty.fasta-out_clean', 'reduced_808_aligned_noempty.fasta-out_clean', 'reduced_1801_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_1815_aligned_noempty.fasta-out_clean', 'reduced_326_aligned_noempty.fasta-out_clean', 'reduced_825_aligned_noempty.fasta-out_clean', 'reduced_182_aligned_noempty.fasta-out_clean', 'reduced_32e_aligned.fasta-out_clean', 'reduced_82_aligned_noempty.fasta-out_clean', 'reduced_1842_aligned_noempty.fasta-out_clean', 'reduced_32s_aligned_noempty.fasta-out_clean', 'reduced_83_aligned_noempty.fasta-out_clean', 'reduced_1854_aligned_noempty.fasta-out_clean', 'reduced_332_aligned_noempty.fasta-out_clean', 'reduced_84_aligned_noempty.fasta-out_clean', 'reduced_1877_aligned_noempty.fasta-out_clean', 'reduced_357_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_1901_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_191_aligned_noempty.fasta-out_clean', 'reduced_362_aligned_noempty.fasta-out_clean', 'reduced_872_aligned_noempty.fasta-out_clean', 'reduced_194_aligned_noempty.fasta-out_clean', 'reduced_363_aligned_noempty.fasta-out_clean', 'reduced_874_aligned_noempty.fasta-out_clean', 'reduced_197_aligned_noempty.fasta-out_clean', 'reduced_369_aligned_noempty.fasta-out_clean', 'reduced_883e_aligned_noempty.fasta-out_clean', 'reduced_1986_aligned_noempty.fasta-out_clean', 'reduced_378e_aligned_noempty.fasta-out_clean', 'reduced_883n_aligned_noempty.fasta-out_clean', 'reduced_201_aligned_noempty.fasta-out_clean', 'reduced_378s_aligned_noempty.fasta-out_clean', 'reduced_886_aligned_noempty.fasta-out_clean', 'reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_38_aligned_noempty.fasta-out_clean', 'reduced_88_aligned_noempty.fasta-out_clean', 'reduced_204s_aligned_noempty.fasta-out_clean', 'reduced_391_aligned_noempty.fasta-out_clean', 'reduced_897_aligned_noempty.fasta-out_clean', 'reduced_2056_aligned_noempty.fasta-out_clean', 'reduced_392_aligned_noempty.fasta-out_clean', 'reduced_89_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_415_aligned_noempty.fasta-out_clean', 'reduced_938_aligned_noempty.fasta-out_clean', 'reduced_215_aligned_noempty.fasta-out_clean', 'reduced_417_aligned_noempty.fasta-out_clean', 'reduced_948_aligned_noempty.fasta-out_clean', 'reduced_2164_aligned_noempty.fasta-out_clean', 'reduced_421_aligned_noempty.fasta-out_clean', 'reduced_94_aligned_noempty.fasta-out_clean', 'reduced_218_aligned_noempty.fasta-out_clean', 'reduced_449_aligned_noempty.fasta-out_clean', 'reduced_950_aligned_noempty.fasta-out_clean', 'reduced_21_aligned_noempty.fasta-out_clean', 'reduced_464_aligned_noempty.fasta-out_clean', 'reduced_958_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_484_aligned_noempty.fasta-out_clean', 'reduced_964_aligned_noempty.fasta-out_clean', 'reduced_225_aligned_noempty.fasta-out_clean', 'reduced_490_aligned_noempty.fasta-out_clean', 'reduced_977_aligned_noempty.fasta-out_clean', 'reduced_226_aligned_noempty.fasta-out_clean', 'reduced_497_aligned_noempty.fasta-out_clean', 'reduced_982_aligned_noempty.fasta-out_clean', 'reduced_2291_aligned_noempty.fasta-out_clean', 'reduced_4_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_231_aligned_noempty.fasta-out_clean', 'reduced_508_aligned_noempty.fasta-out_clean', 'reduced_989_aligned_noempty.fasta-out_clean']



#Running IQTREE for trimmed and aligned files
for i in range(0, len(genes)):
   gwf.target_from_template('Iqtree_'+str(i), iqtree(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/data/Dryad/B_main_analysis/A_alignments/"))

# running AMAS to convert .fasta files to .nex files
for i in range(0, len(genes)):
   gwf.target_from_template('amas_'+str(i), amas_c(genes = genes[i], path_in="/home/kris/dypsidinae/data/nex_files/"))


#running MrBayes
for i in range(0, len(genes)):
    gwf.target_from_template('MrBayes'+str(i), MrBayes(genes = genes[i], path_in ="/home/kris/dypsidinae/data/nex_files/"))

"""

                        
# Running ASTRAL f
gwf.target_from_template('astral_2', astral_2(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/A/3.Astral/",
                                                    gene_tree_file="gene_trees.nex",
                                                    output="astral_tree_probabilities.tre"))
                                                        
# Running ASTRAL f
gwf.target_from_template('astral_gene_supporting', astral_gene_supporting(genes = genes[i],
                                                    path_in = "/home/kris/dypsidinae/A/3.Astral/",
                                                    gene_tree_file="genes_trees.nex",
                                                    output="astral_gene_probability.tre"))

"""
