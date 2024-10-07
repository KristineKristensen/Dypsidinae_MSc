#This script takes the .t files made by MrBayes and then make a new filewith the genetrees  without burnin
import os, argparse
os.chdir("/home/kris/dypsidinae/data/nex_files/")

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("genes")
args = parser.parse_args()
genes=str(args.genes)

genes_3_5=['reduced_280_aligned_noempty.fasta-out_clean', 'reduced_1017_aligned_noempty.fasta-out_clean', 'reduced_252e_aligned_noempty.fasta-out_clean', 'reduced_231_aligned_noempty.fasta-out_clean', 'reduced_563_aligned_noempty.fasta-out_clean', 'reduced_1815_aligned_noempty.fasta-out_clean', 'reduced_1877_aligned_noempty.fasta-out_clean', 'reduced_958_aligned_noempty.fasta-out_clean', 'reduced_89_aligned_noempty.fasta-out_clean', 'reduced_226_aligned_noempty.fasta-out_clean', 'reduced_883e_aligned_noempty.fasta-out_clean', 'reduced_362_aligned_noempty.fasta-out_clean', 'reduced_484_aligned_noempty.fasta-out_clean', 'reduced_1484_aligned_noempty.fasta-out_clean', 'reduced_21_aligned_noempty.fasta-out_clean']
genes_5=['reduced_1020_aligned_noempty.fasta-out_clean', 'reduced_12_aligned_noempty.fasta-out_clean', 'reduced_1007_aligned_noempty.fasta-out_clean', 'reduced_83_aligned_noempty.fasta-out_clean', 'reduced_1168_aligned_noempty.fasta-out_clean', 'reduced_2370_aligned_noempty.fasta-out_clean', 'reduced_391_aligned_noempty.fasta-out_clean', 'reduced_673_aligned_noempty.fasta-out_clean', 'reduced_369_aligned_noempty.fasta-out_clean', 'reduced_1171_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_1201_aligned_noempty.fasta-out_clean', 'reduced_897_aligned_noempty.fasta-out_clean', 'reduced_758_aligned_noempty.fasta-out_clean']
genes_7_5= ['reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean']
not_converged=['reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_191_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean']
list_not_in=[]

if genes in genes_3_5:
    genetrees = open(genes+".fasta-out.nex.run1.t").read()
    to=genetrees[:-5].split("tree gen.1625000")
    to_print="tree gen.1625000"+to[1]
    tree_file= open(genes+"_1_tree_file.nex", 'w')
    tree_file.write(to_print)
    genetrees_2 = open(genes+".fasta-out.nex.run2.t").read()
    to_2=genetrees_2[:-5].split("tree gen.1625000")
    to_print_2="tree gen.1625000"+to_2[1]
    tree_file_2= open(genes+"_2_tree_file.nex", 'w')
    tree_file_2.write(to_print_2)
    nr_3_5=open("nr_3_5", "a")
    nr_3_5.write(genes+"\n")
elif genes in genes_5:
    genetrees = open(genes+".fasta-out.nex.run1.t").read()
    to=genetrees[:-5].split("tree gen.3125000")
    to_print="tree gen.3125000"+to[1]
    tree_file= open(genes+"_1_tree_file.nex", 'w')
    tree_file.write(to_print)
    genetrees_2 = open(genes+".fasta-out.nex.run2.t").read()
    to_2=genetrees_2[:-5].split("tree gen.3125000")
    to_print_2="tree gen.3125000"+to_2[1]
    tree_file_2=open(genes+"_2_tree_file.nex", 'w')
    tree_file_2.write(to_print_2)
    nr_5=open("nr_5", "a")
    nr_5.write(genes+"\n")
elif genes in genes_7_5:
    genetrees = open(genes+".fasta-out.nex.run1.t").read()
    to=genetrees[:-5].split("tree gen.5625000")
    to_print="tree gen.5625000"+to[1]
    tree_file= open(genes+"_1_tree_file.nex", 'w')
    tree_file.write(to_print)
    genetrees_2 = open(genes+".fasta-out.nex.run2.t").read()
    to_2=genetrees_2[:-5].split("tree gen.5625000")
    to_print_2="tree gen.5625000"+to_2[1]
    tree_file_2=open(genes+"_2_tree_file.nex", 'w')
    tree_file_2.write(to_print_2)
    nr_7_5=open("nr_7_5", "a")
    nr_7_5.write(genes+"\n")
elif genes in not_converged:
    nr_nc=open("nr_nc", "a")
    nr_nc.write(genes+"\n")
else:
    genetrees = open(genes+".fasta-out.nex.run1.t").read()
    to=genetrees[:-5].split("tree gen.625000")
    to_print="tree gen.625000"+to[1]
    tree_file= open(genes+"_1_tree_file.nex", 'w')
    tree_file.write(to_print)
    genetrees_2 = open(genes+".fasta-out.nex.run2.t").read()
    to_2=genetrees_2[:-5].split("tree gen.625000")
    to_print_2="tree gen.625000"+to_2[1]
    tree_file_2= open(genes+"_2_tree_file.nex", 'w')
    tree_file_2.write(to_print_2)
    nr_2_5=open("nr_2_5", "a")
    nr_2_5.write(genes+"\n")
