# This script converts the nexus files from MrBayes into Newick files.
# It also saves the tree generated in the last 500 generations into a new file thereby removing the burnin
from Bio import Phylo
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("genes")
args = parser.parse_args()
genes=str(args.genes)

os.chdir("/home/kris/dypsidinae/B_phylopalm/2.Nex_files/")
#converting the files
Phylo.convert(genes+"_no_X-out.nex.run1.t", "nexus", genes+"-out.nex.run1_newick.t", "newick")
Phylo.convert(genes+"_no_X-out.nex.run2.t", "nexus", genes+"-out.nex.run2_newick.t", "newick")

#excluding burnin and creating files
newick_1=open(genes+"-out.nex.run1_newick.t","r").readlines()
run1=newick_1[-500:] #taking the last 500 tree generated by run 1 in MrBayes
newick_2=open(genes+"-out.nex.run2_newick.t","r").readlines()
run2=newick_2[-500:] #taking the last 500 tree generated by run 1 in MrBayes

genetree=open(genes+"_genetree.tre", "a")
for i in run1:
    genetree.write(i)
for t in run2:
    genetree.write(t)
