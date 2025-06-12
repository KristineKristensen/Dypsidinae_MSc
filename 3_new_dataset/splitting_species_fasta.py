# Splitting the gene fasta files into one for each gene and making one file for each with the sequence from each species
import os, argparse

os.chdir("/home/kris/dypsidinae/C_large_dataset/0.new_species/")

# taking the genes argument from the gwf workflow
parser = argparse.ArgumentParser()
parser.add_argument("species_fasta")
args = parser.parse_args()
species_fasta = str(args.species_fasta)

sp= str(species_fasta)[:-14]

# open the file and split into new files
with open(species_fasta, 'r') as sp_file:
     sp_file_lines= sp_file.readlines()

gene_start_index=[]
gene_end_index=[]
for n, lines in enumerate(sp_file_lines):
    if sp+"_gene" in lines:
        gene_start_index.append(n)
        if n> 0:
            gene_end_index.append(n-1)
gene_end_index.append(len(sp_file_lines)-1)

for n,i in enumerate(gene_start_index):
    sp_name= sp_file_lines[i].split("_")[0]
    gene_seq=sp_file_lines[i+1:gene_end_index[n]]
    gene=sp_name+"\n"+"".join(gene_seq)
    gene_name= sp_file_lines[i].split(" ")[0]
    gene_name= gene_name.split("_")[1]
    with open(gene_name+".FNA_output_tapper.fasta", 'a') as f:
        f.write(gene)