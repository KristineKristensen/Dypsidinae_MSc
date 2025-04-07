import os, argparse

os.chdir("/home/kris/dypsidinae/C_large_dataset/0.old_alignments/")

#Species that need to be removed from original MSA
sp_to_remove=["Arec1688", "Arec1539", "Arec1695", "Arec704"]

#Opens the MSA for each gene and finds the line index where the sequence for each species start
# taking the genes argument from the gwf workflow
parser = argparse.ArgumentParser()
parser.add_argument("gene")
args = parser.parse_args()
gene_fasta = str(args.gene)

with open("/home/kris/dypsidinae/data/1.Dypsidiinae/"+gene_fasta, "r") as f:
    file_w_sp=f.readlines()
    gene_indx=[]
    for n,l in enumerate(file_w_sp):
        if ">" in l:
            gene_indx.append(n)

    #Adds the last line in file to list of indexes
    gene_indx.append(len(file_w_sp))

    #Making new .fasta files 
    rm_s=[]
    for i in range(len(gene_indx)-1):
            if not any(s in file_w_sp[gene_indx[i]] for s in sp_to_remove): # checks if any of the unwanted species names is in that line and if it is not in the line, the line and the following sequence is written to a new file       
                with open(gene_fasta, "a") as new_f:
                    new_f.write("".join(file_w_sp[gene_indx[i]:(gene_indx[i+1])]))
            else:
                rm_s.append("".join(file_w_sp[gene_indx[i]:(gene_indx[i+1])]))
               





            
    