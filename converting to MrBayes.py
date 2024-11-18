import os, re
from Bio import Phylo

os.chdir("/home/kris/dypsidinae/3.Astral/")

for n in range(1000):
    with open(str(n)+"_astral_tree_probabilities.tre") as f: #opening the 1000 astral trees
        tree_w_support=f.read()
    bits= re.split("(\))",tree_w_support) #removing the support values
    tree_no_support= ""
    for b in bits:
        if len(b)>1 and ":" != b[0] and bits[0] != b and bits[-1] != b:
            tree_no_support= tree_no_support+b[8:]
        else: 
            tree_no_support=tree_no_support+b
            
    with open("astral_for_mrbayes.tre", "a") as new_tree: # make a new file with all the 1000 trees without support values
        new_tree.write(tree_no_support)


        

#Phylo.convert("astral_for_mrbayes.tre", "newick","astral_for_mrbayes.run1.t","nexus") resulted in a format that MrBayes could not read so ape package in R was used
#After using write(nexus) the following code was used

with open("mrbayes.run1.t") as f:
    file_to_read= f.readlines()

with open('a_2_mb.t', 'a') as file:
    number=0
    for lines in file_to_read:
        if 'TREE * UNTITLED = [&R]' in lines:
            number+=1
            newline=lines.replace('TREE * UNTITLED = [&R]', 'tree gen.'+str(number)+'= [&U]')
            file.write(newline)
        else:
            file.write(lines)


#lavede dummy .nex file as that is required by MrBayes
#Then ran MrBayes