import os, re

os.chdir("/home/kris/dypsidinae/B_phylopalm/4.Astral")

for n in range(1000):
    with open(str(n)+"_astral_tree.tre") as f:#opening one of the 1000 astral tree files
        tree_w_support=f.read()
    #removing the support values
    bits= re.split("(\))",tree_w_support)
    tree_no_support= ""
    for b in bits:
        if len(b)>1 and ":" != b[0] and bits[0] != b and bits[-1] !=b:
            tree_no_support= tree_no_support+b[8:]
        else:
            tree_no_support=tree_no_support+b

 #make a new file with all the 1000 trees without support values
    with open("astraltrees_for_mrbayes.t", "a") as new_tree:
        new_tree.write(tree_no_support)



