#for updating
# get information from best_scheme.nex file generated by IQtree modelfinder and write MrBayes block
# this script is inspired by https://github.com/brettc/partitionfinder/blob/42010508ea382c6450b114c32827ace8727a9d9d/partfinder/model_utils.py#L79 script for converting into MrBayes for partitionfinder
import os, argparse
os.chdir("/home/kris/dypsidinae/2.MrBayes/")

parser = argparse.ArgumentParser()
parser.add_argument("genes")
args = parser.parse_args()
genes = str(args.genes)

charset=open("/home/kris/dypsidinae/1.IQtree_2/"+genes+".part.best_scheme.nex").readlines()

partitions_intron=str(charset[2]).split()
partitions_exon= str(charset[3]).split()
partitions= ""
names=""
bases_intron=[]
bases_exon= []

for i in partitions_intron:
    if "-" in i:
        bases_intron.append(i)
for n in range(len(bases_intron)):
        partitions=partitions +"\t charset intron_"+str(n)+" = "+str(bases_intron[n])+";\n"
        names=names+"intron_"+str(n)+", "
partitions=partitions[:-3]+";\n"

for b in partitions_exon:
    if "-" in b:
        bases_exon.append(b)
for n in range(len(bases_exon)):
        partitions=partitions +"\t charset exon_"+str(n)+" = "+str(bases_exon[n])+";\n"
        names=names+"exon_"+str(n)+", "
partitions=partitions[:-3]+";\n"


begin= "begin mrbayes;\n\n"
autoclose= "\t set autoclose=yes; [makes automatic execution]\n"
make_log= "\t log start filename="+genes+"_2.log replace;\n"
file_to_execute= "\t Execute /home/kris/dypsidinae/data/nex_files/"+ genes+".fasta-out.nex;\n"
defining_partition= "\t partition types ="+str(len(bases_exon)+len(bases_intron))+":"+names[:-2]+";\n \t set partition = types;\n" # make recursive
prset_all="prset applyto=(all) ratepr=variable;\n"
unlink= " unlink revmat=(all) pinvar=(all) shape=(all) statefreq=(all); \n"

block= [begin,autoclose, str(make_log), str(file_to_execute), partitions, str(defining_partition),prset_all, unlink]

#taking model from IQtree and defining the model for the MrBayes block
models=[charset[5],charset[6]]

list_for_iteration= names[:-2].split(", ")
for nr,name in enumerate(list_for_iteration, start=1):
    if "intron" in name:
        i =models[0]
    elif "exon" in name:
         i=models[1]

    if "GTR" in i or "SYM" in i:
        nucleotide_model ="\t lset applyto=("+str(nr)+") nst=6 "
    elif "HKY" in i or "HKY85" in i or "K80" in i or "K2P" in i:
        nucleotide_model= "\t lset applyto=("+str(nr)+") nst=2 "
    elif "F81" in i or "JC69" in i or "JC" in i:
         nucleotide_model= "\t lset applyto=("+str(nr)+") nst=1 "
    else:
         nucleotide_model= "not recognized" 


    if "SYM" in i or "JC" in i or "JC69" in i or "K80" in i or "K2P" in i or "+FQ" in i:
          base_freq= "\t prset applyto=("+str(nr)+") statefreqpr= fixed(equal); \n"
    else:
         base_freq=""

    #rate varation
    #If IQtree gives I
    if "+I" in i and "+G" not in i:
        rate_variation= " rates=propinv ; \n"
    #If IQtree gives G
    elif "+G" in i and "+I" not in i:
        start = "+G"
        end = ":"
        index_s = i.partition(start)
        ngammacat=str(index_s[2]).partition(end)
        rate_variation= " rates=gamma ngammacat="+ngammacat[0]+"; \n"
    #If IQtree gives G+I
    elif "+I" in i and "+G" in i:
        start = "+G"
        end = ":"
        index_s = i.partition(start)
        ngammacat=str(index_s[2]).partition(end)
        rate_variation= " rates=invgamma ngammacat="+ngammacat[0]+"; \n"
    else:
        rate_variation= "; \n"

    block.append(str(nucleotide_model)+str(rate_variation))
    block.append(str(base_freq))


block.append("\t mcmc ngen=3500000 samplefreq=2500 printfreq=2500 diagnfreq=2500 append=yes;\n ")
block.append("\t sump;\n")
block.append("END;")

#Write MrBayes block
MrBayes_block = open(genes+"_MrBayes_block_2.nex", 'w')

for i in block:
      MrBayes_block.write(i)




