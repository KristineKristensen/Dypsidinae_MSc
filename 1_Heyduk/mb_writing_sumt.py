# Writing a MrBayes blok with the information to run sumt on the trees in the posteriror sample
import os, argparse
os.chdir("/home/kris/dypsidinae/2.MrBayes/")

#Lists of genes that converged after the different number of runs
#genes converged after 2.5 milion
genes_v2 = ['reduced_191_aligned_noempty.fasta-out_clean', 'reduced_252e_aligned_noempty.fasta-out_clean', 'reduced_21_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_362_aligned_noempty.fasta-out_clean', 'reduced_883e_aligned_noempty.fasta-out_clean', 'reduced_226_aligned_noempty.fasta-out_clean', 'reduced_563_aligned_noempty.fasta-out_clean', 'reduced_89_aligned_noempty.fasta-out_clean']
#genes that converged after 3.5 milion
genes_v3 =['reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean','reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_1020_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_12_aligned_noempty.fasta-out_clean', 'reduced_1007_aligned_noempty.fasta-out_clean', 'reduced_83_aligned_noempty.fasta-out_clean', 'reduced_1168_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_2370_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean','reduced_191_aligned_noempty.fasta-out_clean', 'reduced_391_aligned_noempty.fasta-out_clean', 'reduced_673_aligned_noempty.fasta-out_clean', 'reduced_369_aligned_noempty.fasta-out_clean', 'reduced_1171_aligned_noempty.fasta-out_clean', 'reduced_207_aligned_noempty.fasta-out_clean', 'reduced_1201_aligned_noempty.fasta-out_clean', 'reduced_897_aligned_noempty.fasta-out_clean', 'reduced_758_aligned_noempty.fasta-out_clean'] 
#genes that converged after 5 milion
genes_v4= ['reduced_204e_aligned_noempty.fasta-out_clean', 'reduced_191_aligned_noempty.fasta-out_clean', 'reduced_1064_aligned_noempty.fasta-out_clean', 'reduced_629_aligned_noempty.fasta-out_clean', 'reduced_282_aligned_noempty.fasta-out_clean', 'reduced_252s_aligned_noempty.fasta-out_clean', 'reduced_757_aligned_noempty.fasta-out_clean', 'reduced_863_aligned_noempty.fasta-out_clean', 'reduced_2238_aligned_noempty.fasta-out_clean', 'reduced_299_aligned_noempty.fasta-out_clean', 'reduced_985_aligned_noempty.fasta-out_clean', 'reduced_2363_aligned_noempty.fasta-out_clean', 'reduced_250_aligned_noempty.fasta-out_clean', 'reduced_360_aligned_noempty.fasta-out_clean', 'reduced_323_aligned_noempty.fasta-out_clean', 'reduced_24_aligned_noempty.fasta-out_clean', 'reduced_822_aligned_noempty.fasta-out_clean', 'reduced_290_aligned_noempty.fasta-out_clean', 'reduced_855_aligned_noempty.fasta-out_clean', 'reduced_136_aligned_noempty.fasta-out_clean', 'reduced_139_aligned_noempty.fasta-out_clean', 'reduced_110_aligned_noempty.fasta-out_clean', 'reduced_150_aligned_noempty.fasta-out_clean', 'reduced_2388_aligned_noempty.fasta-out_clean']
                                                    

# taking the converged genes argument from the gwf workflow
parser = argparse.ArgumentParser()
parser.add_argument("genes")
args = parser.parse_args()
genes = str(args.genes)

# Make the MrBayes block for the genes according to how many generations it ran
if genes not in genes_v2:
        with open(genes+'_mb_sumt.nex', 'a') as file:
            mb_command=["begin mrbayes;\n\n", "\t Execute /home/kris/dypsidinae/data/nex_files/"+genes+"_no_X-out.nex;\n", "\t sumt conformat=simple;\n","END;"]
            for i in mb_command:
                file.write(i)
elif genes in genes_v2 and genes not in genes_v3:
        with open(genes+'_mb_sumt.nex', 'a') as file:
            mb_command=["begin mrbayes;\n\n", "\t Execute /home/kris/dypsidinae/data/nex_files/"+genes+"_no_X-out.nex;\n", "\t sumt relburnin=no burnin=650 conformat=simple;\n","END;"]
            for i in mb_command:
                file.write(i)
elif genes in genes_v3 and genes not in genes_v4:
       with open(genes+'_mb_sumt.nex', 'a') as file:
            mb_command=["begin mrbayes;\n\n", "\t Execute /home/kris/dypsidinae/data/nex_files/"+genes+"_no_X-out.nex;\n", "\t sumt relburnin=no burnin=1250 conformat=simple;\n","END;"]
            for i in mb_command:
                file.write(i) 
elif genes in genes_v4:
       with open(genes+'_mb_sumt.nex', 'a') as file:
            mb_command=["begin mrbayes;\n\n", "\t Execute /home/kris/dypsidinae/data/nex_files/"+genes+"_no_X-out.nex;\n", "\t sumt relburnin=no burnin=2250 conformat=simple;\n","END;"]
            for i in mb_command:
                file.write(i) 



