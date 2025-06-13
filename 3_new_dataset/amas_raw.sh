#!/bin/bash

#Go to folder
cd /home/kris/dypsidinae/C_large_dataset/2.trimAL/

#Activate enviroment
source /home/kris/miniconda3/etc/profile.d/conda.sh
conda activate amas
 

#Calculate AMAS statistics for raw trimal value.
for i in *aligned.fasta; do AMAS.py summary -c 16 -f fasta -d dna -i ${i} -o ${i}_summary.txt; done

#Delete the first line in the summary files
for i in *summary.txt; do sed '1d' ${i} > ${i}_summary_for_real.txt; done

#Save all summary for real in a statistics for each gt
cat *_summary_for_real.txt > summary_0.txt


#Put all information in a new line
sed -e 's/gene/\ngene/g' summary_0.txt > 2.txt

#Remove summary0
rm summary_0.txt

#Rename
mv 2.txt summary_0.txt

#Remove all previous summary
rm *_summary_for_real.txt

