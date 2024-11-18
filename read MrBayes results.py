import os, re, csv

os.chdir("/home/kris/dypsidinae/B.phylopalm/3.Mrbayes/")

parser = argparse.ArgumentParser()
parser.add_argument("genes")
args = parser.parse_args()
genes = str(args.genes)


data = [['gene', 'ASDSF','plot','ESS_min','PSRF_min', 'PSRF_max']]
if os.path.isfile(genes):
    list_with_data=[]
    list_with_data.append(genes)
    with open(genes, 'r') as f:  
        gene_log= f.readlines()
        for n, lines in enumerate(gene_log):
            if 'Average standard deviation of split frequencies:' in lines:
                l_index= n
            if "+------------------------------------------------------------+" in lines:
                index_n=n
            if "Parameter      Mean      Variance     Lower       Upper       Median    min ESS*  avg ESS    PSRF+ " in lines:
                table_index=n
                ESS_index= lines.index("min ESS*")
                PSRF_index= lines.index("PSRF+")
            if "* Convergence diagnostic (ESS = Estimated Sample Size); min and avg values" in lines:
                stop= n
         
        list_with_data.append(gene_log[l_index][-10:-2])

        plots=gene_log[index_n+1:index_n+16]
        if  any("1" in l and "2" in l for l in plots):
            list_with_data.append("good")
        else:
            list_with_data.append("problem")
            
        #check ESS
        table_values= gene_log[table_index+2:stop-1]
        ESS_values=[]
        for v in table_values:
            ESS_values.append(v[ESS_index:ESS_index+6])
        ESS_min=min(ESS_values)
        list_with_data.append(ESS_min)
            
        #Check PSRF
        table_values= gene_log[table_index+2:stop-1]
        PSRF_values=[]
        for v in table_values:
            PSRF_values.append(v[PSRF_index-1:PSRF_index+4])
        PSRF_min=min(PSRF_values)
        PSRF_max=max(PSRF_values)
        list_with_data.append(PSRF_min)
        list_with_data.append(PSRF_max)
        
        data.append(list_with_data)    
else:
    with open('mb_not_made.txt', 'a') as f:
        f.write(genes+"\n")

with open('mb_records.csv', 'a', newline='') as csvfile:
    new_file= csv.writer(csvfile, delimiter =";")
    new_file.writerows(data)


