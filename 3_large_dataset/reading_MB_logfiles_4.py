import os, csv
os.chdir("/home/kris/dypsidinae/C_large_dataset/8.MrBayes/")

# The genes for which convergence is tested
genes = ['gene3388', 'gene4435', 'gene5001', 'gene3087', 'gene1483', 'gene3415', 'gene4668', 'gene4699', 'gene1729', 'gene4898', 'gene2346', 'gene667', 'gene1051', 'gene2318', 'gene2319', 'gene1305', 'gene883', 'gene5367', 'gene2340', 'gene2904', 'gene4264', 'gene1894', 'gene3039', 'gene2422', 'gene2972', 'gene1539', 'gene2954', 'gene1200', 'gene888', 'gene2180', 'gene3700', 'gene5399', 'gene1527', 'gene1945', 'gene4588', 'gene1470', 'gene347', 'gene2285', 'gene331', 'gene2814', 'gene2189', 'gene4235', 'gene1544', 'gene1297', 'gene4444', 'gene4867', 'gene3174', 'gene310', 'gene503', 'gene3944', 'gene106', 'gene3436', 'gene4336', 'gene623', 'gene1652', 'gene4443', 'gene694', 'gene1351', 'gene910', 'gene884', 'gene1523', 'gene1585', 'gene1545', 'gene4007', 'gene4997', 'gene3482', 'gene608', 'gene2391', 'gene4618', 'gene1121', 'gene4619', 'gene1471', 'gene3295', 'gene5165', 'gene330', 'gene5368', 'gene5164', 'gene3732', 'gene1463', 'gene1447', 'gene2436', 'gene807', 'gene2551', 'gene4728', 'gene5232', 'gene2622', 'gene4254', 'gene3201', 'gene2571', 'gene5463', 'gene1346', 'gene4648', 'gene2348', 'gene1338', 'gene4724', 'gene1203', 'gene4470', 'gene4839', 'gene4899', 'gene4686', 'gene4685', 'gene1424', 'gene4253', 'gene2838', 'gene2621', 'gene3438', 'gene5370', 'gene570', 'gene5234', 'gene3858', 'gene504', 'gene4904', 'gene1827', 'gene928', 'gene1040', 'gene2061', 'gene433', 'gene3302', 'gene1457', 'gene512', 'gene3463', 'gene2539', 'gene4804', 'gene5033', 'gene3344', 'gene4445', 'gene3414', 'gene1151', 'gene333', 'gene4624', 'gene2554', 'gene348', 'gene1518', 'gene649', 'gene5259', 'gene3464', 'gene3772', 'gene1142', 'gene1730', 'gene3316', 'gene2950', 'gene4796', 'gene3870', 'gene2790', 'gene5013', 'gene2188', 'gene5257', 'gene3111', 'gene1826', 'gene4198', 'gene3319', 'gene5392', 'gene1307', 'gene2543', 'gene3532', 'gene332', 'gene1404', 'gene3282', 'gene5010', 'gene2392', 'gene988', 'gene2704', 'gene4775', 'gene855', 'gene3097', 'gene5233', 'gene1350', 'gene896', 'gene3125', 'gene2408', 'gene3403', 'gene4183', 'gene1461', 'gene961', 'gene4041', 'gene4614', 'gene2977', 'gene2795', 'gene3896', 'gene4197', 'gene4805', 'gene1446', 'gene5322', 'gene2404', 'gene5011', 'gene607', 'gene3806', 'gene1371', 'gene729', 'gene5305', 'gene2555', 'gene4795', 'gene730', 'gene1733', 'gene69', 'gene904', 'gene4760', 'gene143', 'gene4761', 'gene4213', 'gene2345', 'gene469', 'gene1722', 'gene3031', 'gene86', 'gene911', 'gene88', 'gene2263', 'gene5320', 'gene192', 'gene1279', 'gene2267', 'gene4932', 'gene1797', 'gene1725', 'gene4228', 'gene2788', 'gene3861', 'gene5163', 'gene724', 'gene3015', 'gene3872', 'gene5015', 'gene3926', 'gene2142', 'gene932', 'gene5325', 'gene3811', 'gene4687', 'gene5366', 'gene4706', 'gene972', 'gene939', 'gene1108', 'gene2652', 'gene5361', 'gene572', 'gene1352', 'gene439', 'gene1345', 'gene3301', 'gene2976', 'gene1560', 'gene1090', 'gene4931', 'gene2875', 'gene4181', 'gene1658', 'gene1038', 'gene1342', 'gene1462', 'gene1341', 'gene2164', 'gene4635', 'gene1172', 'gene1776', 'gene2062', 'gene1665', 'gene3100', 'gene2393', 'gene566', 'gene3062', 'gene4419', 'gene1027', 'gene1430', 'gene2295', 'gene4488', 'gene2264', 'gene296', 'gene4214', 'gene1429', 'gene2605', 'gene1417', 'gene2703', 'gene1576', 'gene4946', 'gene4770', 'gene4154', 'gene1439', 'gene2861', 'gene2663', 'gene3281', 'gene3818', 'gene4623', 'gene5098', 'gene152', 'gene4474', 'gene1089', 'gene2825', 'gene3895', 'gene3439', 'gene3779', 'gene3942', 'gene4702']

data = [['gene', 'ASDSF','plot','ESS_min','PSRF_min', 'PSRF_max']]
for gene in genes:
    #Opening the log file from MrBayes
    gene_file=gene+"_output_tapper.fasta_4.log"
    if os.path.isfile(gene_file):
        list_with_data=[]
        list_with_data.append(gene)
        with open(gene_file, 'r') as f: 
            #Reading lines in the log file to find the line where the information is and the index of the line
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
            # Taking ASDF value
            list_with_data.append(gene_log[l_index][-10:-2])
            #Check the plot in the log file that show if the two runs are similar
            plots=gene_log[index_n+1:index_n+16]
            if  any("1" in l and "2" in l for l in plots):
                list_with_data.append("good")
            else:
                list_with_data.append("problem")
                
            #check ESS
            table_values= gene_log[table_index+2:stop-1]
            ESS_values=[]
            for v in table_values:
                ESS_values.append(v[ESS_index-1:ESS_index+6])
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
        #Making a file that lists the genes for which no log file was made by MrBayes
        with open('mb_not_made_4.txt', 'a') as f:
            f.write(gene+"\n")

#Making a CSV file with the convergence information
with open('mb_records_4.csv', 'a', newline='') as csvfile:
    new_file= csv.writer(csvfile, delimiter =";")
    new_file.writerows(data)


