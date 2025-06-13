import os,csv, argparse

#taking the genes argument from the gwf workflow
parser = argparse.ArgumentParser()
parser.add_argument("genes")
args = parser.parse_args()
genes = args.genes.split(',')

data = [['gene', 'mrb_model','freerate','model_test']]


#Finding the models selected
for g in genes:
    list_with_data=[]
    list_with_data.append(g)
    if os.path.exists("/home/kris/dypsidinae/C_large_dataset/10.IQtree_mrb_model/"+g+"_output_tapper.fasta.iqtree"):
        with open("/home/kris/dypsidinae/C_large_dataset/10.IQtree_mrb_model/"+g+"_output_tapper.fasta.iqtree") as f1:
            IQtree_result_mrb_model=f1.readlines()
            for lines in IQtree_result_mrb_model:
                if "Best-fit model according to BIC: " in lines:
                    list_with_data.append(lines)
    if os.path.exists("/home/kris/dypsidinae/C_large_dataset/10.IQtree/"+g+"_output_tapper.fasta.iqtree"):
        with open("/home/kris/dypsidinae/C_large_dataset/10.IQtree/"+g+"_output_tapper.fasta.iqtree") as f2:
            IQtree_result_freerate= f2.readlines()
            for lines2 in IQtree_result_freerate:
                if "Best-fit model according to BIC: " in lines2:
                    list_with_data.append(lines2)
    if os.path.exists("/home/kris/dypsidinae/C_large_dataset/6.IQtree/"+g+"_output_tapper.fasta.iqtree"):
        with open("/home/kris/dypsidinae/C_large_dataset/6.IQtree/"+g+"_output_tapper.fasta.iqtree") as f3:
            IQtree_result_model_test =  f3.readlines()
            for lines3 in IQtree_result_model_test:
                if "Best-fit model according to BIC: " in lines3:
                    list_with_data.append(lines3)
    data.append(list_with_data)

#Making a CSV file with the model information
with open('models_chosen.csv', 'a', newline='') as csvfile:
    new_file= csv.writer(csvfile, delimiter =";")
    new_file.writerows(data)

