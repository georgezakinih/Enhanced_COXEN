import numpy as np
import pandas as pd
from EnhancedCOXEN_Functions import coxen_multi_drug_gene_selection



res = pd.read_csv('../Data/Drug_Response_Data_Of_Set_1.txt', sep='\t', engine='c', na_values=['na', '-', ''],
                  header=0, index_col=None)

data1 = pd.read_csv('../Data/Gene_Expression_Data_Of_Set_1.txt', sep='\t', engine='c', na_values=['na', '-', ''],
                    header=0, index_col=0)

data2 = pd.read_csv('../Data/Gene_Expression_Data_Of_Set_2.txt', sep='\t', engine='c', na_values=['na', '-', ''],
                    header=0, index_col=0)

assert np.sum(data1.columns != data2.columns) == 0

# Use enhanced COXEN method to select genes. First, select 200 predictive genes to form the candidate pool, and then
# select 100 generalizable genes from the candidate pool. The absolute value of Pearson correlation coefficient is
# used as the measure of gene's prediction power.

id = coxen_multi_drug_gene_selection(source_data=data1, target_data=data2, drug_response_data=res,
                                     drug_response_col='Efficacy', tumor_col='Cell_Line', drug_col='Drug',
                                     prediction_power_measure='pearson', num_predictive_gene=200,
                                     num_generalizable_gene=100)

print('Selected genes are:')
print(data1.columns[id])