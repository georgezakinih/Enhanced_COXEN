## Method Description

The enhanced COXEN method is designed for the applications where the drug efficacy data of a set of cancer cases are used to predict the response of another set of cancer cases. Gene expression data of both sets of cancer cases must be available. The enhanced COXEN method identifies the genes predictive of drug response based on the first cancer set. It first ranks the genes according to their prediction power for each individual drug and then takes a union of top predictive genes of all the drugs. For each predictive gene, the method evaluates how well its co-expression pattern with other predictive genes is preserved between the two sets of cancer cases. Predictive genes that best preserve the co-expression patterns will be identified. The output of the enhanced COXEN package are the selected predictive and generalizable genes, based on which users can build drug response prediction models in subsequent analysis.  

## Setup

To set up the Python environment needed to run this algorithm:
1. Install [conda](https://docs.conda.io/en/latest/) package manager.
2. Clone this repository.
3. Enter the directory of Scripts
4. Create the environment as shown below.
    ```
    conda env create -f environment.yml -n Enhanced_COXEN
    conda activate Enhanced_COXEN
    ```
5.  Run the Example_Run.py script for demo. This demo uses the gene expression values for cases in cancer set 1 [Gene_Expression_Data_Of_Set_1.txt](../Data/Gene_Expression_Data_Of_Set_1.txt), the drugs responses fro cases in cancer set 1 [Drug_Response_Data_Of_Set_1.txt](../Data/Drug_Response_Data_Of_Set_1.txt), and the gene experssion values for cases in cancer set 2 [Gene_Expression_Data_Of_Set_2.txt](../Data/Gene_Expression_Data_Of_Set_2.txt) where we would like to predict the durgs responses. The result contains the list of selected genes that has high predictive power of the drugs responses. 

```
python Example_Run.py 
Selected genes are:
Index(['SEC31B', 'PRSS8', 'DCAF12', 'KDM2B', 'MAML1', 'DLG4', 'ACVR2A', 'HNMT',
       'TLCD1', 'EIF4ENIF1', 'RBBP8NL', 'ZMYND19', 'MTA2', 'MAP2K3', 'USP40',
       'DKK4', 'MUTYH', 'ZNF451', 'GRIPAP1', 'RGPD4', 'ZNF33B', 'COA5',
       'CHRNA3', 'FOXJ1', 'UBN2', 'LCT', 'CSK', 'TMX4', 'STRADB', 'AMPD2',
       'SDF2', 'SOWAHB', 'RXRG', 'NOL3', 'PIKFYVE', 'ZFYVE1', 'ETNK2', 'SLTM',
       'KIF2A', 'LIMA1', 'LUC7L2', 'ZNF273', 'AKR1C3', 'FAP', 'NUP160',
       'CYP4F2', 'BRD7', 'CDH2', 'RILPL2', 'ZNF107', 'FGFR3', 'HIPK1',
       'LGALS1', 'EIF6', 'PRF1', 'RHOB', 'GINS1', 'PROS1', 'C1orf210',
       'RSL24D1', 'SLC7A6OS', 'KLK6', 'FGD3', 'PAWR', 'OSM', 'TAF4', 'MYO6',
       'A2ML1', 'CERS5', 'SLC16A10', 'DPH5', 'IREB2', 'IHH', 'PARVG', 'ZNF250',
       'EVI2A', 'PIK3C2G', 'TNFSF12', 'BLM', 'ATP6V1B2', 'PHKA2', 'SLC35D2',
       'NFYA', 'TREH', 'FHDC1', 'MBTD1', 'CD3G', 'P2RX6', 'FRG1', 'TEAD1',
       'C12orf76', 'SDC1', 'RAD23A', 'CMTM3', 'C1orf131', 'ZNF660', 'LMNB1',
       'CDC42EP3', 'ENTPD2', 'CDKN2C'],
      dtype='object')
```


## Use Enhanced COXEN for Gene Selection

The [EnhancedCOXEN_Functions.py](./EnhancedCOXEN_Functions.py) script provides all the functions used by the enhanced COXEN method. coxen_multi_drug_gene_selection is the main function that performs the enhanced COXEN gene selection. Please see its comments explaining the input and output of the function in details. [Example_Run.py](./Example_Run.py) provides an exmaple demonstrating how to use the function. 
