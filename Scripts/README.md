## Method Description

The enhanced COXEN method is designed for the applications where the drug efficacy data of a set of cancer cases are used to predict the response of another set of cancer cases. Gene expression data of both sets of cancer cases must be available. The enhanced COXEN method identifies the genes predictive of drug response based on the first cancer set. It first ranks the genes according to their prediction power for each individual drug and then takes a union of top predictive genes of all the drugs. For each predictive gene, the method evaluates how well its co-expression pattern with other predictive genes is preserved between the two sets of cancer cases. Predictive genes that best preserve the co-expression patterns will be identified. The output of the enhanced COXEN package are the selected predictive and generalizable genes, based on which users can build drug response prediction models in subsequent analysis.  

## Setup

To set up the Python environment needed to train and run this model:
1. Install [conda](https://docs.conda.io/en/latest/) package manager.
2. Clone this repository.
3. Enter the directory of Scripts
3. Create the environment as shown below.
    conda env create -f environment.yml -n Enhanced_COXEN
    conda activate Enhanced_COXEN
4.  Run the Example_Run.py script.

## Use Enhanced COXEN for Gene Selection

coxen_multi_drug_gene_selection is the main function that performs the enhanced COXEN gene selection. It is provided in the EnhancedCOXEN_Functions.py script, with comments explaining the input and output of the function in details. Example_Run.py provides an exmaple demonstrating how to use the function. 