# Enhanced Co-Expression Extrapolation (COXEN) Gene Selection Method for Building Drug Response Prediction Models

## Description

The enhanced co-expression extrapolation COXEN method enhances the original COXEN method to select genes that are predictive of the efficacies of multiple drugs, for the purpose of building general drug response prediction models that are not specific to a particular drug. It is designed for the applications where the drug efficacy data of a set of cancer cases are used to predict the response of another set of cancer cases. A publication describing the enhaced COXEN method is available at https://www.mdpi.com/2073-4425/11/9/1070

## User Community

Primary: Cancer biology data modeling
Secondary: Machine Learning; Bioinformatics; Computational Biology

## Usability

To use this software package for enhanced COXEN analyses, users must possess the basic skills to run Python scripts. Users need to process the gene expression data and drug response data into the data format accepted by the enhanced COXEN package. 

## Uniqueness

The original COXEN method has been successfully used in multiple studies to select genes for predicting the response of tumor cells to a specific drug treatment. The enhanced COXEN method selects genes that are predictive of the efficacies of multiple drugs for building general drug response prediction models that are not specific to a particular drug. It first ranks the genes according to their prediction power for each individual drug and then takes a union of top predictive genes of all the drugs, among which the algorithm further selects genes whose co-expression patterns are well preserved between cancer cases for building prediction models. 

## Components

The package includes two Python scripts. 
1. EnhancedCOXEN_Functions.py provides all the functions used by the enhanced COXEN method
2. Example_Run.py provides example code demonstrating how to use the functions for enhanced COXEN analysis.

The package includes a small dataset for demonstrating its utility, which is composed of three data files including:
1. Gene_Expression_Data_Of_Set_1.txt: the gene expression data of cancer case set 1.
2. Drug_Response_Data_Of_Set_1.txt: the drug response data of cancer cases in set 1.
3. Gene_Expression_Data_Of_Set_2.txt: the gene expression data of cancer case set 2, for which drug response needs to be predicted.

## Technical Details

Refer to this [README](https://github.com/zhuyitan/Enhanced_COXEN/Scripts/README.md).