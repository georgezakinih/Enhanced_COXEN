# Enhanced Co-Expression Extrapolation (COXEN) Gene Selection Method for Building Drug Response Prediction Models

The co-expression extrapolation (COXEN) method has been successfully used in multiple studies to select genes for predicting the response of tumor cells to a specific drug treatment. We enhance the COXEN method to select genes that are predictive of the efficacies of multiple drugs for building general drug response prediction models that are not specific to a particular drug. The enhanced COXEN method first ranks the genes according to their prediction power for each individual drug and then takes a union of top predictive genes of all the drugs, among which the algorithm further selects genes whose co-expression patterns are well preserved between cancer cases for building prediction models. The paper is available at https://www.mdpi.com/2073-4425/11/9/1070

Scripts folder includes two scripts, EnhancedCOXEN_Functions.py and Example_Run.py. EnhancedCOXEN_Functions.py provides all the functions used by the enhanced COXEN algorithm. Comments explaining the input and output of each function are provided in the script. Example_Run.py gives examples showing how to run the enhanced COXEN algorithm for demo purpose.

Data folder includes a small dataset used by the Example_Run.py script for running the demos.
