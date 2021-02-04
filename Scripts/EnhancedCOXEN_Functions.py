import sys
import pandas as pd
import numpy as np
import numpy.linalg as la
import patsy
from scipy import stats
from astropy.stats import median_absolute_deviation
from collections import Counter
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import mutual_info_regression
import matplotlib.pyplot as plt



def select_features_by_variation(data, variation_measure='var', threshold=None, equal_flag=True,
                                 portion=None, draw_histogram=False, bins=100, log=False):
    '''
    This function evaluates the variations of individual features and returns the indices of features with large
    variations. Missing values are ignored in evaluating variation.

    Parameters:
    -----------
    data: numpy array or pandas data frame of numeric values, with a shape of [n_samples, n_features].
    variation_metric: string indicating the metric used for evaluating feature variation. 'var' indicates variance;
        'std' indicates standard deviation; 'mad' indicates median absolute deviation. Default is 'var'.
    threshold: float. Features with a variation larger than threshold will be selected. Default is None.
    equal_flag: boolean, indicating whether features with a variation equal to threshold should be selected.
        Default is True.
    portion: float in the range of [0, 1]. It is the portion of features to be selected based on variation.
        The number of selected features will be the smaller of int(portion * n_features) and the total number of
        features with non-missing variations. Default is None. threshold and portion can not take real values
        and be used simultaneously.
    draw_histogram: boolean, whether to draw a histogram of feature variations. Default is False.
    bins: positive integer, the number of bins in the histogram. Default is the smaller of 50 and the number of
        features with non-missing variations.
    log: boolean, indicating whether the histogram should be drawn on log scale.


    Returns:
    --------
    indices: 1-D numpy array containing the indices of selected features. If both threshold and
        portion are None, indices will be an empty array.
    '''

    if isinstance(data, pd.DataFrame):
        data = data.values
    elif not isinstance(data, np.ndarray):
        print('Input data must be a numpy array or pandas data frame')
        sys.exit(1)

    if variation_measure == 'std':
        v_all = np.nanstd(a=data, axis=0)
    elif variation_measure == 'mad':
        v_all = median_absolute_deviation(data=data, axis=0, ignore_nan=True)
    else:
        v_all = np.nanvar(a=data, axis=0)

    indices = np.where(np.invert(np.isnan(v_all)))[0]
    v = v_all[indices]

    if draw_histogram:
        if len(v) < 50:
            print('There must be at least 50 features with variation measures to draw a histogram')
        else:
            bins = int(min(bins, len(v)))
            _ = plt.hist(v, bins=bins, log=log)
            plt.show()

    if threshold is None and portion is None:
        return np.array([])
    elif threshold is not None and portion is not None:
        print('threshold and portion can not be used simultaneously. Only one of them can take a real value')
        sys.exit(1)

    if threshold is not None:
        if equal_flag:
            indices = indices[np.where(v >= threshold)[0]]
        else:
            indices = indices[np.where(v > threshold)[0]]
    else:
        n_f = int(min(portion * data.shape[1], len(v)))
        indices = indices[np.argsort(-v)[:n_f]]

    indices = np.sort(indices)

    return indices



def calculate_concordance_correlation_coefficient(u, v):
    '''
    This function calculates the concordance correlation coefficient between two input 1-D numpy arrays.

    Parameters:
    -----------
    u: 1-D numpy array of a variable
    v: 1-D numpy array of a variable

    Returns:
    --------
    ccc: a numeric value of concordance correlation coefficient between the two input variables.
    '''
    a = 2 * np.mean((u - np.mean(u)) * (v - np.mean(v)))
    b = np.mean(np.square(u - np.mean(u))) + np.mean(np.square(v - np.mean(v))) + np.square(np.mean(u) - np.mean(v))
    ccc = a/b
    return ccc



def generalization_feature_selection(data1, data2, measure, cutoff):
    '''
    This function uses the Pearson correlation coefficient to select the features that are generalizable
    between data1 and data2.

    Parameters:
    -----------
    data1: 2D numpy array of the first dataset with a size of (n_samples_1, n_features)
    data2: 2D numpy array of the second dataset with a size of (n_samples_2, n_features)
    measure: string. 'pearson' indicates the Pearson correlation coefficient;
        'ccc' indicates the concordance correlation coefficient. Default is 'pearson'.
    cutoff: a positive number for selecting generalizable features. If cutoff < 1, this function selects
        the features with a correlation coefficient >= cutoff. If cutoff >= 1, it must be an
        integer indicating the number of features to be selected based on correlation coefficient.

    Returns:
    --------
    fid: 1-D numpy array containing the indices of selected features.
    '''
    cor1 = np.corrcoef(np.transpose(data1))
    cor2 = np.corrcoef(np.transpose(data2))
    num = data1.shape[1]
    cor = []
    if measure == 'pearson':
        for i in range(num):
            cor.append(np.corrcoef(np.vstack((list(cor1[:i, i]) + list(cor1[(i + 1):, i]),
                       list(cor2[:i, i]) + list(cor2[(i + 1):, i]))))[0, 1])
    elif measure == 'ccc':
        for i in range(num):
            cor.append(calculate_concordance_correlation_coefficient(np.array(list(cor1[:i, i]) + list(cor1[(i + 1):, i])),
                np.array(list(cor2[:i, i]) + list(cor2[(i + 1):, i]))))
    cor = np.array(cor)
    fid = np.argsort(-cor)[:int(cutoff)]
    return fid



def coxen_single_drug_gene_selection(source_data, target_data, drug_response_data, drug_response_col, tumor_col,
                                     prediction_power_measure='pearson', num_predictive_gene=100,
                                     generalization_power_measure='ccc', num_generalizable_gene=50,
                                     multi_drug_mode=False):
    '''
    This function selects genes for drug response prediction using the COXEN approach. The COXEN approach is
    designed for selecting genes to predict the response of tumor cells to a specific drug. This function
    assumes no missing data exist.

    Parameters:
    -----------
    source_data: pandas data frame of gene expressions of tumors, for which drug response is known. Its size is
        [n_source_samples, n_features].
    target_data: pandas data frame of gene expressions of tumors, for which drug response needs to be predicted.
        Its size is [n_target_samples, n_features]. source_data and target_data have the same set
        of features and the orders of features must match.
    drug_response_data: pandas data frame of drug response values for a drug. It must include a column of drug
        response values and a column of tumor IDs.
    drug_response_col: non-negative integer or string. If integer, it is the column index of drug response in
        drug_response_data. If string, it is the column name of drug response.
    tumor_col: non-negative integer or string. If integer, it is the column index of tumor IDs in drug_response_data.
        If string, it is the column name of tumor IDs.
    prediction_power_measure: string. 'pearson' uses the absolute value of Pearson correlation coefficient to
        measure prediction power of gene; 'mutual_info' uses the mutual information to measure prediction power
        of gene. Default is 'pearson'.
    num_predictive_gene: positive integer indicating the number of predictive genes to be selected.
    generalization_power_measure: string. 'pearson' indicates the Pearson correlation coefficient;
        'ccc' indicates the concordance correlation coefficient. Default is 'ccc'.
    num_generalizable_gene: positive integer indicating the number of generalizable genes to be selected.
    multi_drug_mode: boolean, indicating whether the function runs as an auxiliary function of COXEN
        gene selection for multiple drugs. Default is False.

    Returns:
    --------
    indices: 1-D numpy array containing the indices of selected genes, if multi_drug_mode is False;
    1-D numpy array of indices of sorting all genes according to their prediction power, if multi_drug_mode is True.
    '''

    if isinstance(drug_response_col, str):
        drug_response_col = np.where(drug_response_data.columns == drug_response_col)[0][0]

    if isinstance(tumor_col, str):
        tumor_col = np.where(drug_response_data.columns == tumor_col)[0][0]

    drug_response_data = drug_response_data.copy()
    drug_response_data = drug_response_data.iloc[np.where(np.isin(drug_response_data.iloc[:, tumor_col],
                                                                  source_data.index))[0], :]

    source_data = source_data.copy()
    source_data = source_data.iloc[np.where(np.isin(source_data.index, drug_response_data.iloc[:, tumor_col]))[0], :]

    # Remove the genes that do not change over cancer cases
    source_std_id = select_features_by_variation(source_data, variation_measure='std', threshold=0.00000001)
    target_std_id = select_features_by_variation(target_data, variation_measure='std', threshold=0.00000001)
    std_id = np.sort(np.intersect1d(source_std_id, target_std_id))
    source_data = source_data.iloc[:, std_id]
    target_data = target_data.copy()
    target_data = target_data.iloc[:, std_id]

    # Perform the first step of COXEN approach to select predictive genes. To avoid exceeding the memory limit,
    # the prediction power of genes is calculated in batches.
    batchSize = 20000
    numBatch = int(np.ceil(source_data.shape[1]/batchSize))
    prediction_power = np.empty((source_data.shape[1], 1))
    prediction_power.fill(np.nan)
    for i in range(numBatch):
        startIndex = i*batchSize
        endIndex = min((i+1)*batchSize, source_data.shape[1])

        if prediction_power_measure == 'pearson':
            cor_i = np.corrcoef(np.vstack((np.transpose(source_data.iloc[:, startIndex:endIndex].loc[drug_response_data.iloc[:, tumor_col],
                :].values), np.reshape(drug_response_data.iloc[:, drug_response_col].values, (1, drug_response_data.shape[0])))))
            prediction_power[startIndex:endIndex, 0] = abs(cor_i[:-1, -1])

        if prediction_power_measure == 'mutual_info':
            mi = mutual_info_regression(X=source_data.iloc[:, startIndex:endIndex].loc[drug_response_data.iloc[:, tumor_col], :].values,
                                y=drug_response_data.iloc[:, drug_response_col].values)
            prediction_power[startIndex:endIndex, 0] = mi

    if multi_drug_mode:
        indices = np.argsort(-prediction_power[:, 0])
        return std_id[indices]

    num_predictive_gene = int(min(num_predictive_gene, source_data.shape[1]))
    gid1 = np.argsort(-prediction_power[:, 0])[:num_predictive_gene]

    # keep only predictive genes for source and target data
    source_data = source_data.iloc[:, gid1]
    target_data = target_data.iloc[:, gid1]
    num_generalizable_gene = int(min(num_generalizable_gene, len(gid1)))

    # perform the second step of COXEN approach to select generalizable genes among the predictive genes
    if isinstance(generalization_power_measure, list):
        indices = {}
        for g in generalization_power_measure:
            gid2 = generalization_feature_selection(source_data.values, target_data.values, g, num_generalizable_gene)
            indices[g] = np.sort(std_id[gid1[gid2]])
    else:
        gid2 = generalization_feature_selection(source_data.values, target_data.values, generalization_power_measure,
                                                num_generalizable_gene)
        indices = np.sort(std_id[gid1[gid2]])

    return indices



def coxen_multi_drug_gene_selection(source_data, target_data, drug_response_data, drug_response_col, tumor_col, drug_col,
                         prediction_power_measure='pearson', num_predictive_gene=100, num_generalizable_gene=50):
    '''
    This function uses the COXEN approach to select genes for predicting the response of multiple drugs.
    It assumes no missing data exist. For each drug, this functions ranks the genes according to their power
    of predicting the response of the drug. The union of an equal number of predictive genes for every drug
    will be generated, and its size must be at least num_predictive_gene. Then, num_generalizable_gene
    generalizable genes will be selected from the candidate pool.

    Parameters:
    -----------
    source_data: pandas data frame of gene expressions of tumors, for which drug response is known. Its size is
        [n_source_samples, n_features].
    target_data: pandas data frame of gene expressions of tumors, for which drug response needs to be predicted.
        Its size is [n_target_samples, n_features]. source_data and target_data have the same set
        of features and the orders of features must match.
    drug_response_data: pandas data frame of drug response that must include a column of drug response values,
        a column of tumor IDs, and a column of drug IDs.
    drug_response_col: non-negative integer or string. If integer, it is the column index of drug response in
        drug_response_data. If string, it is the column name of drug response.
    tumor_col: non-negative integer or string. If integer, it is the column index of tumor IDs in drug_response_data.
        If string, it is the column name of tumor IDs.
    drug_col: non-negative integer or string. If integer, it is the column index of drugs in drug_response_data.
        If string, it is the column name of drugs.
    prediction_power_measure: 'pearson' or 'mutual_info'. 'pearson' uses the absolute value of Pearson correlation
        coefficient to measure prediction power of a gene; 'mutual_info' uses the mutual information to measure
        prediction power of a gene.
    num_predictive_gene: positive integer indicating the number of predictive genes to be selected.
    num_generalizable_gene: positive integer indicating the number of generalizable genes to be selected.

    Returns:
    --------
    indices: 1-D numpy array containing the indices of selected genes.
    '''

    if isinstance(drug_response_col, str):
        drug_response_col = np.where(drug_response_data.columns == drug_response_col)[0][0]

    if isinstance(tumor_col, str):
        tumor_col = np.where(drug_response_data.columns == tumor_col)[0][0]

    if isinstance(drug_col, str):
        drug_col = np.where(drug_response_data.columns == drug_col)[0][0]

    drug_response_data = drug_response_data.copy()
    drug_response_data = drug_response_data.iloc[np.where(np.isin(drug_response_data.iloc[:, tumor_col],
                                                                  source_data.index))[0], :]
    drugs = np.unique(drug_response_data.iloc[:, drug_col])

    source_data = source_data.copy()
    source_data = source_data.iloc[np.where(np.isin(source_data.index, drug_response_data.iloc[:, tumor_col]))[0], :]

    source_std_id = select_features_by_variation(source_data, variation_measure='std', threshold=0.00000001)
    target_std_id = select_features_by_variation(target_data, variation_measure='std', threshold=0.00000001)
    std_id = np.sort(np.intersect1d(source_std_id, target_std_id))
    source_data = source_data.iloc[:, std_id]
    target_data = target_data.copy()
    target_data = target_data.iloc[:, std_id]

    num_predictive_gene = int(min(num_predictive_gene, source_data.shape[1]))

    gene_rank = np.empty((len(drugs), source_data.shape[1]))
    gene_rank.fill(np.nan)
    gene_rank = pd.DataFrame(gene_rank, index=drugs)
    for d in range(len(drugs)):
        idd = np.where(drug_response_data.iloc[:, drug_col] == drugs[d])[0]
        response_d = drug_response_data.iloc[idd, :]
        temp_rank = coxen_single_drug_gene_selection(source_data, target_data, response_d,
            drug_response_col, tumor_col, prediction_power_measure, num_predictive_gene=None,
            generalization_power_measure=None, num_generalizable_gene=None, multi_drug_mode=True)
        gene_rank.iloc[d, :len(temp_rank)] = temp_rank

    for i in range(int(np.ceil(num_predictive_gene/len(drugs))), source_data.shape[1]+1):
        gid1_temp = np.unique(np.reshape(gene_rank.iloc[:, :i].values, (1, gene_rank.shape[0]*i))[0, :])
        gid1_temp = gid1_temp[np.where(np.invert(np.isnan(gid1_temp)))[0]]
        if len(gid1_temp) >= num_predictive_gene:
            break
    gid1 = gid1_temp.astype(np.int64)

    # keep only predictive genes for source and target data
    temp_source_data = source_data.iloc[:, gid1]
    temp_target_data = target_data.iloc[:, gid1]
    num_generalizable_gene = int(min(num_generalizable_gene, len(gid1)))

    # perform the second step of COXEN approach to select generalizable genes among the predictive genes
    indices = generalization_feature_selection(data1=temp_source_data.values, data2=temp_target_data.values,
        measure='ccc', cutoff=num_generalizable_gene)
    indices = np.sort(std_id[gid1[indices]])

    return indices