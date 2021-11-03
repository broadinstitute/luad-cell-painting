import numpy as np
import scipy
import scipy.stats
import random
from tqdm import tqdm

# Most of the code in this module assumes that the original data comes
# in a Pandas DataFrame.

random.seed(8)

################################################################################
## COMPUTE PEARSON CORRELATIONS
################################################################################

# Compute correlations between samples and save them in a list
def correlation_upper_triangle(sample, first_feature):
    corr_values = []
    for i in range(len(sample)):
        x = sample.iloc[i][first_feature:]
        for j in range(i+1, len(sample)):
            y = sample.iloc[j][first_feature:]
            corr, p_val = scipy.stats.pearsonr(x, y)
            corr_values.append(corr)
    return corr_values

# Compute the full correlation matrix between all samples
# Complexity (N^2)/2
def correlation_matrix(sample, first_feature):
    return np.corrcoef(sample)


################################################################################
## SAMPLE SUB-MATRICES FROM A CORRELATION MATRIX
################################################################################

# Get the upper triangle of a symmetric submatrix
# Takes as input the row indices
# Returns the list of values in the upper triangle
def sample_upper_triangle(sample_rows, cmatrix):
    sub_matrix = cmatrix[sample_rows, sample_rows]
    values = []
    for i in range(len(sample_rows)):
        for j in range(i+1, len(sample_rows)):
            row = sample_rows[i]
            col = sample_rows[j]
            values.append(cmatrix[row, col])
    return values

# Extract a non-symmetric submatrix from the correlation matrix
# Assumes sample1 and sample 2 are disjoint sets of replicates
def sample_rectangular_matrix(sample1, sample2, cmatrix):
    sub_matrix = cmatrix[sample1, :]
    sub_matrix = sub_matrix[:, sample2]
    return sub_matrix


# Compute the mean row of a self-correlation matrix, ignoring the
# diagonal values (assuming these are the maximum values)
def correlation_median_row(cmatrix):
    X = np.sort(cmatrix, axis=0)
    median_idx = int( (X.shape[0]-1)/2 )
    return X[median_idx,:]


# Extract asymmetric matrix of alleles vs controls
# This heavily depends on metadata. For each allele, we want control that are in the same plate
def allele_to_control_matrix(allele_index, metadata, plate_field, ctl_mask, control_samples, corr_matrix):
    control_corr_matrix = []
    for i in allele_index:
        a = metadata.loc[i]
        #WARN: Here a complex indexing is applied using both indices and values
        # of ctl_mask. More info at https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.loc.html
        # "An alignable boolean Series. The index of the key will be aligned
        # before masking.""
        ctl_idx = metadata[metadata[plate_field] == a[plate_field]].loc[ctl_mask].index
        if len(ctl_idx) > control_samples:
            ctl_idx = list(ctl_idx)
            random.shuffle(ctl_idx)
            ctl_idx = ctl_idx[0:control_samples]
        elif len(ctl_idx) < control_samples:
            ctl_idx = list(metadata[ctl_mask].index)
            random.shuffle(ctl_idx)
            ctl_idx = ctl_idx[0:control_samples]
        m = sample_rectangular_matrix([i], ctl_idx, corr_matrix)
        #return
        control_corr_matrix.append(m)
    return np.concatenate(control_corr_matrix)

################################################################################
## COMPUTE STATISTICS ON CORRELATION DISTRIBUTIONS
################################################################################

# Creates groups of size "sample_size" from some rows in the
# correlation matrix "cmatrix". The process is repeated several times
# to create a large enough background distribution
# TODO: choose random groups of samples that do not share the same treatment
def null_distribution(sample_rows, cmatrix, sample_size, repeats=30):
    distribution = []
    for i in tqdm(range(repeats)):
        index = np.asarray([ k for k in sample_rows])
        np.random.shuffle(index)
        for j in range(0, len(index), sample_size):
            median_corr = median_correlation(index[j:j+sample_size], cmatrix)
            distribution.append(median_corr)
    return distribution

# For one particular treatment compute the median of replicates, 
# and also create a null distribution for comparison.
# Return mean of replicate correlations and 95th percentile of the null
def median_correlation(sample_index, cmatrix):
    corr_values = sample_upper_triangle(sample_index, cmatrix)
    median = np.median(corr_values)
    return median

# Compute replicate correlation for all treatments, and report
# the fraction of treatments whose replicate correlation is greater than
# the 95th percentile of the null.
def fraction_strong_test(treated_samples, treatments, corr_matrix, null, treatment_field, sample_size):
    null.sort()
    p95 = null[ int( 0.95*len(null) ) ]
    results = {}
    for t in tqdm(treatments):
        sample = treated_samples[treated_samples[treatment_field] == t]
        replicates = len(sample)
        if replicates >= sample_size:
            sample_data = sample.sample(sample_size)
            median = median_correlation(sample_data.index, corr_matrix)
            results[t] = median
            #replicates -= sample_size
    # Evaluate results:
    fraction = np.sum([m > p95 for m in results.values()])/len(results)
    print("Treatments tested:", len(results))
    print("At 95th percentile of the null")
    print("Fraction strong: {:5.2f}%".format(fraction*100))
    print("Null threshold: {:6.4f}".format(p95))
    return results

