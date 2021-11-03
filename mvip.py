import os
import pandas
import numpy
import scipy
import scipy.stats
import statsmodels.sandbox.stats.multicomp as stats_models
import random
import correlations as corr
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from plotly import express as px
import pandas

H0 = "There is no difference between replicate self correlation and control correlation"
random.seed(8)

CMAP_TYPE = {'EMPTY': '#1f77b4', 'REF':'#2CA02C',
             'VAR': '#FF7F0E', 'VAR_REF': '#8C564B'}


def wilcoxon_test(self_corr_matrix, control_corr_matrix):
    self_corr_median = corr.correlation_median_row(self_corr_matrix)
    control_corr_median = numpy.median(control_corr_matrix, axis=1)
    test_result = scipy.stats.wilcoxon(self_corr_median, control_corr_median)
    return test_result.pvalue


def wilcoxon_ranksums(sample1, sample2):
    test_result = scipy.stats.ranksums(sample1, sample2)
    return test_result.pvalue


def kruskal_wallis_test(self_corr_upper, control_corr_matrix):
    test_result = scipy.stats.kruskal(self_corr_upper, control_corr_matrix.flatten().tolist())
    return test_result.pvalue


def kruskal_wallis_on_medians(self_corr_matrix, control_corr_matrix):
    self_corr_median = corr.correlation_median_row(self_corr_matrix)
    control_corr_median = numpy.median(control_corr_matrix, axis=1)
    test_result = scipy.stats.kruskal(self_corr_median, control_corr_median)
    return test_result.pvalue


## Matrices required for eVIP
## A. Wild type vs Wild type 	= Wild type self-correlation 	= wt_wt
## B. Mutant vs Mutant 		= Mutant self-correlation 	= mut_mut
## C. Wild type vs Control 	= Wild type control correlation = wt_ctl
## D. Mutant vs Control 	= Mutant control correlation 	= mut_ctl
## E. Wild type vs Mutant 	= Wild type mutant correlation 	= wt_mut

class Morphology_VIP(object):

    def __init__(self, metadata, corr_matrix, treatment_samples, control_samples, controls_value="control", perturbation_field="Metadata_x_mutation_status", controls_field="Metadata_broad_sample_type", plate_field="Metadata_Plate"):
        self.metadata = metadata
        self.corr_matrix = corr_matrix
        self.treatment_samples = treatment_samples
        self.control_samples = control_samples
        self.ctl_mask = metadata[controls_field] == controls_value
        self.perturbation_field = perturbation_field
        self.plate_field = plate_field
        self.index = {"name": "genes", "children": []}
        self.test_cols = ["wild_type", "wt_samples", "mutant", "mut_samples", "wt_has_effect", "mut_has_effect", "wt_mut_difference"]


    def evaluate(self, wild_type, mutant, create_images=False, false_positives=False,
                 images_dir='./'):
        results = {}
        results["wild_type"] = wild_type
        results["mutant"] = mutant
        # Get indices of data
        wt_mask = self.metadata[self.perturbation_field] == wild_type
        wt_index = self.metadata[wt_mask].index
        mut_mask = self.metadata[self.perturbation_field] == mutant
        mut_index = self.metadata[mut_mask].index

        if not false_positives: # Regular evaluation
            if len(mut_index) > self.treatment_samples:
                mut_index = list(mut_index)
                random.shuffle(mut_index)
                mut_index = mut_index[0:self.treatment_samples]
            if len(wt_index) > self.treatment_samples:
                wt_index = list(wt_index)
                random.shuffle(wt_index)
                wt_index = wt_index[0:self.treatment_samples]
        else: # False positives evaluation
            tmp_index = list(mut_index)
            random.shuffle(tmp_index)
            mut_index = tmp_index[0:self.treatment_samples]
            wt_index = tmp_index[-self.treatment_samples:]

        results["wt_samples"] = len(wt_index)
        results["mut_samples"] = len(mut_index)
        # Copy matrices
        self.wt_wt = corr.sample_rectangular_matrix(wt_index, wt_index, self.corr_matrix)
        self.mut_mut = corr.sample_rectangular_matrix(mut_index, mut_index, self.corr_matrix)
        self.wt_ctl = corr.allele_to_control_matrix(wt_index, self.metadata, self.plate_field, self.ctl_mask, self.control_samples, self.corr_matrix)
        self.mut_ctl = corr.allele_to_control_matrix(mut_index, self.metadata, self.plate_field, self.ctl_mask, self.control_samples, self.corr_matrix)
        self.wt_mut = corr.sample_rectangular_matrix(wt_index, mut_index, self.corr_matrix)

        # Add json index entry
        # search for wild_type
        wt_idx = [x for x in range(len(self.index["children"])) if self.index["children"][x]["name"] == wild_type ]
        if len(wt_idx) == 0:
            self.index["children"].append({"name": wild_type, "children": []})
            wt_idx = -1
        else:
            wt_idx = wt_idx[0]
        self.index["children"][wt_idx]["children"].append({"name": mutant, "pair": wild_type + "_" + mutant})

        # Run tests
        if create_images:
            self.create_plots(results, images_dir)
        return self.statistical_tests_medians(results)

    def statistical_tests(self, results):
        iu = numpy.triu_indices(self.treatment_samples, 1)
        results["wt_has_effect"] = kruskal_wallis_test(self.wt_wt[iu], self.wt_ctl)
        results["mut_has_effect"] = kruskal_wallis_test(self.mut_mut[iu], self.mut_ctl)
        results["wt_mut_difference"] = kruskal_wallis_test(self.wt_wt[iu], self.wt_mut)
        return results
        

    def statistical_tests_medians(self, results):
        wt_pvalue = wilcoxon_test(self.wt_wt, self.wt_ctl)
        mut_pvalue = wilcoxon_test(self.mut_mut, self.mut_ctl)
        results["wt_has_effect"] = wt_pvalue
        results["mut_has_effect"] = mut_pvalue

        # Fix (see next comment):
        wt_mut_pvalue = wilcoxon_test(self.mut_mut, self.wt_mut)
        results["wt_mut_difference"] = wt_mut_pvalue
        return results

        # The following does not make too much sense. We are comparing mutants against wild types by comparing
        # the distribution of medians in the mutant matrix vs the distribution of median rows AND columns in the cross correlation matrix.
        # Using both, rows and columns is useful if we compare mut_mut VS wt_wt VS wt_mut. i.e., when we compare three distributions.
        # In the current version of the test, we can pair mut_mut medians with wt_mut medians in a Wilconxon test (right?)
        self_corr_median = corr.correlation_median_row(self.mut_mut)
        cross_corr_row = numpy.median(self.wt_mut, axis=0)
        cross_corr_col = numpy.median(self.wt_mut, axis=1)
        test_result = scipy.stats.kruskal(self_corr_median, numpy.concatenate([cross_corr_row, cross_corr_col]))
        wt_mut_pvalue = test_result.pvalue
        results["wt_mut_difference"] = wt_mut_pvalue

        return results


    def create_plots(self, results, images_dir):
        # Dot plots
        row = numpy.median(self.wt_mut, axis=0)
        col = numpy.median(self.wt_mut, axis=1)
        cross_corr_median = numpy.concatenate([row, col])
        dots = pandas.DataFrame()
        # dots = dots.append( [{"Sample":"REF_CTL","Correlation":k} for k in numpy.median(self.wt_ctl, axis=1)] )
        dots = dots.append( [{"Sample":"REF","Correlation":k} for k in corr.correlation_median_row(self.wt_wt)] )
        dots = dots.append( [{"Sample":"VAR_REF","Correlation":k} for k in cross_corr_median] )
        dots = dots.append( [{"Sample":"VAR","Correlation":k} for k in corr.correlation_median_row(self.mut_mut)] )
        # dots = dots.append( [{"Sample":"VAR_CTL","Correlation":k} for k in numpy.median(self.mut_ctl, axis=1)] )

        fig = px.box(dots, x='Sample', y='Correlation', color='Sample',
                     points='all', color_discrete_map=CMAP_TYPE)
        fig.update_layout(showlegend=False,
                          yaxis_range=(-0.2, 1.0),
                          # margin=dict(l=0, r=0, t=0, b=0)
                          )

        output_name = results["mutant"] + "_dots"
        output_name = os.path.join(images_dir, output_name)
        fig.write_image(output_name + '.png')
        with open(output_name + '.json', 'w') as f:
            f.write(fig.to_json(pretty=True))

        # Matrix plots
        matrices = [self.wt_wt, self.wt_mut, self.mut_mut]
        samples = ["REF_REF", "VAR_REF", "VAR_VAR"]

        # zmin, zmax = min(map(numpy.min, matrices)), max(map(numpy.max, matrices))
        # zmin, zmax = dots.Correlation.min(), dots.Correlation.max()
        zmin, zmax = -0.2, 1.0
        zranges = {'zmin': zmin, 'zmax': zmax}
        os.makedirs(images_dir, exist_ok=True)
        fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.05, subplot_titles=samples)
        for i, matrix in enumerate(matrices, 1):
            scaled_matrix = (numpy.clip(matrix, zmin, zmax) - zmin) / (zmax - zmin)
            hmap = go.Heatmap(z=scaled_matrix,
                              colorscale=[(0, "blue"), (0.5, "white"), (1, "red")],
                              **zranges)
            fig.add_trace(hmap, row=1, col=i)
        fig.update_traces(showscale=False)
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False, autorange='reversed')
        fig.update_layout(#margin=dict(l=0, r=0, t=30, b=0),
                          height=340
                          )

        output_name = results["mutant"] + "_matrices"
        output_name = os.path.join(images_dir, output_name)
        fig.write_image(output_name + '.png')
        with open(output_name + '.json', 'w') as f:
            f.write(fig.to_json(pretty=True))


    def search_wild_type(self, mutant, ignore_wt=False):
        # Search the corresponding wild type
        wild_type = mutant.split("_")[0] + "_WT"
        wt_search = self.metadata[self.perturbation_field].str.find(wild_type) != -1
        wt_alternatives = self.metadata[wt_search][self.perturbation_field].unique()
        if len(wt_alternatives) == 1:
            wild_type = wt_alternatives[0]
        elif len(wt_alternatives) > 1:
            wt_alternatives = [alt for alt in wt_alternatives if alt.find(".c") != -1]
            wild_type = wt_alternatives[0]
        elif len(wt_alternatives) == 0:
            if not ignore_wt:
                print("No wild type for:", mutant)
                wild_type = None
            else:
                wild_type = mutant
        return wild_type

 
    def test_allele_set(self, alleles, create_images=False, false_positives=False, null_distribution=None,
            images_dir='./'):
        self.null_dist = null_distribution
        results = pandas.DataFrame()
        for mutant in alleles:
            wild_type = self.search_wild_type(mutant, ignore_wt=false_positives)
            if wild_type is not None:
                r = self.evaluate(wild_type, mutant, create_images=create_images, false_positives=false_positives, images_dir=images_dir)
                results = results.append(r, ignore_index=True)
    
        return results[self.test_cols]


    def adjust_pvalues(self, results, Q=0.05):
        m = len(results)
        test_fields = ["wt_has_effect", "mut_has_effect", "wt_mut_difference"]
        all_tests = []

        ## Adjust p-values of all tests
        for f in test_fields:        
            sig, adj, a, b = stats_models.multipletests(results[f], method="fdr_bh")
            results["is_sig_" + f] = sig
            results["adjusted_" + f] = adj

        # Run the Variant-Impact Phenotyping test
        wt_has_effect = results["is_sig_wt_has_effect"]
        mut_has_effect = results["is_sig_mut_has_effect"]
        wt_mut_diff = results["is_sig_wt_mut_difference"]

        results["prediction"] = "NI"
        results.loc[~ wt_has_effect & mut_has_effect, "prediction"] = "GOF"
        results.loc[wt_has_effect & ~ mut_has_effect, "prediction"] = "LOF"
        results.loc[wt_has_effect & mut_has_effect & wt_mut_diff, "prediction"] = "COF"
        results.loc[wt_has_effect & mut_has_effect & ~ wt_mut_diff, "prediction"] = "NT"
        
        return results


    def eval_pvalues(self, results, threshold=0.05):
        m = len(results)
        test_fields = ["wt_has_effect", "mut_has_effect", "wt_mut_difference"]
        all_tests = []

        for f in test_fields:
            results["is_sig_" + f] = results[f] < threshold

        wt_has_effect = results["is_sig_wt_has_effect"]
        mut_has_effect = results["is_sig_mut_has_effect"]
        wt_mut_diff = results["is_sig_wt_mut_difference"]
        
        results["prediction"] = "NI"
        results.loc[~ wt_has_effect & mut_has_effect, "prediction"] = "GOF"
        results.loc[wt_has_effect & ~ mut_has_effect, "prediction"] = "LOF"
        results.loc[wt_has_effect & mut_has_effect & wt_mut_diff, "prediction"] = "COF"
        results.loc[wt_has_effect & mut_has_effect & ~ wt_mut_diff, "prediction"] = "NT"

        return results


class Morphology_VIP_CNN_Features(Morphology_VIP):


    def statistical_tests_medians(self, results):
        self.test_cols = ['wild_type', 'wt_samples', "mutant", 'mut_samples', 'impact_test', 'strength_test', 'directionality_test', 'power_test']

        wt_self_corr = corr.correlation_median_row(self.wt_wt)
        mut_self_corr = corr.correlation_median_row(self.mut_mut)
        cross_corr_row = numpy.median(self.wt_mut, axis=0)
        cross_corr_col = numpy.median(self.wt_mut, axis=1)
        wt_mut_cross = numpy.concatenate([cross_corr_row, cross_corr_col])
        impact_test = scipy.stats.kruskal(wt_self_corr, mut_self_corr, wt_mut_cross)
        results["impact_test"] = impact_test.pvalue

        strength_pvalue = wilcoxon_ranksums(wt_self_corr, mut_self_corr)
        results["strength_test"] = strength_pvalue

        power_pvalue = wilcoxon_ranksums(wt_mut_cross, self.null_dist)
        results["power_test"] = power_pvalue

        wt_signal = numpy.median(wt_self_corr)
        mut_signal = numpy.median(mut_self_corr)
        results["directionality_test"] =  mut_signal > wt_signal

        return results

    def adjust_pvalues(self, results, Q=0.05):
        m = len(results)
        test_fields = ["impact_test", "strength_test", "power_test"]
        all_tests = []

        ## Adjust p-values of all tests
        for f in test_fields:        
            sig, adj, a, b = stats_models.multipletests(results[f], method="fdr_bh")
            results["is_sig_" + f] = sig
            results["adjusted_" + f] = adj

        # Run the Variant-Impact Phenotyping test
        impact_test = results["is_sig_impact_test"]
        strength_test = results["is_sig_strength_test"]
        power_test = results["is_sig_power_test"]
        direction_test = results["directionality_test"] > 0.0

        results["prediction"] = "NI"
        results.loc[impact_test & strength_test & direction_test , "prediction"] = "GOF"
        results.loc[impact_test & strength_test & ~direction_test, "prediction"] = "LOF"
        results.loc[impact_test & ~ strength_test, "prediction"] = "COF"
        results.loc[~impact_test & power_test, "prediction"] = "NT"

        return results


    def update_index(self, results, field_name):
        for idx, record in results.iterrows():
            # search for wild_type
            wt_idx = [x for x in range(len(self.index["children"])) if self.index["children"][x]["name"] == record["wild_type"] ]
            if len(wt_idx) > 0:
                wt_idx = wt_idx[0]
                # search for mutant
                mut_idx = [x for x in range(len(self.index["children"][wt_idx]["children"])) if self.index["children"][wt_idx]["children"][x]["name"] == record["mutant"]]
                if len(mut_idx) > 0:
                    mut_idx = mut_idx[0]
                    # Add test results to the index:
                    self.index["children"][wt_idx]["children"][mut_idx][field_name] = dict(record)
                else:
                    print(record["mutant"], "missing in index")
                self.index["children"][wt_idx]["children"].sort(key=lambda r: r["name"])
            else:
                print(record["wild_type"], "missing in index")
        self.index["children"].sort(key=lambda r: r["name"])
         

