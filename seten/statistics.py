"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import random
import operator
import numpy
from scipy import stats
from seten.multi_comp import fdr_correction, bonferroni_correction


def compute_gene_level_score(scores, method='highest'):
    """
    Computes a gene level score from list of
    binding scores belonging to the same gene
    Methods: lowest, highest, mean, median
    """
    if method == 'lowest':
        return min(scores)
    elif method == 'highest':
        return max(scores)
    elif method == 'mean':
        return numpy.mean(scores)
    elif method == 'median':
        return numpy.median(scores)
    else:
        raise ValueError('could not find %s in available methods' % (method))


def randomize(scores, overlap_scores, operator=operator.gt, cutoff=0.05):
    """
    Random sample from data and compare the random set
    with the overlap set and return the p-value if
    it meets the conditions
    """
    random_scores = random.sample(scores, len(overlap_scores))
    try:
        if operator(numpy.median(overlap_scores), numpy.median(random_scores)):
            statistic, pvalue = stats.mannwhitneyu(overlap_scores, random_scores)
            # two-tailed test fix
            pvalue *= 2
            if pvalue < cutoff:
                return pvalue
    except ValueError:
        pass
    return None


def compute_gse_pvalue(scores, overlap_scores,
        operator=operator.gt, cutoff=0.05, times=1000):
    """
    Computes a p-value by random sampling from data
    and doing Mann-Whitney U test and median comparison
    """
    # serial execution
    pvalues = [randomize(scores, overlap_scores, operator=operator,
                          cutoff=cutoff) for _ in xrange(times)]
    # count the p-values
    n = times - pvalues.count(None)
    return max(1 - (n / float(times)), 1 / float(times))


def compute_fe_pvalue(a, b, c, d, method='fishers'):
    """
    Returns p-value for functional enrichment
    Only tests for over-enrichment
    a: Overlap size
    b: Gene set size
    c: Data size
    d: Collections size
    """
    if method == 'fishers':
        odds_ratio, pvalue = stats.fisher_exact(
            [[a, b], [c-a, d-b]], alternative='greater')
        return pvalue
    else:
        raise ValueError('could not find %s in available methods' % (method))


def _correct_pvalues(pvalues, alpha=0.05, method='bh'):
    """
    Correct for multiple p-values
    """
    if method == 'bh':
        return fdr_correction(pvalues, alpha=alpha, method='indep')
    elif method == 'by':
        return fdr_correction(pvalues, alpha=alpha, method='negcorr')
    elif method == 'bon':
        return bonferroni_correction(pvalues, alpha=alpha)
    raise ValueError('could not find %s in available methods' % (method))


def correct_and_combine_pvalues(results,
        alpha=0.05, correction_method='bh', combination_method='fisher'):
    """
    Corrects p-values in from functional enrichment
    """
    results_final = []
    # raw p-values from functional enrichment results
    fe_pvalues = [result['fe_pvalue'] for result in results]
    # correct for p-values
    reject, fe_pvalues_cor = _correct_pvalues(
        fe_pvalues, alpha=alpha, method=correction_method
        )
    # combine corrected funtional enrichment p-values with
    # gene set enrichment p-values
    for i, result in enumerate(results):
        # add corrected pvalue to the results
        result['fe_pvalue_corrected'] = fe_pvalues_cor[i]
        # combine the two p-values
        statistic, combined_pvalue = stats.combine_pvalues([
            result['gse_pvalue'],
            result['fe_pvalue_corrected']
            ], method=combination_method)
        # add combined pvalue to the results
        result['combined_pvalue'] = combined_pvalue
        results_final.append(result)
    return results_final
