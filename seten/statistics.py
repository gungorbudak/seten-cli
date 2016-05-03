"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import random
import numpy as np
from scipy import stats
from seten.multi_comp import fdr_correction
from seten.multi_comp import bonferroni_correction


def compute_gene_level_score(scores, method='max'):
    """
    Computes a gene level score from list of
    binding scores belonging to the same gene
    Methods: min, max, mean, median
    """
    if method == 'min':
        return min(scores)
    elif method == 'max':
        return max(scores)
    elif method == 'mean':
        return np.mean(scores)
    elif method == 'median':
        return np.median(scores)
    elif method == 'sum':
        return np.sum(scores)
    else:
        raise ValueError('could not find %s in available methods' % (method))


def _gse_test(scores, overlap_scores):
    """
    Random sample from data and compare the random set
    with the overlap set and return the p-value if
    it meets the conditions
    """
    random_scores = random.sample(scores, len(overlap_scores))
    # Mann Whitney U test for testing if overlap scores
    # have significantly larger values than random scores
    statistic, pvalue = stats.mannwhitneyu(overlap_scores, random_scores)
    return pvalue


def compute_gse_pvalue(scores, overlap_scores, significance_cutoff=0.05, iters=1000):
    """
    Computes a p-value by random sampling from data
    and doing Mann-Whitney U test and median comparison
    """
    pvalues = [_gse_test(scores, overlap_scores) for _ in xrange(iters)]
    # count the p-values
    n = np.sum(np.array(pvalues, dtype=np.float64) < significance_cutoff)
    return max(1 - (n / float(iters)), 1 / float(iters))


def compute_fe_pvalue(a, b, c, d, method='fishers'):
    """
    Returns p-value for functional enrichment
    a: Overlap size
    b: Gene set size
    c: Data size
    d: Collections size
    """
    if method == 'fishers':
        odds_ratio, pvalue = stats.fisher_exact([[a, b], [c-a, d-b]])
        return pvalue
    else:
        raise ValueError('could not find %s in available methods' % (method))


def _correct_pvalues(pvalues, alpha=0.05, method='fdr'):
    """
    Correct for multiple p-values
    """
    if method == 'bh' or method == 'fdr':
        return fdr_correction(pvalues, alpha=alpha, method='indep')
    elif method == 'by':
        return fdr_correction(pvalues, alpha=alpha, method='negcorr')
    elif method == 'bon':
        return bonferroni_correction(pvalues, alpha=alpha)
    raise ValueError('could not find %s in available methods' % (method))


def correct_pvalues(results, alpha=0.05, method='fdr'):
    """
    Corrects p-values in from functional enrichment
    """
    pvalues = [result['fe_pvalue'] for result in results]
    # correct for p-values
    reject, pvalues_corr = _correct_pvalues(
        pvalues, alpha=alpha, method=method
        )
    for i in xrange(len(results)):
        results[i]['fe_pvalue_corr'] = pvalues_corr[i]
    return results
