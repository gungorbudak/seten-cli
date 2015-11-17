"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten/blob/master/LICENSE
for full license details.
"""
import random
import operator
import numpy
from scipy import stats
from seten.multi_comp import fdr_correction, bonferroni_correction


def gene_level_score(scores, method='highest'):
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
        z, p = stats.mannwhitneyu(overlap_scores, random_scores)
        # two-tailed test fix
        p *= 2
        # check medians and cutoff
        if operator(numpy.median(overlap_scores), numpy.median(random_scores)) and p < cutoff:
            return p
    except ValueError:
        pass
    return None

def gene_set_p_value(scores, overlap_scores, operator=operator.gt, cutoff=0.05, times=1000):
    """
    Computes a p-value by random sampling from data
    and doing Mann-Whitney U test and median comparison
    """
    # serial execution
    p_values = [randomize(scores, overlap_scores, operator=operator, cutoff=cutoff) for _ in xrange(times)]
    # count the p values
    n = len([_ for _ in p_values if _ != None])
    return max(1 - (n / float(times)), 1 / float(times))

def functional_p_value(gs_size, ov_size, gc_size, dt_size):
    """
    Does a Fisher's exact test
    """
    odds_ratio, p_value = stats.fisher_exact([[gs_size, ov_size], [gc_size, dt_size]])
    return p_value

def p_values_correction(p_values, alpha=0.05, method='bh'):
    """
    Correct for multiple p-values
    """
    if method == 'bh':
        return fdr_correction(p_values, alpha=alpha, method='indep')
    elif method == 'by':
        return fdr_correction(p_values, alpha=alpha, method='negcorr')
    elif method == 'bon':
        return bonferroni_correction(p_values, alpha=alpha)
    raise ValueError('could not find %s in available methods' % (method))
