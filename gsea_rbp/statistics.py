import random
import numpy
import scipy


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

def randomize(scores, overlap_scores, operator=operator.gt, cutoff=0.01):
    """
    Random sample from data and compare the random set
    with the overlap set and return the p-value if
    it meets the conditions
    """
    random_scores = random.sample(scores, len(overlap_scores))
    try:
        z, p = scipy.stats.mannwhitneyu(overlap_scores, random_scores)
        # two-tailed test fix
        p *= 2
        # check medians and cutoff
        if operator(numpy.median(overlap_scores), numpy.median(random_scores)) and p < cutoff:
            return p
    except ValueError:
        pass
    return None

def gene_set_p_value(scores, overlap_scores, operator=operator.gt, cutoff=0.01, times=1000):
    """
    Computes a p-value by random sampling from data
    and doing Mann-Whitney U test and median comparison
    """
    # serial execution
    p_values = [randomize(scores, overlap_scores, operator=operator, cutoff=cutoff) for _ in xrange(times)]

    # count the p values
    n = len([_ for _ in p_values if _ != None])

    return max(1 - (n / float(times)), 1 / float(times))
