import random
import numpy
import scipy


def gene_level_score(scores, method='highest'):
    """
    Computes a gene level score from list of
    binding scores belonging to the same gene
    Methods: lowest, highest, mean
    """
    if method == 'lowest':
        return min(scores)
    elif method == 'highest':
        return max(scores)
    elif method == 'mean':
        return float(sum(scores)) / max(len(scores), 1)
    else:
        raise ValueError('could not find %s in available methods' % (method))

def randomize(scores, original_scores, operator=operator.gt, cutoff=0.01):
    """
    Random sample from data and compare the random set
    with the original set and return the p-value if
    it meets the conditions
    """
    random_scores = [scores[k] for k in random.sample(scores, len(original_scores))]
    try:
        z, p = scipy.stats.mannwhitneyu(original_scores, random_scores)
        # two-tailed test fix
        p *= 2
        # median of original set should be lower and cutoff achieved
        if operator(numpy.median(original_scores), numpy.median(random_scores)) and p < cutoff:
            return p
    except ValueError:
        pass
    return None

def gene_set_p_value(scores, original_scores, operator=operator.gt, cutoff=0.01, times=1000):
    """
    Computes a p-value by random sampling from data
    and doing Mann-Whitney U test and median comparison
    """
    # serial execution
    p_values = [randomize(scores, original_scores, operator=operator, cutoff=cutoff) for _ in xrange(times)]
    # count the p values
    n = len([_ for _ in p_values if _ != None])
    return max(1 - (n / float(times)), 1 / float(times))
