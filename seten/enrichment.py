"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten/blob/master/LICENSE
for full license details.
"""
import operator
from itertools import chain
from collections import defaultdict
from seten.statistics import gene_level_score
from seten.statistics import gene_set_p_value
from seten.statistics import functional_p_value
from seten.statistics import p_values_correction
from seten.mapping import search


def _overlapping_genes(genes, gs_genes):
    """
    Returns overlapping genes in a gene set
    """
    return list(set.intersection(set(genes), gs_genes))


def collect_scores(path='', mapping=None, index=4, method='highest'):
    """
    Collect scores from a BED file
    """
    scores = defaultdict(list)
    with open(path, 'r') as rows:
        for row in rows:
            cols = row.strip().split()
            # value available?
            if cols[index]:
                genes = search(cols[0], cols[1], cols[2], mapping)
                # any genes available?
                if genes:
                    # collect values for each gene
                    for gene in genes:
                        # cast value as float
                        scores[gene].append(float(cols[index]))

    # compute a gene level score from list of binding scores belonging to the
    # same gene
    scores = dict((g, gene_level_score(s, method=method))
                  for g, s in scores.iteritems())
    return scores


def collect_gene_sets(path='', resources_dir='resources'):
    """
    Collects gene sets as a list
    """
    gene_sets = []
    try:
        with open(path, 'r') as rows:
            for row in rows:
                cols = row.strip().split('\t')
                gene_sets.append(dict(name=cols[0],
                                      genes=list(set(cols[2:])),
                                      size=len(set(cols[2:]))))
    except IOError:
        pass
    gc_size = len(set(chain.from_iterable(
        [gene_set['genes'] for gene_set in gene_sets])))
    return gene_sets, gc_size


def gene_set_enrichment(scores, gene_set, operator=operator.gt, cutoff=0.05, count=3):
    """
    Applies gene set enrichment on the data
    """
    overlap = _overlapping_genes(scores.keys(), gene_set['genes'])
    if len(overlap) >= count:
        overlap_scores = [scores[str(g)] for g in overlap]
        p_value = gene_set_p_value(scores.values(), overlap_scores, operator)
        # there is no cutoff anymore
        return dict(name=gene_set['name'],
                    genes=len(overlap),
                    size=gene_set['size'],
                    p_value=p_value)
    return None


def functional_enrichment(genes, gene_set, gc_size, cutoff=0.05, count=3):
    """
    Applies functional enrichment on the data
    """
    overlap = _overlapping_genes(genes, gene_set['genes'])
    if len(overlap) >= count:
        p_value = functional_p_value(
            len(gene_set['genes']), len(overlap), gc_size, len(genes))
        return dict(name=gene_set['name'],
                    genes=len(overlap),
                    size=gene_set['size'],
                    p_value=p_value)
    return None


def integrated_enrichment(scores, gene_set, gc_size, operator=operator.gt, cutoff=0.05, count=3):
    """
    Applies integrated enrichment on the data
    """
    overlap = _overlapping_genes(scores.keys(), gene_set['genes'])
    if len(overlap) >= count:
        overlap_scores = [scores[str(g)] for g in overlap]
        gse_p_value = gene_set_p_value(
            scores.values(), overlap_scores, operator)
        fe_p_value = functional_p_value(
            len(gene_set['genes']), len(overlap), gc_size, len(scores.keys()))
        return dict(name=gene_set['name'],
                    genes=len(overlap),
                    size=gene_set['size'],
                    gse_p_value=gse_p_value,
                    fe_p_value=fe_p_value)
    return None


def correct_results(results, alpha=0.05, method='bh'):
    """
    Corrects p-values in integrated and functional enrichment
    """
    results_cor = []
    fe_p_values = [result['fe_p_value'] for result in results]
    reject, fe_p_values_cor = p_values_correction(
        fe_p_values, alpha=alpha, method=method)
    for i, result in enumerate(results):
        result['fe_p_value_corrected'] = fe_p_values_cor[i]
        result['integrated_p_value'] = 1 - \
            ((1 - result['gse_p_value']) *
             (1 - result['fe_p_value_corrected']))
        results_cor.append(result)
    return results_cor
