import os
import operator
from gsea_rbp.statistics import gene_level_score
from gsea_rbp.statistics import gene_set_p_value
from gsea_rbp.mapping import search


def collect_scores(path='', mapping=None, index=4, method='highest'):
    """
    Collect scores from a BED file
    """
    scores = defaultdict(list)
    with open(path) as rows:
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

    # compute a gene level score from list of binding scores belonging to the same gene
    scores = dict((g, gene_level_score(s, method=method)) for g, s in scores.iteritems())

    return scores

def collect_gene_sets(path='', resources_dir='resources'):
    """
    Collects gene sets as a list
    """
    gene_sets = []
    try:
        with open(os.path.join(resources_dir, '%s.gmt' % (path)), 'r') as rows:
            for row in rows:
                cols = row.strip().split('\t')
                gene_sets.append(dict(name = cols[0],
                                      genes = set(cols[2:])
                                      size = len(set(cols[2:]))))
    except IOError:
        pass

    return gene_sets

def gene_set_enrichment(scores, gene_set_collection='', resources_dir='resources', operator=operator.gt, cutoff=0.01, count=3):
    """
    Applies gene set enrichment on the data
    """
    genes = set(scores.keys())
    gene_sets = collect_gene_sets(gene_set_collection, resources_dir)
    results = []
    for n, gene_set in enumerate(gene_sets):
        overlap = list(set.intersection(genes, gene_set['genes']))
        if len(overlap) >= count:
            overlap_scores = dict((str(g), scores[str(g)]) for g in overlap)
            p_value = gene_set_p_value(scores.values(), overlap_scores.values(), operator)
            if p_value < cutoff:
                results.append(dict(name = gene_set['name'],
                                    genes = overlap_scores.keys(),
                                    size = gene_set['size'],
                                    p_value = p_value))

    return results
