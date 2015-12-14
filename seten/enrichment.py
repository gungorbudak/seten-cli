"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import os
import operator
import multiprocessing
from itertools import chain
from collections import defaultdict
from seten.statistics import compute_gene_level_score
from seten.statistics import compute_gse_pvalue
from seten.statistics import compute_fe_pvalue
from seten.statistics import correct_and_combine_pvalues
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
    # compute a gene level score from list of binding scores
    # belonging to the same gene
    scores = dict((g, compute_gene_level_score(s, method=method))
                  for g, s in scores.iteritems())
    return scores


def collect_collections(collections_dir='resources/collections'):
    """
    Collects collections as a dict
    """
    collections = defaultdict(list)
    try:
        for collection in os.listdir(collections_dir):
            # proceed if it's a GMT file
            if collection.endswith('.gmt'):
                col_path = os.path.join(collections_dir, collection)
                # get collection ID from the file name
                # e.g. c2.cp.{kegg}.v5.0.symbols.gmt
                col_name = os.path.basename(col_path).split('.')[-5]
                with open(col_path, 'r') as rows:
                    for row in rows:
                        cols = row.strip().split('\t')
                        collections[col_name].append(dict(
                                name=cols[0],
                                genes=list(set(cols[2:])),
                                size=len(set(cols[2:]))
                                ))
    except IOError:
        raise Exception('cannot find collections')
    # count the number of unique genes in all collections
    collections_size = len(set(chain.from_iterable(
        [_set['genes'] for sets in collections.values() for _set in sets]
        )))
    return collections, collections_size


def _enrichment_worker(job):
    # some local variables for job
    scores = job['scores']
    gene_set = job['gene_set']
    collections_size = job['collections_size']
    operator = job['operator']
    cutoff = job['cutoff']
    count = job['count']
    # obtain overlapping genes between the data and the gene set
    overlap = _overlapping_genes(scores.keys(), gene_set['genes'])
    if len(overlap) >= count:
        # obtain scores of overlapping genes for testing
        overlap_scores = [scores[str(g)] for g in overlap]
        # test for gene set enrichment
        gse_pvalue = compute_gse_pvalue(
            scores.values(), overlap_scores, operator)
        # test for functional enrichment
        fe_pvalue = compute_fe_pvalue(
            len(overlap), gene_set['size'],
            len(scores.keys()), collections_size)
        # add to the results
        return dict(
            name=gene_set['name'],
            genes=overlap,
            overlap_size=len(overlap),
            gene_set_size=gene_set['size'],
            percent=(len(overlap)/float(gene_set['size']))*100,
            gse_pvalue=gse_pvalue,
            fe_pvalue=fe_pvalue
            )
    return None


def integrated_enrichment(scores, collection, collections_size,
        operator=operator.gt, cutoff=0.05, count=5, processes=4):
    """
    Applies integrated enrichment on the data
    """
    jobs = [
        dict(
            scores = scores,
            gene_set = gene_set,
            collections_size = collections_size,
            operator = operator,
            cutoff = cutoff,
            count = count
            )
        for gene_set in collection]
    # limit CPU count to the machines CPU count if a higher count given
    cpu_count = min(processes, multiprocessing.cpu_count())
    # pool of processes for parallel execution
    pool = multiprocessing.Pool(cpu_count)
    # submit jobs to worker
    results = pool.map(_enrichment_worker, jobs)
    # close not to cause high memory use and
    # join to wait for collecting results
    pool.close()
    pool.join()
    # remove None results
    results = [result for result in results if result != None]
    # correct for multiple p-values from functional enrichment
    # and combined corrected functional enrichment p-values
    # with gene set enrichment p-values
    return correct_and_combine_pvalues(results)
