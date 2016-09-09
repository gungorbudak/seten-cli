"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import os
import json
import operator
import multiprocessing
from itertools import chain
from collections import defaultdict
from seten.statistics import compute_gene_level_score
from seten.statistics import compute_gse_pvalue
from seten.statistics import compute_fe_pvalue
from seten.statistics import correct_pvalues
from seten.mapping import search
from seten.utils import get_resources_dir


def _read_gmt_file(coll_file):
    coll = {
        'collectionId': 'given',
        'collection': 'Given Collection',
        'geneSets': [],
        'organismId': 'NA',
        'organism': 'NA',
        'identifier': 'geneName',
        'size': 0,
    }

    with open(coll_file) as rows:
        for row in rows:
            cols = row.strip().split('\t')
            genes = cols[2:]
            coll['geneSets'].append({
                'id': 'NA',
                'name': cols[0],
                'genes': genes,
                'size': len(genes),
            })
            coll['size'] += 1

    return coll


def collect_collections(organism, organism_file, given_colls, coll_file, all_colls):
    """
    Collects all available collections as a list
    and returns together with number of unique genes in all
    """
    colls = []
    genes = []

    if organism_file is None:
        colls_dir = os.path.join(
            get_resources_dir(),
            'collections',
            organism.split('_')[0]
            )
        for coll in all_colls:
            coll_path = os.path.join(colls_dir, coll + '.json')
            if os.path.exists(coll_path):
                with open(coll_path, 'r') as f: _coll = json.load(f)
                genes.extend([gene_set['genes'] for gene_set in _coll['geneSets']])
                if coll in given_colls:
                    colls.append(_coll)

    if coll_file is not None:
        _coll = _read_gmt_file(coll_file)
        genes.extend([gene_set['genes'] for gene_set in _coll['geneSets']])
        colls.append(_coll)

    # count the number of unique genes in all collections
    colls_size = len(set(chain.from_iterable(genes)))
    return colls, colls_size


def collect_scores(path, mapping, scr_method='max', index=4):
    """
    Collect scores from a BED or two-column gene - score file
    """
    scores = defaultdict(list)
    with open(path, 'r') as rows:
        for row in rows:
            if row.startswith(('#', 'browser', 'track')):
                continue
            cols = row.strip().split()
            # is it two-column gene - score file?
            if len(cols) == 2:
                scores[cols[0]].append(float(cols[1]))
            else:
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
    scores = dict((g, compute_gene_level_score(s, method=scr_method))
                  for g, s in scores.iteritems())

    return scores


def _overlapping_genes(genes, gs_genes):
    """
    Returns overlapping genes in a gene set
    """
    return list(set.intersection(set(genes), gs_genes))


def _gene_set_enrichment(job):

    # some local variables for job
    scores = job['scores']
    gene_set = job['gene_set']
    overlap_cutoff = job['overlap_cutoff']
    sign_cutoff = job['significance_cutoff']
    iters = job['iters']
    # obtain overlapping genes between the data and the gene set
    overlap = _overlapping_genes(scores.keys(), gene_set['genes'])

    if overlap_cutoff <= len(overlap):
        # obtain scores of overlapping genes for testing
        overlap_scores = [scores[str(g)] for g in overlap]
        # test for gene set enrichment
        gse_pvalue = compute_gse_pvalue(scores.values(),
            overlap_scores, sign_cutoff, iters)
        # add to the results
        return dict(
            id = gene_set['id'],
            name=gene_set['name'],
            genes=overlap,
            overlap_size=len(overlap),
            gene_set_size=gene_set['size'],
            percent=(len(overlap)/float(gene_set['size']))*100,
            gse_pvalue=gse_pvalue
            )

    return None


def _functional_enrichment(scores, gene_sets, colls_size,
        gene_set_cutoff=350, overlap_cutoff=5):

    results = []

    for gene_set in gene_sets:
        overlap = _overlapping_genes(scores.keys(), gene_set['genes'])
        if overlap_cutoff <= len(overlap) <= gene_set['size'] <= gene_set_cutoff:
            # test for functional enrichment
            fe_pvalue = compute_fe_pvalue(
                len(overlap), gene_set['size'],
                len(scores.keys()), colls_size)
            results.append(dict(
                id = gene_set['id'],
                name=gene_set['name'],
                genes=overlap,
                overlap_size=len(overlap),
                gene_set_size=gene_set['size'],
                percent=(len(overlap)/float(gene_set['size']))*100,
                fe_pvalue=fe_pvalue
                ))

    return results


def enrichment_handler(scores, collection, collections_size, gene_set_cutoff=350,
        overlap_cutoff=5, significance_cutoff=0.05, iters=1000, corr_method='fdr',
        enr_method='gse', processes=4):
    """
    Applies gene set enrichment or functional enrichment on the data
    """

    if enr_method == 'gse':
        jobs = [
            dict(
                scores = scores,
                gene_set = gene_set,
                overlap_cutoff = overlap_cutoff,
                significance_cutoff = significance_cutoff,
                iters = iters
                )
            for gene_set in collection['geneSets']
            if gene_set['size'] <= gene_set_cutoff
            ]
        # limit CPU count to the machines CPU count if a higher count given
        cpu_count = min(processes, multiprocessing.cpu_count())
        # pool of processes for parallel execution
        pool = multiprocessing.Pool(cpu_count)
        # submit jobs to worker
        results = pool.map(_gene_set_enrichment, jobs)
        # close not to cause high memory use
        pool.close()
        # join to wait for collecting results
        pool.join()
        # remove None results
        results = [result for result in results if result != None]

    elif enr_method == 'fe':
        results = _functional_enrichment(scores, collection['geneSets'],
            collections_size, gene_set_cutoff, overlap_cutoff)
        results = correct_pvalues(results, method=corr_method)

    return results
