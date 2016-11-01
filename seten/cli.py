"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import os
import argparse
from time import time
from seten.mapping import generate
from seten.enrichment import collect_collections
from seten.enrichment import collect_scores
from seten.enrichment import enrichment_handler
from seten.utils import output_results


def main():
    # parse terminal arguments
    parser = argparse.ArgumentParser(
        description='Gene set enrichment on \
                     CLIP-seq RNA-binding protein binding signals datasets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # all available gene set collections
    COLLS_CHOICES = [
        'biocarta', 'kegg', 'reactome', # pathways
        'gobp', 'gomf', 'gocc', # gene ontology
        'hpo', # phenotype
        'malacards' # disease
    ]
    # all available organisms
    ORG_CHOICES = [
        'hsa_hg19', 'hsa_hg38', # human builds
        'mmu_mm10', # mouse
        'rno_rn6', # rat
        'dme_dme3', # fruit fly
        'cel_cel235', # worm
        'sce_r6411' # yeast
    ]
    # enrichment, scoring and correction choices
    ENR_CHOICES = ['gse', 'fe']
    SCR_CHOICES = ['min', 'max', 'mean', 'median', 'sum']
    CORR_CHOICES = ['fdr', 'bh', 'by', 'bon']

    parser.add_argument(
        'data',
        help='path to an input file or a directory of input files, \
              input files can be UCSC BED formatted text files, or \
              two-column gene - score pairs'
        )
    parser.add_argument(
        '--colls',
        metavar='LIST',
        default=['kegg', 'gobp'],
        choices=COLLS_CHOICES,
        nargs='+',
        help='gene set collections to do enrichment analysis on'
    )
    parser.add_argument(
        '--coll-file',
        metavar='FILE',
        default=None,
        help='GMT-formatted gene set collection file \
              to do enrichment analysis on'
    )
    org_group = parser.add_mutually_exclusive_group()
    org_group.add_argument(
        '--org',
        metavar='ORG',
        default='hsa_hg19',
        choices=ORG_CHOICES,
        help='organism'
    )
    org_group.add_argument(
        '--org-file',
        metavar='FILE',
        default=None,
        help='organism file, Tab-separated rows of \
              chromosomal location and gene name'
    )
    parser.add_argument(
        '--enr-mtd',
        metavar='MTD',
        default='gse',
        choices=ENR_CHOICES,
        help='enrichment method, gene set enrichment (gse) or \
              functional enrichment (fe) using Fisher\'s exact test'
    )
    parser.add_argument(
        '--scr-mtd',
        metavar='MTD',
        default='max',
        choices=SCR_CHOICES,
        help='method to compute a gene level score from \
              multiple binding scores for the same gene'
    )
    parser.add_argument(
        '--corr-mtd',
        metavar='MTD',
        default='fdr',
        choices=CORR_CHOICES,
        help='correction method after Fisher\'s exact test for \
              functional enrichment, Benjamini & Hochberg (fdr or bh), \
              Benjamini & Yekutieli (by) and Bonferroni (bon)'
    )
    parser.add_argument(
        '--pc',
        metavar='PVAL',
        default=0.05,
        type=float,
        help='p-value cutoff for significant gene set enrichment \
              or corrected functional enrichment results'
    )
    parser.add_argument(
        '--gsc',
        metavar='NUM',
        default=350,
        type=int,
        help='gene set cutoff, maximum number of genes in \
              gene sets in selected gene set collections'
    )
    parser.add_argument(
        '--oc',
        metavar='NUM',
        default=5,
        type=int,
        help='overlap cutoff, minimum number of overlapping genes \
              between the dataset and each gene set'
    )
    parser.add_argument(
        '--sc',
        metavar='PVAL',
        default=0.05,
        type=float,
        help='significance cutoff for significant Mann-Whitney U test \
              result in gene set enrichment iterations'
    )
    parser.add_argument(
        '--iter',
        metavar='NUM',
        default=1000,
        type=int,
        help='number of iterations for gene set enrichment analysis'
    )
    parser.add_argument(
        '--proc',
        metavar='NUM',
        default=4,
        type=int,
        help='number of processes to use for analyses'
    )
    parser.add_argument(
        '--out',
        metavar='DIR',
        default='output',
        help='path to the output directory for storing results'
    )
    args = parser.parse_args()

    # timer starts
    start_time = time()

    # collect paths to data files here
    data_paths = []
    # is the data a directory?
    if os.path.isdir(args.data):
        for data_file in os.listdir(args.data):
            data_paths.append(os.path.join(args.data, data_file))
    else:
        data_paths.append(args.data)

    # generate the mapping
    mapping = generate(args.org, args.org_file)
    # collect gene sets and collection size
    colls, colls_size = collect_collections(
        args.org, args.org_file, args.colls, args.coll_file, COLLS_CHOICES)
    print '[#]', colls_size, 'unique genes found in gene set collections'

    # run for each data file
    for data_path in data_paths:
        # collect scores
        scores = collect_scores(
            data_path, mapping=mapping, scr_method=args.scr_mtd)
        print '[#]', len(scores.keys()), 'unique genes found in', data_path

        # start analyses for each collection
        for coll in colls:
            # collection timer
            start_collection_time = time()

            # collect results
            results = enrichment_handler(scores, coll, colls_size,
                gene_set_cutoff=args.gsc, overlap_cutoff=args.oc,
                significance_cutoff=args.sc, iters=args.iter,
                corr_method=args.corr_mtd, enr_method=args.enr_mtd,
                processes=args.proc)

            # output results
            output_results(results, data_path, args.out,
                coll['collectionId'], args.enr_mtd, args.pc)

            # collection timer ends
            print ' '.join([
                '[#] Completed in',
                str(round(time() - start_collection_time, 2)),
                'seconds'
            ])

    # timer ends
    print '[#] Took', round(time() - start_time, 2), 'seconds in total'


if __name__ == '__main__':
    main()
