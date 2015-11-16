import os
import operator
import argparse
from time import time
from seten.mapping import generate
from seten.enrichment import *
from seten.utils import filename, output_results


def main():
    # timer starts
    start_time = time()

    # parse terminal arguments
    parser = argparse.ArgumentParser(
                description='Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('data', help='can be a path to a BED file or a directory of BED files')
    parser.add_argument('collection', help='is a path to a gene set collection GMT file containing gene sets and associated genes')
    parser.add_argument('-r', default='resources', help='is a path to the directory that can store resources such as mapping and gene set collections')
    parser.add_argument('-o', default='output', help='is a path to the output directory that will store results')
    parser.add_argument('-i', default=4, type=int, help='is the zero-based index of the score column in a BED file')
    parser.add_argument('-m', default='highest', help='is the method to compute a gene level score from multiple binding scores for the same gene')
    parser.add_argument('-p', default='gt', help='relates the operator for comparing the median of overlap set and every randomly sampled sets')
    args = parser.parse_args()

    # collect paths to data files here
    data_paths = []
    # is the data a directory?
    if os.path.isdir(args.data):
        for data_file in os.listdir(args.data):
            data_paths.append(os.path.join(args.data, data_file))
    else:
        data_paths.append(args.data)

    # relate operators
    ops = {'gt': operator.gt,
           'lt': operator.lt}
    op = ops[args.p]

    # generate the mapping
    mapping = generate('grch37', 'hsapiens_gene_ensembl', args.r)

    # run for each data file
    for data_path in data_paths:
        # collect scores
        scores = collect_scores(data_path, mapping=mapping, index=args.i, method=args.m)

        # collect gene sets and gene collection size
        gene_sets, gc_size = collect_gene_sets(args.collection)

        # collect results
        results = [integrated_enrichment(scores, gene_set=gene_set, gc_size=gc_size, operator=op) for gene_set in gene_sets]
        results = [result for result in results if result != None]

        # correct p-values in results
        results_cor = correct_results(results)

        # output results
        output = os.path.join(
                        args.o,
                        filename(data_path) + '_' + filename(args.collection) + '.tsv'
                        )
        output_results(output, results_cor)

    # timer ends
    print 'Took', round(time() - start_time, 2), 'seconds'

if __name__ == '__main__':
    main()
