"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import os
import operator
import argparse
from time import time
from seten.mapping import generate
from seten.enrichment import *
from seten.utils import get_filename, output_results


def main():
    # timer starts
    start_time = time()

    # parse terminal arguments
    parser = argparse.ArgumentParser(
        description='Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'data', help='can be a path to a BED file or a directory of BED files')
    parser.add_argument('-c', default='resources/collections',
                        help='is a path to the collections directory that stores GMT formatted gene collection files')
    parser.add_argument('-o', default='output',
                        help='is a path to the output directory that stores results')
    parser.add_argument('-i', default=4, type=int,
                        help='is the zero-based index of the score column in a BED file')
    parser.add_argument('-m', default='highest',
                        help='is the method to compute a gene level score from multiple binding scores for the same gene')
    parser.add_argument(
        '-r', default='gt', help='relates the operator for comparing the median of overlap set and every randomly sampled sets')
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
    op = ops[args.r]

    # generate the mapping
    mapping = generate('grch37', 'hsapiens_gene_ensembl', 'resources')

    # collect gene sets and collection size
    collections, collections_size = collect_collections(args.c)
    print '[#]', collections_size, 'unique genes found in collections'

    # run for each data file
    for data_path in data_paths:
        # collect scores
        scores = collect_scores(
            data_path, mapping=mapping, index=args.i, method=args.m)
        print '[#]', len(scores.keys()), 'unique genes found in', data_path

        # collect results
        for name, collection in collections.iteritems():
            results = integrated_enrichment(
                scores,
                collection=collection,
                collections_size=collections_size,
                operator=op)
            print '[#]', len(results), 'gene set enrichment analysis results from', name

            # output results
            output = os.path.join(
                args.o,
                get_filename(data_path) + '-' + name + '.tsv'
            )
            output_results(output, results)

    # timer ends
    print '[#] Took', round(time() - start_time, 2), 'seconds'

if __name__ == '__main__':
    main()
