import os
import csv
import operator
import argparse
from time import time
from itertools import chain
from seten.mapping import generate
from seten.enrichment import *


def output_results(output, results):
    """
    Outputs gene set enrichment results with their gene counts and p-values
    """
    # output directory exists?
    if not os.path.exists(os.path.dirname(output)):
        os.mkdir(os.path.dirname(output))

    # any result available
    if results:
        with open(output, 'w') as f:
            header = results[0].keys()
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            for result in results:
                if result:
                    writer.writerow(result)

def main():
    # timer starts
    start_time = time()

    # parse terminal arguments
    parser = argparse.ArgumentParser(
                description='Gene set enrichment on CHIP-seq RBA-binding protein binding signals datasets',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('data', help='can be a path to a BED file or a directory of BED files')
    parser.add_argument('-r', default='resources', help='path to the directory that stores/will store resources such as mapping')
    parser.add_argument('-o', default='output', help='path to the output directory that will store results')
    parser.add_argument('-i', default=4, type=int, help='index of the score column in a BED file')
    parser.add_argument('-m', default='highest', help='method to compute a gene level score from multiple binding scores for the same gene')
    parser.add_argument('-p', default='gt', help='relate operator for comparing the median of overlap set and every randomly sampled sets')
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

        # define gene set collections matching the name of the file in resources
        gene_set_collections = [
            # 'c2.cp.biocarta.v5.0.symbols',
            'c2.cp.kegg.v5.0.symbols',
            # 'c2.cp.reactome.v5.0.symbols',
            # 'c5.bp.v5.0.symbols',
            # 'c5.cc.v5.0.symbols',
            # 'c5.mf.v5.0.symbols',
            # 'cx.hpo.v5.0.symbols',
            # 'cx.malacard.v5.0.symbols'
        ]

        # collect gene sets, gene collection sizes tuples
        gene_sets = []
        for gene_set_collection in gene_set_collections:
            gc_sets = collect_gene_sets(gene_set_collection, args.r)
            gc_size = len(set(chain.from_iterable([gc_set['genes'] for gc_set in gc_sets])))
            gene_sets.append((gc_sets, gc_size))

        # collect results
        results = []
        for gc_sets, gc_size in gene_sets:
            for gene_set in gc_sets:
                results.append(integrated_enrichment(scores, gene_set=gene_set, gc_size=gc_size, operator=op))

        # output results
        output = os.path.join(args.o, os.path.basename(data_path))
        output_results(output, results)

    # timer ends
    print 'Took', round(time() - start_time, 2), 'seconds'

if __name__ == '__main__':
    main()
