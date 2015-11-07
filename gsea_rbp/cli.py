import os
import operator
import argparse
from gsea_rbp.mapping import generate
from gsea_rbp.enrichment import collect_scores
from gsea_rbp.enrichment import gene_set_enrichment


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
            for result in results:
                f.write(result['name'] + '\t' + str(len(result['genes'])) + '\t' + str(result['size']) + '\t' + str(result['p_value']) + '\n')

def main():
    # parse terminal arguments
    parser = argparse.ArgumentParser(
                description='Inferring associations between RNA-binding proteins (RBPs) and gene sets based on binding signals.',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('data', help='can be a path to a BED file or a directory of BED files')
    parser.add_argument('-r', default='resources', help='path to directory that stores/will store resources such as mapping')
    parser.add_argument('-o', default='output', help='path output directory that will store results')
    parser.add_argument('-i', default=4, type=int, help='index of the score column in a BED file')
    parser.add_argument('-m', default='highest', help='method to compute a gene level score from multiple binding scores for the same gene')
    parser.add_argument('-p', default='gt', help='relate operator for comparing the median of overlap set and randomly sampled sets')
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
    ops = {
           'gt': operator.gt,
           'lt': operator.lt,
           'ge': operator.ge,
           'le': operator.le,
           'eq': operator.eq
          }
    op = ops[args.p]

    # generate the mapping
    mapping = generate('grch37', 'hsapiens_gene_ensembl', 'resources')

    # run for each data file
    for data_path in data_paths:
        # collect scores
        scores = collect_scores(data_path, mapping=mapping, index=args.i, method=args.m)

        # define gene set collections matching the name of the file in resources
        gene_set_collections = [
            'c2.cp.biocarta.v5.0.symbols',
            'c2.cp.kegg.v5.0.symbols',
            'c2.cp.reactome.v5.0.symbols',
            'c5.bp.v5.0.symbols',
            'c5.cc.v5.0.symbols',
            'c5.mf.v5.0.symbols',
            'hpo',
            'malacard'
        ]

        # collect results
        results = []
        for gene_set_collection in gene_set_collections:
            results.extend(gene_set_enrichment(scores, gene_set_collection=gene_set_collection, resources_dir=args.r, operator=op))

        # output results
        output = os.path.join(args.o, os.path.basename(data_path))
        output_results(output, results)

if __name__ == '__main__':
    main()