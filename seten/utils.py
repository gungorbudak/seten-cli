"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
from os import path, mkdir


def get_filename(filepath):
    """
    Obtain filename from the filepath
    """
    return path.splitext(path.basename(filepath))[0]


def output_results(output, results):
    """
    Outputs gene set enrichment results with their gene counts and p-values
    """
    # output directory exists?
    if not path.exists(path.dirname(output)):
        mkdir(path.dirname(output))

    # any result available
    if results:
        with open(output, 'w') as f:
            # construct header
            header = [
                'name',
                'genes',
                'overlap_size',
                'gene_set_size',
                'percent',
                'fe_pvalue',
                'fe_pvalue_corrected',
                'gse_pvalue',
                'combined_pvalue'
            ]
            f.write('\t'.join(header) + '\n')
            for result in results:
                if result:
                    # convert gene list to comma separated string
                    result['genes'] = ', '.join(result['genes'])
                    # construct row same order as header
                    row = [
                        result['name'],
                        result['genes'],
                        result['overlap_size'],
                        result['gene_set_size'],
                        result['percent'],
                        result['fe_pvalue'],
                        result['fe_pvalue_corrected'],
                        result['gse_pvalue'],
                        result['combined_pvalue']
                    ]
                    f.write('\t'.join(map(str, row)) + '\n')
