"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
from os import path, mkdir
import seten


def get_filename(filepath):
    """
    Obtain filename from the filepath
    """
    return path.splitext(path.basename(filepath))[0]


def get_resources_dir():
    """
    Obtain path to resources directory in package directory
    """
    return path.join(path.dirname(seten.__file__), 'resources')


def output_results(results, data_path, out_dir, coll_id, enr_method, pval_cutoff):
    """
    Outputs gene set enrichment results with their gene counts and p-values
    """
    output = path.join(
        out_dir,
        get_filename(data_path) + '-' + coll_id + '-' + enr_method + '.tsv'
        )
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
                'percent'
            ]
            if enr_method == 'gse':
                header.append('gse_pvalue')
            elif enr_method == 'fe':
                header.extend(['fe_pvalue', 'fe_pvalue_corr'])
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
                    ]
                    if enr_method == 'gse':
                        row.append(result['gse_pvalue'])
                        is_significant = result['gse_pvalue'] < pval_cutoff
                    elif enr_method == 'fe':
                        row.extend([result['fe_pvalue'], result['fe_pvalue_corr']])
                        is_significant = result['fe_pvalue_corr'] < pval_cutoff
                    if is_significant:
                        f.write('\t'.join(map(str, row)) + '\n')
