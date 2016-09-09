"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import json
from os import path
from intervaltree import IntervalTree
from seten.utils import get_resources_dir


def generate(organism, organism_file):
    """
    Generates an interval tree of mapping from downloaded file
    """
    # parse either JSON file or given Tab-separated file as mapping
    mapping = []

    if organism_file is None:
        # Read mapping to a list
        mapping_path = path.join(get_resources_dir(), 'mappings', organism + '.json')
        if path.exists(mapping_path):
            # Collect chromosomes in a dictionary as an interval tree
            with open(mapping_path, 'r') as f: mapping = json.load(f)

    else:
        with open(organism_file, 'r') as rows:
            for row in rows:
                cols = row.strip().split('\t')
                if len(cols) == 4:
                    mapping.append({
                        'chrName': cols[0],
                        'start': cols[1],
                        'end': cols[2],
                        'symbol': cols[3],
                    })

    # generate mapping interval tree for fast lookup
    mappingTree = {}
    for each in mapping:
        # have we added that chromosome name before?
        if not mappingTree.has_key(each['chrName']):
            mappingTree[each['chrName']] = IntervalTree()
        mappingTree[each['chrName']].addi(
            int(each['start']), int(each['end']), each['symbol'])

    return mappingTree


def search(chr_name, start, end, mapping=None):
    """
    Searches for chromosomal locations on the interval tree
    """
    # TODO handle chromosome names better to cover more options
    chr_name = str(chr_name).replace('chr', '') if str(
        chr_name).startswith('chr') else str(chr_name)
    # fixes the name incorrectly given as M instead of MT
    chr_name = 'MT' if chr_name == 'M' else chr_name
    genes = mapping[chr_name].search(int(start), int(end))
    if genes:
        return [gene.data for gene in genes]
    return []
