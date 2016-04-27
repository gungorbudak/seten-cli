"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten-cli/blob/master/LICENSE
for full license details.
"""
import json
from os import path
from intervaltree import IntervalTree
from seten.utils import get_resources_dir


def generate(organism):
    """
    Generates an interval tree of mapping from downloaded file
    """
    mappingTree = {}
    # Read mapping to a list
    mapping_path = path.join(get_resources_dir(), 'mappings', organism + '.json')
    if path.exists(mapping_path):
        with open(mapping_path, 'r') as f: mapping = json.load(f)
        # Collect chromosomes in a dictionary as an interval tree
        for each in mapping:
            # have we added that chromosome name before?
            if not mappingTree.has_key(each['chrName']):
                mappingTree[each['chrName']] = IntervalTree()
            mappingTree[each['chrName']].addi(
                int(each['start']), int(each['end']), each['hgncSymbol'])
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
