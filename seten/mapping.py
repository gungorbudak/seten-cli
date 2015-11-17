"""
This file is part of Seten which is released under the MIT License (MIT).
See file LICENSE or go to https://github.com/gungorbudak/seten/blob/master/LICENSE
for full license details.
"""
from os import path, mkdir
import requests
import csv
from intervaltree import IntervalTree


def _path(version='grch37', dataset='hsapiens_gene_ensembl', resources_dir='resources'):
    """
    Construct path to mapping
    """
    return path.join(resources_dir, '_'.join([version, dataset, 'mapping.csv']))

def _download(version='grch37', dataset='hsapiens_gene_ensembl', resources_dir='resources'):
    """
    Downloads mapping of chromosomal locations to HGNC symbols
    """
    # TODO Handle chromosome names for different species
    url_template = \
        '''http://%s.ensembl.org/biomart/martservice?query=''' \
        '''<?xml version="1.0" encoding="UTF-8" ?>''' \
        '''<!DOCTYPE Query>''' \
        '''<Query virtualSchemaName="default" formatter="CSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">''' \
        '''<Dataset name="%s" interface="default">''' \
        '''<Filter name="chromosome_name" value="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"/>''' \
        '''<Filter name="with_hgnc" excluded="0"/>''' \
        '''<Attribute name="chromosome_name"/>''' \
        '''<Attribute name="start_position"/>''' \
        '''<Attribute name="end_position"/>''' \
        '''<Attribute name="hgnc_symbol"/>''' \
        '''</Dataset>''' \
        '''</Query>'''
    url = url_template % (version, dataset)

    request = requests.get(url, stream=True)

    if not path.exists(resources_dir):
        mkdir(resources_dir)

    mapping_path = _path(version, dataset, resources_dir)
    with open(mapping_path, 'w') as f:
        f.write(','.join(['chromosome_name', 'start_position', 'end_position', 'hgnc_symbol']))
        f.write('\n')
        for line in request.iter_lines():
            f.write(line.strip())
            f.write('\n')
    return mapping_path

def generate(version='grch37', dataset='hsapiens_gene_ensembl', resources_dir='resources'):
    """
    Generates an interval tree of mapping from downloaded file
    """
    mapping_path = _path(version, dataset, resources_dir)
    if not path.exists(mapping_path):
        mapping_path = _download(version, dataset, resources_dir)
    # Collect chromosomes in a dictionary as an interval tree
    mapping = {}
    with open(mapping_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # have we added that chromosome name before?
            if not mapping.has_key(row['chromosome_name']):
                mapping[row['chromosome_name']] = IntervalTree()
            mapping[row['chromosome_name']].addi(int(row['start_position']), int(row['end_position']), row['hgnc_symbol'])
    return mapping

def search(chromosome_name, start_position, end_position, mapping=None):
    """
    Searches for chromosomal locations on the interval tree
    """
    # TODO handle chromosome names better to cover more options
    chromosome_name = str(chromosome_name).replace('chr', '') if str(chromosome_name).startswith('chr') else str(chromosome_name)
    # fixes the name incorrectly given as M instead of MT
    chromosome_name = 'MT' if chromosome_name == 'M' else chromosome_name
    # TODO better handle this exception
    try:
        return [gene.data for gene in mapping[chromosome_name].search(int(start_position), int(end_position))]
    except:
        return []
