# Seten (command line interface)

Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets. Given a peak called BED file, Seten will result in gene set enrichment analysis results for following gene collections:

* Pathways: REACTOME
* Pathways: BIOCARTA
* Pathways: KEGG
* GO: Biological Process
* GO: Molecular Function
* GO: Cellular Compartment
* Human Phenotype Ontology
* MalaCards Disease Ontology

## Installation

### Using package manager

    pip install git+https://github.com/gungorbudak/seten-cli.git

### From source

    git clone https://github.com/gungorbudak/seten-cli.git
    cd seten-cli/
    python setup.py install

## Usage

    seten -h
    usage: seten [-h] [-o O] [-i I] [-m M] [-r R] [-p P] data

    Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets

    positional arguments:
      data        can be a path to a BED file or a directory of BED files

    optional arguments:
      -h, --help  show this help message and exit
      -o O        is a path to the output directory that stores results (default:
                  output)
      -i I        is the zero-based index of the score column in a BED file
                  (default: 4)
      -m M        is the method to compute a gene level score from multiple
                  binding scores for the same gene (default: highest)
      -r R        relates the operator for comparing the median of overlap set and
                  every randomly sampled sets (default: gt)
      -p P        is the number of processes to use for analyses (default: 4)

## Disclamer about the collections

REACTOME, BIOCARTA, KEGG, GO biological process, GO molecular function and GO cellular compartment collections have been obtained from Broad Institute's Molecular Signature Database version 5.0 (Subramanian, Tamayo, et al. 2005). Human phenotype ontology has been obtained from The Human Phenotype Ontology project  (Köhler et al. 2014). MalaCards disease ontology has been obtained from MalaCards: The human disease database (Rappaport et al. 2014). The use of these collections are for academic purpose and for the use of any purposes the developers and maintainers of the collections should be contacted.
