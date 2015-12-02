# Seten

Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets

## Installation

### Using package manager

    pip install git+https://github.com/gungorbudak/seten-cli.git

### From source

    git clone https://github.com/gungorbudak/seten-cli.git
    cd seten-cli/
    python setup.py install

## Usage

    seten -h
    usage: seten [-h] [-r R] [-o O] [-i I] [-m M] [-p P] data collection

    Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets

    positional arguments:
      data        can be a path to a BED file or a directory of BED files
      collection  is a path to a gene set collection GMT file containing gene sets
                  and associated genes

    optional arguments:
      -h, --help  show this help message and exit
      -r R        is a path to the directory that can store resources such as
                  mapping and gene set collections (default: resources)
      -o O        is a path to the output directory that will store results
                  (default: output)
      -i I        is the index of the score column in a BED file (default: 4)
      -m M        is the method to compute a gene level score from multiple
                  binding scores for the same gene (default: highest)
      -p P        relates the operator for comparing the median of overlap set and
                  every randomly sampled sets (default: gt)
