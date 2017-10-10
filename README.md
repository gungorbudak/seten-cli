# Seten (command line interface)

Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets. Given a peak-detected BED file or a two-column gene - score file, Seten does a gene set enrichment analysis.

Supports following organisms:

* Fruit fly (dme6 build)
* Human (hg19 build and hg38 build)
* Mouse (mm10 build)
* Rat (rn6 build)
* Worm (cel235 build)
* Yeast (r64-1-1 build)

Supports following gene set collections for the above organisms:

* Pathways: BIOCARTA*
* Pathways: KEGG
* Pathways: REACTOME
* GO: Biological Process
* GO: Molecular Function
* GO: Cellular Compartment
* Human Phenotype Ontology**
* MalaCards Disease Ontology**

\* Only available for human and mouse \** Only available for human

If the organism or gene set collection of your interest is not listed here, you can extend Seten for them. Please go to [the guide for extending Seten to new organisms and/or gene set collections](#extending-to-new-organisms-andor-gene-set-collections).

## Installation

### Using package manager

    pip install git+https://github.com/gungorbudak/seten-cli.git

### Using source code

    git clone https://github.com/gungorbudak/seten-cli.git
    cd seten-cli/
    python setup.py install

### Locally

If you don't have superuser rights to install Seten, you can add `--user` to the either installation command to make it available locally. In this case, depending on your system you might have to add the local binaries directory to the PATH environment variable (e.g. `/home/user/.local/bin` in Linux).

## Usage

Use below command and its help option (`-h`or `--help`) to get usage of Seten.

    $ seten -h

### Extending to new organisms and/or gene set collections

#### Generating a mapping file for an organism

A mapping file is a Tab-separated table including columns chromosome name, gene start position, gene end position and corresponding gene name for all genes available for that organism. Such a table might be obtained from Ensembl BioMart. Please make sure for genes, you select gene name field for corresponding organism (i.e. HGNC symbols for human) and do not include header row in the tables.

#### Generating a gene set collection file (in GMT format)

GMT (Gene Matrix Transposed) format is a tab delimited file format that describes gene sets in a gene set collection. In the GMT format, each row represents a gene set; in the GMX format, each column represents a gene set. Go to [Broad Institute's official page for organizing a gene set collection in GMT format](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). Please make sure you are using correct identifier for the genes which should be corresponding given gene names. For example, if you are generating a gene set collection file for fruit fly, the correct gene identifier should be FlyBase gene name.

## Disclamer about the collections

The use of gene set collections are for academic purpose and for the use of any purposes the developers and maintainers of the collections should be contacted.
