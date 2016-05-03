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

### Using source code

    git clone https://github.com/gungorbudak/seten-cli.git
    cd seten-cli/
    python setup.py install

### Locally

If you don't have superuser rights to install Seten, you can add `--user` to installation commands to make it available locally. In this case, depending on your system you might have to add the local binaries directory to the PATH environment variable (e.g. `/home/user/.local/bin` in Linux).

## Usage

    $ seten -h
    usage: seten [-h]
                 [--colls {biocarta,kegg,reactome,gobp,gomf,gocc,hpo,malacards} [{biocarta,kegg,reactome,gobp,gomf,gocc,hpo,malacards} ...]]
                 [--org {hsa_hg19,mmu_mm10,rno_rn6,dme_bdgp6}]
                 [--enr-mtd {gse,fe}] [--scr-mtd {min,max,mean,median,sum}]
                 [--corr-mtd {fdr,bh,by,bon}] [--pc PVAL] [--gsc NUM] [--oc NUM]
                 [--sc PVAL] [--iter NUM] [--proc NUM] [--out DIR]
                 data

    Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets

    positional arguments:
      data                  path to an input file or a directory of input files,
                            input files can be UCSC BED formatted text files, or
                            two-column gene - score pairs

    optional arguments:
      -h, --help            show this help message and exit
      --colls {biocarta,kegg,reactome,gobp,gomf,gocc,hpo,malacards} [{biocarta,kegg,reactome,gobp,gomf,gocc,hpo,malacards} ...]
                            gene set collections to do enrichment analysis on
                            (default: ['biocarta', 'kegg', 'reactome', 'gobp',
                            'gomf', 'gocc', 'hpo', 'malacards'])
      --org {hsa_hg19,mmu_mm10,rno_rn6,dme_bdgp6}
                            organism (default: hsa_hg19)
      --enr-mtd {gse,fe}    enrichment method, gene set enrichment (gse) or
                            functional enrichment (fe) using Fisher's exact test
                            (default: gse)
      --scr-mtd {min,max,mean,median,sum}
                            method to compute a gene level score from multiple
                            binding scores for the same gene (default: max)
      --corr-mtd {fdr,bh,by,bon}
                            correction method after Fisher's exact test for
                            functional enrichment, Benjamini & Hochberg (fdr or
                            bh), Benjamini & Yekutieli (by) and Bonferroni (bon)
                            (default: fdr)
      --pc PVAL             p-value cutoff for significant gene set enrichment or
                            corrected functional enrichment results (default:
                            0.05)
      --gsc NUM             gene set cutoff, maximum number of genes in gene sets
                            in selected gene set collections (default: 350)
      --oc NUM              overlap cutoff, minimum number of overlapping genes
                            between the dataset and each gene set (default: 5)
      --sc PVAL             significance cutoff for significant Mann-Whitney U
                            test result in gene set enrichment iterations
                            (default: 0.05)
      --iter NUM            number of iterations for gene set enrichment analysis
                            (default: 1000)
      --proc NUM            number of processes to use for analyses (default: 4)
      --out DIR             path to the output directory for storing results
                            (default: output)

## Disclamer about the collections

REACTOME, BIOCARTA, KEGG, GO biological process, GO molecular function and GO cellular compartment collections have been obtained from Broad Institute's Molecular Signature Database version 5.0 (Subramanian, Tamayo, et al. 2005). Human phenotype ontology has been obtained from The Human Phenotype Ontology project  (Köhler et al. 2014). MalaCards disease ontology has been obtained from MalaCards: The human disease database (Rappaport et al. 2014). The use of these collections are for academic purpose and for the use of any purposes the developers and maintainers of the collections should be contacted.
