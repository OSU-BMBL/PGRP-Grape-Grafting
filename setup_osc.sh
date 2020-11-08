#/bin/bash

# The bash script is used for Ohio Super Computer(OSC) environment.
conda create -n pgrp-grape -y
source activate pgrp-grape
conda install -c bioconda snakemake -y
conda install -c bioconda fastqc -y
conda install -c bioconda bbmap -y
conda install -c bioconda star -y
conda install -c bioconda subread -y
conda install -c conda-forge lzstring -y
conda install -c conda-forge spectra -y
conda install -c bioconda multiqc -y
