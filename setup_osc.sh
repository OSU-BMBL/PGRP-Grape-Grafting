#/bin/bash

# The bash script is used for Ohio Super Computer(OSC) environment.
module load python/3.6
conda create -n pgrp-grape
conda activate pgrp-grape
conda install -c bioconda snakemake -y
conda install -c bioconda fastqc -y
conda install -c bioconda multiqc -y
conda install -c bioconda bbmap -y
conda install -c bioconda star -y
conda install -c bioconda subread -y