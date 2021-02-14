conda create -n pgrp-grape
conda activate pgrp-grape
conda install -c bioconda -c conda-forge snakemake -y
# conda install -c bioconda agat -y
conda install -c bioconda fastqc -y
conda install -c bioconda multiqc -y
conda install -c bioconda trimmomatic -y
conda install -c bioconda bwa -y
conda install -c bioconda samtools -y
conda install -c bioconda subread -y
conda install -c bioconda gffread -y 