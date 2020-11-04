#!/bin/bash

#SBATCH --output=star_test.out
#SBATCH --job-name=star_test
#SBATCH -p defq

# queue options = [defq, chem, genomic]

module load conda/3.6
source activate /xfs2/millerlab/tools/conda_envs/rna_env

###############################################################################

for file in /xfs2/millerlab/angelawu/2017/trim/*

do
        filename=`echo $file | egrep -o 'A[^.]*'`

        echo $file
        echo $filename
        STAR --genomeDir /xfs2/millerlab/tools/genomes/PN40024/STAR_indices --runThreadN 8 --readFilesIn $file --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $filename
        echo "$filename done"

done
