#!/bin/bash

#SBATCH -p defq
#SBATCH --job-name=htseq
#SBATCH --output=htseq.out

module load conda/3.6
source activate /xfs2/millerlab/tools/conda_envs/rna_env

###############################################################################

for file in /xfs2/millerlab/angelawu/2017/alignment_STAR/PN40024/*.bam

do
        filename=`echo $file | egrep -o 'A1Y1_[^.]*'`

        htseq-count -s no -f bam -r pos -t gene -i ID $file /xfs2/millerlab/tools/genomes/PN40024/V3_annotation_500extend_edited_chr_geneonly.gff3 > $filename.counts

done
