#!/bin/bash

#SBATCH --output=bbduk.out
#SBATCH --job-name=bbduk
#SBATCH -p defq

# queue options = [defq, medmem, himem]

module load conda/3.6
source activate /xfs2/millerlab/tools/conda_envs/rna_env

###############################################################################

for file in /xfs2/millerlab/angelawu/2017/raw/*;

do
        filename=`echo $file | egrep -o 'A1.*' | sed 's/\.fastq.gz//g'`

        bbduk.sh in=$file out=/xfs2/millerlab/angelawu/2017/trim/trim_$filename.fq.gz ref=/xfs2/millerlab/angelawu/2017/adapters.fa qtrim=rl trimq=10

        echo $filename

done 2> 2017_bbduk.stderr
~                          
