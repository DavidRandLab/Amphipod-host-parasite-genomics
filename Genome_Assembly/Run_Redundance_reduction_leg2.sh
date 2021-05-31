#!/bin/bash

#SBATCH -J redundans_amphy_Leg2
#SBATCH -n 12
#SBATCH --mem=60GB
#SBATCH -t 60:00:00

module load redundans/1.0
module load python/2.7.12

F=~/data/shared_folder/For_AS_Spierer/amphi_2_leg_q20_l5_r140.F.SGAfiltered.fastq
R=~/data/shared_folder/For_AS_Spierer/amphi_2_leg_q20_l5_r140.R.SGAfiltered.fastq
CONTIGS=~/data/shared_folder/For_AS_Spierer/amphi_genome_Leg2.51.12.fasta
CPU=12
MQ=20

redundans.py -v -i $F $R -f $CONTIGS -t $CPU -q $MQ

