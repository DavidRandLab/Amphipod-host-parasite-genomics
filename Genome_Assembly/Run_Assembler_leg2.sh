#!/bin/bash

#SBATCH -J sparse_amphi_Leg2
#SBATCH -n 8
#SBATCH --mem=120GB
#SBATCH -t 72:00:00

module load dbg2olc/Sep2018 
module load quast
module load boost/1.63.0

g=12
GS=3000000000
ProjectID=amphi_genome

F=~/data/shared_folder/For_AS_Spierer/amphi_2_leg_q20_l5_r140.F.SGAfiltered.fastq
R=~/data/shared_folder/For_AS_Spierer/amphi_2_leg_q20_l5_r140.R.SGAfiltered.fastq
S=~/data/shared_folder/For_AS_Spierer/amphi_2_leg_q20_l5_r140.Singletons.SGAfiltered.fastq

for K in 51 
do

SparseAssembler LD 0 k $K g $g GS $GS  i1 $F i2 $R f $S NodeCovTh 0 EdgeCovTh 1 ResolveBranchesPE 1 LinkCovTh 1
mv Contigs.txt $ProjectID.$K.$g.fasta


rm Asm_HT_content ContigsLog.txt Asm_HT_idx.txt Cov_HistContigs.txt Assembly_Log.txt CovHist.txt BasesCount.txt EdgeCovHist.txt ContigGraph.txt InsertSizeEst.txt Contigs_Cov.txt Contigs_HP.txt MergeHT_idx.txt Contigs_info.txt nLB_Contigs.txt Contigs_Len.txt nRB_Contigs.txt MergeHT_content
#done

# Run K-mer optimization analysis
quast.py $ProjectID.$K.$g.fasta



