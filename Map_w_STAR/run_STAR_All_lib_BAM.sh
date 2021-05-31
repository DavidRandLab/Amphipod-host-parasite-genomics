#!/bin/bash

#SBATCH -J Run_STAR_for_gene_prediction
#SBATCH -n 
#SBATCH --mem=
#SBATCH -t 80:00:00

######################
# Load Oscar Modules #
######################

module load star/2.6.1b 
module load samtools/1.9
 
##########################
# Declare User variables #
##########################

CPU=

STAR_index=~/data/O_grillus_genomic_resources/Genomes/STAR_index

ProjectName=AllRNA_libs

###############
# Run Command #
###############

F=$ProjectName.F.fq.gz
R=$ProjectName.R.fq.gz

echo "Running STAR"
#STAR portion
STAR --runThreadN $CPU \
--genomeDir $STAR_index \
--readFilesIn $F $R  \
--readFilesCommand gunzip -c \
--outFileNamePrefix $ProjectName.

#samtools portion
samtools view -q $QUAL -F 0x0004 --threads $CPU -b $ProjectName.sam > flt.$ProjectName.bam

# Remove non-mapping reads, fix paired SAM flags added by BWA
samtools flagstat --threads $CPU flt.$ProjectName.bam > flagstats_filter_$ProjectName.txt

# Sort with picard(DNA/RNA)
java -Xmx$JAVAMEM -jar $PICARD SortSam I=flt.$ProjectName.bam O=flt.srt.$ProjectName.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

#flagstats
samtools flagstat --threads $CPU flt.srt.$ProjectName.bam > flagstats_final_$ProjectName.txt

#index
samtools index flt.srt.$ProjectName.bam

#Quality Check
qualimap bamqc -bam flt.srt.$ProjectName.bam  -outdir ./Qualimap_$ProjectName --java-mem-size=$JAVAMEM

rm $ProjectName.sam 
rm $ProjectName.bam
rm flt.$ProjectName.bam

echo "Done!!"