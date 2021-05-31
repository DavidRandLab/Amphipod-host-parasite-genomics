#!/bin/bash

#SBATCH -J SGA_filtering_on_amphi_CYST
#SBATCH -n 6
#SBATCH --mem=120GB
#SBATCH -t 72:00:00


################################ USER DEFINED VARIABLES ################################

#Load Raw Reads as forward and reverse. fq.gz files are OK!
forward=~/igert/amphipods/raw_data/GENEWIZ_JULY2016/DNAseq/DR02-Cyst_S32_L005_R1_001.fastq.gz
#Reverse Reads of the PE set. These should be in FastQ format
reverse=~/igert/amphipods/raw_data/GENEWIZ_JULY2016/DNAseq/DR02-Cyst_S32_L005_R2_001.fastq.gz
#File indentifier. This is the name given to the output files
#NOTE=avoid "." at the end of the workid.
workid=amphi_CYST_q20_l5_r140

#System Variables # Number of CPUs for processing 
CPU=6

#Program Variables
#SGA minumun k-mer coverage for considering a true K-mer
COV_FILTER=2 
#Java Memory
JAVAMEM=110g

#BBduk locaiton
#bbtools=~/data/Jcbn/Software/bbmap  ##BBTools is located in OSCAR in: /users/drand/data/Jcbn/Software/bbmap
preqc=~/data/S_balanoides_genomics_resources/Misc_resources/sga_preqc_report.py

######### PIPELINE ######################################## PIPELINE ####################
#Load softwares
module load fastqc/0.11.5
module load sga/0.10.15
module load bbmap/38.23 

#Starting clean up 
sga preprocess -v --pe-mode 1 --pe-orphans $workid.orphans.fastq -q20 -m 80 -o $workid.interleaved.SGA.fastq $forward $reverse # Convert reads into an interleaved format for SGA to work with it. 

#Build an FM-Index
sga index -a ropebwt -t $CPU $workid.interleaved.SGA.fastq # Build an FM index for the data. SGA internal 

#complexity filtering via kmer filter
# For BBMap to work the file extension in the data must be traditional. E.g. ".fasta" or ".fastq". SGA output ".fa" thus we must modify the extension to ".fastq"

sga filter -x $COV_FILTER -t $CPU -o $workid.interleaved.SGA.filter.pass.fastq --homopolymer-check --low-complexity-check $workid.interleaved.SGA.fastq # filter the data for low complexity and homopolimeric reads 

###################################
#RUN PREQC PIPELINE               #
###################################

# Do preqc on raw files
sga preqc -t $CPU $workid.interleaved.SGA.fastq > $workid.interleaved.interleaved.SGA.preqc

# Do preqc on SGA filtered files
sga index -a ropebwt --no-reverse -t 8 $workid.interleaved.SGA.filter.pass.trim.fastq
sga preqc -t $CPU $workid.interleaved.SGA.filter.pass.trim.fastq > $workid.interleaved.SGA.filter.pass.trim.preqc


#Generate report for both pre and post filtered
$preqc $workid.interleaved.SGA.filter.pass.trim.preqc $workid.interleaved.interleaved.SGA.preqc

###################################
#RUN FASTQ PIPELINE               #
###################################

# Compare read qualities post and pre filtering
fastqc -t $CPU $workid.interleaved.SGA.filter.pass.trim.fastq $workid.interleaved.SGA.fastq


#################################################
#Trim the edges of reads for polymerase effects #
#################################################

# trimming the right side of the read to 115

bbduk.sh in=$workid.interleaved.SGA.filter.pass.trim.fastq out=$workid.interleaved.SGA.filter.pass.trim.2r115.fastq  -Xmx$JAVAMEM ftr=115

fastqc -t $CPU $workid.interleaved.SGA.filter.pass.trim.2r115.fastq

###################################
#Restore PE Structure using BBMap #
###################################

#Reformat reads as PE reads
reformat.sh in=$workid.interleaved.SGA.filter.pass.trim.2r115.fastq verifypairing -Xmx$JAVAMEM

repair.sh in=$workid.interleaved.SGA.filter.pass.trim.2r115.fastq out1=$workid.F.SGAfiltered.fastq out2=$workid.R.SGAfiltered.fastq outsingle=$workid.Singletons.SGAfiltered.fastq -Xmx$JAVAMEM #Restore the pair end nature of the reads

##################################################
#Perform final QC on filtered and Unfiltered data#
##################################################

fastqc -t $CPU $workid.F.SGAfiltered.fastq $workid.R.SGAfiltered.fastq $workid.Singletons.SGAfiltered.fastq $workid.orphans.fastq

fastqc -t $CPU $workid.interleaved.SGA.filter.pass.discard.fa

############################################################
#Perform final gzipping of reads that did not passed filter#
############################################################

gzip $workid.F.SGAfiltered.fastq
gzip $workid.R.SGAfiltered.fastq
gzip $workid.Singletons.SGAfiltered.fastq

#This is crucial to ensure proper space optimizations
gzip $workid.interleaved.SGA.filter.pass.discard.fa
gzip $workid.orphans.fastq


echo "Filtering pipeline Completed"