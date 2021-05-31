#!/bin/bash

#SBATCH -J Local_blast
#SBATCH -n 6
#SBATCH --mem=60GB
#SBATCH -t 60:00:00


### or

#interact -t 120:00 -n 6 -m 60Gb

##########################
# User defined variables #
##########################

getAnnoFasta=~/data/S_balanoides_genomics_resources/Misc_resources/getAnnoFast.pl
GFF=~/data/O_grillus_genomic_resources/Gene_models/Hints_iterated_fly_Augutus/Ogril1_1k_masked_RNAevidence.iteratedHints.gff3

CPU=6  #Number of cores for the analysis
#DB=~/data/LOCAL_BLAST/SWISS_PROT_2019/swissprot
DB=~/data/LOCAL_BLAST/DMELANOGASTER_2019/DMEL_GCF_000001215.4_Release_6_plus_ISO1_MT_protein.fasta
QUERY=Ogril1_1k_masked_RNAevidence.iteratedHints3.aa.fasta
OUT=Ogril1__RNAevidence_DMEL_BLAST.txt

###################
# get fasta files #
###################

#$getAnnoFasta $GFF

################################
# make the swiss prot database #
################################

#-#cd  ~/data/LOCAL_BLAST/SWISS_PROT_2019
#-#curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz > local_swissprot.tar.gz
#-#tar -xvf local_swissprot.tar.gz

###################
# run local blast #
###################

module load blast/2.7.1+ 

blastp -db $DB -evalue 1E-10 -query $QUERY -out $OUT -outfmt 6 -max_target_seqs 1  -num_threads $CPU