# Repo: Parasite manipulation of host phenotypes inferred from transcriptional analyses in a trematode-amphipod system

The files in this GitHub site support the analyses in the manuscript “Parasite manipulation of host phenotypes inferred from transcriptional analyses in a trematode-amphipod system”.

The genomic data consist of an assembly of the amphipod (Orchestia grillus) nuclear genome from paired end DNA Illumina reads. The transcriptome data consist of 48 individual paired end Illumina RNAseq reads from 24 amphipods infected with the trematode Levinseniella byrdi and 24 uninfected individuals. 

The data for this paper can be downloaded from NCBI using the following accession IDs:

The RNA samples are in BioProject: PRJNA557566; SRAs: SRR9870971-SRR9870986

The genome is in BioProject: PRJNA557538; GeneBank: VOMN00000000

The genome assembly code is in the Genome_Assembly folder. 
The read mapping code is in the Map_w_STAR folder. 
Gene models for the transcriptional analysis were developed from homology to Drosophila melanogaster to take advantage of gene ontology (GO) annotations. 
The gene models are contained in the gene_models folder. 

The differential expression analyses were conducted with three different R packages: DEseq2, EBseq and edgeR. The R scripts for each of those analyses are in their respective folders. Each script uses the same read-count file and the same master annotation file specifying the gene models: 
Ogril1_clean_count.txt 
Ogril1_RNAevidence_MASTER_ANNOTATION_FILE.txt 
These files are located in the main GitHub folder and are called from each R script as the URL for the raw text contained in each file. 

The Functions folder contains code for several functions that are called from other scripts in other directories. 

Older versions of these scripts are in the DE_Analysis folder.

The Figures folder contains code for generating some of the figures; other figures are generated in the individual R scripts for differential expression.
