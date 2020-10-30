##################################################################################
## This directory holds processed data for NASA GLDS-72                         ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/                 ##
##                                                                              ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
## Based on GL-DPPD-XXXX for 454 or IonTorrent amplicon data                    ##
## Annotated code is in "processing_info.tar"                                   ##
##################################################################################

Summary of contents:

    - README.txt                        - this file

    - FastQC_Outputs/                   - multiQC summary reports of FastQC runs

    - Trimmed_Sequence_Data/            - primer-trimmed sequences and trimmed-read count table
        - *.fastq.gz                    - primer-trimmed sequence fastq files
        - trimmed-read-counts.tsv       - starting and trimmed-read counts per sample
        - cutadapt.log                  - log file from cutadapt

    - Filtered_Sequence_Data/           - quality-filtered sequence fastq files and filtered-read count table
        - *.fastq.gz                    - quality-filtered reads
        - filtered-read-counts.tsv      - starting and filtered-read counts per sample
        - bbduk.log                     - log file from bbduk

    - Final_Outputs/                    - primary output files
        - OTUs.fasta                    - fasta file of generated OTUs
        - counts.tsv                    - count table of OTUs across samples
        - taxonomy.tsv                  - taxonomy of OTUs
        - taxonomy-and-counts.biom      - biom-formatted output of taxonomy and counts
        - taxonomy-and-counts.tsv       - combined table of counts and taxonomy across samples
        - read-count-tracking.tsv       - read counts at each processing step

    - processing_info.tar               - a compressed directory holding info related to processing
        - unique-sample-IDs.txt         - single-column file of unique sample identifiers
        - Snakefile                     - snakemake workflow file
        - full-R-processing.R           - all R processing code used
        - primers.fa                    - primer sequences
        - snakemake-run.log             - log file from workflow execution
        - R-processing.log              - log file from R processing
        - vsearch.log                   - log file from vsearch 
