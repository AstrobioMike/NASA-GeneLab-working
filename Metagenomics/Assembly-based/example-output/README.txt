##################################################################################
## This directory holds processed data for NASA GLDS-286 (example)              ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-286/                ##
##                                                                              ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
## Based on GL-DPPD-XXXX for Illumina metagenomics data (assembly-based)        ##
## Annotated code is in "processing_info.tar"                                   ##
##################################################################################

Summary of contents:

    - README.txt                                - this file

    - FastQC_Outputs/                           - multiQC summary reports of FastQC runs

    - Filtered_Reads/                           - quality-filtered sequence fastq files
        - *.fastq.gz                            - quality-filtered reads
        - *.log                           - log file of stdout/stderr from bbduk

    - Assemblies/                               - assembly files
        - *-assembly.fasta                      - fasta files of individual sample assemblies
        - assembly-summaries.tsv                - table of all assemblies' summary statistics
        - *.log                                 - log files of assembly runs
        
    - Predicted_Genes/                          - per-sample predicted gene files
        - *.faa                                 - gene amino-acid sequences
        - *.fasta                               - gene nucleotide sequences
        - *.gff                                 - predicted genes in general feature format
    
    - Annotations_and_Taxonomy/                 - Kegg Orthology (KO) annotations, taxonomy, and coverages 
        - *-gene-coverage-annotation-tax.tsv    - tables with gene coverage, annotation, and taxonomy info
        - *-contig-coverage-and-tax.tsv         - tables with contig coverage and taxonomy info
    
    - Mapping_Files/                            - bam and coverage files
        - *.bam                                 - bam files
        - *.log                                 - mapping log files

    - Combined_Outputs/                         - summary outputs with all samples combined
        - GLDS-286-KO-function-coverages.tsv    - table of combined, normalized KO coverages
        - GLDS-286-gene-taxonomy-coverages.tsv  - table of combined, normalized taxonomy coverages

    - processing_info.tar                       - a tar'd directory holding processing code and files
        - Snakefile                             - Snakemake workflow file
        - unique-sample-IDs.txt                 - single-column file of unique sample identifiers
