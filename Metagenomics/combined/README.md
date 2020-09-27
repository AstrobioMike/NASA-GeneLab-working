
# GeneLab bioinformatics processing protocol for Illumina metagenomics data – in review

**Date:** September 26, 2020  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

**Submitted by:**  
Michael D. Lee

---

# Details for reviewing

* The primary GeneLab protocol document for review, holding main steps/programs and example code, is ['GL-DPPD-XXXX.md' here](GL-DPPD-XXXX.md). 

* An example was done with 2 samples from [GLDS-286](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-286/). 
  * The ['example-output' directory](example-output) in this repo holds the output files for your perusing pleasure, or the whole thing can be downloaded to your computer as a zip file by [clicking here **not present yet**]() (it is ~63MB compressed and ~243MB uncompressed).

  * The processing is implemented as a Snakemake workflow. The primary snakefile is [here](example-output/processing_info/Snakefile).

  * The `README.txt` in the ['example-output' directory](example-output) describes the output contents and can be found [here](example-output/README.txt), and is also pasted below.
  
<br>
Thanks! 🙂


---

Contents of ['example-output' directory](example-output):

```
##################################################################################
## This directory holds processed data for NASA GLDS-286 (example)              ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-286/                ##
##                                                                              ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
## Based on GL-DPPD-XXXX for Illumina metagenomics data (assembly-based)        ##
## Annotated code is in "processing_info/"                                      ##
##################################################################################

Summary of contents:

    - README.txt                                    - this file

    - Raw_Data/                                     - initial fastq files

    - FastQC_Outputs/                               - multiQC summary reports of FastQC runs

    - Filtered_Sequence_Data/                       - quality-filtered sequence fastq files
        - *.fastq.gz                                - quality-filtered reads
        - *.log                                     - log file of stdout/stderr from bbduk

    - Assembly_based_processing/                    - files generated from an assembly-based approach

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
            - GLDS-286-taxonomy-coverages.tsv       - table of combined, normalized gene-level taxonomy coverages

    - Read_based_processing/                        - files generated from a read-based approach

        - FILL IN


    - processing_info                           - processing code and files
        - Snakefile                             - Snakemake workflow file
        - unique-sample-IDs.txt                 - single-column file of unique sample identifiers
```
