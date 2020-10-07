
# GeneLab bioinformatics processing protocol for Illumina metagenomics data – in review

**Date:** September 26, 2020  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

**Submitted by:**  
Michael D. Lee

---

<p align="center">
<a href="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GL-illumiina-metagenomics-overview.pdf"><img src="https://github.com/AstrobioMike/AstrobioMike.github.io/blob/master/images/GL-illumiina-metagenomics-overview.png"></a>
</p>

--- 

# Details for reviewing
* The primary GeneLab protocol document for review, holding main steps/programs and example code, is ['GL-DPPD-XXXX-Illumina-metagenomics.md'](GL-DPPD-XXXX-Illumina-metagenomics.md). 

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
## Based on GL-DPPD-XXXX for Illumina metagenomics data                         ##
## Annotated code is in "processing_info/"                                      ##
##################################################################################

Summary of contents:

    - README.txt                                             - this file

    - processing_info                                        - processing code and files
        - unique-sample-IDs.txt                              - single-column file of unique sample identifiers
        - Snakefile                                          - Snakemake workflow file
        - snakemake-run.log                                  - Snakemake log file
        - environment.yml                                    - conda environment file

    - Raw_Data/                                              - initial fastq files

    - FastQC_Outputs/                                        - multiQC summary reports of FastQC runs

    - Filtered_Sequence_Data/                                - quality-filtered sequence fastq files
        - *.fastq.gz                                         - quality-filtered reads
        - *.log                                              - log file of stdout/stderr from bbduk

    - Assembly-based_processing/                             - results generated from an assembly-based approach

        - Assemblies/                                        - per-sample assembly files and info
            - *-assembly.fasta                               - fasta files of individual sample assemblies
            - *.log                                          - log files of assembly runs
            - assembly-summaries.tsv                         - table of all assemblies' summary statistics

        - Predicted_Genes/                                   - per-sample predicted gene files
            - *.faa                                          - gene amino-acid sequences
            - *.fasta                                        - gene nucleotide sequences
            - *.gff                                          - predicted genes in general feature format

        - Annotations_and_Taxonomy/                          - per-sample Kegg Orthology (KO) annotations, taxonomy, and coverages
            - *-gene-coverage-annotation-tax.tsv             - tables with gene coverage, annotation, and taxonomy info
            - *-contig-coverage-and-tax.tsv                  - tables with contig coverage and taxonomy info

        - Mapping_Files/                                     - per-sample bam and coverage files
            - *.bam                                          - bam files
            - *.log                                          - mapping log files

        - Combined_Outputs/                                  - summary outputs with all samples combined
            - GLDS-286-KO-function-coverages.tsv             - table of combined, normalized KO coverages
            - GLDS-286-taxonomy-coverages.tsv                - table of combined, normalized gene-level taxonomy coverages

    - Read-based_processing/                                 - results generated from a read-based approach

        - GLDS-286-gene-families.tsv                         - gene-family abundances
        - GLDS-286-gene-families-grouped-by-taxa.tsv         - gene-family abundances grouped by taxa
        - GLDS-286-gene-families-cpm.tsv                     - gene-family abundances normalized to copies-per-million
        - GLDS-286-gene-families-KO-cpm.tsv                  - KO term abundances normalized to copies-per-million
        - GLDS-286-pathway-abundances.tsv                    - pathway abundances
        - GLDS-286-pathway-abundances-grouped-by-taxa.tsv    - pathway abundances grouped by taxa
        - GLDS-286-pathway-abundances-cpm.tsv                - pathway abundances normalized to copies-per-million
        - GLDS-286-pathway-coverages.tsv                     - pathway coverages
        - GLDS-286-pathway-coverages-grouped-by-taxa.tsv     - pathway coverages grouped by taxa
        - GLDS-286-metaphlan-taxonomy.tsv                    - metaphlan estimated taxonomic relative abundances
```
