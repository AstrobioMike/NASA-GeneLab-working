#!/bin/bash

########## GENERAL INFO ##########
# Example of GL filtering of human reads on a paired-end dataset
# Put together by Homer Fogle and Michael D. Lee (Mike.Lee@nasa.gov)

# Initial database construction done with kraken2 v2.1.1 on 29-Nov-2020
  # details of building can be found on the corresponding README.md, or on this page: https://hackmd.io/@astrobiomike/GL-kraken2-human-db-setup
  # that database build can be downloaded and unpacked as follows:
# curl -L -o kraken2-human-db.tar.gz https://ndownloader.figshare.com/files/25627058
# tar -xzvf kraken2-human-db.tar.gz

# tiny example reads with 1 human read, 1 phiX read, and 1 E. coli read can be downloaded as follows for testing:
# curl -L -o sample-A-R1.fq.gz https://ndownloader.figshare.com/files/25627406
# curl -L -o sample-A-R2.fq.gz https://ndownloader.figshare.com/files/25627406

# As written, this will grab all files in the 'input_fastq_dir' that end with the suffixes specified in the below variables.
#############################################

########## VARIABLES TO ADJUST ##########

# directory holding the kraken2 human database
kraken2_db_dir="/path/to/kraken2-human-db"

# directory holding input fastq files
input_fastq_dir="/path/to/input-fastq-files"
# directory for output files
outdir="/path/to/output-directory"

# R1 and R2 suffixes following unique portions of sample names (as written, this script finds all "samples" based on them existing in the 'input_fastq_dir' with the following provided 'R1_suffix')
R1_suffix="-R1.fastq.gz"
R2_suffix="-R2.fastq.gz"

# number of threads to pass to kraken2 call
threads="10"
#############################################

# log file path
log=${outdir}/kraken_run.log

printf "kraken2 human-read removal log file - " > $log
printf "Starting human-read removal - "
date | tee -a $log


rm -rf building.tmp

for R1_file in ${input_fastq_dir}/*${R1_suffix}
do

    # getting base path of forward read input file
    base=$(basename ${R1_file})

    # getting unique portion of sample name
    sample=${base%${R1_suffix}}

    # building reverse read file path
    R2_file=${R1_file%${R1_suffix}}${R2_suffix}

    printf "\n  Working on sample: ${sample}\n" | tee -a $log

    kraken2 --db ${kraken2_db_dir} --gzip-compressed --threads ${threads} --output ${outdir}/kraken_${sample}.txt --use-names --report ${outdir}/report_${sample}.txt --unclassified-out ${outdir}/${sample}_R#.fastq --paired ${R1_file} ${R2_file} >> $log 2>&1

    # renaming outputs
    mv ${outdir}/${sample}_R_1.fastq ${outdir}/${sample}${R1_suffix%.gz}
    mv ${outdir}/${sample}_R_2.fastq ${outdir}/${sample}${R2_suffix%.gz}

    # gzipping outputs
    gzip ${outdir}/${sample}${R1_suffix%.gz}
    gzip ${outdir}/${sample}${R2_suffix%.gz}

    # building summary table
    total_fragments=$(wc -l ${outdir}/kraken_${sample}.txt | cut -f 1 -d " ")
    fragments_retained=$(grep -w -m 1 "unclassified" ${outdir}/report_${sample}.txt | cut -f 2)
    perc_removed=$(printf "%.2f\n" $(echo "scale=4; 100 - ${fragments_retained} / ${total_fragments} * 100" | bc -l))

    printf "${sample}\t${total_fragments}\t${fragments_retained}\t${perc_removed}\n" >> building.tmp

done

# adding header and combining all sample summaries into one table
cat <( printf "sample\tTotal_fragments_before\tTotal_fragments_after\tPercent_human_reads_removed\n" ) building.tmp > ${outdir}/Human-read-removal-summary.tsv && rm building.tmp

printf "\nFinished - " | tee -a $log
date | tee -a $log
