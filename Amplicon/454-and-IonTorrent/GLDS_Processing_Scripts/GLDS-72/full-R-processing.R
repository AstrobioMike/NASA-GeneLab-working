##################################################################################
## R processing script for 16S data for GLDS-72                                 ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/                 ##
##                                                                              ##
## This code as written expects to be run within the Processing_Info/ directory ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
##################################################################################

# this is called by the corresponding Snakefile

  # loading libraries
library(DECIPHER); packageVersion("DECIPHER") # 2.12.0
library(biomformat); packageVersion("biomformat") # 1.12.0


#write.table(count_summary_tab, "../Final_Outputs/read-count-tracking.tsv", sep = "\t", quote=F, row.names=F)

    ### assigning taxonomy ###
  # reading OTUs into a DNAStringSet object
dna <- readDNAStringSet("../Final_Outputs/OTUs.fasta")

  # downloading reference R taxonomy object
download.file("http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", "SILVA_SSU_r138_2019.RData")
  # loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")
  # removing downloaded file
file.remove("SILVA_SSU_r138_2019.RData")

# assignig taxonomy
tax_info <- IdTaxa(dna, trainingSet, strand="both", processors=NULL)

  # making and writing out a taxonomy table:
    # creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(tax_tab) <- ranks
row.names(tax_tab) <- NULL
otu_ids <- names(tax_info)
tax_tab <- data.frame("OTU_ID"=otu_ids, tax_tab, check.names=FALSE)

write.table(tax_tab, "../Final_Outputs/taxonomy.tsv", sep = "\t", quote=F, row.names=FALSE)

    # reading in counts table to generate other outputs
otu_tab <- read.table("../Final_Outputs/counts.tsv", sep="\t", header=TRUE, check.names=FALSE)

    # generating and writing out biom file format
biom_object <- make_biom(data=otu_tab, observation_metadata=tax_tab)
write_biom(biom_object, "../Final_Outputs/taxonomy-and-counts.biom")

    # making a tsv of combined tax and counts
tax_and_count_tab <- merge(tax_tab, otu_tab)
write.table(tax_and_count_tab, "../Final_Outputs/taxonomy-and-counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# making final count summary table
cutadapt_tab <- read.table("../Trimmed_Sequence_Data/trimmed-read-counts.tsv", sep="\t", header=TRUE)
bbduk_tab <- read.table("../Filtered_Sequence_Data/filtered-read-counts.tsv", sep="\t", header=TRUE)[,c(1,3)]
otu_tab <- read.table("../Final_Outputs/counts.tsv", sep="\t", header=TRUE, check.names=FALSE, row.names=1)
mapped_sums <- colSums(otu_tab)
mapped_tab <- data.frame(sample=names(mapped_sums), mapped_to_OTUs=mapped_sums, row.names=NULL)

t1 <- merge(cutadapt_tab, bbduk_tab)
count_summary_tab <- merge(t1, mapped_tab)

write.table(count_summary_tab, "../Final_Outputs/read-count-tracking.tsv", sep="\t", quote=FALSE, row.names=FALSE)
