# This script download GenBank sequences using a set of accession numbers read from a txt file.
  # Please cite www.github.com/melisaolave/

######## Requiered ############################################################################
  # libraries: ape
  # a txt file with accession numbers separeted by line breaks (\n)

######## Settings (check and change) ##########################################################
wd <- "Users/my_path"                   # change with working directory path
accessionList <- "my_accessions.txt";   # change for txt file with accession name list
format <- "fasta"                       # change for desired output format
outputName <- "GenBank-"                #change for prefix in output matrix

######## don't change anything from here ######################################################
library(ape);

setwd(wd)
accessions <- scan(accessionList, what="character", sep="\n");
seq <- read.GenBank(access.nb = accessions, seq.names = accessions, species.names = TRUE,
                    gene.names = T, as.character = FALSE);
species <- attr(seq, "species"); # get species names
names(seq) <- paste(species, accessions, sep="_");
gene <- attr(seq, "gene"); # get gene names
write.dna(seq, paste(outputName, accessionList, sep=""), format = "fasta");
table <- data.frame(species, accessions, gene);
colnames(table) <- c("species", "accessions", "gene");
if(any(list.files() %in% "GenBank-SummaryTable.txt")){
  old.table <- read.table("GenBank-SummaryTable.txt", sep="\t", header=T, quote=NULL);
  write.table(rbind(old.table, table), "GenBank-SummaryTable.txt", sep="\t", quote=F, row.names=F);
}else{
  write.table(table, "GenBank-SummaryTable.txt", sep="\t", quote=F, row.names=F);
}

