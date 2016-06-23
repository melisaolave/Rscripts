#Please cite github.com/melisaolave/Rscripts
  ## This script will download sequences form GenBank, EMBL, SWISSPROT and/or NBRF, using sp names and gene names
  ## will get: a fasta matrix per gene and a table detailing accession numbers

###################### Required ####################################################
# 1. you will need seqinr library installed
# 2. a file listing species names, separated by breaks:
    #Liolaemus boulengeri
    #Liolaemus rothi
    #Homo sapiens

##################### Settings (check) #############################################
wd <- "/Users/melisaolave/Desktop/GitHub/Rscripts/" # change path
bank <- c("genbank","embl", "swissprot", "nbrf"); # Specify which source want to check
tableName <- "table.txt"; # Provide a table with a list of species you want to check
genNames <- c("cytb", "ND2") # add as many genes you want

####################### don't change anything from here ############################
library(seqinr)
setwd(wd);
table <- scan(tableName, what="character", sep="\n")
genes <- genNames;
new.table <- table;
species <- table;
for(k in bank){
  choosebank(k);
  for(i in 1:length(genes)){
    matrices <- NULL;
    new.table <- data.frame(new.table, rep(0, length(species)));
    matrix <- character();
    accession <- character();
    sp.vec <- character();
    seq <- list();
    for(j in 1:length(species)){
      seq$nelem <- 0;
      if(genes[i] == "12S" | genes[i] == "16S"){
        try(seq <- query(genes[i], query=paste("K=", genes[i]," ribosomal RNA AND sp=", species[j], sep="")));
      }else{
        try(seq <- query(genes[i], query=paste("K=", genes[i]," AND sp=", species[j], sep="")));
      }
      if(seq$nelem == 0){
        cat("NOT found: Gene:", genes[i], "species:", species[j], "\n");
      }else{
        sequence <- getSequence(seq$req[1]);
        accession <- getName(seq)[1];
        if(genes[i] %in% matrices){
          fastaMatrix <- scan(paste(genes[i], "-", k, ".fas", sep=""), what="character", sep="\n");
          final.matrix <- paste(">", species[j], "\n", paste(sequence[[1]], collapse=""), sep="");
          final.fasta <- c(fastaMatrix, final.matrix);
          write(final.fasta, paste(genes[i], "-", k, ".fas", sep=""));
        }else{
          final.matrix <- paste(">", species[j], "\n", paste(sequence[[1]], collapse=""), sep="");
          write(final.matrix, paste(genes[i], "-", k ,".fas", sep=""));
          matrices[length(matrices)+1] <- genes[i];
        }
        colnames(new.table) <- c("species", genes[1:i]);
        new.table[new.table$species == species[j], i+1] <- accession;
        write.table(new.table, "output.txt", sep="\t", quote=F, row.names=F, col.names=T);
        cat("Gene:", genes[i], "species:", species[j], "found = ", accession,"\n");
      }
    }
  }
}

