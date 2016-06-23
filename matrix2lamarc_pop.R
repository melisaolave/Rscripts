# Please cite GitHub/melisaolave/Rscripts
# This script will convert matrices (format: just indiviualNames tab sequences) to Lamarc population format:
  # ind1  ATGCTGC
  # ind2  GTCGTCC
  # ind3  TTTTGGC

# Required:
  # You will also need an Imap File with the species/populations and individuals references, exactly as required by BEAST program

#################### Settings (please check) ####################################

wd <- "Users/myfolder" # change path
pattern <- ".nxs" # change using an id pattern that is shared for all your matrix files (and is not in any other file in the folder), for example ".fas"
ImapFile <- "Imap.txt" # change for Imap name
outputName <- "2lamarc.txt" # name added to output

################ Running script ##################################

setwd(wd);
dir <- getwd();
fileNames <- list.files(path = dir, pattern = pattern, all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE); # except for pattern, other commands are defaults
cat("A total of", length(fileNames), pattern, "files were found in", dir, ":\n", fileNames,"\n\t Ready to read Imap\n");

Imap  <- read.table(ImapFile, header = T, sep="\t");
sp <- as.character(unique(Imap$species));

for (i in fileNames){
  matrix <- read.table(i, sep="\t", header=F, col.names=c("ind", "seq"));
  vseq.final <- character();
  vind.final <- character();
  vsp.final <- character();
  for (j in sp){
    subset <- Imap[Imap$species == j,];
    ind <- as.character(subset$traits);
    vseq <- character();
    vind <- character();
    vsp <- character();
    for(k in ind){
      submatrix <- matrix[as.character(matrix$ind) == k,];
      if(nrow(submatrix) >= 1){
        seq <- as.character(submatrix$seq);
        vseq <- c(vseq, seq);
        new.name <- paste(unlist(strsplit(k, split=""))[1:9], collapse="");
        vind <- c(vind, new.name);
        vsp <- c(vsp, j);
      }
    }
    vseq.final <- c(vseq.final, j, vseq);
    vind.final <- c(vind.final, as.character(length(vind)), vind);
    vsp.final <- c(vsp.final, unique(vsp));
  }
  seq.length <- character();
  seq.length <- as.character(length(unlist(strsplit(seq, split=""))));
  vind.final <- c(seq.length, vind.final);
  vseq.final <- c("", vseq.final);
  
  new.matrix <- data.frame(vind.final, vseq.final);
  colnames(new.matrix) <- c(length(unique(vsp.final)), "1");
  write.table(new.matrix, paste(i, outputName, sep = ""), sep=" ", row.names=F, col.names=T, quote=F);
}




