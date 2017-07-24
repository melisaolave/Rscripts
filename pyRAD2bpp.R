# pyRAD sequence format to bpp
  # when using please cite www.github.com/melisaolave/ (pyRAD2bpp.R)

# sequence matrix format example:
  #spA_1       TGATCTTCAATCTGCATTTT
  #spA_2       TGATCTTCAATCTGCATTTT
  #spB_1       TGATCTTCAATCTGCATTTT
  #//             **     -- - * 
  #spA_1       CGTTCTTCAATCTGCAGGGG
  #spB_1       CGTTCTTCAATCTGCAGGGG
  #spB_2       CGTTCTTCAATCTGCAGGGG

##### Please check and modify setting below ######
wd <- "/Users/mypath" # change here for the path to the dataset
myseq <- "lusc_50.loci" # change here the sequence matrix name
mypattern <- "//" # this is a pattern that indicates the end of an individual gene
sp.ind.pattern <- "_"; # symbol that separates species of individual names (e.g.: in spA_1 the species name is separated by a _ of the individual number)

##### don't change anything from here ######
setwd(wd)
matrix <- scan(myseq, what="character", sep="\n");
end <- grep(pattern=mypattern, matrix) -1;
start <- c(1, end[1:length(end)-1]+2);
new.matrix <- NULL;
taxVec <- NULL;
for(i in 1:length(start)){
  submatrix <- matrix[start[i]:end[i]];
  taxNames <- gsub(pattern="^(.*) +.*", replacement="\\1", submatrix);
  taxNames <- gsub(pattern=" ", replacement="", taxNames);
  taxVec <- c(taxVec, taxNames);
  n.ind <- length(submatrix);
  bp <- gsub(pattern="^.* +(\\w+)", replacement="\\1", submatrix[1]);
  bp <- length(unlist(strsplit(bp,"")));
  header <- paste(n.ind, bp)
  new.matrix <- c(new.matrix, header, submatrix, "\n");
  cat(i, "/", length(start), "done\n")
}
taxa <- sort(unique(taxVec));
species <- gsub(pattern=paste("(\\w+)", sp.ind.pattern, "\\w+", sep=""), replacement="\\1", taxa);
Imap <- data.frame(1:length(taxa), species);
write.table(Imap, "Imap.txt", row.names=F, col.names=F, quote=F)
for(i in 1:nrow(Imap)){
  new.matrix <- gsub(pattern=taxa[i], replacement=paste(species[i],"^",i,sep=""),new.matrix);
  cat(i, "/", nrow(Imap), "done\n");
}
write(new.matrix, "bpp.seqInput.txt");
cat("program done!\n");

taxa

