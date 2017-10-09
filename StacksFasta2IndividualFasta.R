## Script for converting a fasta matrix (Stacks output) to spearate alignment per locus in fasta or phylip format, and including single SNP or full sequences
# by Melisa Olave, when using please cite www.github/melisaolave
	# required: libraries foreach and doMC

#### Check settings below #####
wd <- "/Users/myworkingdirectory"; # provide working directory path
fasta.name <- "batch_1-full.fa" # name of fasta input file 
format <- "fasta"; # output format: fasta or phyilip;
single.snp <- F; # if TRUE, then only the first snp per locus is taken, if FALSE, then full sequences are included
adapter <- "TGCAG"; # usually first positions are adapters (usually 5), add here the number of bp that want to remove from fasta matrix. If don't want to remove any, then 0.
cores <- 8; # to run in parallel, specify number of cores

#### don't change anything from here #######
list.of.packages <- c("foreach", "doMC");
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])];
if(length(new.packages)){
 install.packages(new.packages);
}
library(foreach);
library(doMC);
registerDoMC(cores);

setwd(wd);
cat("working directory:", wd, "\nreading file:", fasta.name, "\n");
fasta <- scan(fasta.name, what="character", sep="\n");
taxonLines <- grep(pattern="^>", fasta);
all.loci <- gsub(pattern="^(.*_Locus_\\d+)_Allele_.*$", replacement="\\1", fasta[taxonLines]);
loci <- unique(gsub(pattern="^>(.*)_Sample_\\d+_.*", replacement="\\1", all.loci));
cat(length(loci), "loci found in the fasta matrix\n");

foreach(l=1:length(loci)) %dopar%{
  cat("working on locus", loci[l],"\n");
  locusLines <- grep(pattern=paste("^>",loci[l], "_", sep=""), fasta);
  taxa <- gsub(pattern="^>.*(Allele_\\d) \\[(.*); .*\\]$", replacement="\\2_\\1", fasta[locusLines]);
  cat("number of taxa found", length(taxa), "\n");
  if(length(taxa) != length(unique(taxa))){
    stop("duplicated taxa name\n" );
  }
  seq <- fasta[locusLines+1];
  cat("number of sequences before removing adapters:", length(seq), "\n");
  seq <- gsub(pattern=paste("^", adapter, "(.*$)",sep=""), replacement="\\1", seq);
  if(length(taxa) != length(seq)){
    stop("Error line 28: number of taxa not identical to number of sequences, taxa: ", length(taxa), " seq: ",length(seq));
  }
  if(single.snp){
    table <- data.frame();
    table <- do.call(rbind, strsplit(seq, ""));
    i <- TRUE;
    count <- 1
    while(i){
      bp <- as.character(table[,count]);
      if(length(unique(bp)) > 1){
        i <- F;
      }
      count <- count+1;
    }
    seq <- bp;
  }
  add.taxon <- NULL;
  add.seq <- NULL;
  count <- 1;
  for(j in 1:length(taxa)){ # Allele_1 will only appear in the fasta input if the indidivual is heterozygous, so here I duplicated phases in cases of homozygous to also have a identical Allele_1
    taxon <- gsub(pattern="^(.*)_Allele_\\d$", replacement="\\1", taxa[j]);
    if(any(grepl(pattern=paste(taxon, "_Allele_1", sep=""), taxa)) == FALSE){
     add.taxon <-  c(add.taxon,paste(taxon, "_Allele_1", sep=""));
     add.seq <- c(add.seq, seq[count]);
    }
    count <- count+1;
  }
  taxa <- c(taxa, add.taxon);
  seq <- c(seq, add.seq);
  if(length(taxa) != length(seq)){
    stop("Error line 57: number of taxa not identical to number of sequences, taxa: ", length(taxa), " seq: ",length(seq));
  }
  if(format == "fasta"){
    matrix <- paste(">", taxa, "\n",seq, sep="");
    outputName <- paste(loci[l], ".fa", sep="");
  }else if(format == "phylip"){
    seqLength <- length(unlist(strsplit(seq[1],"")));
    ntaxa <- length(taxa);
    header <- as.character(paste(ntaxa, seqLength));
    matrix <- paste(taxa, "     ", seq, sep="");
    matrix <- c(header, matrix);
    outputName <- paste(loci[l], ".phy", sep="");
  }else{
    stop("Error lines 60-70: unknown format type. Only fasta or phylip are possible");
  }
  write(matrix, outputName);
  cat(loci[l], "saved!", l, "/", length(loci), "\n");
  cat(" number of lines in written file:", length(matrix)*2,"\n")
    # preventing error, removed information below
  taxa <- NULL; 
  seq <- NULL; 
  matrix <- NULL;
  outputName <- NULL; 
}