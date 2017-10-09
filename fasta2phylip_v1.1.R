### Script for taking fasta matrices to phylip v1.1.
		### when using, please cite www.github.com/melisaolave
  ## what's new: write single snp added, run in parallel
  ## requires libraries: foreach and doMC
  

wd <- "/Users/myworkingdirectory"; 
outputDir <- "/Users/myworkingdirectory"
FastaInputPattern <- "\\.fa$"
deleteFullMissing <- F; # True if want to delete those individuals with full missing data sequences
aa <- F; # if input format contains amino acids, set TRUE here
single.snp <- T; # if TRUE, then only the first snp per locus is taken, if FALSE, then full sequence is included
cores <- 1; # number of cores to run in parallel

###### don't change anything from here #######
list.of.packages <- c("foreach", "doMC");
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])];
if(length(new.packages)){
  install.packages(new.packages);
}
library(foreach);
library(doMC);
setwd(wd);
files <- list.files(pattern=FastaInputPattern);
registerDoMC(cores);

foreach(i=1:length(files)) %dopar%{
  cat("reading matrix", files[i], "\n");
  matrix <- scan(files[i], what="character", sep="\n")
  start <- grep(pattern=">", matrix);
  end <- c(start[2:length(start)]-1, length(matrix));
  ind <- gsub(pattern = ">(.*)$", replacement="\\1", matrix[start]);
  ind <- gsub(pattern=" ", replacement="_", ind);
  seqVec <- NULL;
  for(j in 1:length(ind)){
    seqVec[j] <- paste(matrix[(start[j]+1):end[j]], collapse="");
  }
  if(single.snp){
    table <- data.frame();
    table <- do.call(rbind, strsplit(seqVec, ""));
    k <- TRUE;
    count <- 1
    while(k){
      bp <- as.character(table[,count]);
      if(length(unique(bp)) > 1){
        k <- F;
      }
      count <- count+1;
    }
    seqVec <- bp;
    seqLength <- 1;
  }else{
    seqLength <- length(unlist(strsplit(seqVec[1], split="")));
  }
  maxCharLength <- NULL;
  for(l in 1:length(ind)){ # get number of spaces from name to seq
    maxCharLength <- max(c(maxCharLength, length(strsplit(ind, "")[[l]])));
  }
  spacesVec <- NULL;
  for(k in 1:length(ind)){
    spaces <- maxCharLength - length(strsplit(ind, "")[[k]]) +1;
    spacesVec[k] <- paste(rep(" ",times=spaces), collapse="");
  }
  matrix <- paste(ind, spacesVec, seqVec); 
  header <- paste(length(ind), seqLength);
  matrix <- c(header, matrix);
  outputName <- paste(files[i], ".phy", sep="");
  setwd(outputDir);
  write(matrix, outputName);
  if(deleteFullMissing){
    matrix.txt <- scan(outputName, what="character", sep="\n");
    if(aa){
      del <- grep(pattern="^.* +-+$", matrix.txt);
      del <- c(del, grep(pattern="^.* +\\?+$", matrix.txt)); #?????
      del <- c(del, grep(pattern="^.* +-+\\?+$", matrix.txt)); # ----????
      del <- c(del, grep(pattern="^.* +\\?+-+$", matrix.txt)); # ???----
    }else{
      del <- grep(pattern="^.* +-+$", matrix.txt); # getting those rows with just ambiguities, below all combination of ? n -
      del <- c(del, grep(pattern="^.* +n+-+$", matrix.txt)); # n----
      del <- c(del, grep(pattern="^.* +-+n+$", matrix.txt)); #----n
      del <- c(del, grep(pattern="^.* +n+$", matrix.txt)); #nnnnnn
      del <- c(del, grep(pattern="^.* +\\?+$", matrix.txt)); #?????
      del <- c(del, grep(pattern="^.* +n+\\?+$", matrix.txt)); #nnn????
      del <- c(del, grep(pattern="^.* +-+\\?+$", matrix.txt)); # ----????
      del <- c(del, grep(pattern="^.* +\\?+-+$", matrix.txt)); # ???----
    }
    keep <- setdiff(1:length(matrix.txt), del);
    old.n.taxa <- as.integer(gsub(pattern="^(\\d+) .*", replacement="\\1",matrix.txt[1]));
    new.n.taxa <- old.n.taxa - length(del);
    matrix.txt[1] <- gsub(pattern="^(\\d+)", replacement=new.n.taxa,matrix.txt[1])
    write(matrix.txt[keep],  outputName);
  }
  cat("Matrix", outputName, "saved!", i, "/", length(files), "\n");
  setwd(wd);
}

