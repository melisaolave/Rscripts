# Please cite GitHub/melisolve/Rscripts
# Code to create IMa2 input
      ## by Melisa Olave (11/21/2012)

############################## Please README ############################################
# Before running:
    # 1. Check settings below
    # 2. You will need this following files contained at the working directory:
        # I. File with individual-species associations in exactly the same format requiered by *BEAST (and separeted by tab)
        # II. Matrix sequences with names locusName_whatever (Ej: cytb_boulengeri) and with JUST the taxaName "space" sequences. Please, avoid using uppercase letter in taxa name. Ej:
            #3160boul ATGCTAAATTT
            #3091rothi ATCCTAAATTG
            #5530tel ATGCTAATTTT
          # Remember that:
            # IMa2 doesn't allow ambiguities or missing data. This script will replace them with gaps (-), and IMa2 will ignore those possitions.
            # Max taxa name length = 10 characters.
        # III. A newick species tree with same names used in the (I) and branch length.
        # IV. Library ape is required.
    # 3. Every locus is going to be load as nuclear autosome (Inheritance scalar = 1). If you have mit or sex liked loci, make sure change that by hand later.
############################### Settings ####################################################

wd <- "/Users/myfolder/"; # change for the path of the folder cointaing files
pattern <- ".tnt";                 # change using an id pattern that is shared for all your matrix files (and is not in any other file in the folder), for example ".fas"
ImapFile <- "Imap.txt";            # change for Imap name (as used by *BEAST) for individual-species associations. IMa2 input will cointain taxa specified here. This mean you can have more sequences in your fasta matrices, but only this taxa will be used.
newickTree <- "sp_tree.newick";       # change for newick species tree file name
outputName <- "IMa2_infile.txt";   # change name of output file 
model <- "H";                      # Ej. HKY

#################### Running script (don't change anything from here) #######################
library(ape);
setwd(wd);
dir <- getwd();
fileNames <- list.files(path = dir, pattern = pattern, all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE); # except for pattern, other commands are defaults
cat("A total of", length(fileNames), pattern, "files were found in", dir, ":\n", fileNames,"\n\t Ready to read Imap\n");

Imap  <- read.table(ImapFile, header = T, sep="\t");
cat("A total of", length(unique(Imap$species)), "species/populations will be written\n");
tree <- read.tree(newickTree);

################################## to fix tree #############################################
keepSp <- as.character(unique(Imap$species));
dropSp <- setdiff(tree$tip.label, keepSp);
newTree <- drop.tip(tree, tip=dropSp);
taxaNames <- newTree$tip.label;
ntaxa <- length(newTree$tip.label)
newTree$tip.label <- as.character(1:length(newTree$tip.label));
nodes <- newTree$Nnode;
fixedTree <- write.tree(newTree);

fixedTree <- gsub(pattern="(\\d+):\\d+\\.\\d+", x = fixedTree, replacement="\\1");

for(i in 1:nodes){
  ntaxa <- ntaxa+1;
  fixedTree <- sub("\\d+\\.\\d+", x = fixedTree, replacement=as.character(ntaxa));
}

################################ locus information ###########################################
count2 <- 1;
locusInfo <- character();
new.matrix <- character();
for (i in 1:length(fileNames)){
  locus <- unlist(strsplit(fileNames[i], split="_"));
  locus <- locus[1];
  matrix <- scan(fileNames[i], what="character", sep="\n");
  count <- 1;
  ind <- character();
  indLength <- character();
  for (j in taxaNames){
    subImap <- Imap[Imap$species == j,];
    ind <- as.character(subImap$traits);
    nTaxa <- integer();
    count3 <- 1;
    orderedSeq <- character();
    for (k in ind){
      if(length(grep(pattern=k, matrix)) >0){
        name <- k;
        nTaxa[count3] <- grep(pattern=k, matrix); # detecting taxon coincident in matrix
        orderedSeq[count3] <- grep(pattern=k, matrix, value=T);
        orderedSeq[count3] <- gsub(pattern="[WSMYNRKVBDH?]", replacement="-", x=orderedSeq[count3]); #replace ambiguities or missing data for "-"
        sequence <- sub(pattern="^.* +([ATGC-]+)", replacement="\\1", orderedSeq[count3]); #taking seq
        nameLength <- length(unlist(strsplit(k, split="")));
        if (nameLength > 10) {
          cat("You taxa name", k," in ", fileNames[i], " file is longer than 10 characters!\n");
          name <- paste(unlist(strsplit(k, split=""))[1:9], collapse="");
          spaces <- " ";
        }else{
          nspaces <- 10-nameLength;
          spaces <- " ";
          for (l in 2:nspaces){
            spaces <- paste(spaces, " ", sep="");
          }
        }
        orderedSeq[count3] <- paste(name, spaces, sequence, sep="");
        count3 <- count3 + 1;
      }
    }
    indLength[count] <- as.character(length(nTaxa));
    new.matrix[count] <- paste(orderedSeq, collapse="\n");
    count <- count+1;
  }
  seq <- sub(pattern="^.* (.*)", x = matrix[1],replacement = "\\1");
  Length <- as.character(length(unlist(strsplit(seq, split=""))));
  locusInfo[count2] <- paste(locus, paste(indLength, collapse=" "), Length, model, "1\n", sep=" ");
  locusInfo[count2] <- paste(locusInfo[count2], paste(new.matrix, collapse="\n"),"\n", sep="");
  count2 <- count2 +1;
}

######################### pasting and saving all info ########################################
nGroups <- length(taxaNames);
IMa.input <- paste("IMa2.0 - input by Melisa Olave (11-21-2012)\n", # 1st line
                   nGroups, "\n", # number of species/pop in the matrix
                   paste(taxaNames, collapse=" "), "\n", #species/pop names
                   fixedTree, "\n", # tree
                   length(fileNames), "\n", #number of loci
                   paste(locusInfo, collapse=""), "\n", sep="");

write(file = outputName, x = IMa.input);
cat("Program done!\n\tIMa2 input file was written successfully under name", outputName, "\n");
