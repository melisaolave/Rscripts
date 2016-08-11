# Script to use ms and seq-gen to simulate gene trees and mutations
  # by Melisa Olave, please cite http://www.github.com/melisaolave/

############################################ README ########################################################
#Before running, you will need:
# 1. a file with one or more species tree in newick format

########################################### Settings ########################################################
#rm(list=ls(all=TRUE)); #don't change
wd <- "Users/mypath" # change for your working directory

# requiered input files
# to simulate genetrees and sequences
sp.treeName <- "my_sptrees.txt" # change for your species tree file in newick format
fixed.nsamples <- 5;                      # if nsamples.from.matrices = FALSE; then provide a fixed number of individuals per lineage you want to simulate

# simulation settings and parameters
nloci <- 1000;                           # change for number of independent loci you want to simulate
geneflow <- FALSE;                      # TRUE if you want to simulate gene flow between species
gf.code <- NULL;                        # if geneflow is TRUE, then provide "time sp1 sp2 geneflowParameter". Original taxa names are permited (this script will change them later)
popSize <- FALSE;                       # TRUE if you want to model population size changes
popSize.code <- NULL;           # if popSize is TRUE, then provide "time sp NeParameter". Original taxa names are permited (this script will change them later)

# simulate sequences?
seq.gen <- F;                       # TRUE if you want to simulate DNA sequences, using ms gene trees
seqLength <- 500;                        # if seq.gen is TRUE change for sequence length
model <- "HKY";                           # if seq.gen is TRUE change for evolutionary model to be simulated
theta <- 0.01;                            # if seq.gen is TRUE change for theta value
seq.outputName <- "seq";                  # if seq.gen is TRUEchange sequence output file name

# output names
outputName <- "genetrees";              # change output file name
seq.outputName <- "seq"

############################# Running the program (don't change anything from here) ###########################

############################################### getting started  ##################################################
library(ape);
setwd(wd);
source(file.path(wd, "functions.R"));
dir <- getwd();
################################################## go!  ########################################################

# reading files (species tree)
  # sp tree
tree <- read.tree(sp.treeName);
tree.newick <- scan(sp.treeName, what="character", sep="\n"); #reading tree as txt file
if(length(tree.newick) == 1){
  species <- sort(tree$tip.label);
}else{
  species <- sort(tree[[1]]$tip.label); # if nsamples.from.matrices is TRUE, we'll only use the species that are in the Imap file
}
nspecies <- length(species);

# getting tau info from species tree file
topology <- get.topo2ms(tree.newick, species, brlength.correction = NULL);

# geneflow and popSize parameters
count <- 1;
if(geneflow){
  for(i in species){ # to replace species names by numbers, then ms can take them
    gf.code <- gsub(pattern=paste("([[:space:]])",i, "([[:space:]])", sep=""), # the space is just to be sure it's matching taxa names, even when taxa names are numbers.
                    replacement=paste("\\1",count,"\\2", sep=""), 
                    x=gf.code);
    count <- count+1;
  }
  gf.code <- paste("-em", gf.code);  
}else{
  gf.code <- NULL;
}
if(popSize){
  count <- 1;
  for(i in species){ # to replace species names by numbers, then ms can take them
    popSize.code <- gsub(pattern=paste("([[:space:]])",i, "([[:space:]])", sep=""), # the space is just to be sure it's matching taxa names, even when taxa names are numbers.
                         replacement=paste("\\1",count,"\\2", sep=""), 
                         x=popSize.code);
    count <- count+1;
  }
  popSize.code <- paste("-en", popSize.code);
}else{
  popSize.code <- NULL;
}

taxaVec <- integer();
# simulating gene trees for each matrix and comparing gene trees and sp tree
 if (length(fixed.nsamples) == 1){  # if same number of individual per every sp
    for (z in 1:length(species)){
      taxaVec <- c(taxaVec, fixed.nsamples);
    }
  }else if (length(fixed.nsamples) == length(species)){ # if want to simulate sp with different number of individuals
    taxaVec <- fixed.nsamples;
  }else{
    stop(paste("Error: number of species (=", length(species), 
               ") is not equal to number of fixed.nsamples (=", 
               length(fixed.nsamples), 
               ") you specified\n", sep=""));
  }
  totalTaxa <- sum(taxaVec);
  taxaVec <- paste(taxaVec, collapse=" ");

for (j in 1:length(tree.newick)){ # in case you have multiple sp trees, add a prefix to each simulated file
  if(length(tree.newick) > 1){
    ST <- paste("_ST", j, sep="");
  }else{
    ST <- NULL;
  }
  new.outputName <- paste(outputName, ST,".txt", sep="");
  command <- paste("./ms", 
                   totalTaxa, 
                   nloci, 
                   "-T -I", 
                   length(species), 
                   taxaVec, 
                   topology,
                   gf.code,
                   popSize.code,
                   "| tail +4 | grep -v // >", 
                   new.outputName);
  system(command);
  
  # simulating sequencies from genetrees
  if(seq.gen){
    newSeq.outputName <- paste(seq.outputName, ST,".txt", sep="");
    command <- paste("./seq-gen", 
                     " -m", model, 
                     " -l", seqLength,
                     " -n1",
                     " -s", theta,
                     " -a0.8 -g4 -f0.3 0.2 0.3 0.2 -t3.0 -op",
                     " <", new.outputName,
                     " >", newSeq.outputName,
                     sep="");
    system(command);
    read.dna(seq.outputName, format="sequential");
  }
}


