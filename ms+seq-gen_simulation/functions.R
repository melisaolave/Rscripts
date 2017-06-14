# Many functions to do many things
  ### by Melisa Olave, please cite http://www.github.com/melisaolave/

library(ape);

# get topological info from tree and codify to ms

get.topo2ms <- function(tree.newick, species, brlength.correction = NULL){
  count <- 1;
  new.tree.newick <- tree.newick;
  for(i in species){ # to replace species names by numbers, then ms can take them
    new.tree.newick <- gsub(pattern=i, replacement=count, x=new.tree.newick);
    count <- count+1;
  }
  topology <- character();
  for(j in 1:length(tree.newick)){ # in case you provided more than one sp tree topology
    if(length(tree.newick) > 1){
      print(paste("Working on tree number", j,"\n"));
    }
    if(grepl(pattern=";", x=new.tree.newick[j]) == F){ # if your sp tree doesn't have a ";" at the end, then add it
      new.tree.newick[j] <- gsub(pattern="(.*)$", replacement="\\1;", x=new.tree.newick[j]);
    }
    taxa1.Vec <- character();
    taxa2.Vec <- character();
    blength1.Vec <- numeric();
    blength2.Vec <- numeric();
    i <- T;
    while(i){
      if(length(unlist(strsplit(gsub(pattern="[^,]", replacement="", new.tree.newick), ","))) == 1){
        clade <- gsub(pattern="^.*(\\([[:alnum:]]+:\\d+\\.*\\d*,[[:alnum:]]+:\\d+\\.*\\d*\\):\\d+\\.*\\d*).*", 
                      replacement="\\1;", x=new.tree.newick[j]);
      }else{
        clade <- gsub(pattern="^.*(\\([[:alnum:]]+:\\d+\\.*\\d*,[[:alnum:]]+:\\d+\\.*\\d*\\):\\d+\\.*\\d*).*", 
                      replacement="\\1", x=new.tree.newick[j]);
      }
      if(length(unlist(strsplit(gsub(pattern="[^,]", replacement="", clade), ","))) > 1){ # if there is more than one "," in clade, means a polytomy
        #clade <- gsub(pattern="")
        cladeBackup <- clade;
        while(length(unlist(strsplit(gsub(pattern="[^,]", replacement="", clade), ","))) > 1){
          taxa1 <- gsub(pattern="^\\(([[:alnum:]]+):.*", replacement="\\1", x=clade);
          taxa2 <- gsub(pattern="^\\(.*,([[:alnum:]]+):.*", replacement="\\1", x=clade);
          blength1 <- as.numeric(gsub(pattern=paste("^\\(.*,[[:alnum:]]+:.*", taxa2, ":(\\d+\\.*\\d*).*", sep=""), replacement="\\1", x=clade));
          
          taxa1.Vec <- c(taxa1.Vec, taxa1);
          taxa2.Vec <- c(taxa2.Vec, taxa2);
          blength1.Vec <- c(blength1.Vec, blength1);
          clade <- gsub(pattern=paste("^(\\(", taxa1, ":.*),", taxa2, ":.*(\\).*)", sep=""),#"^.*(\\([[:alnum:]]+:\\d+\\.*\\d*,[[:alnum:]]+:\\d+\\.*\\d*\\):\\d+\\.*\\d*).*", 
                        replacement="\\1\\2", x=clade);     
        }
        new.clade <- clade;
        cladeBackup <- gsub(pattern="\\(", replacement="\\\\(", cladeBackup);
        cladeBackup <- gsub(pattern="\\)", replacement="\\\\)", cladeBackup);
        new.tree.newick[j] <- gsub(pattern=paste("(^.*)(", cladeBackup, ")(.*)", sep=""),            
                                   replacement=paste("\\1", new.clade, "\\3", sep=""), x=new.tree.newick[j]);
      }else if(grepl(pattern=";", x=clade) == F){
        taxa1 <- gsub(pattern="^\\(([[:alnum:]]+):.*", replacement="\\1", x=clade);
        taxa2 <- gsub(pattern="^\\(.*,([[:alnum:]]+):.*", replacement="\\1", x=clade);
        blength1 <- as.numeric(gsub(pattern="^\\([[:alnum:]]+:(.*),.*", replacement="\\1", x=clade));
        taxa1.Vec <- c(taxa1.Vec, taxa1);
        taxa2.Vec <- c(taxa2.Vec, taxa2);
        blength1.Vec <- c(blength1.Vec, blength1);
        blength2 <- as.numeric(gsub(pattern="^\\(.*,[[:alnum:]]+:.*:(\\d+\\.*\\d*)", replacement="\\1", x=clade));
        blength2.Vec <- c(blength2.Vec, blength2);
        new.blength <- blength1 + blength2;
        new.clade <- paste(taxa1, ":", new.blength, sep="");
        new.tree.newick[j] <- gsub(pattern="(^.*)(\\([[:alnum:]]+:\\d+\\.*\\d*,[[:alnum:]]+:\\d+\\.*\\d*\\):\\d+\\.*\\d*)(.*)",            
                                   replacement=paste("\\1", new.clade, "\\3", sep=""), x=new.tree.newick[j]);
      }else{
        taxa1 <- gsub(pattern="^\\(([[:alnum:]]+):.*", replacement="\\1", x=clade);
        taxa2 <- gsub(pattern="^\\(.*,([[:alnum:]]+):.*", replacement="\\1", x=clade);
        blength1 <- as.numeric(gsub(pattern="^\\([[:alnum:]]+:(.*),.*", replacement="\\1", x=clade));
        taxa1.Vec <- c(taxa1.Vec, taxa1);
        taxa2.Vec <- c(taxa2.Vec, taxa2);
        blength1.Vec <- c(blength1.Vec, blength1);
        new.blength <- blength1;
        new.clade <- paste(taxa1, ":", new.blength, sep="");
        new.tree.newick[j] <- gsub(pattern="(^.*)(\\([[:alnum:]]+:\\d+\\.*\\d*,[[:alnum:]]+:\\d+\\.*\\d*\\):\\d+\\.*\\d*)(.*)",            
                                   replacement=paste("\\1", new.clade, "\\3", sep=""), x=new.tree.newick[j]);
        i <- F;
      }
    }
    if(length(brlength.correction) != 0){ # this is to correct branch length unites if sp tree was estimated by a phylo inference program
      blength1.Vec <- blength1.Vec/brlength.correction;
    }
    blength1.Vec <- round(blength1.Vec, 4);
    topology[j] <- paste("-ej", blength1.Vec, taxa2.Vec, taxa1.Vec, collapse = " ");
  }
  return(topology)
}


# gene flow and pop size changes to ms

gf.code2ms <- function(gf.code, species){
  # -em t i j x
    # t = time; i = species i; j = species j; x = migration parameter (4Nm)
  if(gf.code == "-"){
    gf.code <- NULL;
  }else{
    count <- 1;
    gf.code <- gsub(pattern="(\\w+)$", replacement="\\1 ", x=gf.code); # to add an space at the end of gf.code provided
    for(i in species){ # to replace species names by numbers, then ms can take them
      gf.code <- gsub(pattern=paste("([[:space:]])",i, "([[:space:]])", sep=""), # the space is just to be sure it's matching taxa names, even when taxa names are numbers.
                      replacement=paste("\\1",count,"\\2", sep=""), 
                      x=gf.code);
      count <- count+1;
    }
    gf.code <- gsub(pattern="(\\w+)[[:space:]]$", replacement="\\1", x=gf.code); # to delete space at the end of code
    gf.code <- paste("-em", gf.code);
  }
  return(gf.code);
}


popSize.code2ms <- function(popSize, species){
  count <- 1;
  popSize.code <- gsub(pattern="(\\w+)$", replacement="\\1 ", x=popSize.code); # to add an space at the end of popSize.code provided
  for(i in species){ # to replace species names by numbers, then ms can take them
    popSize.code <- gsub(pattern=paste("([[:space:]])",i, "([[:space:]])", sep=""), # the space is just to be sure it's matching taxa names, even when taxa names are numbers.
                         replacement=paste("\\1",count,"\\2", sep=""), 
                         x=popSize.code);
    count <- count+1;
  }
  popSize.code <- gsub(pattern="(\\w+)[[:space:]]$", replacement="\\1", x=popSize.code); # to delete space at the end of code
  popSize.code <- paste("-en", popSize.code);
  return(popSize.code)
}

# get individual number from a fasta matrix and codify ms code

getTaxa2ms <- function(fileName, species, Imap){
  # read fasta matrix
  matrix <- scan(fileName, what="character", sep="\n");
  seq.starts <- grep(pattern=">", x=matrix);
  seq <- character();
  taxa <- character();
  taxaVec <- integer();
  for(j in 1:length(seq.starts)){ # to get taxa names vector
    if(j != max(length(seq.starts))){
      sequence <- matrix[(seq.starts[j]+1):(seq.starts[j+1]-1)];
    }else{
      sequence <- matrix[(seq.starts[j]+1):length(matrix)];
    }
    taxa[j] <- matrix[seq.starts[j]];
  }
  taxa <- gsub(pattern=">", replacement="", x=taxa); # deleting ">"
  
  for (k in species){ # to recognize number of individuals to be simulated per species
    subset <- Imap[Imap$species == k,];
    ind <- as.character(subset$traits);
    matched.taxa <- sum(ind %in% taxa);
    taxaVec <- c(taxaVec, matched.taxa); # vector of number of individual per species from seq matrices
  }
  return <- taxaVec
}


