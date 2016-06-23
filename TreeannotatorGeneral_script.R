#### Please cite GitHub/melisaolve/Rscripts #####

# This script will:
  ##1. take all .trees files within a folder
  ##2. run them in Treeannotar 1.6 
  ##3. write the final tree.

# you'll need ape library

########## settings ####################
wd <- "Users/mydir"; # path to working directory;
burnin <- 1000; # set burnin
treeannotator_path <- "Users/BEAST.1.6.1/bin/" #set directory to TreeAnnotator executable
pattern <- "\\.trees"; # set the pattern you want to match in all files.
newick <- TRUE; # do you want to also save a tree in newick format? TRUE - yes; FALSE - no

########## running Treeannotator # don't change anything from here ####################
setwd(wd);
dir <- getwd();
files <- list.files(pattern=pattern);
library(ape);

for(i_files in files){
  fileName <- gsub(pattern=pattern, replacement="", x=i_files);
  output_name <- paste(dir,"/", fileName, "tree", sep="");  
  input <- paste(dir, "/",i_files, sep="");
  command <- paste("./treeannotator","-burnin", burnin, input, output_name);
  command <- paste (treeannotator_path, command, sep="");
  system(command);
    
  if(newick){
    myTree <- read.nexus(paste(fileName, "tree", sep="")); 
    write.tree(myTree, file=paste(fileName, "newick", sep="")); 
  }
}  

cat("Program done!\n");

