## Please cite GitHub/melisaolave/Rscripts
## by Melisa Olave, 10th Octuber 2012;

################################################################################################
#1. Description
## This script will list a set of DNA matrices within a folder and run them into jModeltest 2.2.1
## The program will search among 24 models.
  ### JC, HKY, GTR,
  ### plus I and G
## The program will compute AIC and BIC.
## An output per matrix with results will be written in the same folder than the files are.

################################################################################################
# 2. Settings

    # You should check and specify the following settings

################################################################################################

files_dir <- "Users/myfiles/"; # change for your matrices directory
pattern <- "\\.fas" # change using an id pattern that is shared for all your matrix files (and is not in any other file in the folder), for example ".fas"
wd <- "Users/jModeltest/" # set jModeltest 2.2.1 directory

################################################################################################
# 3. Running the program (don't change anything from here!!)

################################################################################################
setwd(wd); # change for your jModeltest directory
dir <- getwd();
fileNames <- list.files(path = files_dir, pattern = pattern, all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE); # except for pattern, other commands are defaults
cat("A total of", length(fileNames), pattern, "files were found in", files_dir, ":\n", fileNames,"\n\t Ready to run jModeltest\n");

for (i in fileNames){
  cat("\tRunning", i,"file in jModeltest\n");
  command <- paste("java -jar jModelTest.jar -d ", files_dir,i, " -g 4 -i -f -AICc -BIC -a", sep="");
  # For example: java -jar jModelTest.jar -d example-data/aP6.fas -g 4 -i -f -AICc -BIC -a
  output <- system(command, 
                   intern = T,
                   ignore.stdout = F);
  output_name <- paste(i, "_jModeltest-result.txt", sep="");
  write.table(output, file.path(files_dir, output_name), quote = F, col.names=F,row.names=F, sep="\t");
  cat("\tOutput succesffully written under name:", output_name, "\n");
}

cat("Program done!\n");