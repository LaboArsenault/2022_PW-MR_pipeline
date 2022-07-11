#!/usr/bin/env Rscript



#### Bind MR results if FROM et TO were specified, especially if the 'decodetxt.loop.sh' script was used.
library(data.table)
library(optparse)
{
  option_list = list(
    ## GLOBAL OPTIONS
    make_option("--wd", action="store", default="./results/MRres/", type='character',
                help=" MR results directory (defaults to ./results/MRres/ but should also preferentially include the study name folder if using default folders."),
    make_option("--pattern", action="store", default=NA, type='character',
                help=" filename pattern "),
    make_option("--out", action="store", default=NA, type='character',
                help=" outname with path. Best to keep 'NA' for the rest of the pipeline to run with default MR results input argument."),
    make_option("--remove_files", action="store", default=FALSE, type='logical',
                help=" remove individual files after they are concatenated ? (default FALSE) ")
    )
  
  opt = parse_args(OptionParser(option_list=option_list))
}

files <- paste0(opt$wd, list.files(path = opt$wd))
if(!is.na(opt$pattern)){
  files <- subset(x = files, subset = grepl(pattern = opt$pattern, x = files))
}
if(is.na(opt$out)){
  opt$out <- paste0(opt$wd, "MR_results_sensitivity.txt")
}

total_rows <- 0
for(file in files){
  res_file <- fread(file = file, header = TRUE, stringsAsFactors = FALSE)
  if(file == files[1]){
    res <- res_file
  } else {
    res <- rbind(res, res_file)
  }
  total_rows <- total_rows + nrow(res_file)
}
if((nrow(res) == 0) & (total_rows != 0)){
  stop(" Concatenated results file is empty, but individual results files are not. 
       Before writing the concatenated results file and deleting the individual files, make sure this is what is intented.
       To go on with the script, set '--force TRUE'")
}

fwrite(x = res, file = opt$out, append = FALSE, sep = "\t")

# Removing files here to make sure script got to the end and concatenated results before deleting them.
if(opt$remove_files){
  for(file in files){
    unlink(file)
  }
}

