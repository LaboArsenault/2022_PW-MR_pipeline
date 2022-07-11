#!/bin/env Rscript

library(data.table)
library(optparse)
{
  option_list = list(
    make_option("--wd", action="store", default="./results/MRres/", type='character',
                help=" working directory in which to save results (defaults to './results/MRres/')"),
    make_option("--mr_res1", action="store", default=NA, type='character',
                help=" Path to MR results file for analysis 1 "),
    make_option("--mr_res2", action="store", default=NA, type='character',
                help=" Path to MR results file for analysis 2 "),
    make_option("--mr_res3", action="store", default=NA, type='character',
                help=" Path to MR results file for analysis 3 "),
    make_option("--study1", action="store", default=NA, type='character',
                help=" Study 1 name "),
    make_option("--study2", action="store", default=NA, type='character',
                help=" Study 2 name "),
    make_option("--study3", action="store", default=NA, type='character',
                help=" Study 3 name "),
    make_option("--method", action="store", default="fdr", type='character',
                help=" Method to select significant proteins (see p.adjust() for methods) "),
    make_option("--or", action="store", default=TRUE, type='character',
                help=" Output proteins that are significant in at least one dataset, whether it is present or not in the other. "),
    make_option("--and", action="store", default=FALSE, type='character',
                help=" Output proteins that are significant (and present) in both datasets. "),
    # make_option("--outcome_name", action="store", default="outcome", type='character',
    #             help=" outcome name (Default 'outcome')"),
    make_option("--out", action="store", default='all.proteins.signif.fdr.test.', type='character',
                help=" output name (Default to 'all.proteins.signif.fdr.[study1]_[OR/AND]_[study2][_[OR/AND]_[study3]].txt')
                For this argument, only specify the name pattern, not the studies, otherwise if you run 'or' and 'and', 'and' results will overwrite 'or' results.")
  )
  opt = parse_args(OptionParser(option_list=option_list))
  
  # opt <- data.frame(wd = "/mnt/sda/boujer01/Pancreatite/New_Analyses/results/MRres/",
  #                   mr_res1 = "/mnt/sda/boujer01/Pancreatite/New_Analyses/results/MRres/ARIC/MR_results_sensitivity_pval.5e-08_r2.0.1.txt",
  #                   mr_res2 = "/mnt/sda/boujer01/Pancreatite/New_Analyses/results/MRres/deCODE/MR_results_sensitivity_pval.5e-08_r2.0.1.txt",
  #                   mr_res3 = NA,
  #                   study1 = "ARIC",
  #                   study2 = "deCODE",
  #                   study3 = NA,
  #                   method = 'fdr',
  #                   or = TRUE,
  #                   and = TRUE,
  #                   out = 'all.proteins.signif.fdr.')
  
}
setwd(opt$wd)


# n_mr_res <- sum(!is.na(opt$mr_res1), !is.na(opt$mr_res2), !is.na(opt$mr_res3))
n_mr_res <- data.frame(s1=!is.na(opt$mr_res1), s2=!is.na(opt$mr_res2), s3=!is.na(opt$mr_res3))
n_mr_res <- data.frame(value = t(n_mr_res))
n_mr_res <- as.numeric(sapply(strsplit(x = rownames(subset(x = n_mr_res, subset = value)) , split = "s"), `[[`, 2))
if((3 %in% n_mr_res) & (!(1 %in% n_mr_res) | !(2 %in% n_mr_res))){
  stop("You specified MR results for study 3, but not for study 1 or 2. If you use only 2 studies, remove 'mr_res3' argument.")
}
# n_studies <- sum(!is.na(opt$study1), !is.na(opt$study2), !is.na(opt$study3))
n_studies <- data.frame(s1=!is.na(opt$study1), s2=!is.na(opt$study2), s3=!is.na(opt$study3))
n_studies <- subset(x = data.frame(value = t(n_studies)), subset = !is.na(value))
n_studies <- as.numeric(sapply(strsplit(x = rownames(subset(x = n_studies, subset = value)) , split = "s"), `[[`, 2))
if((3 %in% n_studies) & (!(1 %in% n_studies) | !(2 %in% n_studies))){
  stop("You specified a study name for study 3, but not for study 1 or 2. If you use only 2 studies, remove 'study3' argument.")
}

if((length(n_mr_res) < 2) | (length(n_studies) < 2)){
  stop("You have to specify at least 2 studies.")
}
if(length(n_mr_res) != length(n_studies)){
  stop(paste0("MR results provided for ", length(n_mr_res), " studies, but names for provided for ", length(n_studies), " studies."))
}
# bad_studies <- c()
# for(i in 1:max(n_studies, n_mr_res)){
#   eval(parse(text = paste0("if(xor(is.na(opt$study", i, "), is.na(opt$mr_res", i, "))){bad_studies <- c(bad_studies, ", i, ")} ")))
# }
# if(length(bad_studies) !=0){
#   stop(paste0("\nMissing argument for study : ", bad_studies))
# }

x_or <- opt$or
x_and <- opt$and

if(x_or){
  ## SIGNIF IN ONE OF ALL OR ALL
  subset_command <- c()
  outname_command <- c()
  for( i in 1:length(n_studies)){
    file <- eval(parse(text = paste0("fread(file = opt$mr_res", i, ", header = TRUE, stringsAsFactors=FALSE)")))
    file$padj <- p.adjust(p = file$pval, method = opt$method)
    file <- subset(x = file, select = c("exposure", "b", "se", "pval", "padj"))
    file <- file[with(file, order(padj, decreasing=FALSE)),]
    file <- subset(x = file, subset = !duplicated(exposure))
    file.clean <- file
    file.clean$`Beta(se)` <- paste0(round(x = file.clean$b, digits = 4), "(", round(x = file.clean$se, digits = 4), ")")
    file.clean$pval <- ifelse(file.clean$pval < 0.01, formatC(file.clean$pval, format = "e", digits = 2), round(file.clean$pval, digits = 3))
    file.clean$padj <- ifelse(file.clean$padj < 0.01, formatC(file.clean$padj, format = "e", digits = 2), round(file.clean$padj, digits = 3))
    file.clean$`P-value(padj)` <- paste0(file.clean$pval, "(", file.clean$padj, ")")
    file.clean <- subset(x = file.clean, select = c("exposure", "padj", "Beta(se)", "P-value(padj)"))
    eval(parse(text = paste0("colnames(file.clean) <- c('exposure', 'padj.", i, "', 'Beta(se).", i, "', 'P-value(padj).", i, "')")))
    file.clean <- as.data.frame(file.clean)
    if(i == 1){
      merged <- file.clean
    } else {
      merged <- merge(x = merged, y = file.clean, by = "exposure", all = TRUE)
    }
    # merged.clean <- merged
    eval(parse(text = paste0("merged$padj.", i, " <- ifelse(test = is.na(merged$padj.", i, "), yes = 99, no = merged$padj.", i, ")")))
    eval(parse(text = paste0("merged$padj.", i, " <- as.numeric(merged$padj.", i, ")")))
    subset_command <- c(subset_command, paste0("padj.", i))
    eval(parse(text = paste0("outname_command <- c(outname_command, opt$study",i,")")))
    if(i == length(n_studies)){
      for(j in 1:length(subset_command)){
        if(j == 1){
          subset_command.clean <- paste0('(', subset_command[j], ' < 0.05)')
          outname_command.clean <- paste0(outname_command[j])
        } else {
          subset_command.clean <- paste(subset_command.clean, paste0('(', subset_command[j], ' < 0.05)'), sep = " | ")
          outname_command.clean <- paste(outname_command.clean, outname_command[j], sep = "_OR_")
        }
      }
      out <- paste0(opt$wd, opt$out, outname_command.clean)
      eval(parse(text = paste0("merged.clean <- subset(x = merged, subset = ", subset_command.clean, ")")))
      eval(parse(text = paste0("merged.clean <- merged.clean[with(merged.clean, order(padj.", i, ", decreasing = FALSE)),]")))
      merged.clean <- merged.clean[, -c((n_studies*3) - 1)]
      fwrite(x = merged.clean, file = paste0(out, ".txt"), sep = "\t")
      message(paste0("Wrote 'OR' results to ", paste0(out, ".txt")))
    }
  }
}


if(x_and){
  ## ONLY SIGNIF IN ALL 
  outname_command <- c()
  for( i in 1:length(n_studies)){
    file <- eval(parse(text = paste0("fread(file = opt$mr_res", i, ", header = TRUE, stringsAsFactors=FALSE)")))
    file$padj <- p.adjust(p = file$pval, method = opt$method)
    file <- subset(x = file, subset = (padj < 0.05), select = c("exposure", "b", "se", "pval", "padj"))
    file <- file[with(file, order(padj, decreasing=FALSE)),]
    file <- subset(x = file, subset = !duplicated(exposure))
    file.clean <- file
    file.clean$`Beta(se)` <- paste0(round(x = file.clean$b, digits = 4), "(", round(x = file.clean$se, digits = 4), ")")
    file.clean$pval <- ifelse(file.clean$pval < 0.01, formatC(file.clean$pval, format = "e", digits = 2), round(file.clean$pval, digits = 3))
    file.clean$padj <- ifelse(file.clean$padj < 0.01, formatC(file.clean$padj, format = "e", digits = 2), round(file.clean$padj, digits = 3))
    file.clean$`P-value(padj)` <- paste0(file.clean$pval, "(", file.clean$padj, ")")
    file.clean <- subset(x = file.clean, select = c("exposure", "padj", "Beta(se)", "P-value(padj)"))
    eval(parse(text = paste0("colnames(file.clean) <- c('exposure', 'padj.", i, "', 'Beta(se).", i, "', 'P-value(padj).", i, "')")))
    file.clean <- as.data.frame(file.clean)
    if(i == 1){
      merged <- file.clean
    } else {
      merged <- merge(x = merged, y = file.clean, by = "exposure", all = FALSE)
    }
    eval(parse(text = paste0("merged$padj.", i, " <- as.numeric(merged$padj.", i, ")")))
    eval(parse(text = paste0("outname_command <- c(outname_command, opt$study",i,")")))
    if(i == length(n_studies)){
      for(j in 1:length(n_studies)){
        if(j == 1){
          outname_command.clean <- paste0(outname_command[j])
        } else {
          outname_command.clean <- paste(outname_command.clean, outname_command[j], sep = "_AND_")
        }
      }
      out <- paste0(opt$wd, opt$out, outname_command.clean)
      eval(parse(text = paste0("merged.clean <- merged[with(merged, order(padj.", i, ", decreasing = FALSE)),]")))
      merged.clean <- merged.clean[, -c((n_studies*3) - 1)]
      fwrite(x = merged.clean, file = paste0(out, ".txt"), sep = "\t")
      message(paste0("Wrote 'AND' results to ", paste0(out, ".txt")))
    }
  }
}
