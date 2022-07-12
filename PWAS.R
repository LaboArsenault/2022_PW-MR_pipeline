#!/usr/bin/env Rscript


############################################################################### #
#### i. Script options ####
library(optparse)
{
  option_list = list(
    ## GLOBAL OPTIONS
    make_option("--wd", action="store", default="./", type='character',
                help=" working directory from which to access files"),
    make_option("--window", action="store", default=5e05, type='numeric',
                help=" inclusion and clumping window (default 500 kb)"),
    make_option("--pval", action="store", default=5e-08, type='numeric',
                help=" pval threshold (default 5e-05)"),
    make_option("--r2", action="store", default=0.1, type='numeric',
                help=" clumping threshold (default 0.01)"),
    make_option("--from", action="store", default=NA, type='numeric',
                help=" If you want to run this script in multiple sequences, select which min protein number to include"),
    make_option("--to", action="store", default=NA, type='numeric',
                help=" If you want to run this script in multiple sequences, select which max protein number to include"),
    make_option("--study", action="store", default="deCODE", type='character',
                help=" which study as exposure ? (INTERVAL, deCODE (default) or ARIC)"),
    make_option("--study_dir", action="store", default="/home/couchr02/Mendel_Commun/dbs_Web/deCODE/", type='character',
                help=" study directory where GWASes are stored (deCODE dir (default))"),
    make_option("--genes_ref_file", action="store", default="./data/GRCh38_chr_pos_genes.txt", type='character',
                help=" path to genes reference file (name, pos, chr) ( HGNC Build 38 (default))"),
    make_option("--snps_ref_file", action="store", default="./data/ukb_rsid_v3.txt", type='character',
                help=" (Optional) path to SNPs reference file (rsid, pos, chr). This is only required if SNPs positions/chr or RSIDs are required."),
    make_option("--plink_binaries", action="store", default="./data/EUR_rs", type='character',
                help=" (Optional) path to plink binary files for local clumping "),
    make_option("--exp_info_file", action="store", default="./data/decode_infos.txt", type='character',
                help=" (Optional) path to protein info file. IF study is not INTERVAL, deCODE or ARIC, this file must absolutely contain prot_name, id, trait and file columns. 'study_dir' and 'file' concatenated must correspond to the absolute path to the file. "),
    make_option("--cojo_ldref", action="store", default="./cojo/ldref/chr1_22.b37/all_phase3", type='character',
                help=" (Optional) path to cojo binary PLINK files including everything except extension. 
                File name must therefore be in that exact format (COJO_LDREF_FILE.{ext}), so not splitted by chromosome. 
                Make sure you have the right genome build (must be the same between exposure and LD reference"),
    make_option("--r2_corr", action="store", default=0.6, type='numeric',
                help=" (Optional) R2 threshold for analyses including LD-correction (default 0.6) "),
    make_option("--cojo_ldref_rsids", action="store", default=1, type='numeric',
                help=" (Optional) Does LD reference file for CoJo contains SNPs reference as rsids (1) instead of chr:pos:a1:a2 (0) ? (default 0) "),
    make_option("--reverse", action="store", default=FALSE, type='logical',
                    help=" (Optional) Are you performing a reverse MR ? (default FALSE)
                    In that case, keep the options as they are, exposure will automatically by reversed as outcome, and vice-versa.
                    If performing MR reverse, only "),
    
    ## EXPOSURE OPTIONS
    make_option("--exp_beta_col", action="store", default="Beta", type='character',
                help=" beta col"),
    make_option("--exp_se_col", action="store", default="SE", type='character',
                help=" SE col"),
    make_option("--exp_chr_col", action="store", default=NA, type='character',
                help=" chromosome col"),
    make_option("--exp_pos_col", action="store", default=NA, type='character',
                help=" position col"),
    make_option("--exp_pval_col", action="store", default="Pval", type='character',
                help=" P-value col"),
    make_option("--exp_ea_col", action="store", default="effectAllele", type='character',
                help=" effect allele col"),
    make_option("--exp_oa_col", action="store", default="otherAllele", type='character',
                help=" other allele col"),
    make_option("--exp_eaf_col", action="store", default=NA, type='character',
                help=" effect allele frequency col"),
    make_option("--exp_ss_col", action="store", default=NA, type='character',
                help=" sample size col"),
    make_option("--exp_snp_col", action="store", default="SNP", type='character',
                help=" SNP rsids col"),
    make_option("--exp_ss_total", action="store", default=NA, type='numeric',
                help=" Total samplesize for exposure"),
    make_option("--exp_add_rsids", action="store", default=FALSE, type='logical',
                help=" (Optional) should you add rsids to your exposure file"),
    make_option("--exp_p_is_log", action="store", default=FALSE, type='logical',
                help=" (Optional) Is exposure pvalue in log10 ? "),
    make_option("--exp_check_alleles", action="store", default=FALSE, type='logical',
                help=" (Optional) If dataset used is ARIC, do alleles need to be checked to make sur ALT is always A1 (default FALSE) ? "),
    make_option("--exp_is_vcf", action="store", default=FALSE, type='logical',
                help=" is exposure file a vcf file ? "),
    make_option("--exp_vcf_skip", action="store", default=1, type='numeric',
                help=" number of header lines to skip in vcf file. Do not specify anything if file is not vcf, otherwise CoJo will crash. "),
    make_option("--exp_snp_col_num", action="store", default=3, type='numeric',
                help=" exposure SNP column number "),
    make_option("--exp_chr_col_num", action="store", default=1, type='numeric',
                help=" exposure chromosome column number "),
    make_option("--exp_pos_col_num", action="store", default=2, type='numeric',
                help=" exposure position column number "),
    make_option("--exp_a1_col_num", action="store", default=4, type='numeric',
                help=" exposure allele 1 column number "),
    make_option("--exp_a2_col_num", action="store", default=5, type='numeric',
                help=" exposure allele 2 column number "),
    make_option("--exp_form_col_num", action="store", default=12, type='numeric',
                help=" exposure beta column number "),
    ## IF FILE IS VCF, the next five columns number should be specified 
    ## according to their position in the FORMAT column, not in 
    make_option("--exp_maf_col_num", action="store", default=4, type='numeric',
                help=" exposure maf column number "),
    make_option("--exp_beta_col_num", action="store", default=5, type='numeric',
                help=" exposure beta column number "),
    make_option("--exp_se_col_num", action="store", default=6, type='numeric',
                help=" exposure se column number "),
    make_option("--exp_p_col_num", action="store", default=7, type='numeric',
                help=" exposure pval column number "),
    make_option("--exp_n_col_num", action="store", default=8, type='numeric',
                help=" exposure N column number "),
    make_option("--exp_ncase", action="store", default=NA, type='numeric',
                help=" exposure ncases for coloc only (if NA, coloc will use 'type = quant' instead of 'type = cc'). Note that 'exp_ss_total' is mandatory if this argument is specified. "),
    make_option("--exp_cojo_snp_col", action="store", default="Name", type='character',
                help=" exposure SNP column for CoJo (SNPs in format chr:pos:a1:a2 ) "),
    make_option("--exp_cojo_snp_col_num", action="store", default=3, type='numeric',
                help=" exposure SNP column number for CoJo (SNPs in format chr:pos:a1:a2 ) "),
    
    ## OUTCOME OPTIONS 
    make_option("--out_gwas_file", action="store", default="./data/meta_clean_Acute_Pancreatitis.txt.gz", type='character',
                help=" GWAS file for outcome. MUST NOT CONTAIN ANY OTHER HEADER THAN COLUMN NAMES (if vcf, for example). Other methods to read outcome GWAS are not implemented yet."),
    make_option("--out_pheno", action="store", default="Acute_Pancreatitis", type='character',
                help=" Phenotype for outcome"),
    make_option("--out_beta_col", action="store", default="Beta", type='character',
                help=" beta col"),
    make_option("--out_se_col", action="store", default="SE", type='character',
                help=" SE col"),
    make_option("--out_chr_col", action="store", default=NA, type='character',
                help=" chromosome col"),
    make_option("--out_pos_col", action="store", default=NA, type='character',
                help=" position col"),
    make_option("--out_pval_col", action="store", default="Pval", type='character',
                help=" P-value col"),
    make_option("--out_ea_col", action="store", default="effectAllele", type='character',
                help=" effect allele col"),
    make_option("--out_oa_col", action="store", default="otherAllele", type='character',
                help=" other allele col"),
    make_option("--out_eaf_col", action="store", default=NA, type='character',
                help=" effect allele frequency col"),
    make_option("--out_ss_col", action="store", default=NA, type='character',
                help=" sample size col"),
    make_option("--out_snp_col", action="store", default="rsids", type='character',
                help=" SNP rsids col"),
    make_option("--out_ss_total", action="store", default=NA, type='numeric',
                help=" Total samplesize for outcome"),
    make_option("--out_add_rsids", action="store", default=FALSE, type='logical',
                help=" (Optional) should you add rsids to your outcome file"),
    make_option("--out_p_is_log", action="store", default=FALSE, type='logical',
                help=" (Optional) Is outcome pvalue in log10 ? "),
    make_option("--out_ncase", action="store", default=NA, type='numeric',
                help=" outcome ncases for coloc only (if NA, coloc will use 'type = quant' instead of 'type = cc'). Note that 'out_ss_total' is mandatory if this argument is specified. "),
    
    ## OUTPUT OPTIONS
    make_option("--save_every", action="store", default=1, type='numeric',
                help=" save results every X analysis (default 1)"),
    make_option("--outname_res", action="store", default=NA, type='character',
                help=" output name of results"),
    make_option("--outname_dat", action="store", default=NA, type='character',
                help=" output name of MR data (harmonized / clumped)"),
    make_option("--outname_coloc", action="store", default=NA, type='character',
                help=" output name of coloc results"),
    make_option("--save_dat", action="store", default=TRUE, type='logical',
                help=" should you save MR data for each protein? (default TRUE) "),
    make_option("--nthreads", action="store", default=1, type='numeric',
                help=" How many threads should you use to perform analyses ? (default 1) ")
  )
  opt = parse_args(OptionParser(option_list=option_list))
}
############################################################################### #


############################################################################### #
#### ii. Loading libraries ####
library(data.table)
library(TwoSampleMR)
library(R2BGLiMS)
library(ieugwasr)
library(MendelianRandomization)
library(coloc)
setwd(opt$wd)
setDTthreads(threads = opt$nthreads, throttle = 1024*opt$nthreads)
############################################################################### #


############################################################################### #
#### 1. Parameters ####
{
  for(i in 1:length(opt)){
    if(!is.na(opt[[i]])){
      if(opt[[i]] == "NA"){
        opt[[i]] <- NA
      }
    }
  }
  from_na <- !is.na(opt$from)
  to_na <- !is.na(opt$to)
  use_from_to <- ifelse(test = (from_na+to_na)!=1,
                        yes = TRUE,
                        no = FALSE)
  if(!use_from_to){
    stop(" You specified only one of the 'from' and 'to' arguments. Specify either none or both of them. ")
  }
  # if(!(opt$study %in% c("INTERVAL", "deCODE", "ARIC"))){
  if(!(grepl(pattern = "deCODE", x = opt$study) | !grepl(pattern = "ARIC", x = opt$study) | !grepl(pattern = "INTERVAL", x = opt$study))){
    message(" Study name provided is not in INTERVAL, deCODE or ARIC. This script has not been tested yet on other studies and may be unstable for this analysis. ")
  }
  if(is.na(opt$exp_ss_total)){
    stop("Total samplesize for exposure is required.")
  }
  # from <- as.numeric(opt$from)
  # to <- as.numeric(opt$to)
  window <- as.numeric(opt$window)
  pval_threshold <- as.numeric(opt$pval)
  clump_r2 <- as.numeric(opt$r2)
  if(is.na(opt$outname_res)){
    if(use_from_to & (from_na+to_na==2)){
      outname <- paste0("./results/MRres/", opt$study, "/MR_results_sensitivity_pval.", pval_threshold, "_r2.", clump_r2, "_from.", opt$from, "_to.", opt$to, ".txt")
    } else {
      outname <- paste0("./results/MRres/", opt$study, "/MR_results_sensitivity_pval.", pval_threshold, "_r2.", clump_r2, ".txt")
    }
  } else {
    outname <- opt$outname_res
  }
  if(opt$save_dat){
    if(is.na(opt$outname_dat)){
      outname_dat <- paste0("./results/MRdat/", opt$study, "/MR_dat_")
    } else {
      outname_dat <- opt$outname_dat
    }
  }
  if(is.na(opt$outname_coloc)){
    outname_coloc <- paste0("./results/coloc/", opt$study, "/coloc_res_")
  } else {
    outname_coloc <- opt$outname_coloc
  }
  unprocessed_proteins <- data.frame(prot = character(), reason = character())
  if(!dir.exists(paste0("./results/unprocessed/", opt$study, "/"))){
    dir.create(paste0("./results/unprocessed/", opt$study, "/"))
  }
  if(use_from_to & (from_na+to_na==2)){
    unprocessed_proteins_file <- paste0("./results/unprocessed/", opt$study, "/", opt$study, "_from.", opt$from, "_to.", opt$to, ".unprocessed.txt")
  } else {
    unprocessed_proteins_file <- paste0("./results/unprocessed/", opt$study, "/", opt$study, ".unprocessed.txt")
  }
  if((is.na(opt$exp_ss_total) & !is.na(opt$exp_ncase)) | (is.na(opt$out_ss_total) & !is.na(opt$out_ncase))){
    stop("'exp/out_ncase' argument is provided, but 'exp/out_ss_total' is not. This is required to calculate the ratio of cases to controls for COLOC.")
  }
  for(directory in c("MRres", "MRdat", "coloc")){
    if(!dir.exists(paste0("./results/", directory, "/", opt$study))){
      dir.create(path = paste0("./results/", directory, "/", opt$study))
    }
  }
}
############################################################################### #


############################################################################### #
#### 2. Proteins directory ####
{
  prots_wd <- opt$study_dir
  if(!opt$exp_is_vcf){
    prots_dir <- list.files(path = opt$study_dir)
    # if(is.na(from) | is.na(to)){
    #   from = 1
    #   to = length(prots_dir)
    # }
  } else {
    prots_dir <- prots_wd
  }
}
############################################################################### #


############################################################################### #
#### 3. Proteins files ####
{
  if(opt$exp_is_vcf){
    vcf_info <- fread(opt$exp_info_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads)
    # vcf_info <- subset(x = vcf_info, subset = grepl(pattern = opt$study, x = consortium))
    vcf_info$path <- paste0(prots_wd, "/", vcf_info$id, "/", vcf_info$id, ".vcf.gz")
    prot_dataset <- vcf_info
    prot_dataset <- subset(x = prot_dataset, subset = id %in% list.files(prots_wd))
    prot_dataset$prot_name <- prot_dataset$note
  } else {
    if(grepl(pattern = "INTERVAL", x = opt$study)){
      # SOMAMER_ID,UniProt,TargetFullName,Target
      prot_dataset <- as.data.frame(fread(file = opt$exp_info_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
      prot_dataset$file <- paste0(prot_dataset$id, "_chrom_")
      prot_dataset$path <- paste0(opt$study_dir, prot_dataset$id, "/", prot_dataset$file)
      prot_dataset <- subset(x = prot_dataset, subset = path %in% prots_dir)
    } 
    if(grepl(pattern = "deCODE", x = opt$study)){
      prot_dataset <- data.frame(file = sapply(strsplit(x = list.files(prots_wd), split = "[.]"), `[[`, 1))
      prot_dataset <- subset(x = prot_dataset, subset = paste0(file, ".txt.gz") %in% prots_dir)
      # prot_dataset$prot_name <- paste0(sapply(strsplit(x = prot_dataset$file, split = "_"), `[[`, 3))
      prot_dataset$path <- paste0(prots_wd, prots_dir)
    }
    if(grepl(pattern = "ARIC", x = opt$study)){
      prot_dataset <- as.data.frame(fread(file = opt$exp_info_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
      # prot_dataset$prot_name <- prot_dataset$entrezgenesymbol
      prot_dataset$file <- paste0(prot_dataset$file, ".PHENO1.glm.linear")
      prot_dataset$path <- paste0(opt$study_dir, prot_dataset$file)
      prot_dataset <- subset(x = prot_dataset, subset = file %in% prots_dir)
    } else {
      prot_dataset <- as.data.frame(fread(file = opt$exp_info_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
      prot_dataset$path <- paste0(opt$study_dir, prot_dataset$file)
      prot_dataset <- subset(x = prot_dataset, subset = file %in% prots_dir)
    }
    rm(prots_dir)
  }
}
############################################################################### #


############################################################################### #
#### 4. Genes reference files ####
{
  if(opt$study == "INTERVAL"){
    # INTERVAL is b37 and requires a different reference file
    genes <- as.data.frame(fread(file = opt$genes_ref_file, header = FALSE, stringsAsFactors = FALSE, nThread = opt$nthreads))
    colnames(genes) <- c("CHR", "start", "end", "GeneName")
  } else {
    genes <- as.data.frame(fread(file = opt$genes_ref_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
    genes <- subset(x = genes, subset = `Chromosome/scaffold name` %in% c(1:22), nThread = opt$nthreads)
    genes$`Chromosome/scaffold name` <- as.numeric(genes$`Chromosome/scaffold name`)
    genes <- genes[,c(2:5)]
    colnames(genes) <- c("GeneName", "CHR", "start", "end")
  }
  genes <- subset(x = genes, subset = CHR <= 22)
  genes$CHR <- as.numeric(genes$CHR)
  genes$start <- as.numeric(genes$start)
  genes$end <- as.numeric(genes$end)
  # Genes that are common with proteins
  genes_dataset <- subset(x = genes, subset = GeneName %in% prot_dataset$prot_name)
  genes_dataset <- subset(x = genes_dataset, subset = !duplicated(GeneName))
}
############################################################################### #


############################################################################### #
#### 5. Merge genes reference and proteins files ####
{
  gene_prot_file <- merge(x = prot_dataset, y = genes_dataset, by.x = "prot_name", by.y = "GeneName", all.x = FALSE, all.y = FALSE)
  # if(use_from_to  & (from_na+to_na==2)){
  #   gene_prot_file <- gene_prot_file[opt$from:opt$to,]
  # }
}
############################################################################### #


############################################################################### #
#### 6. SNPs rsids and position reference file (*b38*) ####
{
  if(opt$exp_add_rsids | opt$out_add_rsids){
    snps_ref <- as.data.frame(fread(file = opt$snps_ref_file, header = FALSE, stringsAsFactors = FALSE, nThread = opt$nthreads))
  }
  chr_lengths <- as.data.frame(fread(file = "./data/chromosomes_lengths.txt", header = TRUE, stringsAsFactors = FALSE))
  chr_lengths <- subset(x = chr_lengths, subset = Chr %in% c(1:22))
  chr_lengths$Chr <- as.numeric(chr_lengths$Chr)
}
############################################################################### #


############################################################################### #
#### 7. Outcome GWAS ####
{
  out_gwas <- as.data.frame(fread(file = opt$out_gwas_file, header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
  if(opt$out_add_rsids){
    if(is.na(opt$out_chr_col) | is.na(opt$out_pos_col)){
      stop(" *** No chromosome and/or position columns were supplied for outcome GWAS. Cannot add rsids to the file. *** \n Please set the 'add_rsids' option to FALSE or specify the 'out_chr_col' and 'out_pos_col' arguments. ")
    } else {
      ## Add rsids
      colnames(snps_ref) <- c("SNP", "pos", "chr")
      out_gwas <- merge(x = out_gwas, y = snps_ref, by.x = c(opt$out_chr_col, opt$out_pos_col), by.y = c("chr", "pos"))
    }
  }
  out_gwas$phenotype <- opt$out_pheno
  out_gwas$samplesize <- opt$out_ss_total
  
  # Format AP GWAS 
  command_format <- paste0("outcome_format <- TwoSampleMR::format_data(dat = out_gwas, type = 'outcome'",
                           ", snp_col = '", opt$out_snp_col,
                           "', beta_col = '", opt$out_beta_col,
                           "', se_col = '", opt$out_se_col,
                           "', effect_allele_col = '", opt$out_ea_col,
                           "', other_allele_col = '", opt$out_oa_col,
                           "', pval_col = '", opt$out_pval_col,
                           "', phenotype_col = 'phenotype'")
  if(!is.na(opt$out_pos_col)){
    command_format <- paste0(command_format, ", pos_col = '", opt$out_pos_col, "'")
  }
  if(!is.na(opt$out_chr_col)){
    command_format <- paste0(command_format, ", chr_col = '", opt$out_chr_col, "'")
  }
  if(!is.na(opt$out_eaf_col)){
    command_format <- paste0(command_format, ", eaf_col = '", opt$out_eaf_col, "'")
  }
  if(!is.na(opt$out_ss_col)){
    command_format <- paste0(command_format, ", samplesize_col = '", opt$out_ss_col, "'")
  }
  command_format <- paste0(command_format, ")")
  eval(parse(text = command_format))
  rm(out_gwas)
  if(is.na(opt$out_ss_col)){
    outcome_format$samplesize.outcome <- opt$out_ss_total
  }
  if(opt$out_p_is_log){
    outcome_format$pval.outcome <- 10**(outcome_format$pval.outcome)
  }
  # If performing reverse MR, keeping only SNPs with p < threshold in outcome.
  if(opt$reverse){
    outcome_all.clean <- outcome_format
    outcome_format.clean <- subset(x = outcome_format, subset = pval.outcome <= pval_threshold)
  }
}
############################################################################### #


############################################################################### #
#### 8. Preparing results ####
{
  mr_res <- data.frame(matrix(nrow = 0, ncol = 48))
  colnames(mr_res) <- c("outcome", "exposure", "method", "nsnp", 
                        "b", "b_CIlow", "b_CIhigh", "se", "pval", 
                        "b_egger", "b_egger_CIlow", "b_egger_CIhigh","se_egger", "pval_egger",
                        "nsnp_egger_i", "b_egger_i", "se_egger_i", "pval_egger_i",
                        "b_weighted_median", "b_weightedmed_CIlow", "b_weightedmed_CIhigh", "se_weighted_median", "pval_weighted_median",
                        "Q_Egger", "P_Q_Egger", "Q_IVW", "P_Q_IVW",
                        "jammr.b", "jammr.se", "jammr.nsnps", "pca.b", "pca.se", "pca.nsnps", "pca.npc",
                        "cojo.b.slct", "cojo.se.slct", "cojo.p.slct", "cojo.nsnps.slct", "cojo.b.cond", "cojo.se.cond", "cojo.p.cond", "cojo.nsnps.cond",
                        "coloc.nsnps", "coloc.PPH0", "coloc.PPH1", "coloc.PPH2", "coloc.PPH3", "coloc.PPH4")
}
save_res <- FALSE
############################################################################### #


############################################################################### #
#### 9. Begin loop on each protein ####
{
  start_all <- Sys.time()
  loop_command <- ifelse(test = (use_from_to  & (from_na+to_na==2)),
                         yes = paste0(opt$from, ":", opt$to),
                         no = "1:nrow(gene_prot_file)")
  # for(i in 1:nrow(gene_prot_file)){
  for(i in eval(parse(text = loop_command))){
    message(paste0("\n\tProcessing protein ", gene_prot_file$prot_name[i], " (", gene_prot_file$file[i], ") (num :", i,") ...\n"))
    start <- Sys.time()
    if((i %% opt$save_every) == 0){
      save_res <- TRUE
    }
    window <- as.numeric(opt$window)
    
    
    ############################################################## #
    #### __9.1. Reading and preparing exposure files ####
    if(opt$exp_add_rsids){
      ## Add rsids
      colnames(snps_ref) <- c("SNP", "pos", "chr")
      snps_ref_clean <- subset(x = snps_ref, subset = chr == gene_prot_file$CHR[i])
    }
    prot_chr <- as.numeric(gene_prot_file[i, "CHR"][[1]])
    prot_chr_bash <- ifelse(test = opt$study == "deCODE",
                            yes = paste0("chr", prot_chr),
                            no = prot_chr)
    if(grepl(pattern = "deCODE", x = opt$study)){
      # if(!(prot_chr_bash %in% unique(outcome_format$chr.outcome))){
      if(!(gsub(pattern = 'chr', replacement = '', x = prot_chr_bash) %in% unique(outcome_format$chr.outcome))){
        message(paste0("Protein located on chromosome ", prot_chr, ", which is absent from the outcome (exposure if reverse) GWAS. Skipping protein..."))
        next
      }
    } else {
      if(!(prot_chr %in% unique(outcome_format$chr.outcome))){
        message(paste0("Protein located on chromosome ", prot_chr, ", which is absent from the outcome (exposure if reverse) GWAS. Skipping protein..."))
        next
      }
    }
    # prot_chr <- gene_prot_file[i, "CHR"]
    if((gene_prot_file$start[i]-window) <= 0){
      window_start <- (gene_prot_file$start[i]) - 1
    } else {
      window_start <- window
    }
    if((gene_prot_file$end[i]+window) > chr_lengths[which(chr_lengths$Chr == prot_chr), "Total_length_db"]){
      window_end <- (chr_lengths[which(chr_lengths$Chr == prot_chr), "Total_length_db"]) - (gene_prot_file$end[i])
    } else {
      window_end <- window
    }
    if(opt$exp_is_vcf){
      gwas.vcf <- gwasvcf::query_gwas(vcf = paste0(gene_prot_file$path[i]), chrompos = paste0(prot_chr, ":", gene_prot_file$start[i]-window_start, "-", gene_prot_file$end[i]+window_end))
      # eval(parse(text = paste0("pwas_file$", opt$exp_chr_col, " <- sapply(strsplit(x = pwas_file$", opt$exp_chr_col, ", split = 'chr'), `[[`, 2)")))
      pwas_gene <- try(gwasglue::gwasvcf_to_TwoSampleMR(vcf = gwas.vcf, type = "exposure"))
      if(inherits(pwas_gene, "try-error")){
        #rm(list = c("pwas_gene"))
        message("No SNP in vcf file for specified region on specified chromosome. Skipping protein...")
        unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character(attr(x = pwas_gene, which = "condition"))))
        fwrite(x = unprocessed_proteins,
               file = unprocessed_proteins_file,
               sep = "\t", append = FALSE)
        next
      }
      rm(gwas.vcf)
      if(opt$exp_add_rsids){
        pwas_gene <- merge(x = pwas_gene, y = snps_ref_clean, by.x = c("chr.exposure", "pos.exposure"), by.y = c("chr","pos"), all = FALSE)
      }
      exposure <- subset(x = pwas_gene, subset = !is.na(SNP))
      exposure <- subset(x = exposure, subset = !grepl(pattern = ",", x = SNP))
      exposure_all <- exposure
      # Only if NOT performing reverse MR, keeping only SNPs with p < threshold in exposure. If performing reverse MR, keeping all SNPs.
      if(!opt$reverse){
        exposure <- subset(x = exposure, subset = pval.exposure <= pval_threshold)
      }
    } else {
      if(grepl(pattern = " ", x = opt$exp_pval_col) | grepl(pattern = "[(]", x = opt$exp_pval_col)){
        opt$exp_pval_col <- paste0("`", opt$exp_pval_col, "`")
      }
      if((opt$study == "INTERVAL") | grepl(pattern = "INTERVAL", x = opt$study)){
        pwas_file <- as.data.frame(fread(file = paste0(gene_prot_file$path[i], gene_prot_file$CHR[i],"_meta_final_v1.tsv"), header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
        pwas_gene <- subset(x = pwas_file, subset = (eval(parse(text = opt$exp_pos_col)) >= ((gene_prot_file$start[i]) - window_start)) & (eval(parse(text = opt$exp_pos_col)) <= ((gene_prot_file$end[i]) + window_end)))
        if(nrow(pwas_gene) == 0){
          message(" No SNP found for specified region in protein GWAS file. Skipping protein...")
          rm(list = c("pwas_file", "pwas_gene", "exposure"))
          unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("Exposure : No SNP found for specified region in GWAS file")))
          fwrite(x = unprocessed_proteins,
                 file = unprocessed_proteins_file,
                 sep = "\t", append = FALSE)
          next
        }
        pwas_gene <- subset(x = pwas_gene, select = c(opt$exp_chr_col, opt$exp_pos_col, opt$exp_ea_col, opt$exp_oa_col, opt$exp_beta_col, opt$exp_se_col, opt$exp_pval_col))
        if(opt$exp_add_rsids){
          pwas_gene <- merge(x = pwas_gene, y = snps_ref_clean, by.x = c(opt$exp_chr_col, opt$exp_pos_col), by.y = c("chr","pos"), all = FALSE)
        }
      } else {
        exp_samplesize <- opt$exp_ss_total
        if((opt$study == "ARIC") | grepl(pattern = "ARIC", x = opt$study)){
          pwas_file <- as.data.frame(fread(file = paste0(gene_prot_file$path[i]), header = TRUE, stringsAsFactors = FALSE, nThread = opt$nthreads))
        } else {
          if(!dir.exists(paste0("./temp_files/", opt$study))){
            dir.create(path = paste0("./temp_files/", opt$study))
          }
          word1 <- paste0("bash")
          args1 <- paste0(" ./extract_chr_pos.sh", 
                          " -p ", gene_prot_file$prot_name[i], 
                          " -f ", gene_prot_file$path[i], 
                          " -c ", prot_chr_bash, 
                          " -s ", gene_prot_file$start[i]-window_start, 
                          " -e ", gene_prot_file$end[i]+window_end,
                          " -t ", opt$exp_is_vcf,
                          " -n ", opt$study,
                          " -A ", opt$exp_snp_col_num,
                          " -B ", opt$exp_chr_col_num,
                          " -C ", opt$exp_pos_col_num, 
                          " -D ", opt$exp_vcf_skip)
          system2(command = word1, args = args1)
          pwas_file <- fread(file = paste0("./temp_files/", opt$study, "/", gene_prot_file$prot_name[i], ".region.txt"), header = TRUE, stringsAsFactors = FALSE)
          unlink(paste0("./temp_files/", opt$study, "/", gene_prot_file$prot_name[i], ".region.txt"))
        }
        if((opt$study == "ARIC") | grepl(pattern = "ARIC", x = opt$study)){
          pwas_gene <- subset(x = pwas_file, subset = (eval(parse(text = paste0("`", opt$exp_chr_col, "`"))) == prot_chr) & (eval(parse(text = opt$exp_pos_col)) >= ((gene_prot_file$start[i]) - window_start)) & (eval(parse(text = opt$exp_pos_col)) <= ((gene_prot_file$end[i]) + window_end)))
          if(nrow(pwas_gene) == 0){
            message(" No SNP found for specified region in protein GWAS file. Skipping protein...")
            rm(list = c("pwas_file", "pwas_gene", "exposure"))
            unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("Exposure : No SNP found for specified region in GWAS file")))
            fwrite(x = unprocessed_proteins,
                   file = unprocessed_proteins_file,
                   sep = "\t", append = FALSE)
            next
          }
          pwas_gene <- subset(x = pwas_gene, select = c(opt$exp_chr_col, opt$exp_pos_col, opt$exp_ea_col, opt$exp_oa_col, opt$exp_beta_col, opt$exp_se_col, opt$exp_pval_col, opt$exp_snp_col, opt$exp_eaf_col, opt$exp_ss_col))
          eval(parse(text = paste0("pwas_gene <- ", "subset(x = pwas_gene, subset = `", opt$exp_chr_col, "` != 'X')")))
          eval(parse(text = paste0("pwas_gene$`", opt$exp_chr_col, "`<- as.numeric(pwas_gene$`", opt$exp_chr_col, "`)")))
        } else {
          pwas_gene <- subset(x = pwas_file, subset = (eval(parse(text = opt$exp_chr_col)) == prot_chr_bash) & (eval(parse(text = opt$exp_pos_col)) >= ((gene_prot_file$start[i]) - window_start)) & (eval(parse(text = opt$exp_pos_col)) <= ((gene_prot_file$end[i]) + window_end)))
          if(nrow(pwas_gene) == 0){
            message(" No SNP found for specified region in protein GWAS file. Skipping protein...")
            rm(list = c("pwas_file", "pwas_gene", "exposure"))
            unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("Exposure : No SNP found for specified region in GWAS file")))
            fwrite(x = unprocessed_proteins,
                   file = unprocessed_proteins_file,
                   sep = "\t", append = FALSE)
            next
          }
          pwas_gene <- subset(x = pwas_gene, select = c(opt$exp_chr_col, opt$exp_pos_col, opt$exp_ea_col, opt$exp_oa_col, opt$exp_beta_col, opt$exp_se_col, opt$exp_pval_col, opt$exp_snp_col, opt$exp_eaf_col, opt$exp_ss_col))
          pwas_gene <- as.data.frame(pwas_gene)
          if(grepl(pattern = "chr", x = prot_chr_bash)){
            eval(parse(text = paste0("pwas_gene$", opt$exp_chr_col, " <- as.numeric(sapply(strsplit(x = pwas_gene$", opt$exp_chr_col, ", split = 'chr'), `[[`, 2))")))
          }
          eval(parse(text = paste0("pwas_gene <- ", "subset(x = pwas_gene, subset = ", opt$exp_chr_col, " != 'X')")))
          eval(parse(text = paste0("pwas_gene$", opt$exp_chr_col, "<- as.numeric(pwas_gene$", opt$exp_chr_col, ")")))
        }
        if(opt$exp_add_rsids){
          pwas_gene <- merge(x = pwas_gene, y = snps_ref_clean, by.x = c(opt$exp_chr_col, opt$exp_pos_col), by.y = c("chr","pos"), all = FALSE)
          opt$exp_snp_col <- "SNP"
        }
        if((opt$study == "ARIC") | grepl(pattern = "ARIC", x = opt$study) & (opt$exp_check_alleles)){
          # In ARIC, BETA and FREQ_A1 refer to A1, but A1 can be REF or ALT. Making sure the right reference allele is used.
          pwas_gene$temp.ref <- ifelse(test = eval(parse(text = paste0("pwas_gene$", opt$exp_ea_col))) == pwas_gene$A1,
                                       yes = eval(parse(text = paste0("pwas_gene$", opt$exp_ea_col))),
                                       no = eval(parse(text = paste0("pwas_gene$", opt$exp_oa_col))))
          pwas_gene$temp.alt <- ifelse(test = eval(parse(text = paste0("pwas_gene$", opt$exp_ea_col))) != pwas_gene$A1,
                                       yes = eval(parse(text = paste0("pwas_gene$", opt$exp_ea_col))),
                                       no = eval(parse(text = paste0("pwas_gene$", opt$exp_oa_col))))
          eval(parse(text = paste0("pwas_gene$", opt$exp_ea_col, " <- pwas_gene$temp.ref")))
          eval(parse(text = paste0("pwas_gene$", opt$exp_oa_col, " <- pwas_gene$temp.alt")))
        }
      }
      exposure <- subset(x = pwas_gene, subset = !is.na(eval(parse(text = opt$exp_snp_col))))
      exposure <- subset(x = exposure, subset = !grepl(pattern = ",", x = eval(parse(text = opt$exp_snp_col))))
      if(opt$exp_p_is_log){
        eval(parse(text = paste0("exposure$", opt$exp_pval_col, "<- 10**(exposure$", opt$exp_pval_col, ")")))
      }
      exposure_all <- exposure
      # Only if NOT performing reverse MR, keeping only SNPs with p < threshold in exposure. If performing reverse MR, keeping all SNPs.
      if(!opt$reverse){
        exposure <- subset(x = exposure, subset = eval(parse(text = opt$exp_pval_col)) <= pval_threshold)
      }
    }
    if(nrow(exposure) == 0){
      message(paste0("\n No SNPs remaining after pval filtering for ", gene_prot_file$prot_name[i], ". Skipping protein...\n"))
      rm(list = c("pwas_file", "pwas_gene", "exposure"))
      unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("Exposure : No SNPs remaining after pval filtering")))
      fwrite(x = unprocessed_proteins,
             file = unprocessed_proteins_file,
             sep = "\t", append = FALSE)
      next
    }
    exposure$phenotype <- gene_prot_file$prot_name[i]
    # exposure$samplesize <- unique(eval(parse(text = paste0("exposure$", opt$exp_ss_col))))
    ############################################################## #
    
    
    ############################################################## #
    ## __9.2. Format exposure ####
    if(opt$exp_is_vcf){
      exposure_format <- exposure
      exp_samplesize <- round(mean(exposure_format$samplesize.exposure), digits = 0)
      exposure_format$exposure <- exposure_format$phenotype
    } else {
      command_format <- paste0("exposure_format <- TwoSampleMR::format_data(dat = exposure, type = 'exposure'",
                               ", snp_col = '", opt$exp_snp_col,
                               "', beta_col = '", opt$exp_beta_col,
                               "', se_col = '", opt$exp_se_col,
                               "', effect_allele_col = '", opt$exp_ea_col,
                               "', other_allele_col = '", opt$exp_oa_col,
                               "', pval_col = '", opt$exp_pval_col,
                               "', phenotype_col = 'phenotype'")
      if(!is.na(opt$exp_pos_col)){
        command_format <- paste0(command_format, ", pos_col = '", opt$exp_pos_col, "'")
      }
      if(!is.na(opt$exp_chr_col)){
        command_format <- paste0(command_format, ", chr_col = '", opt$exp_chr_col, "'")
      }
      if(!is.na(opt$exp_eaf_col)){
        command_format <- paste0(command_format, ", eaf_col = '", opt$exp_eaf_col, "'")
      }
      if(!is.na(opt$exp_ss_col)){
        command_format <- paste0(command_format, ", samplesize_col = '", opt$exp_ss_col, "'")
      }
      command_format <- paste0(command_format, ")")
      eval(parse(text = command_format))
      rm(list = c("pwas_file"))
    }
    # If performing reverse MR, changing exposure to outcome and vice-versa.
    if(opt$reverse){
      exp_all.temp <- outcome_all.clean
      colnames(exp_all.temp) <- gsub(pattern = 'outcome', replacement = 'exposure', colnames(exp_all.temp))
      exp_format.temp <- outcome_format.clean
      colnames(exp_format.temp) <- gsub(pattern = 'outcome', replacement = 'exposure', colnames(exp_format.temp))
      exp.temp <- outcome_all.clean
      colnames(exp.temp) <- gsub(pattern = 'outcome', replacement = 'exposure', colnames(exp.temp))
      
      out_all.temp <- exposure_format
      colnames(out_all.temp) <- gsub(pattern = 'exposure', replacement = 'outcome', colnames(out_all.temp))
      out_format.temp <- exposure_format
      colnames(out_format.temp) <- gsub(pattern = 'exposure', replacement = 'outcome', colnames(out_format.temp))
      
      exposure_all <- exp_all.temp
      exposure_format <- exp_format.temp
      exposure <- exp.temp
      
      outcome_all <- out_all.temp
      outcome_format <- out_format.temp
      
      rm(list = c("exp_all.temp", "exp_format.temp", "exp.temp", "out_all.temp", "out_format.temp"))
    }
    
    ## Keeping only SNPs with maf >= 1%
    exposure_format <- subset(x = exposure_format, subset = eaf.exposure >= 0.01)
    ############################################################## #
    
    
    ############################################################## #
    ## __9.3. Harmonizing ####
    message("\n\t Performing harmonization...\n")
    dat_all <- TwoSampleMR::harmonise_data(exposure_dat = exposure_format,
                                           outcome_dat = outcome_format)
    if(!opt$reverse){
      dat_all$pos.outcome <- dat_all$pos.exposure
    }
    else {
      dat_all$pos.exposure <- dat_all$pos.outcome
    }
    
    
    if(nrow(dat_all)==0) {
      message("\t\t *** WARNING : no SNPs in the european reference panel. \n\t\t NO CLUMPING WILL BE MADE ***")
      message("No SNP found in reference panel")
      rm(list = c("pwas_gene", "exposure"))
      unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("Harmonized Data : No SNPs in the european reference panel")))
      fwrite(x = unprocessed_proteins,
             file = unprocessed_proteins_file,
             sep = "\t", append = FALSE)
      next
    }
    ############################################################## #
    
    
    ############################################################## #
    ## __9.4. Clumping ####
    message("\n\t Performing clumping...\n")
    out_clump = try(ld_clump(data.frame(rsid=dat_all$SNP, pval=dat_all$pval.exposure),
                             clump_kb=(window/1000),
                             clump_r2=clump_r2,
                             plink_bin=genetics.binaRies::get_plink_binary(),
                             bfile=opt$plink_binaries))
    if(inherits(out_clump, "try-error")){
      out_clump = try(TwoSampleMR::clump_data(dat = dat_all, clump_kb = (window/1000), clump_r2 = clump_r2, pop = "EUR"))
      if(inherits(out_clump, "try-error")){
        rm(list = c("pwas_gene", "exposure"))
        unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("Clumped Data : No SNPs remaining after clumping")))
        fwrite(x = unprocessed_proteins,
               file = unprocessed_proteins_file,
               sep = "\t", append = FALSE)
        next
      }
      dat = dat_all[which(dat_all$SNP %in% out_clump$SNP),]
    } else {
      dat = dat_all[which(dat_all$SNP %in% out_clump$rsid),]
    }
    message("\n\tClumping done!\n")
    
    #### ____9.4.1. Saving MR dat (optional) ####
    if(opt$save_dat){
      fwrite(x = dat, file = paste0(outname_dat, gene_prot_file$prot_name[i], ".", gene_prot_file$id[i], "_SNPs_p.", pval_threshold, "_r2.", clump_r2,".txt"), sep = "\t")
    }
    ############################################################## #
    
    
    ############################################################## #
    ## __9.5. MR ####
    if(nrow(dat) >= 3){
      # mr_prot <- mr(dat = dat, method_list = c("mr_ivw"))
      mr_prot <- mr(dat = dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
      if(nrow(mr_prot) == 0){
        unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("MR : No SNPs remaining for MR")))
        message("No SNPs available for MR for ", gene_prot_file$prot_name[i], "\n")
        fwrite(x = unprocessed_proteins,
               file = unprocessed_proteins_file,
               sep = "\t", append = FALSE)
        next
      }
      CIs <- generate_odds_ratios(mr_prot)
      if(sum(grepl(pattern = "Inverse variance weighted", x = CIs$method)) >= 1){
        mr_prot$b_CIlow <- subset(x = CIs, subset = method == "Inverse variance weighted", select = "lo_ci")[[1]]
        mr_prot$b_CIhigh <- subset(x = CIs, subset = method == "Inverse variance weighted", select = "up_ci")[[1]]
      } else {
        mr_prot$b_CIlow <- NA
        mr_prot$b_CIhigh <- NA
      }
      if(sum(grepl(pattern = "Weighted median", x = CIs$method)) >= 1){
        mr_prot$b_weightedmed_CIlow <- subset(x = CIs, subset = method == "Weighted median", select = "lo_ci")[[1]]
        mr_prot$b_weightedmed_CIhigh <- subset(x = CIs, subset = method == "Weighted median", select = "up_ci")[[1]]
      } else {
        mr_prot$b_weightedmed_CIlow <- NA
        mr_prot$b_weightedmed_CIhigh <- NA
      }
      if(sum(grepl(pattern = "Weighted median", x = CIs$method)) >= 1){
        mr_prot$b_egger_CIlow <- subset(x = CIs, subset = method == "MR Egger", select = "lo_ci")[[1]]
        mr_prot$b_egger_CIhigh <- subset(x = CIs, subset = method == "MR Egger", select = "up_ci")[[1]]
      } else {
        mr_prot$b_egger_CIlow <- NA
        mr_prot$b_egger_CIhigh <- NA
      }
      mr_egg <- mr_egger_regression(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome)
      if(length(mr_egg) == 0){
        mr_prot$nsnp_egger_i <- NA
        mr_prot$b_egger_i <- NA
        mr_prot$se_egger_i <- NA
        mr_prot$pval_egger_i <- NA
      } else {
        egger_inter <- data.frame(nsnp = mr_egg$nsnp, b = mr_egg$b_i, se = mr_egg$se_i, pval = mr_egg$pval_i)
        mr_prot$nsnp_egger_i <- egger_inter$nsnp
        mr_prot$b_egger_i <- egger_inter$b
        mr_prot$se_egger_i <- egger_inter$se
        mr_prot$pval_egger_i <- egger_inter$pval
      }
      mr_het <- mr_heterogeneity(dat)
      if(nrow(mr_het) == 0){
        mr_prot$Q_IVW <- NA
        mr_prot$P_Q_IVW <- NA
        mr_prot$Q_Egger <- NA
        mr_prot$P_Q_Egger <- NA
      } else {
        if(sum(grepl(pattern = "Inverse variance weighted", x = mr_het$method)) >= 1){
          mr_prot$Q_IVW <- subset(x = mr_het, subset = method == "Inverse variance weighted", select = "Q")[[1]]
          mr_prot$P_Q_IVW <- subset(x = mr_het, subset = method == "Inverse variance weighted", select = "Q_pval")[[1]]
        } else {
          mr_prot$Q_IVW <- NA
          mr_prot$P_Q_IVW <- NA
        }
        if(sum(grepl(pattern = "MR Egger", x = mr_het$method)) >= 1){
          mr_prot$Q_Egger <- subset(x = mr_het, subset = method == "MR Egger", select = "Q")[[1]]
          mr_prot$P_Q_Egger <- subset(x = mr_het, subset = method == "MR Egger", select = "Q_pval")[[1]]
        } else {
          mr_prot$Q_Egger <- NA
          mr_prot$P_Q_Egger <- NA
        }
      }
      mr_prot.clean <- mr_prot[1,]
      if(sum(grepl(pattern = "MR Egger", x = mr_prot$method)) >= 1){
        mr_prot.clean$b_egger <- subset(x = mr_prot, subset = method == "MR Egger", select = "b")[[1]]
        mr_prot.clean$se_egger <- subset(x = mr_prot, subset = method == "MR Egger", select = "se")[[1]]
        mr_prot.clean$pval_egger <- subset(x = mr_prot, subset = method == "MR Egger", select = "pval")[[1]]
      } else {
        mr_prot.clean$b_egger <- NA
        mr_prot.clean$se_egger <- NA
        mr_prot.clean$pval_egger <-NA
      }
      if(sum(grepl(pattern = "Weighted median", x = mr_prot$method)) >= 1){
        mr_prot.clean$b_weighted_median <- subset(x = mr_prot, subset = method == "Weighted median", select = "b")[[1]]
        mr_prot.clean$se_weighted_median <- subset(x = mr_prot, subset = method == "Weighted median", select = "se")[[1]]
        mr_prot.clean$pval_weighted_median <- subset(x = mr_prot, subset = method == "Weighted median", select = "pval")[[1]]
      } else {
        mr_prot.clean$b_weighted_median <- NA
        mr_prot.clean$se_weighted_median <- NA
        mr_prot.clean$pval_weighted_median <- NA
      }
      mr_prot <- mr_prot.clean
    }
    if(nrow(dat) == 2){
      mr_prot <- mr(dat = dat)
      if(nrow(mr_prot) == 0){
        message("No SNPs available for MR for ", gene_prot_file$prot_name[i], "\n")
        unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("MR : No SNPs remaining for MR")))
        fwrite(x = unprocessed_proteins,
               file = unprocessed_proteins_file,
               sep = "\t", append = FALSE)
        next
      }
      CIs <- generate_odds_ratios(mr_prot)
      mr_prot$b_CIlow <- subset(x = CIs, subset = method %in% c("Inverse variance weighted", "Wald ratio"), select = "lo_ci")[[1]]
      mr_prot$b_CIhigh <- subset(x = CIs, subset = method %in% c("Inverse variance weighted", "Wald ratio"), select = "up_ci")[[1]]
      mr_het <- mr_heterogeneity(dat)
      if(nrow(mr_het) == 0){
        mr_prot$Q_Egger <- NA
        mr_prot$P_Q_Egger <- NA
        mr_prot$Q_IVW <- NA
        mr_prot$P_Q_IVW <- NA
      } else {
        mr_prot$Q_Egger <- NA
        mr_prot$P_Q_Egger <- NA
        mr_prot$Q_IVW <- mr_het$Q[[1]]
        mr_prot$P_Q_IVW <- mr_het$Q_pval[[1]]
      }
    }
    if(nrow(dat) == 1){
      mr_wald <- mr_wald_ratio(b_exp = dat$beta.exposure, b_out = dat$beta.outcome, se_exp = dat$se.exposure, se_out = dat$se.outcome)
      if(length(mr_wald) == 0){
        message("No SNPs available for MR for ", gene_prot_file$prot_name[i], "\n")
        unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("MR : No SNPs remaining for MR")))
        fwrite(x = unprocessed_proteins,
               file = unprocessed_proteins_file,
               sep = "\t", append = FALSE)
        next
      }
      CIs <- generate_odds_ratios(mr_wald)
      mr_prot <- data.frame(outcome=outcome_format$outcome[1], exposure=gene_prot_file$prot_name[i], method="Wald ratio", nsnp = mr_wald$nsnp, b = mr_wald$b, b_CIlow = CIs$lo_ci, b_CIhigh = CIs$up_ci, se = mr_wald$se, pval = mr_wald$pval)
      rm(mr_wald)
      mr_prot$Q_Egger <- NA
      mr_prot$P_Q_Egger <- NA
      mr_prot$Q_IVW <- NA
      mr_prot$P_Q_IVW <- NA
    }
    if(nrow(dat) == 0){
      message("\n\tNo SNPs remaining after MR. Skipping protein...\n")
      rm(list = c("pwas_file", "pwas_gene", "merged_file", "exposure"))
      unprocessed_proteins <- rbind(unprocessed_proteins, data.frame(prot = gene_prot_file$prot_name[i], reason = as.character("MR : No SNPs remaining for MR")))
      fwrite(x = unprocessed_proteins,
             file = unprocessed_proteins_file,
             sep = "\t", append = FALSE)
      next
    }
    
    if((nrow(dat) == 2) | (nrow(dat) == 1)){
      mr_prot$b_egger <- NA
      mr_prot$se_egger <- NA
      mr_prot$pval_egger <- NA
      
      mr_prot$nsnp_egger_i <- NA
      mr_prot$b_egger_i <- NA
      mr_prot$se_egger_i <- NA
      mr_prot$pval_egger_i <- NA
      
      mr_prot$b_weighted_median <- NA
      mr_prot$se_weighted_median <- NA
      mr_prot$pval_weighted_median <- NA
      
      mr_prot$b_egger_CIlow <- NA
      mr_prot$b_egger_CIhigh <- NA
      mr_prot$b_weightedmed_CIlow <- NA
      mr_prot$b_weightedmed_CIhigh <- NA
    }
    ############################################################## #
    
    
    ############################################################## #
    #### __9.6. Analyses including correlation correction (JAM, PCA, CoJo) ####
    jammr_analyses_okay <- FALSE
    pca_analyses_okay <- FALSE
    cojo_analyses_okay <- FALSE
    # If IVW p-value for LD-pruning is not below 0.05, these analyses will be skipped to avoid useless waste of computing resources.
    if(mr_prot$pval < 0.05){
      out_clump_corr = try(ld_clump(data.frame(rsid=dat_all$SNP, pval=dat_all$pval.exposure),
                                    clump_kb=(window/1000),
                                    clump_r2=opt$r2_corr,
                                    plink_bin=genetics.binaRies::get_plink_binary(),
                                    bfile=opt$plink_binaries))
      if(inherits(out_clump_corr, "try-error")){
        out_clump_corr = try(TwoSampleMR::clump_data(dat = dat_all, clump_kb = (window/1000), clump_r2 = opt$r2_corr, pop = "EUR"))
        if(inherits(out_clump_corr, "try-error")){
          rm(list = c("exposure"))
        }
      } else {
        dat_corr <- subset(x = dat_all, subset = SNP %in% out_clump_corr$rsid)
        
        #### ____9.6.1. JAM ####
        
        jammr.results <- try(JAMMR(bx = dat_corr$beta.exposure, sx = dat_corr$se.exposure, by = dat_corr$beta.outcome, sy = dat_corr$se.outcome,
                                   N1 = exp_samplesize, eafs = dat_corr$eaf.exposure, w = c(0, exp_samplesize*10^c(seq(from = -2, to = 2, by = 1))), jam.seed = 4189))
        # jammr.results <- try(JAMMR(bx = dat_corr$beta.exposure, sx = dat_corr$se.exposure, by = dat_corr$beta.outcome, sy = dat_corr$se.outcome,
        #                            N1 = exp_samplesize, eafs = dat_corr$eaf.exposure, n.grid = 6, grid.limits = c(100, 100*exp_samplesize), jam.seed = 4189))
        if(inherits(jammr.results, "try-error")){
          message("JAMMR could not be performed.")
          jammr_analyses_okay <- FALSE
        } else {
          jammr_analyses_okay <- TRUE
        }
        
        #### ____9.6.2. PCA ####
        ## IVW estimate (accounting for correlation) using principal components
        rho <- try(ieugwasr::ld_matrix(variants = dat_corr$SNP, 
                                       with_alleles = FALSE, 
                                       plink_bin = genetics.binaRies::get_plink_binary(),
                                       bfile = opt$plink_binaries))
        if(inherits(rho, "try-error")){
          message("PCA could not be performed.")
          pca_analyses_okay <- FALSE
        } else {
          pca_analyses_okay <- TRUE
          if(nrow(rho) != nrow(dat_corr)){
            pca_analyses_okay <- FALSE
          } else {
            Phi <- (dat_corr$beta.exposure/dat_corr$se.outcome)%o%(dat_corr$beta.exposure/dat_corr$se.outcome)*rho
            # summary(prcomp(Phi, scale=FALSE))
            # K is number of principal components to include in analysis. Here, we include PC's that explain > 99% of the variance
            K <- which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]
            if(!is.na(K)){
              betaXG0 = as.numeric(dat_corr$beta.exposure%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
              betaYG0 = as.numeric(dat_corr$beta.outcome%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
              Omega = dat_corr$se.outcome%o%dat_corr$se.outcome*rho
              pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]
              beta_IVWcorrel.pc = solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
              se_IVWcorrel.fixed.pc = sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0))
            } else {
              pca_analyses_okay <- FALSE
            }
          }
        }
        ## IVW estimate (accounting for correlation) using the standard linear regression command after weighting the data by the Cholesky decomposition
        # Omega <- dat_corr$se.outcome%o%dat_corr$se.outcome*rho
        # c_betaXG = solve(t(chol(Omega)))%*%dat_corr$beta.exposure
        # c_betaYG = solve(t(chol(Omega)))%*%dat_corr$beta.outcome
        # beta_IVWcorrel <- lm(c_betaYG~c_betaXG-1)$coef[1]
        # se_IVWcorrel.fixed  = sqrt(1/(t(dat_corr$beta.exposure)%*%solve(Omega)%*%dat_corr$beta.exposure))
        # se_IVWcorrel.random = sqrt(1/(t(dat_corr$beta.exposure)%*%solve(Omega)%*%dat_corr$beta.exposure))*max(summary(lm(c_betaYG~c_betaXG-1))$sigma,1)
        
        #### ____9.6.3. CoJo ####
        ## IVW estimate (accounting for correlation) using conditional joint analysis (CoJo)
        for(cojo_directory in c("input", "output", "temp_files")){
          if(!dir.exists(paste0("./cojo/", cojo_directory, "/", opt$study))){
            dir.create(paste0("./cojo/", cojo_directory, "/", opt$study))
          }
        }
        
        args_cojo = paste0(" ./cojo/cojo_command.sh ", 
                           " -p ", gene_prot_file$prot_name[i], 
                           " -f ", gene_prot_file$path[i], 
                           " -c ", prot_chr_bash, 
                           " -s ", gene_prot_file$start[i]-window_start, 
                           " -e ", gene_prot_file$end[i]+window_end,
                           " -t ", opt$exp_is_vcf,
                           " -l ", opt$exp_p_is_log,
                           " -n ", opt$study,
                           " -i ", opt$cojo_ldref_rsids,
                           " -A ", opt$exp_cojo_snp_col_num,
                           " -B ", opt$exp_a1_col_num,
                           " -C ", opt$exp_a2_col_num,
                           " -D ", opt$exp_maf_col_num,
                           " -E ", opt$exp_beta_col_num,
                           " -F ", opt$exp_se_col_num,
                           " -G ", opt$exp_p_col_num,
                           " -H ", opt$exp_n_col_num,
                           " -I ", opt$exp_vcf_skip,
                           " -J ", opt$exp_form_col_num,
                           " -K ", opt$cojo_ldref,
                           " -M ", opt$exp_chr_col_num,
                           " -N ", opt$exp_pos_col_num)
        
        system2(command = "bash", 
                args = args_cojo)
        if(file.exists(paste0("./cojo/output/", opt$study, "/", gene_prot_file$prot_name[i], ".jma.cojo"))){
          cojo_analyses_okay <- TRUE
        } else {
          if(file.exists(paste0("./cojo/input/", opt$study, "/", gene_prot_file$prot_name[i],".txt"))){
            unlink(paste0("./cojo/input/", opt$study, "/", gene_prot_file$prot_name[i],".txt"))
          }
          if(opt$cojo_ldref_rsids == 0){
            system2(command = "bash", 
                    args = paste0(args_cojo, " -r 1"))
          }
          if(file.exists(paste0("./cojo/output/", opt$study, "/", gene_prot_file$prot_name[i], ".jma.cojo"))){
            cojo_analyses_okay <- TRUE
          } else {
            message("CoJo was performed, but JMA files were not created either due to a crash or a lack of SNPs. Check that alleles in GWAS file SNP names (chr:pos:a1:a2) are in the same order as in the CoJo LD reference file.")
            cojo_analyses_okay <- FALSE
          }
        }
      }
    }
    ############################################################## #
    
    
    ############################################################## #
    #### __9.7. Finalizing results according to IVW (LD-pruning) and correlation correction analyses (JAM, PCA, CoJo) outputs ####
    if(jammr_analyses_okay){
      mr_prot$jammr.b <- jammr.results$causal
      mr_prot$jammr.se <- jammr.results$se
      mr_prot$jammr.nsnps <- length(jammr.results$snp.probs)
    } else {
      mr_prot$jammr.b <- NA
      mr_prot$jammr.se <- NA
      mr_prot$jammr.nsnps <- NA
    }
    if(pca_analyses_okay){
      mr_prot$pca.b <- beta_IVWcorrel.pc
      mr_prot$pca.se <- se_IVWcorrel.fixed.pc
      mr_prot$pca.nsnps <- nrow(rho)
      mr_prot$pca.npc <- K
    } else {
      mr_prot$pca.b <- NA
      mr_prot$pca.se <- NA
      mr_prot$pca.nsnps <- NA
      mr_prot$pca.npc <- NA
    }
    if(cojo_analyses_okay){
      for(cojo_type in c("slct", "cond")){
        cojo_file <- ifelse(test = cojo_type == "slct",
                            yes = ".jma.cojo",
                            no = ".cma.cojo")
        if(file.exists(paste0("./cojo/output/", opt$study, "/", gene_prot_file$prot_name[i], cojo_file))){
          cojo_snps_list <- fread(file = paste0("./cojo/output/", opt$study, "/", gene_prot_file$prot_name[i], cojo_file), header = TRUE, stringsAsFactors = FALSE)
          dat_cojo <- dat_all
          if(sum(grepl(pattern = "rs", x = cojo_snps_list$SNP)) != 0){
            dat_cojo <- subset(x = dat_cojo, subset = SNP %in% cojo_snps_list$SNP)
          } else {
            dat_cojo$snp.name <- paste0(dat_cojo$chr.exposure, ":", dat_cojo$pos.exposure)
            cojo_snps_list$snp.name <- paste0(cojo_snps_list$Chr, ":", cojo_snps_list$bp)
            dat_cojo <- subset(x = dat_cojo, subset = snp.name %in% cojo_snps_list$snp.name)
          }
          if(nrow(dat_cojo) > 0){
            cojo_res <- try(mr(dat = dat_cojo))
            if(inherits(cojo_res, "try-error")){
              message("IVW on CoJo could not be performed.")
              eval(parse(text = paste0("mr_prot$cojo.b.", cojo_type, " <- NA")))
              eval(parse(text = paste0("mr_prot$cojo.se.", cojo_type, " <- NA")))
              eval(parse(text = paste0("mr_prot$cojo.p.", cojo_type, " <- NA")))
              eval(parse(text = paste0("mr_prot$cojo.nsnps.", cojo_type, " <- NA")))
            } else {
              if(sum(grepl(pattern = "Inverse variance weighted", x = cojo_res$method)) > 0){
                eval(parse(text = paste0("mr_prot$cojo.b.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Inverse variance weighted', select = 'b')[[1]]")))
                eval(parse(text = paste0("mr_prot$cojo.se.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Inverse variance weighted', select = 'se')[[1]]")))
                eval(parse(text = paste0("mr_prot$cojo.p.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Inverse variance weighted', select = 'pval')[[1]]")))
                eval(parse(text = paste0("mr_prot$cojo.nsnps.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Inverse variance weighted', select = 'nsnp')[[1]]")))
              } else {
                if(sum(grepl(pattern = "Wald ratio", x = cojo_res$method)) > 0){
                  eval(parse(text = paste0("mr_prot$cojo.b.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Wald ratio', select = 'b')[[1]]")))
                  eval(parse(text = paste0("mr_prot$cojo.se.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Wald ratio', select = 'se')[[1]]")))
                  eval(parse(text = paste0("mr_prot$cojo.p.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Wald ratio', select = 'pval')[[1]]")))
                  eval(parse(text = paste0("mr_prot$cojo.nsnps.", cojo_type, " <- subset(x = cojo_res, subset = method == 'Wald ratio', select = 'nsnp')[[1]]")))
                } else {
                  eval(parse(text = paste0("mr_prot$cojo.b.", cojo_type, " <- NA")))
                  eval(parse(text = paste0("mr_prot$cojo.se.", cojo_type, " <- NA")))
                  eval(parse(text = paste0("mr_prot$cojo.p.", cojo_type, " <- NA")))
                  eval(parse(text = paste0("mr_prot$cojo.nsnps.", cojo_type, " <- NA")))
                }
              }
            }
          } else {
            eval(parse(text = paste0("mr_prot$cojo.b.", cojo_type, " <- NA")))
            eval(parse(text = paste0("mr_prot$cojo.se.", cojo_type, " <- NA")))
            eval(parse(text = paste0("mr_prot$cojo.p.", cojo_type, " <- NA")))
            eval(parse(text = paste0("mr_prot$cojo.nsnps.", cojo_type, " <- NA")))
          }
        } else {
          eval(parse(text = paste0("mr_prot$cojo.b.", cojo_type, " <- NA")))
          eval(parse(text = paste0("mr_prot$cojo.se.", cojo_type, " <- NA")))
          eval(parse(text = paste0("mr_prot$cojo.p.", cojo_type, " <- NA")))
          eval(parse(text = paste0("mr_prot$cojo.nsnps.", cojo_type, " <- NA")))
        }
      }
      for(file_cojo in paste0("./cojo/output/", opt$study, "/", gene_prot_file$prot_name[i], c(".log", ".cma.cojo", ".jma.cojo", ".ldr.cojo", ".freq.badsnps", ".badsnps"))){
        if(file.exists(file_cojo)){
          unlink(file_cojo)
        }
      }
    } else {
      for(file_cojo in paste0("./cojo/output/", opt$study, "/", gene_prot_file$prot_name[i], c(".log", ".badsnps", ".freq.badsnps"))){
        if(file.exists(file_cojo)){
          unlink(file_cojo)
        }
      }
      mr_prot$cojo.b.slct <- NA
      mr_prot$cojo.se.slct <- NA
      mr_prot$cojo.p.slct <- NA
      mr_prot$cojo.nsnps.slct <- NA
      mr_prot$cojo.b.cond <- NA
      mr_prot$cojo.se.cond <- NA
      mr_prot$cojo.p.cond <- NA
      mr_prot$cojo.nsnps.cond <- NA
    }
    message("\n\tMR done!\n")
    mr_prot <- subset(x = mr_prot, select = c("outcome", "exposure", "method", "nsnp", 
                                              "b", "b_CIlow", "b_CIhigh", "se", "pval", 
                                              "b_egger", "b_egger_CIlow", "b_egger_CIhigh","se_egger", "pval_egger",
                                              "nsnp_egger_i", "b_egger_i", "se_egger_i", "pval_egger_i",
                                              "b_weighted_median", "b_weightedmed_CIlow", "b_weightedmed_CIhigh", "se_weighted_median", "pval_weighted_median",
                                              "Q_Egger", "P_Q_Egger", "Q_IVW", "P_Q_IVW", 
                                              "cojo.b.slct", "cojo.se.slct", "cojo.p.slct", "cojo.nsnps.slct", "cojo.b.cond", "cojo.se.cond", "cojo.p.cond", "cojo.nsnps.cond",
                                              "jammr.b", "jammr.se", "jammr.nsnps", "pca.b", "pca.se", "pca.nsnps", "pca.npc"
    ))
    
    mr_prot$exposure_exact <- gene_prot_file$id[i]
    mr_prot <- mr_prot[,c(1:2,43,3:42)]
    ############################################################## #
    
    
    ############################################################## #
    ## __9.8. Coloc ####
    exposure_coloc <- pwas_gene
    if(opt$exp_is_vcf){
      exposure_coloc <- subset(x = exposure_coloc, select = c("pval.exposure", "samplesize.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "SNP"))
    } else {
      exposure_coloc <- subset(x = exposure_coloc, select = c(opt$exp_pval_col, opt$exp_ss_col, opt$exp_eaf_col, opt$exp_beta_col, opt$exp_se_col, opt$exp_snp_col))
    }
    if(is.na(opt$exp_ss_col)){
      exposure_coloc$samplesize.exposure <- opt$exp_ss_total
    } else {
      exposure_coloc$samplesize.exposure <- eval(parse(text = paste0("exposure_coloc$", opt$exp_ss_col)))
    }
    outcome_coloc <- outcome_format
    outcome_coloc$samplesize.outcome <- ifelse(test = !opt$reverse, 
                                               yes = opt$out_ss_total,
                                               no = opt$exp_ss_total)
    outcome_coloc <- subset(x = outcome_coloc, select = c("pval.outcome", "samplesize.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "SNP"))
    if(opt$exp_is_vcf){
      exposure_coloc <- subset(x = exposure_coloc, subset = !is.na(SNP))
      exposure_coloc <- subset(x = exposure_coloc, subset = !grepl(pattern = ",", x = SNP))
      exposure_coloc <- subset(x = exposure_coloc, subset = !duplicated(SNP))
      exposure_coloc <- subset(x = exposure_coloc, subset = SNP %in% outcome_format$SNP)
    } else {
      exposure_coloc <- subset(x = exposure_coloc, subset = !is.na(eval(parse(text = opt$exp_snp_col))))
      exposure_coloc <- subset(x = exposure_coloc, subset = !grepl(pattern = ",", x = eval(parse(text = opt$exp_snp_col))))
      exposure_coloc <- subset(x = exposure_coloc, subset = !duplicated(eval(parse(text = opt$exp_snp_col))))
      exposure_coloc <- subset(x = exposure_coloc, subset = eval(parse(text = opt$exp_snp_col)) %in% outcome_format$SNP)
      eval(parse(text = paste0("exposure_coloc$eaf.exposure <- exposure_coloc$", opt$exp_eaf_col)))
      eval(parse(text = paste0("exposure_coloc$pval.exposure <- exposure_coloc$", opt$exp_pval_col)))
      eval(parse(text = paste0("exposure_coloc$beta.exposure <- exposure_coloc$", opt$exp_beta_col)))
      eval(parse(text = paste0("exposure_coloc$se.exposure <- exposure_coloc$", opt$exp_se_col)))
      eval(parse(text = paste0("exposure_coloc$SNP <- exposure_coloc$", opt$exp_snp_col)))
    }
    
    # exposure_coloc$chr.exposure <- as.numeric(exposure_coloc$chr.exposure)
    outcome_coloc <- subset(x = outcome_coloc, subset = SNP %in% exposure_coloc$SNP)
    if((nrow(exposure_coloc) == 0) | (nrow(outcome_coloc) == 0)){
      message(paste0("\n No SNPs remaining in exposure or outcome to perform coloc on ", gene_prot_file$prot_name[i], ". Coloc will not be performed...\n"))
      mr_prot$coloc.nsnps <- NA
      mr_prot$coloc.PPH0 <- NA
      mr_prot$coloc.PPH1 <- NA
      mr_prot$coloc.PPH2 <- NA
      mr_prot$coloc.PPH3 <- NA
      mr_prot$coloc.PPH4 <- NA
    } else {
      coloc_type_exp <- ifelse(test = !is.na(opt$exp_ncase),
                               yes = paste0("type = 'cc', s = (opt$exp_ncase / opt$exp_ss_total),"),
                               no = paste0("type = 'quant',"))
      coloc_type_out <- ifelse(test = !is.na(opt$out_ncase),
                               yes = paste0("type = 'cc', s = (opt$out_ncase / opt$out_ss_total),"),
                               no = paste0("type = 'quant',"))
      if(!opt$reverse){
        type_out.temp <- coloc_type_exp
        type_exp.temp <- coloc_type_out
        coloc_type_exp <- type_exp.temp
        coloc_type_out <- type_out.temp
        
        exposure_coloc <- exposure_format
        
        rm(list = c("type_out.temp", "type_exp.temp"))
      }
      eval(parse(text = paste0("dataset1 <- data.frame(pvalues = exposure_coloc$pval.exposure, 
                             N = exposure_coloc$samplesize.exposure, 
                             MAF = exposure_coloc$eaf.exposure, 
                             beta = exposure_coloc$beta.exposure, 
                             varbeta = exposure_coloc$se.exposure,",
                               coloc_type_exp,
                               "snp = exposure_coloc$SNP)")))
      eval(parse(text = paste0("dataset2 <- data.frame(pvalues = outcome_coloc$pval.outcome,
                             N = outcome_coloc$samplesize.outcome, 
                             MAF = outcome_coloc$eaf.outcome, 
                             beta = outcome_coloc$beta.outcome, 
                             varbeta = outcome_coloc$se.outcome, ", 
                               coloc_type_out, 
                               "snp = outcome_coloc$SNP)")))
      coloc_res = try(coloc::coloc.abf(dataset1 = dataset1, dataset2 = dataset2))
      if(inherits(coloc_res, "try-error")){
        message("Coloc did not work for protein.")
        mr_prot$coloc.nsnps <- NA
        mr_prot$coloc.PPH0 <- NA
        mr_prot$coloc.PPH1 <- NA
        mr_prot$coloc.PPH2 <- NA
        mr_prot$coloc.PPH3 <- NA
        mr_prot$coloc.PPH4 <- NA
      } else {
        mr_prot$coloc.nsnps <- coloc_res$summary[[1]]
        mr_prot$coloc.PPH0 <- coloc_res$summary[[2]]
        mr_prot$coloc.PPH1 <- coloc_res$summary[[3]]
        mr_prot$coloc.PPH2 <- coloc_res$summary[[4]]
        mr_prot$coloc.PPH3 <- coloc_res$summary[[5]]
        mr_prot$coloc.PPH4 <- coloc_res$summary[[6]]
        fwrite(x = coloc_res$results, file = paste0(outname_coloc, gene_prot_file$prot_name[i], ".txt"))
      }
    }
    ############################################################## #
    
    
    ############################################################## #
    ## __9.9. Saving results ####
    mr_res <- rbind(mr_res, mr_prot)
    if(((i %% opt$save_every) == 0) | save_res){
      fwrite(x = mr_res,
             file = outname,
             sep = "\t", append = FALSE)
      save_res <- FALSE
    }
    rm(list = c("pwas_gene", "exposure", "mr_prot", "dat", "exposure_format"))
    end <- Sys.time()
    print(end - start)
  }
  fwrite(x = mr_res,
         file = outname,
         sep = "\t", append = FALSE)
  fwrite(x = unprocessed_proteins,
         file = unprocessed_proteins_file,
         sep = "\t", append = FALSE)
  end_all <- Sys.time()
  message("Total time to run loop on all specified proteins : \n")
  print(end_all - start_all)
  rm(list=ls())
  message("==================== DONE ==================== \n")
  
}




