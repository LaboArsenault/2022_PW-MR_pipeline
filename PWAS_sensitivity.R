#!/usr/bin/env Rscript

library(hyprcoloc)
library(TwoSampleMR)
library(MRPRESSO)
library(MendelianRandomization)
library(data.table)
library(optparse)
{
  option_list = list(
    make_option("--wd", action="store", default="./", type='character',
                help=" working directory from which to access files"),
    make_option("--proteins", action="store", default=NA, type='character',
                help=" if specified, this script will use a protein names file to extract proteins. 
                * This file must contain the following protein information with these exact column names : 'name', 'id'. 
                * If running HyPrColoc, this file must contain the mandatory 2 additional columns with the exact same names : 'trait', 'path'
                Where path contains the absolute path to the tissue file. 
                      * You can include many 'path' columns (i.e. multiple traits), but they must all be in the same column, semicolon-separated (if using data from same MR results, write path to MR data)
                      * If you do include many traits, you must do the same thing with 'trait' column : all traits in the same column, semicolon-separated. 
                        If using MR data from same MR results, just add 'type@' in the trait name, where 'type' corresponds to type of data for the protein (either exposure or outcome) in MR data. 
                      * In that case, the additional GWASes included must only contain the following data (columns) with these column names : 'SNP', 'beta' 'se'
                * If you want to make a volcano plot, you should also include the 'beta' and 'pval' columns to use.
                
                This is to make sure to get the right protein if two or more proteins have the same name.
                Specifying this option makes the script ignore the 'pval', 'pval_col', 'use_adjp', 'adjp_method' and 'prot_col' options. 
                "),
    make_option("--mr_res", action="store", default=NA, type='character',
                help=" if specified, this script will use this MR results file to extract proteins. "),
    make_option("--pval", action="store", default=5e-08, type='numeric',
                help=" if 'proteins' argument is not provided, analyze proteins with p < pval (Default 5e-08)"),
    make_option("--pval_col", action="store", default="pval", type='character',
                help=" if 'proteins' argument is not provided, pval column to extract proteins (no matter the method) (Default 'pval'). 
                Note that this can be modified to use adjusted p-values instead. 
                In that case, the 'pval' option should also be modified to FDR of Bonferroni threshold (usually 0.05)."),
    make_option("--use_adjp", action="store", default=FALSE, type='logical',
                help=" Should you calculate adjusted p-values to select proteins on which to perform the analysis ? (Default FALSE).
                If this option is selected, make sure the 'pval' option is set to adjusted p-value threshold. "),
    make_option("--adjp_method", action="store", default="bonferroni", type='character',
                help=" If calculating adjusted p-values (use_adjp is TRUE), what method should you use ? (Default 'Bonferroni').
                For information on all methods, see 'p.adjust.methods' from the 'stats' package. "),
    make_option("--prot_col", action="store", default="exposure", type='character',
                help=" if 'proteins' argument is not provided, protein name column name (Default 'exposure')"),
    make_option("--id_col", action="store", default="exposure_exact", type='character',
                help=" if 'proteins' argument is not provided, protein id column name (Default 'exposure_exact')"),
    make_option("--beta_col", action="store", default="beta.exposure", type='character',
                help=" if 'proteins' argument is not provided, beta column name (Default 'b')"),
    make_option("--study", action="store", default="deCODE", type='character',
                help=" which study as exposure ? (INTERVAL, deCODE (default) or ARIC)"),
    make_option("--dat_dir", action="store", default="./results/MRdat/", type='character',
                help=" study directory where MR data are stored (Default to this script results folder : ./results/MRdat/)"),
    
    make_option("--hyprcoloc", action="store", default=FALSE, type='logical',
                help=" should you run HyPrColoc ? (default FALSE).
                This only works if you supply a protein list to the 'proteins' option. 
                If TRUE, see 'proteins' option requirements."),
    make_option("--volcano_plot_adjp", action="store", default=FALSE, type='logical',
                help=" Should the volcano plot y-axis correspond to adjusted p-values ? (default FALSE) "),
    
    make_option("--outcome_name", action="store", default="outcome", type='character',
                help=" outcome name (Default 'outcome')"),
    make_option("--volcano", action="store", default=TRUE, type='character',
                help=" should you produce a volcano plot of results ? (default TRUE)"),
    make_option("--tissue", action="store", default=TRUE, type='character',
                help=" should you produce a tissue enrichment plot of results ? (default TRUE)"),
    make_option("--out", action="store", default=NA, type='character',
                help=" output path and name (Default to this script results folder : ./results/sensitivity/*study*/sensitivity_analysis.*study*)")
  )
  opt = parse_args(OptionParser(option_list=option_list))
}
setwd(opt$wd)

############################################################################### #
#### 1. Parameters and selecting protein files ####
{
  hypercoloc <- opt$hyprcoloc
  volcano <- opt$volcano
  pval_col_noadj <- opt$pval_col
  if(is.na(opt$proteins)){
    # If proteins not provided, use MR results file and p-value threshold
    proteins <- as.data.frame(fread(file = opt$mr_res))
    proteins <- subset(x = proteins, select = c(opt$prot_col, opt$id_col, opt$beta_col, opt$pval_col))
    if(opt$use_adjp){
      proteins$padj <- p.adjust(p = eval(parse(text = paste0("proteins$", opt$pval_col))), method = opt$adjp_method)
      # eval(parse(text = paste0("proteins$", opt$pval_col, " <- proteins$padj")))
      opt$pval_col <- "padj"
      colnames(proteins) <- c("name", "id", "beta", "pval", "padj")
      proteins_all <- proteins
      proteins <- as.data.frame(subset(x = proteins, subset = padj < opt$pval))
    } else {
      colnames(proteins) <- c("name", "id", "beta", "pval")
      proteins_all <- proteins
      proteins <- as.data.frame(subset(x = proteins, subset = pval < opt$pval))
    }
    #### MAKE SURE THAT PVAL / PVAL_COL / PROT_COL / MR_RES file ARE ALSO SPECIFIED
  } else {
    proteins <- as.data.frame(fread(file = opt$proteins, header = TRUE, stringsAsFactors = FALSE))
    proteins_all <- proteins
    if((grepl(pattern = "trait", x = colnames(proteins), ignore.case = FALSE)) & (grepl(pattern = "path", x = colnames(proteins), ignore.case = FALSE))){
      hypercoloc <- TRUE
    } else {
      hypercoloc <- FALSE
    }
    if((grepl(pattern = "beta", x = colnames(proteins), ignore.case = FALSE)) & (grepl(pattern = "pval", x = colnames(proteins), ignore.case = FALSE))){
      volcano <- TRUE
    } else {
      volcano <- FALSE
    }
  }
  outname <- ifelse(test = is.na(opt$out),
                    yes = paste0("./results/sensitivity/", opt$study, "/sensitivity_analysis.", opt$study, ".txt"),
                    no = opt$out)
  outname_tissue_plot <- ifelse(test = is.na(opt$out),
                         yes = paste0("./results/sensitivity/", opt$study, "/sensitivity_analysis.tissue.", opt$study, ".png"),
                         no = opt$out)
  outname_volcano_plot <- ifelse(test = is.na(opt$out),
                                yes = paste0("./results/sensitivity/", opt$study, "/sensitivity_analysis.volcano.", opt$study, ".png"),
                                no = opt$out)
  if((!dir.exists(paste0("./results/sensitivity/", opt$study))) & ((is.na(opt$out)))){
    dir.create(path = paste0("./results/sensitivity/", opt$study))
  }
}

############################################################################### #


############################################################################### #
#### 2. Begin loop on each protein ####
{
  ## Preparing results file
  mr_res <- data.frame(matrix(nrow = 0, ncol = 18))
  colnames(mr_res) <- c("outcome", "exposure", 
                        "nsnp_conmix", "b_conmix", "se_conmix", "pval_conmix",
                        "nsnp_PRESSO", "b_PRESSO", "se_PRESSO", "pval_PRESSO", 
                        "nsnp_PRESSO_corr", "b_PRESSO_corr", "se_PRESSO_corr", "pval_PRESSO_corr",
                        "hyprcoloc_traits", "hyprcoloc_posterior_prob", "hyprcoloc_candidate_snp", "hyprcoloc_posterior_explained_by_snp")
  start <- Sys.time()
  for(i in 1:nrow(proteins)){
    message(paste0("\n\tProcessing protein ", proteins$name[i], " (", proteins$id[i], ") (num :", i,") ...\n"))
    
    
    ############################################################## #
    #### __2.1. Reading and preparing dat files ####
    {
      dat_dir_files <- list.files(path = opt$dat_dir)
      dat_file <- dat_dir_files[grepl(pattern = proteins$name[i], x = dat_dir_files)]
      dat_file <- dat_file[grepl(pattern = proteins$id[i], x = dat_file)]
      dat_file <- as.data.frame(fread(paste0(opt$dat_dir, "/", dat_file)))
      mr_prot <- data.frame(matrix(nrow = 1, ncol = 0))
      mr_prot$outcome <- opt$outcome_name
      mr_prot$exposure <- paste0(proteins$name[i], ".", proteins$id[i])
    }
    ############################################################## #
    
    
    ############################################################## #
    #### __2.2. Contamination Mixture ####
    {
      dat <- dat_file
      conmix_dat <- MendelianRandomization::mr_input(bx = dat$beta.exposure,
                                                     bxse = dat$se.exposure,
                                                     by = dat$beta.outcome,
                                                     byse = dat$se.outcome,
                                                     exposure = dat$exposure[[1]],
                                                     outcome = dat$outcome[[1]],
                                                     effect_allele = dat$effect_allele.exposure,
                                                     other_allele = dat$other_allele.exposure,
                                                     eaf = dat$eaf.exposure,
                                                     snps = dat$SNP)
      conmix_res <- try(MendelianRandomization::mr_conmix(conmix_dat, alpha = 0.05))
      if(inherits(conmix_res, "try-error")){
        message("\nError while running Contamination Mixture. Setting NA values for ConMix for this protein.\n")
        mr_prot$nsnp_conmix <- NA
        mr_prot$b_conmix <- NA
        mr_prot$se_conmix <- NA
        mr_prot$pval_conmix <- NA
      } else {
        conmix <- data.frame(nsnp = length(conmix_res@ValidSNPs), b = conmix_res@Estimate, se = (conmix_res@CIUpper[1] - conmix_res@CILower[1])/3.92, pval = conmix_res@Pvalue)
        mr_prot$nsnp_conmix <- conmix$nsnp
        mr_prot$b_conmix <- conmix$b
        mr_prot$se_conmix <- conmix$se
        mr_prot$pval_conmix <- conmix$pval
        # conmix_res_clean <- subset(dat, subset = dat$SNP %in% conmix_res@ValidSNPs)
        # conmix_res_clean <- subset(conmix_res_clean, subset = conmix_res_clean$pval.exposure <= opt$pval)
      }
    }
    ############################################################## #
    
    
    ############################################################## #
    #### __2.3. MR-PRESSO ####
    {
      dat <- dat_file
      if(nrow(dat) >= 4){
        presso_res <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                          BetaExposure = "beta.exposure",
                                          SdOutcome = "se.outcome",
                                          SdExposure = "se.exposure",
                                          OUTLIERtest = TRUE,
                                          DISTORTIONtest = TRUE,
                                          data = dat,
                                          NbDistribution = 10000,
                                          SignifThreshold = 0.05)
        out_presso <- data.frame(nsnp = ifelse(test = is.null(presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                                               yes = nrow(dat),
                                               no = (nrow(dat) - length(presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))),
                                 b = presso_res$`Main MR results`$`Causal Estimate`,
                                 se = presso_res$`Main MR results`$Sd,
                                 pval = presso_res$`Main MR results`$`P-value`)
        mr_prot$nsnp_PRESSO <- out_presso$nsnp[1]
        mr_prot$b_PRESSO <- out_presso$b[1]
        mr_prot$se_PRESSO <- out_presso$se[1]
        mr_prot$pval_PRESSO <- out_presso$pval[1]
        if(!is.null(presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)){
          mr_prot$nsnp_PRESSO_corr <- (out_presso$nsnp[1]) - (length(presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))
          fwrite(x = as.list(dat[presso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, "SNP"]), file = paste0(gsub(pattern = ".txt", replacement = "", x = outname), ".", proteins$name[i], ".pleiotropic.txt"))
        } else {
          mr_prot$nsnp_PRESSO_corr <- NA
        }
        mr_prot$b_PRESSO_corr <- out_presso$b[2]
        mr_prot$se_PRESSO_corr <- out_presso$se[2]
        mr_prot$pval_PRESSO_corr <- out_presso$pval[2]
      } else {
        mr_prot$nsnp_PRESSO <- NA
        mr_prot$b_PRESSO <- NA
        mr_prot$se_PRESSO <- NA
        mr_prot$pval_PRESSO <- NA
        mr_prot$nsnp_PRESSO_corr <- NA
        mr_prot$b_PRESSO_corr <- NA
        mr_prot$se_PRESSO_corr <- NA
        mr_prot$pval_PRESSO_corr <- NA
      }
    }
    ############################################################## #
    
    
    ############################################################## #
    #### __2.4. HyPrColoc ####
    {
      if(hypercoloc){
        dat_hyprcoloc <- dat
        dat_hyprcoloc <- subset(x = dat_hyprcoloc, select = c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome"))
        for(k in c("beta", "se")){
          for(l in c("exposure", "outcome")){
            if(l == 'outcome'){eval(parse(text = paste0("colnames(dat_hyprcoloc)[which(colnames(dat_hyprcoloc) == '", k, ".", l, "')] <- paste0('", k, ".', dat$outcome[[1]])")))} else {
              eval(parse(text = paste0("colnames(dat_hyprcoloc)[which(colnames(dat_hyprcoloc) == '", k, ".", l, "')] <- paste0('", k, ".', proteins$name[i])")))
            }
          }
        }
        if(!is.na(opt$proteins)){
          # n_traits <- length(colnames(proteins)[grepl(pattern = "trait", x = colnames(proteins))])
          traits <- strsplit(x = proteins$trait[[i]], split = ';')[[1]]
          n_traits <- length(traits)
          paths <- strsplit(x = proteins$path[[i]], split = ';')[[1]]
          for(j in 1:n_traits){
            trait_data <- data.table::fread(file = paths[j], header = TRUE, stringsAsFactors = FALSE)
            if(grepl(pattern = "@", x = traits[j])){
              type <- sapply(strsplit(x = traits[j], split = "@"), `[[`, 1)
              trait_data <- subset(x = trait_data, select = c("SNP", paste0("beta.", type), paste0("se.", type)))
              colnames(trait_data) <- c("SNP", "beta", "se")
              traits[j] <- sapply(strsplit(x = traits[j], split = "@"), `[[`, 2)
            }
            trait_name <- traits[j]
            dat_hyprcoloc_trait <- merge(x = dat_hyprcoloc, y = trait_data, by = "SNP", all = FALSE)
            rm(trait_data)
            if(nrow(dat_hyprcoloc_trait) == 0){
              message("Trait ", traits[j], " has no SNPs in common with other traits data so far. This trait will not be accounted for.")
              next
            }
            dat_hyprcoloc <- dat_hyprcoloc_trait
            for(k in c("beta", "se")){
              colnames(dat_hyprcoloc)[which(colnames(dat_hyprcoloc) == "beta")] <- paste0("beta.", traits[j])
            }
          }
        }
        for(k in c("beta", "se")){
          eval(parse(text = paste0(k, "s <- dat_hyprcoloc[grepl('", k, "', colnames(dat_hyprcoloc))]")))
          eval(parse(text = paste0("rownames(", k, "s) <- dat_hyprcoloc$SNP")))
          eval(parse(text = paste0("colnames(", k, "s) <- gsub(pattern = '", k, ".', replacement = '', x = colnames(dat_hyprcoloc)[grepl(pattern = '", k, "', x = colnames(dat_hyprcoloc))])")))
        }
        traits <- unique(gsub(pattern = "beta.", replacement = "", x = colnames(dat_hyprcoloc)[grepl(pattern = "beta", x = colnames(dat_hyprcoloc))]))
        rsid <- dat_hyprcoloc$SNP
        
        # Colocalization analysis
        hyprcoloc_res <- hyprcoloc(as.matrix(betas), as.matrix(ses), trait.names=traits, snp.id=rsid)[["results"]]
        hyprcoloc_res <- subset(x = hyprcoloc_res, select = c("traits", "posterior_prob", "candidate_snp", "posterior_explained_by_snp"))
        colnames(hyprcoloc_res) <- paste0("hyprcoloc_", colnames(hyprcoloc_res))
        
        traits <- c()
        post_probs <- c()
        candidates <- c()
        post_exps <- c()
        for(m in 1:nrow(hyprcoloc_res)){
          traits <- paste0(traits, hyprcoloc_res[i,"hyprcoloc_traits"], ";")
          post_probs <- paste0(post_probs, hyprcoloc_res[i,"hyprcoloc_posterior_prob"], ";")
          candidates <- paste0(candidates, hyprcoloc_res[i,"hyprcoloc_candidate_snp"], ";")
          post_exps <- paste0(post_exps, hyprcoloc_res[i,"hyprcoloc_posterior_explained_by_snp"], ";")
        }
        mr_prot$hyprcoloc_traits <- traits
        mr_prot$hyprcoloc_posterior_prob <- post_probs
        mr_prot$hyprcoloc_candidate_snp <- candidates
        mr_prot$hyprcoloc_posterior_explained_by_snp <- post_exps
        # mr_prot$hyprcoloc[i] <- hyprcoloc_res[[1]]
        # for(m in 2:4){
        #   mr_prot$hyprcoloc[i] <- paste(mr_prot$hyprcoloc[i], hyprcoloc_res[[m]], sep = ";")
        # }
        
      }
    }
    
    
    ############################################################## #
    
    
    mr_res <- rbind(mr_res, mr_prot)
    fwrite(x = mr_res,
           file = outname,
           sep = "\t", append = FALSE)
  }
  end <- Sys.time()
  print(end - start)
  fwrite(x = mr_res,
         file = outname,
         sep = "\t", append = FALSE)
}
############################################################################### #


############################################################################### #
#### 3. Volcano Plot ####
{
  if(volcano){
    library(ggplot2)
    library(ggrepel)
    library(magrittr)
    library(dplyr)
    plot_pval_col <- ifelse(test = sum(grepl(pattern = "padj", x = colnames(proteins_all)) != 0),
                            yes = "padj",
                            no = "pval")
    pval_col_volcano <- ifelse(test = opt$volcano_plot_adjp,
                               yes = plot_pval_col,
                               no = pval_col_noadj)
    proteins_all$col <- ifelse(test = eval(parse(text = paste0("proteins_all$", plot_pval_col))) < opt$pval,
                               yes = ifelse(test = proteins_all$beta < 0,
                                            yes = "green",
                                            no = "red"),
                               no = "black")
    lim <- ifelse(test = abs(min(proteins_all$beta)) > abs(max(proteins_all$beta)),
                  yes = abs(min(proteins_all$beta)),
                  no = abs(max(proteins_all$beta)))
    lim_y <- -log10(min(eval(parse(text = paste0("proteins_all$", pval_col_volcano)))))
    breaks_x_to <- ifelse(test = (ceiling(lim/0.25)*0.25)%%0.25 == 0,
                          yes = (ceiling(lim/0.25)*0.25)+0.25,
                          no = (ceiling(lim/0.25)*0.25))
    breaks_y_to <- ifelse(test = (ceiling(lim_y/2.5)*2.5)%%2.5 == 0,
                          yes = (ceiling(lim_y/2.5)*2.5) + 2.5,
                          no = (ceiling(lim_y/2.5)*2.5))
    # Plotting the volcano plot
    v_plot <- ggplot(data=proteins_all, aes(x=beta, y=-log10(eval(parse(text = pval_col_volcano))), col=col)) +
      geom_point(size =1) + 
      theme_minimal() +
      theme(legend.position="none") +
      scale_color_manual(values=c("black", "green3", "red3")) +
      geom_text_repel(data = . %>% mutate(label = ifelse(test = eval(parse(text = paste0("proteins_all$", plot_pval_col))) < opt$pval,
                                                         yes = proteins_all$name,
                                                         no = "")),
                      aes(label = label), 
                      box.padding = 0.2,
                      show.legend = FALSE, 
                      max.overlaps = 300,
                      size = 3.5,
                      segment.size = 0.5,
                      segment.alpha = 0.3,
                      segment.linetype = "solid",
                      min.segment.length = 0,
                      force = 5,
                      force_pull = 1, 
                      nudge_y = -log(eval(parse(text = paste0("proteins_all$", pval_col_volcano))))/10, 
                      nudge_x = -log2(abs(proteins_all$beta)) *proteins_all$beta/2
      ) +
      xlim(-lim, lim) +
      xlab("Effect (Beta)") +
      ylab("-log10(p)") +
      ggtitle(label = paste0("PW-MR (", opt$study, ")")) +
      scale_x_continuous(limits=c(-lim, lim), breaks = c(seq(from = -1, to = 1, by = 0.25))) +
      scale_y_continuous(limits=c(0, lim_y),breaks = c(seq(from = 0, to = breaks_y_to, by = 2.5)))
    v_plot
    ggsave(filename = outname_volcano_plot, plot = v_plot, device = "png", bg = "white", width = 4000, height = 2400, units = 'px')
  }
}
############################################################################### #


############################################################################### #
#### 4. Tissue Enrichment analysis ####
{
  if(opt$tissue){
    ## TissueEnrich R (https://bioconductor.org/packages/release/bioc/html/TissueEnrich.html)
    library(dplyr)
    library(ggplot2)
    library(TissueEnrich)
    
    inputGenes <- proteins$name
    inputGenes <- inputGenes[!duplicated(inputGenes)]
    gs<-GeneSet(geneIds=inputGenes,organism='Homo Sapiens',
                geneIdType=SymbolIdentifier())
    teOutput<-teEnrichment(gs)
    seEnrichmentOutput<-teOutput[[1]]
    enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
                                          row.names = rowData(seEnrichmentOutput)[,1]),
                               colData(seEnrichmentOutput)[,1])
    enrichmentOutput$Tissue<-row.names(enrichmentOutput)
    #Plotting the P-Values
    #png(filename = outname_tissue_plot, width = 1000, height = 800, type = "cairo")
    t_plot <- ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
                                label = Tissue.Specific.Genes,fill = Tissue))+
      geom_bar(stat = 'identity')+
      labs(x='', y = '-log10(P-value)')+
      theme_bw()+
      theme(legend.position='none')+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
              element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.grid.major= element_blank(),panel.grid.minor = element_blank())
    ggsave(filename = outname_tissue_plot, plot = t_plot, device = "png", bg = "white")
    #dev.off()
  }
}


