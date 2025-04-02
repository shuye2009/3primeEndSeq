## Nujhat ZNF281 3'end seq 
# 
options(error=recover)
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1] # wd <- "C:/GREENBLATT/Nujhat/3endseq/Jun04_2024"
factor <- args[2] # "gU2AF1", "gU170K", "gSF1"
option <- args[3]  # global, peak, read, other

if (Sys.info()[["sysname"]] == "Windows") {
   resources <- setupProjects(genome="hg19", 
                              baseDir="C:/GREENBLATT",
                              overwrite = FALSE)
            
   scriptdir <- "C:/GREENBLATT/Rscripts/3primeEndSeq/R"
   #scriptdir <- dirname(rstudioapi::getSourceEditorContext()$path)
}else{
   source("~/R_script/misc_genomics_functions.R") 
   resources <- setupProjects(genome="hg19", 
                              baseDir="~",
                              overwrite = FALSE) 
   scriptdir <- "~/R_script/3endseq"
}

source(file.path(scriptdir,"misc_genomics_functions.R"))
source(file.path(scriptdir,"3endseq_lib.R"))
source(file.path(scriptdir,"3endseq.R"))

  
## Process single nucleotide resolution cleavage (CPA) sites and internal priming sites,
## cluster CPA sites into consolidated CPA peaks,
## Run Dexseq analysis to define shorteing and lengthening genes 
if(option == "global"){
   outdir <- file.path(wd, "global")
   if(!dir.exists(outdir)) dir.create(outdir)
   
   if(file.exists(file.path(outdir, "res.RDS"))){
      res <- readRDS(file.path("res.RDS"))
   }else{
      res <- process_global(indir, outdir, resources=resources)
   }
      
   CScluster <- res[res$status=="kept"]
   
   outdir <- file.path(wd, factor)
   if(!dir.exists(outdir)) dir.create(outdir)
   sample_info <- paste0("sample_info_", factor, ".txt")
   
   if(file.exists(file.path(outdir, "cluster_res.RDS"))){
      cluster_res <- readRDS(file.path(outdir, "cluster_res.RDS"))
   }else{
      cluster_res <- process_factor(indir=wd, outdir=outdir, factor=factor, metaData=sample_info, 
                  resources=resources, CScluster=CScluster)
   }

   if(!file.exists(file.path(outdir, "3UTR", "DEXSeq_results_annotated.tab")))
      DEXSeq_analysis(cluster_res=cluster_res, outdir=outdir, feature="3UTR", factor=factor)
   if(!file.exists(file.path(outdir, "Intron", "DEXSeq_results_annotated.tab")))
      DEXSeq_analysis(cluster_res=cluster_res, outdir=outdir, feature="Intron", factor=factor)
   if(!file.exists(file.path(outdir, "Transcript", "DEXSeq_results_annotated.tab")))
      DEXSeq_analysis(cluster_res=cluster_res, outdir=outdir, feature="Transcript", factor=factor)
   
   
   if(!file.exists(file.path(outdir, "Transcript", "premature_CPA_annotated.tab")))
   intronic_CPA(file.path(outdir, "Transcript"), "DEXSeq_results_annotated.tab")
   
}
  
## Plot iCLIP peaks around proximal and distal cleavage sites
if(option == "peak"){
   clipDir <- "C:/GREENBLATT/Nabeel/CLIP_peaks_Jun2024"
   # "U2AF1", "U1C", "SF3A1", "SF3B4",
   for (clipFactor in c("PTBP1", "PRPF31", "MATR3")){
      integrate_3endseq_with_iCLIP(clipDir, wd, clipFactor=clipFactor, 
                                   cleavageFactor=factor, feature="3UTR")
   }
}

## Plot iCLIP reads around proximal and distal cleavage sites
if(option == "read"){
   clipDir <- "C:/GREENBLATT/Nabeel/CLIP_bam_merged"
   # "U2AF1", "U1C", "SF3B4", 
   if(0){
   for (clipFactor in c("U2AF1", "U1C", "SF3B4", "PTBP1", "NUDT21", "SRSF2")){
       lh_list <- list(proximal = list(shortening = list(first = c(0, 100), second = c(300, 500)),
                                       lengthening = list(first = c(100, 300), second = c(300, 500))),
                       distal = list(shortening = list(first = c(-200, 0), second = c(0, 200)),
                                     lengthening = list(first = c(-200, 0), second = c(0, 200))))
       ext_list <- list(proximal = c(-100, 500), distal = c(-400, 200))
      integrate_3endseq_with_iCLIP_bam(clipDir, wd, clipFactor=clipFactor, 
                                   cleavageFactor=factor, feature="3UTR", lh_list, ext_list)
   }
   }
   
   #For SF3B4, U2AF1, SRSF2 and NUDT21: 
        #stats for 0 to +200 around proximal of shortening genes in all 3 conditions (gU170K, gSF1 and gU2AF1)
   #For U1C only:
        #Stats for 0 to +200 around proximal of shortening in gU170K
        #Stats for +100 to +400 around proximal of shortening in gSF1 and gU2AF1
   
   #For SF3B4 in gU2AF1 condition only:
        #Stats for +200 to +400 around proximal of shortening
   if(0){     
   for (clipFactor in c("U2AF1", "SF3B4", "NUDT21", "SRSF2", "U1C")){
       
       stat1 <- c(0, 200)
       stat2 <- NULL
       if(clipFactor == "SF3B4" && factor == "gU2AF1") stat2 <- c(200, 400)
       if(clipFactor == "U1C" && factor %in% c("gU2AF1", "gSF1")) stat1 <- c(100, 400)
           
       lh_list <- list(proximal = list(shortening = list(first = stat1, second = stat2),
                                       lengthening = list(first = NULL, second = NULL)),
                       distal = list(shortening = list(first = NULL, second = NULL),
                                     lengthening = list(first = NULL, second = NULL)))
       ext_list <- list(proximal = c(-100, 500), distal = c(-400, 200))
       integrate_3endseq_with_iCLIP_bam(clipDir, wd, clipFactor=clipFactor, 
                                        cleavageFactor=factor, feature="3UTR", lh_list, ext_list)
   }
   }
   
   #For the premature CPA shortening category only (where proximal is in intron and distal is in 3'UTR), could you please plot: 
   #CLIP-over-input of U2AF1 around the proximal and distal sites of gU2AF1-shortening and "no  change" genes  
   #CLIP-over-input of SF3B4 around the proximal and distal sites of gSF1-shortening and "no change" genes
   #CLIP-over-input of U1C around the proximal and distal sites of gU170K-shortening and "no change" genes
   #We can try 2 different windows around the cleavage site (+/-1kb and +/-500bp).
   #As done for previous CLIP-APA correlation plots, please match the number and expression of the "no change" genes in each case to the "shortening" genes under wild-type conditions.

   if(0){     
      for (clipFactor in c("U2AF1", "SF3B4", "U1C")){
         
         
         if((clipFactor == "U2AF1" && factor == "gU2AF1") | 
            (clipFactor == "SF3B4" && factor %in% c("gSF1")) |
            (clipFactor == "U1C" && factor %in% c("gU170K"))){
            stat1 <- c(0, 200)
            stat2 <- c(-200, 0)
         }else{
            stat1 <- NULL
            stat2 <- NULL
         }
         
         lh_list <- list(proximal = list(shortening = list(first = stat1, second = NULL),
                                         lengthening = list(first = stat1, second = NULL)),
                         distal = list(shortening = list(first = stat2, second = NULL),
                                       lengthening = list(first = stat2, second = NULL)))
         ext_list <- list(proximal = c(-1000, 1000), distal = c(-1000, 1000))
         integrate_3endseq_with_iCLIP_bam(clipDir, wd, clipFactor=clipFactor, 
                                          cleavageFactor=factor, feature="Transcript", lh_list, ext_list)
         
         ext_list <- list(proximal = c(-500, 500), distal = c(-500, 500))
         integrate_3endseq_with_iCLIP_bam(clipDir, wd, clipFactor=clipFactor, 
                                          cleavageFactor=factor, feature="Transcript", lh_list, ext_list)
      }
   }
   
   # Dec19, 2024
   #For U1C with gU170K premature CPA: Plot window +/-500bp around both cleavage sites and statistical tests for -300 to -100. 
   #For SF3B4 with gSF1 premature CPA: Plot window +/-500bp around both cleavage sites and statistical tests for -100 to +100.
   #Could you also generate plots for U1C CLIP around the gSF1 premature CPA events, and SF3B4 around the gU170K premature events. Window can be +/-500 bp for both proximal and distal. Region for stats can be arbitrary for now.
   #
   #
   if(0){     
       for (clipFactor in c("SF3B4", "U1C")){
           if(factor %in% c("gU170K")){
               stat1 <- c(-300, -100)
               stat2 <- NULL
           }else if(factor %in% c("gSF1")){
               stat1 <- c(-100, 100)
               stat2 <- NULL
           }else{
               stat1 <- NULL
               stat2 <- NULL
           }
           
           lh_list <- list(proximal = list(shortening = list(first = stat1, second = NULL),
                                           lengthening = list(first = stat1, second = NULL)),
                           distal = list(shortening = list(first = stat1, second = NULL),
                                         lengthening = list(first = stat1, second = NULL)))
           
           ext_list <- list(proximal = c(-500, 500), distal = c(-500, 500))
           integrate_3endseq_with_iCLIP_bam(clipDir, wd, clipFactor=clipFactor, 
                                            cleavageFactor=factor, feature="Transcript", lh_list, ext_list)
       }     
   }
   
   # Jan 06, 2025
   #Would you be able to give me the stats for these specific regions?
   #For U1C CLIP around gSF1 plots: -300 to -100
   #For SF3B4 around gU170K plots: -400 to -200 and 0 to 200
   
   #Could you also generate similar plots for U1C around gU2AF1 premature CPA events, as well as SF3B4 around gU2AF1 premature events? Plot window can be +/-500bp for all.
   #Region for stats can be -300 to -100 and -100 to +100 for these.
   #
   if(1){     
       for (clipFactor in c("SF3B4", "U1C", "U2AF1")){
           if(clipFactor %in% c("SF3B4") && factor %in% c("gU170K")){
               stat1 <- c(-400, -200)
               stat2 <- c(0, 200)
           }else if(clipFactor %in% c("U1C") && factor %in% c("gSF1")){
               stat1 <- c(-300, -100)
               stat2 <- NULL
           }else if(clipFactor %in% c("U1C", "SF3B4") && factor %in% c("gU2AF1")){
               stat1 <- c(-300, -100)
               stat2 <- c(-100, 100)
           }else {
               stat1 <- NULL
               stat2 <- NULL
           }
           
           lh_list <- list(proximal = list(shortening = list(first = stat1, second = stat2),
                                           lengthening = list(first = stat1, second = stat2)),
                           distal = list(shortening = list(first = stat1, second = stat2),
                                         lengthening = list(first = stat1, second = stat2)))
           
           ext_list <- list(proximal = c(-500, 500), distal = c(-500, 500))
           integrate_3endseq_with_iCLIP_bam(clipDir, wd, clipFactor=clipFactor, 
                                            cleavageFactor=factor, feature="Transcript", lh_list, ext_list)
       }     
   }
}

## Plot iCLIP peaks with respect to genomic features
if(option == "other"){
   clipDir <- "C:/GREENBLATT/Nabeel/CLIP_peaks_Jun2024"
   outDir = file.path(wd, "genomicplot")
   if(!dir.exists(outDir)) dir.create(outDir)
   #"U1C", "SF3B4", "PEG", "PCIF1", "NUDT16L1"
   for(factorName in c("U2AF1")){
      run_GenomicPlot_Intron(clipDir, outDir, clipFactors=factorName, 
                            txdb=resources$txdb)
      run_GenomicPlot_metagene(clipDir, outDir, clipFactors=factorName,
                               metaFeature=resources$metaFeatures)
      run_GenomicPlot_peak_annotation(clipDir, outDir, clipFactors=factorName,
                                      gtfFile=resources$gtffile)
   }
   
}

## Plot distance between proximal and distal CPA sites for lengthening and shortenig genes
if(option == "distance"){
    for(nc in c("noChange_", "noChangeS_", "noChangeL_")){
        distance_between_proximal_distal(wd, cleavageFactor=factor, 
                                         feature="3UTR", nc = nc)    
    }
    
}


# run in bash terminal
# Rscript C:/GREENBLATT/Rscripts/3primeEndSeq/R/run_3endseq.R C:/GREENBLATT/Nujhat/3endseq/Jun04_2024 gU2AF1 read