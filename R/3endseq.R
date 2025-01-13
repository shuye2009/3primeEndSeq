

process_factor <- function(indir, outdir, factor, metaData, resources, CScluster){
   print("[process_factor] started ...")
   setwd(outdir)
   
   if(file.exists("cluster_res.RDS")){
      cluster_res <- readRDS("cluster_res.RDS")
   }else{
      force(CScluster)
      gtffile <- resources$gtffile
      txdb <- resources$txdb
      metaFeatures <- resources$metaFeatures
      geneFeatures <- resources$geneFeatures
      
      PAS_strong <- c("AAUAAA", "AUUAAA")
      PAS_weak <- c("AGUAAA", "UAUAAA", "CAUAAA", "GAUAAA", "AAUAUA", "AAUACA", "AAUAGA", "ACUAAA", "AAGAAA", "AAUGAA")
      
      sampleData <- read.delim(metaData) %>% 
         dplyr::mutate_if(is.character, as.factor)
      rownames(sampleData) <- sampleData$sample
      sampleData <- sampleData[, -1]
      
      ### metagene plot ####
      print("[process_factor] metagene plot ...")
      CS_files <- match_samples(path=indir, pattern="_CS.bed", 
                                  samples=rownames(sampleData))
      
      importParams <- setImportParams(offset = 0, fix_width=0, 
                                      fix_point="center", norm=FALSE, 
                                      useScore=FALSE,
                                      outRle=TRUE, useSizeFactor=FALSE, 
                                      genome="hg19")
      
      op <- paste0(factor, "_3primeseq_CS_sites_5parts_metagene")
      if(!file.exists(paste0(op, ".pdf")))
      df <- plot_5parts_metagene(CS_files, gFeatures_list=list("meta"=metaFeatures), heatRange=NULL, inputFiles=NULL, scale=FALSE, verbose=TRUE, transform=NA, smooth=FALSE, stranded=TRUE, outPrefix=op, importParams=importParams, heatmap=TRUE, rmOutlier=0, nc=5)
      
      IP_files <- match_samples(path=indir, pattern="_IP.bed", 
                                samples=rownames(sampleData))
   
      
      op <- paste0(factor, "_3primeseq_IP_sites_5parts_metagene")
      if(!file.exists(paste0(op, ".pdf")))
      df <- plot_5parts_metagene(IP_files, gFeatures_list=list("meta"=metaFeatures), heatRange=NULL, inputFiles=NULL, scale=FALSE, verbose=TRUE, transform=NA,  smooth=FALSE, stranded=TRUE, outPrefix=op, importParams=importParams, heatmap=TRUE, rmOutlier=0, nc=5)
      
      
      ### count polyT length ####
      print("[process_factor] counting polyT length ...")
      fas <- match_samples(path=indir, pattern="_filtered_flankL20.fa",
                        samples=rownames(sampleData))
      
      noT <- match_samples(path=indir, pattern="_noT.fa",
                           samples=rownames(sampleData))
      
      if(!file.exists("Occurance_of_G_upstream20.pdf")){
         for(ch in c("A", "C", "G")){
            pdf(paste0("Occurance_of_",ch,"_upstream20.pdf"))
            opar <- par(mfrow=c(2,2))
            
            for (i in 1:length(fas)){
               fa <- fas[i]
               character_occurance_fasta(fa, ch)
            }
            dev.off()
            par(opar)
            
            
            pdf(paste0("Occurance_of_no",ch,"_upstream20.pdf"))
            opar <- par(mfrow=c(2,2))
            
            for (i in 1:length(noT)){
               fa <- noT[i]
               character_occurance_fasta(fa, ch)
            }
            dev.off()
            par(opar)
         }
      }
      
      
      ### overlap merged with individuals CS.bed ####
      print("[process_factor] overlap CS with CS clusters ...")
      importParams <- setImportParams(offset = 0, fix_width=0, 
                                      fix_point="center", norm=FALSE, 
                                      useScore=FALSE,
                                      outRle=FALSE, useSizeFactor=FALSE, 
                                      genome="hg19")
      if(!file.exists("CS_gr.RDS")){
         cl <- start_parallel(nc=length(CS_files))
         clusterExport(cl, varlist=c("handle_bed", "seqlevelsStyle"))
         CS_list <- parLapply(cl, names(CS_files), function(x) handle_bed(CS_files[x], importParams = importParams))
         CS_gr_list <- lapply(CS_list, function(x)x$query)
         stop_parallel(cl)
         
         saveRDS(CS_gr_list, "CS_gr.RDS")
      }else{
         CS_gr_list <- readRDS("CS_gr.RDS")
      }
      
      if(!file.exists("overlap_gr.RDS")){
         cl <- start_parallel(nc=length(CS_gr_list))
         clusterExport(cl, varlist=c("CScluster"))
         overlap_list <- parLapply(cl, CS_gr_list, function(x) {
            GenomicRanges::countOverlaps(subject=x, query=CScluster)})
         stop_parallel(cl)
         
         saveRDS(overlap_list, "overlap_gr.RDS")
      }else{
         overlap_list <- readRDS("overlap_gr.RDS")
      }
      
      names(overlap_list) <- names(CS_files)
      overlap_mat <- bind_cols(overlap_list)
      
      
      pdf("sample_correlation_plot.pdf")
      #DescTools::PlotPairs(log10(overlap_mat+1))
      limma::plotMDS(log10(overlap_mat+1))
      print(pheatmap(cor(overlap_mat)))
      dev.off()
      
      
      mcols(CScluster) <- cbind(mcols(CScluster), overlap_mat)
      CScluster <- detect_PAS_signal(BSgenome.Hsapiens.UCSC.hg19, CScluster, upstream=50)
      
      pdf("cluster_width_readCounts_plot.pdf")
      print(plot(x=width(CScluster), y=CScluster$score, xlab="Cluster_width", ylab="Read_count", log="y"))
      print(hist(log10(CScluster$score)))
      print(hist(width(CScluster)))
      
      dev.off()
      
      cluster_out <- gr2df(CScluster)
      colnames(cluster_out)[1] <- "#chr"
      write.table(cluster_out, "cluster_read_counts.bed", 
                  sep="\t", row.names=FALSE, quote=FALSE)
      
      
      cluster_file <- "cluster_read_counts.bed"
      names(cluster_file) <- "cluster"
      
      if(!file.exists("cluster_annotation.pdf"))
      annot <- plot_peak_annotation(cluster_file, gtffile, importParams=importParams, fiveP=0, threeP=5000, outPrefix = "cluster_annotation", verbose=TRUE)
      
      ## associate cluster with genes
      cluster_out <- read.table(cluster_file, header=TRUE,  comment.char = "",)
      colnames(cluster_out)[1] <- "chr"
      cluster_out <- cluster_out %>% arrange(chr, start, strand)
      
      annotation_file <- "cluster_peak_annotations.tab"
      cluster_tx_association <- read.table(annotation_file, header=TRUE)
      
      ## if a peak is associated with multiple genes (transcripts), select only one arbitrarily
      cluster_tx_association <- cluster_tx_association[!duplicated(cluster_tx_association$idPeak),] 
      
      gene_annotation_file <- "cluster_targeted_annotated_gene.tab"
      tx_gene_association <- read.table(gene_annotation_file, header=TRUE)
      
      pdf("histogram_of_number_of_CS_in_3UTRIntron.pdf", width=10, height=6)
      
      p1 <- ggplot(tx_gene_association, aes(x=X3.UTR)) +
         geom_histogram(position="identity", bins=max(tx_gene_association$X3.UTR), fill="red4", color="black") +
         ggtitle("Cleavage sites in 3'UTR") +
         xlab("Number of cleavage sites") +
         ylab("Number of genes") +
         theme_classic()
      p2 <- ggplot(tx_gene_association, aes(x=TTS)) +
         geom_histogram(position="identity", bins=max(tx_gene_association$TTS), fill="red4", color="black") +
         ggtitle("Cleavage sites in TTS") +
         xlab("Number of cleavage sites") +
         ylab("Number of genes") +
         theme_classic()
      p3 <- ggplot(tx_gene_association, aes(x=Intron)) +
         geom_histogram(position="identity", bins=max(tx_gene_association$Intron), fill="red", color="black") +
         ggtitle("Cleavage sites in Intron") +
         xlab("Number of cleavage sites") +
         ylab("Number of genes") +
         theme_classic()
      print(plot_grid(p1,p2,p3, nrow=1))
      
      dev.off()
      
      
      cluster_tx_association <- cluster_tx_association %>%
         select(c(chrPeak, startPeak, endPeak, strandPeak, tx_name, feature_name)) %>%
         rename(chr=chrPeak, start=startPeak, end=endPeak, strand=strandPeak) %>%
         distinct() %>%
         arrange(chr, start, strand) %>%
         mutate(start=start+1)
      
         
      tx_cluster_count <- inner_join(cluster_out, cluster_tx_association)
      
      cluster_res <- list(count=tx_cluster_count, association=tx_gene_association,
                          meta=sampleData)
      saveRDS(cluster_res, "cluster_res.RDS")
   }
   
   return(cluster_res)
}
 
DEXSeq_analysis <- function(cluster_res, outdir, feature, factor){   
   tx_cluster_count <- cluster_res$count
   tx_gene_association <- cluster_res$association
   sampleData <- cluster_res$meta
   
   ctl_col <- rownames(sampleData %>% filter(condition=="control"))
   exp_col <- rownames(sampleData %>% filter(condition=="treat"))
   
   if(feature == "3UTR"){
      featureNames <- c("3'UTR", "TTS")
      result_feature <- tx_cluster_count %>%
         filter(feature_name %in% featureNames)
   }else if(feature == "Intron"){
      featureNames <- c("Intron")
      result_feature <- tx_cluster_count %>%
         filter(feature_name %in% featureNames)
   }else if(feature == "Transcript"){
      result_feature <- tx_cluster_count 
   }else {
      stop("feature", feature, "is not supported")
   }
   
   message("[DEXSeq analysis] for feature ", feature, " started ...")
   
   subfolder <- file.path(outdir, feature)
   dir.create(subfolder)
   setwd(subfolder)
   
   message("[DEXSeq analysis] associate cleavage-sites cluster with gene")
   result_out <- left_join(result_feature, tx_gene_association, by=c("tx_name"), 
                           suffix=c("", "y")) %>%
      select(c(colnames(result_feature), gene_name)) %>%
      mutate(gene_name=as.factor(gene_name), name=as.factor(name)) %>%
      arrange(gene_name, start)
   
   write.table(result_out, paste0("cluster_read_counts_", feature, ".bed"), 
               sep="\t", row.names=FALSE, quote=FALSE)
   result_out <- read.table(paste0("cluster_read_counts_", feature, ".bed"),  
                            header=TRUE,  comment.char = "",)
   ### APA analysis ####
   message("[DEXSeq analysis] Fisher exact test-based APA analysis")
   if(!file.exists("apa_res.RDS")){
      apa_res <- cal_RCE(result_out, ctl_col, exp_col, 0.02)
      apa_res <- APA_assign(apa_res)
      apa_res <- UTR_change(apa_res)
      
      saveRDS(apa_res, "apa_res.RDS")
   }else{
      apa_res <- readRDS("apa_res.RDS")
   }
   
   write.table(apa_res, "APA_results.tab", col.names=TRUE, 
               row.names=FALSE, sep="\t", quote=FALSE)
   gene_change_apa <- unique(apa_res[, c("gene_name", "change")])
   colnames(gene_change_apa) <- c("gene_name", "change_apa")
   
   ### DEXseq analysis ####
   message("[DEXSeq analysis] DEXSeq analysis ...")
   if(!file.exists("dxr.RDS")){
      countData <- as.matrix(result_out[, rownames(sampleData)])
      
      featureID <- as.character(result_out$name)
      groupID <- result_out$gene_name
      featureRanges <- makeGRangesFromDataFrame(result_out)
      dxd <- DEXSeqDataSet( countData, sampleData,
                            design= ~ sample + exon + condition:exon,
                            featureID, groupID, featureRanges,
                            transcripts=NULL, alternativeCountData=NULL)
      
      print(colData(dxd))
      print(head(featureCounts(dxd), 5))
      print(head(counts(dxd), 5))
      
      if (Sys.info()[["sysname"]] == "Windows") {
         BPPARAM = SnowParam(workers=nrow(sampleData), type="SOCK")
      }else{
         BPPARAM = MulticoreParam(workers=nrow(sampleData))
      }
      message("estimating size factors")
      dxd = estimateSizeFactors(dxd)
      message("estimating dispersions")
      dxd = estimateDispersions(dxd, BPPARAM=BPPARAM)
      message("testing for DEU")
      dxd = testForDEU(dxd, BPPARAM=BPPARAM)
      message("estimating exon fold changes")
      dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
      message("preparing results dxr")
      dxr = DEXSeqResults(dxd)
      
      saveRDS(dxr, "dxr.RDS")
   } else {
      dxr <- readRDS("dxr.RDS")
   }
   
   dxr_df <- as.data.frame(dxr)
   write.table(dxr_df, "DEXSeq_results.tab", col.names=TRUE, 
               row.names=FALSE, sep="\t", quote=FALSE)
   
   dxr_df <- dxr_df[!is.na(dxr_df$padj), ]
   dxr_df <- APA_assign(dxr_df, gene_col="groupID", start_col="genomicData.start", strand_col="genomicData.strand")
   dxr_df <- UTR_change(dxr_df, gene_col="groupID", start_col="genomicData.start", strand_col="genomicData.strand", sorting_col="exonBaseMean", padj_col="padj", direction_col="log2fold_treat_control", padj_cutoff=0.1, norm=0)
   
   dxr_gr <- makeGRangesFromDataFrame(dxr_df)
   dxr_gr <- detect_PAS_signal(BSgenome.Hsapiens.UCSC.hg19, dxr_gr, upstream=50)
   dxr_df$status <- dxr_gr$status
   
   # add feature_name annotation
   cluster_feature <- tx_cluster_count %>%
      dplyr::select(name, feature_name) %>%
      dplyr::mutate(name = as.character(name))
   
   dxr_df <- left_join(dxr_df, cluster_feature, by = join_by(featureID == name))
   write.table(dxr_df, "DEXSeq_results_annotated.tab", 
               col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
   
   
   ## compare dxr with apa ####
   message("[DEXSeq analysis] compare dxr with apa ...")
   gene_change_dxr <- unique(dxr_df[, c("groupID", "change")])
   colnames(gene_change_dxr) <- c("gene_name", "change_dxr")
   combined_change <- inner_join(gene_change_dxr, gene_change_apa)
   confusion_mat <- as.matrix(table(combined_change[, c("change_dxr", "change_apa")]))
   colnames(confusion_mat) <- paste0(colnames(confusion_mat), "_apa")
   rownames(confusion_mat) <- paste0(rownames(confusion_mat), "_dxr")
   write.table(confusion_mat, "confusion_matrix.tab", 
               col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
   
   combined_res <- full_join(apa_res, dxr_df, by=c("name"="featureID", "gene_name"="groupID")) %>%
      arrange(gene_name)
   
   write.table(combined_res, "combined_results.tab", 
               col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
   
   dim(combined_res)
   common_res <- combined_res %>%
      tidyr::drop_na()
   discrepancy_res <- common_res %>%
      filter((change.x=="shortening" & change.y=="lengthening") | (change.y=="shortening" & change.x=="lengthening"))
   write.table(discrepancy_res, "discrepancy.tab", 
               col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
   
   pdf("FANCI_example_plot.pdf", width=10, height=8)
   print(plotDEXSeq(dxr, "FANCI", norCounts=TRUE, expression=FALSE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 ))
   dev.off()
   
   pdf("PAS_readCounts_plot.pdf", width=10, height=8)
   p1 <- draw_boxplot_by_factor(apa_res, xc="status", yc="score", Ylab="Read count", comp=combn(c(1,2,3), 2, simplify = FALSE), nf=1)
   p2 <- draw_boxplot_by_factor(apa_res, xc="status", yc="score", Ylab="Read count", comp=combn(c(1,2,3), 2, simplify = FALSE), nf=1)
   p3 <- draw_rank_plot(apa_res, xc="status", yc="score", Ylab="Read count")
   p4 <- draw_mean_se_barplot(apa_res, xc="status", yc="score", comp=combn(c(1,2,3), 2, simplify = FALSE), Ylab="Read count")
   
   print(plot_grid(p1, p2, p3, p4, ncol=2))
   dev.off()
   
   
   sink("Tally_table.tab")
   cat("\nNumber of genes with APA changes\n")
   print(table(gene_change_dxr$change_dxr))
   cat("\nNumber of genes with CS in 3'UTR\n")
   print(table(tx_gene_association$X3.UTR))
   cat("\nNumber of genes with CS in TTS\n")
   print(table(tx_gene_association$TTS))
   cat("\nNumber of genes with CS in Intron\n")
   print(table(tx_gene_association$Intron))
   cat("\nNumber of genes with CS in 5'UTR\n")
   print(table(tx_gene_association$X5.UTR))
   sink()
   
   ## output bed files ####
   message("[DEXseq analysis] output bed files ...")
   beddir <- file.path(subfolder, "bed_output")
   if(!dir.exists(beddir)) dir.create(beddir)
   results_from_dexseq_relative(factor = factor, dxrFile="DEXSeq_results_annotated.tab", outdir = beddir)
   
   message("[DEXseq analysis] finished ...")
}


process_global <- function(indir, outdir, resources=resources){
   setwd(outdir)
   
   if(file.exists("res.RDS")){
      res <- readRDS("res.RDS")
   }else{
      gtffile <- resources$gtffile
      txdb <- resources$txdb
      metaFeatures <- resources$metaFeatures
      geneFeatures <- resources$geneFeatures
   
      
      CS_file <- file.path(indir, "merged_CSs.bed")
      names(CS_file) <- "SC_cluster"
      IP_file <- file.path(indir, "merged_IPs.bed")
      names(IP_file) <- "IP_cluster"
      queryfiles <- c(CS_file, IP_file)
      
      importParams <- setImportParams(offset=0, fix_width=0, fix_point="center", 
                                      norm=FALSE, useScore=TRUE, outRle=TRUE,
                                      useSizeFactor=FALSE, genome="hg19")
      
      op <- "3primeseq_merged_sites_5parts_metagene"
      if(!file.exists(paste0(op, ".pdf")))
      df <- plot_5parts_metagene(queryfiles, gFeatures_list=list("meta"=metaFeatures),
                                 heatRange=NULL, inputFiles=NULL, scale=FALSE,
                                 verbose=TRUE, transform=NA, smooth=TRUE, stranded=TRUE,
                                 outPrefix=op, importParams=importParams, heatmap=TRUE,
                                 rmOutlier=0, nc=5)
      op <- "cluster_utr3_plots"
      if(!file.exists(paste0(op, ".pdf")))
      plot_start_end(queryFiles=queryfiles, inputFiles=NULL, centerFiles="utr3", 
                     txdb=txdb, importParams = importParams, binSize=10, insert=0,
                     verbose=FALSE, ext=c(-500, 200, -200, 500), 
                     hl=c(-50, 50, -50, 50), stranded=TRUE, scale=FALSE, 
                     smooth=FALSE, rmOutlier=0, outPrefix=op, 
                     transform=NA, shade=TRUE, nc=5)
      
      importParams <- setImportParams(offset=0, fix_width=0, fix_point="center", 
                                      norm=FALSE, useScore=TRUE, outRle=FALSE,
                                      useSizeFactor=FALSE, genome="hg19")
      
      CSs <- handle_bed(CS_file, importParams = importParams)$query
      IPs <- handle_bed(IP_file, importParams = importParams)$query
      #IPs_in_3UTR <- filter_bed_by_genomicFeature(IP_file, "utr3", txdb, importParams)$query
      
      merged <- c(CSs)
      merged <- sort(merged, decreasing=TRUE, by=~score)
      merged$name <- seq_along(merged) ## merged is sorted by score in descending order on RC-cluster
      
      res <- CS2peak(merged, separatingDist=75, addonDist=5)
      saveRDS(res, "res.RDS")
     
      CScluster <- res[res$status=="kept"]
      eliminated <- res[res$status=="eliminated"]
      added <- res[res$status=="added"]
      res_df <- gr2df(res) %>%
         mutate(score = log10(score))
      
      pdf("Read_counts_distribution.pdf", width=10, height=6)
      p1 <- ggplot(res_df, aes(x=score, fill=status)) +
         geom_density(position="identity", color="black", alpha=0.5) +
         ggtitle("Density plot of read counts") +
         xlab("log10 (read counts)") +
         scale_fill_npg() +
         ylab("Density") +
         theme_classic()
      
      p2 <- draw_boxplot_by_factor(res_df, xc="status", yc="score", Ylab="log10 (read counts)", comp=list(c(1,2), c(1,3), c(2,3)), nf=1)
      print(plot_grid(p1,p2, nrow=1, align="h", axis="tb"))
      dev.off()
      
   }
   
   return(res)
}

integrate_3endseq_with_iCLIP <- function(clipDir, wd, clipFactor = "U2AF1", 
                                         cleavageFactor="gU2AF1", feature="3UTR"){
   clipbed <- match_samples(path=clipDir, pattern="\\.bed$", samples=clipFactor)
   names(clipbed) <- paste0(clipFactor, "_clip_peak")
   
   exp <- gsub("^g","", cleavageFactor)
   setwd(file.path(wd, cleavageFactor, feature, "bed_output"))
   
   bedfiles <- list.files(pattern="\\.bed$")
   short <- gsub(paste0(paste0(cleavageFactor,"_"),"|\\.bed"), "", bedfiles, fixed = FALSE)
   names(bedfiles) <- short
   
   for(position in c("proximal", "distal")){
      selected <- short[grep(position , short)]
      selected <- selected[grep("shortening|lengthening|noChange", selected)]
   
      op <- paste0(clipFactor, "_clip_peaks_around_", position, "_cleavage_sites_", exp)
      if(position == "proximal"){
         ext <- c(-200, 500)
      }else{
         ext <- c(-500, 200)
      }
      plot_locus(
         queryFiles = clipbed,
         centerFiles = bedfiles[selected],
         txdb = NULL,
         ext = ext,
         hl = c(-100, 100),
         shade = TRUE,
         smooth = TRUE,
         importParams = setImportParams(),
         verbose = FALSE,
         binSize = 10,
         refPoint = "center",
         Xlab = position,
         Ylab = "Coverage/base/site",
         inputFiles = NULL,
         stranded = TRUE,
         heatmap = TRUE,
         scale = FALSE,
         outPrefix = op,
         rmOutlier = 0,
         transform = NA,
         statsMethod = "wilcox.test",
         heatRange = NULL,
         hw = c(8, 8),
         nc = 5
      )
   }
}

integrate_3endseq_with_iCLIP_bam <- function(clipDir, wd, clipFactor = "U2AF1", 
                                         cleavageFactor="gU2AF1", feature="3UTR",
                                         lh_list, ext_list){
   clipbam <- match_samples(path=clipDir, pattern="\\.bam$", samples=clipFactor, start = TRUE)
   names(clipbam) <- paste0(clipFactor, "_clip_read")
   
   inputbam <- file.path(clipDir, "INP_U1C_merged_noHotspot.thUni.bam")
   if(clipFactor == "U2AF1"){
      inputbam <- file.path(clipDir, "Input_NUDT21_merged_noHotspot.thUni.bam")
   }else if(clipFactor == "PTBP1"){
      inputbam <- file.path(clipDir, "Input_PTBP1_merged_noHotspot.thUni.bam")
   }else if(clipFactor == "NUDT21"){
      inputbam <- file.path(clipDir, "Input_NUDT21_merged_noHotspot.thUni.bam")
   }
   names(inputbam) <- "Input"
   
   print(clipbam)
   print(inputbam)
   
   
   exp <- gsub("^g","", cleavageFactor)
   setwd(file.path(wd, cleavageFactor, feature, "bed_output"))
   
   bedfiles <- list.files(pattern="\\.bed$")
   short <- gsub(paste0(paste0(cleavageFactor,"_"),"|\\.bed"), "", bedfiles, fixed = FALSE)
   names(bedfiles) <- short
   
   if(feature == "Transcript"){
       CPAs <- c("premature") #, "UTR3_only", "intron_only")
   }else{
       CPAs <- NA
   }
   
   
   string_list <- list(shortening = "shortening|noChangeS", lengthening = "lengthening|noChangeL")
   
   for(position in c("proximal", "distal")){
      ext <- ext_list[[position]]
      
      for(apa in c("shortening", "lengthening")){
         selected <- short[grep(position , short)]
         selected <- selected[grep(string_list[[apa]], selected)]
         
         for(CPA in CPAs){
             print(CPA)
             if(!is.na(CPA)){
                 CPA_selected <- selected[grep(CPA, selected)]
             }else{
                 CPA_selected <- selected
             }
             lhs <- lh_list[[position]][[apa]]
             
             for(n in names(lhs)){
                 lh <- lhs[[n]]
                 if(!is.null(lh)){
                     if(!is.na(CPA)){
                         op <- paste0(clipFactor, "_clip_reads_around_", apa, "_", position, "_", CPA,
                                      "_", exp, "_", paste(lh, collapse = "-"), "_", paste(ext, collapse = "-"))
                     }else{
                         op <- paste0(clipFactor, "_clip_reads_around_", apa, "_", position, 
                                      "_", exp, "_", paste(lh, collapse = "-"), "_", paste(ext, collapse = "-"))
                     }
                     
                     plot_locus(
                         queryFiles = clipbam,
                         centerFiles = bedfiles[CPA_selected],
                         txdb = NULL,
                         ext = ext,
                         hl = lh,
                         shade = TRUE,
                         smooth = TRUE,
                         importParams = setImportParams(norm = TRUE),
                         verbose = FALSE,
                         binSize = 10,
                         refPoint = "center",
                         Xlab = position,
                         Ylab = "Coverage/base/site",
                         inputFiles = inputbam,
                         stranded = TRUE,
                         heatmap = TRUE,
                         scale = FALSE,
                         outPrefix = op,
                         rmOutlier = 0,
                         transform = NA,
                         statsMethod = "wilcox.test",
                         heatRange = c(0,1,2),
                         hw = c(8, 8),
                         nc = 5
                     )
                 }
             }
          }
       }
    }
}
run_GenomicPlot_Intron <- function(clipDir, outDir, clipFactors, txdb){
   peakfiles <- match_samples(path=clipDir, pattern="\\.bed$", 
                              samples=clipFactors, start = FALSE)
   print(peakfiles)
   print(clipFactors)
   names(peakfiles) <- clipFactors
   
   setwd(outDir)
   op <- paste0("CLIP_peaks_at_splice_junction_", paste0(clipFactors, collapse = "_"))
   plot_start_end(
      queryFiles = peakfiles,
      inputFiles = NULL,
      centerFiles = "intron",
      txdb = txdb,
      importParams = setImportParams(),
      binSize = 10,
      insert = 100,
      verbose = TRUE,
      ext = c(-100, 100, -100, 100),
      hl = c(-50, 50, -50, 50),
      stranded = TRUE,
      scale = FALSE,
      smooth = TRUE,
      rmOutlier = 0,
      outPrefix = file.path(outDir, op),
      transform = NA,
      shade = TRUE,
      Ylab = "Coverage/base/gene",
      hw = c(10, 10),
      nc = 5
   )
   
}

run_GenomicPlot_metagene <- function(clipDir, outDir, clipFactors, metaFeature){
   peakfiles <- match_samples(path=clipDir, pattern="\\.bed$", samples=clipFactors)
   names(peakfiles) <- clipFactors
   
   setwd(outDir)
   op <- paste0("CLIP_peaks_metagene_plot_", paste0(clipFactors, collapse = "_"))
   
   plot_5parts_metagene(
      queryFiles = peakfiles,
      gFeatures_list = list(meta=metaFeature),
      inputFiles = NULL,
      importParams = setImportParams(),
      verbose = TRUE,
      transform = NA,
      smooth = TRUE,
      scale = FALSE,
      stranded = TRUE,
      outPrefix = file.path(outDir, op),
      heatmap = TRUE,
      heatRange = NULL,
      rmOutlier = 0,
      Ylab = "Coverage/base/gene",
      hw = c(10, 10),
      nc = 5
   )
}

run_GenomicPlot_peak_annotation <- function(clipDir, outDir, clipFactors, gtfFile){
   peakfiles <- match_samples(path=clipDir, pattern="\\.bed$", samples=clipFactors)
   names(peakfiles) <- clipFactors
   
   setwd(outDir)
   op <- paste0("CLIP_peaks_", paste0(clipFactors, collapse = "_"))
   
   plot_peak_annotation(
      peakFile=peakfiles,
      gtfFile=gtfFile,
      importParams = setImportParams(outRle = FALSE),
      fiveP = -1000,
      dsTSS = 300,
      threeP = 1000,
      simple = FALSE,
      outPrefix = file.path(outDir, op),
      verbose = TRUE,
      hw = c(8, 8),
      nc = 5
   )
}


distance_between_proximal_distal <- function(wd, cleavageFactor="gU2AF1", 
                                             feature="3UTR", nc="noChange_"){
    setwd(file.path(wd, cleavageFactor, feature, "bed_output"))
    
    bedfiles <- list.files(pattern="\\.bed$")
    short <- gsub(paste0(paste0(cleavageFactor,"_"),"|\\.bed"), "", bedfiles, fixed = FALSE)
    names(bedfiles) <- short
    
    distance_list <- list()
        
    for(apa in c("shortening_", "lengthening_", nc)){
        selected <- short[grep(apa , short)]
        proximal <- selected[grep("proximal", selected)]
        distal <-selected[grep("distal", selected)]
        
        message(proximal, " ", distal)
        
        proximal_df <- read.delim(bedfiles[proximal], header = FALSE)
        distal_df <- read.delim(bedfiles[distal], header = FALSE)
        
        merged_df <- merge(proximal_df, distal_df, by=c("V1", "V4"))
        
        Distance <- apply(merged_df, 1, function(x){
            abs(pd_distance(as.numeric(x[3:4]), as.numeric(x[7:8])))
        })
        
        APA <- rep(gsub("_", "", apa), length(Distance))
        
        df <- data.frame(cbind(APA, Distance))
        
        distance_list[[apa]] <- df
        
    }
    
    distance_df <- dplyr::bind_rows(distance_list) %>%
        mutate(Distance = as.numeric(Distance)) %>%
        mutate("Log10_Distance" = log10(Distance))
        
    
    p1 <- ggplot(distance_df, aes(x = log10(Distance), color = APA)) +
        geom_density() +
        theme_classic()
    
    p2 <- draw_combo_plot(distance_df,
                          xc = "APA", 
                          yc = "Log10_Distance",
                          Xlab = "APA",
                          Ylab = "Log10(Distance)", 
                          comp = list(c(1, 2), c(1, 3), c(2, 3)), 
                          nf = 1)
    
    pdf(paste0(nc,"distance_between_proximal_distal_cleavae_sites.pdf"), height = 8, width = 8)
    print(p1)
    print(p2)
    dev.off()
}

intronic_CPA <- function(wd, dxr_transcript_file){
    #print(getwd())
    #setwd(wd)
   
   if(!file.exists(file.path(wd,dxr_transcript_file))) stop("dxr file not found!")
   dxr_df <- read.table(file.path(wd, dxr_transcript_file), header=TRUE,  comment.char = "",)
   
   intronic_proximal <- dxr_df %>%
      dplyr::filter(position == "proximal", !feature_name %in% c("3'UTR", "TTS"))
   intronic_df <- dxr_df %>%
      dplyr::filter(groupID %in% intronic_proximal$groupID)
   intronic_distal <- dxr_df %>%
      dplyr::filter(position == "distal", !feature_name %in% c("3'UTR", "TTS"))
   
   UTR3_proximal <- dxr_df %>%
      dplyr::filter(position == "proximal", feature_name %in% c("3'UTR", "TTS"))
   UTR3_df <- dxr_df %>%
      dplyr::filter(groupID %in% UTR3_proximal$groupID)
   UTR3_distal <- dxr_df %>%
      dplyr::filter(position == "distal", feature_name %in% c("3'UTR", "TTS"))
   
   premature_df <- intronic_df %>%
      dplyr::filter(groupID %in% UTR3_distal$groupID)
   
   intron_only_df <- intronic_df %>%
      dplyr::filter(groupID %in% intronic_distal$groupID)
   
   UTR3_only_df <- UTR3_df %>%
      dplyr::filter(groupID %in% UTR3_distal$groupID)
   
   
   length(unique(premature_df$groupID))
   length(unique(intron_only_df$groupID))
   length(unique(UTR3_only_df$groupID))
   
   df_list <- list("premature"=premature_df, 
                   "intron_only"=intron_only_df,
                   "UTR3_only"=UTR3_only_df)
   
   bedDir <- file.path(wd, "bed_output")
   bedfiles <- list.files(path=bedDir, pattern="\\.bed$")
   
   # produce premature bed files
   filter_bedfile <- function(bedDir, bed, df, aname){
       
       shortName <- gsub("\\.bed", "", bed, fixed = FALSE)
       bedfile <- paste0(shortName, "_", aname, ".bed")
       
       if(!file.exists(file.path(bedDir, bedfile))){
           bed_df <- read.delim(file.path(bedDir, bed), header= FALSE)
           
           colnames(bed_df) <- c("chr", "start", "end", "name", "score", "strand")
           print(head(bed_df))
           
           bedout <- inner_join(bed_df, df, by=join_by(chr == genomicData.seqnames,
                                                      start == genomicData.start,
                                                      end == genomicData.end,
                                                      name == groupID,
                                                      strand == genomicData.strand))
           bedout <- bedout[,1:6]
           
           write.table(bedout, file.path(bedDir, bedfile), 
                       col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
       }
   }
   
   for(aname in names(df_list)){
       df <- df_list[[aname]]
       stats <- df %>%
           group_by(change) %>%
           summarize(gene_count = n_distinct(groupID))
       
       write.table(df, file.path(wd, paste0(aname, "_CPA_annotated.tab")), 
                   col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
       write.table(stats, file.path(wd, paste0(aname, "_CPA_stats.tab")), 
                   col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
       
       for (i in seq_along(bedfiles)){
           bed <- bedfiles[i]
           filter_bedfile(bedDir, bed, df, aname)
       }
   }
 
   # make noChange for shortening and lengthening
   setwd(bedDir)
   bedfiles <- list.files(path=bedDir, pattern="\\.bed$")
   
   for(aname in names(df_list)){
      for (p in c("proximal", "distal")){
         
         shortening_file <- bedfiles[grep(paste("shortening", p, aname, sep=".*"), bedfiles, fixed = FALSE)] 
         lengthening_file <- bedfiles[grep(paste("lengthening", p, aname, sep=".*"), bedfiles, fixed = FALSE)] 
         noChange_file <- bedfiles[grep(paste("noChange_", p, aname, sep=".*"), bedfiles, fixed = FALSE)]
         
         lens <- sapply(c(shortening_file, lengthening_file, noChange_file), R.utils::countLines)
         nc_df <- read.delim(noChange_file, header = FALSE)
         nc_shortening <- nc_df[sample(lens[3], lens[1]), ]
         nc_lengthening <- nc_df[sample(lens[3], lens[2]), ]
         
         noChangeS_file <- gsub("noChange", "noChangeS", noChange_file)
         noChangeL_file <- gsub("noChange", "noChangeL", noChange_file)
         
         write.table(nc_shortening, noChangeS_file, col.names = F, row.names = F, sep="\t", quote = F)
         write.table(nc_lengthening, noChangeL_file, col.names = F, row.names = F, sep="\t", quote = F)
      }
   }
}
