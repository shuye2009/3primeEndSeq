
# script.dir <- dirname(sys.frame(1)$ofile)
# source(file.path(script.dir, "GenomicPlot.R"))

setupPackages <- function(install = FALSE) {
  print("Prepare for required packages")
  list.of.packages <- c(
   # "BiocManager",
    "matrixStats",
    "parallel",
    "data.table",
    "MatrixGenerics",
    "tidyr",
    "dplyr",
    "DESeq2",
    "magrittr",
    "cowplot",
    "rtracklayer",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "BSgenome.Hsapiens.UCSC.hg19",
    "genomation",
    "GenomicRanges",
    "plyranges",
    "GenomicFeatures",
    "GenomicPlot",
    "GenomicAlignments",
    "AnnotationHub",
    "gridExtra",
    "ggplot2",
    "ggsignif",
    "ggpubr",
    "R.utils",
    "hrbrthemes",
    "pheatmap",
    "scales",
    "RMariaDB",
    "RCAS",
    "tictoc",
    "BiocStyle",
    "factR",
    "chromstaR",
    "ggpval",
    "ggsci",
    "Repitools",
    "VennDiagram",
    "ggplotify",
    "ComplexHeatmap",
    "forcats",
    "circlize",
    "viridis",
    "DEXSeq",
    "DescTools"
  )
  if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")


  for (apackage in list.of.packages) {
    print(apackage)
    if (apackage %in% installed.packages()) {
      print(paste(apackage, "installed"))
      library(apackage, character.only = T)
    } else if (install) {
      print(paste("trying to install", apackage))
      try(install.packages(apackage))

      if (!apackage %in% installed.packages()) {
        BiocManager::install(apackage)
      }

      if (apackage %in% installed.packages()) {
        library(apackage, character.only = T)
        print(paste(apackage, "installed"))
      } else {
        warning(paste("could not install", apackage))
      }
    }


    # suppressPackageStartupMessages(
    # library(package.i, character.only = TRUE)
  }
}


setupPackages(install=FALSE)

setupProjects <- function(genome = "hg19", baseDir = ".", overwrite = FALSE) {

  if(genome == "hg19"){
     dir <- file.path(baseDir, "genomic_feature")
     gtffile <- file.path(dir, "gencode.v19.annotation.gtf")
  }else if(genome == "hg38"){
     dir <- file.path(baseDir, "genomic_feature")
     gtffile <- file.path(dir, "gencode.v45.primary_assembly.basic.annotation.gtf")
   }else if(genome == "mm9"){
     dir <- file.path(baseDir, "genomic_feature_mouse")
     gtffile <- file.path(dir, "gencode.mouse.v1.annotation.gtf")
  }else if(genome == "mm10"){
     dir <- file.path(baseDir, "genomic_feature_mouse")
     gtffile <- file.path(dir, "gencode.vM22.annotation.gtf")
  }else {
    stop("genome is not supported")
  }


  if (file.exists(file.path(dir, paste0(genome, "_txdb.sql"))) && !overwrite) {
    txdb <- AnnotationDbi::loadDb(file.path(dir, paste0(genome, "_txdb.sql")))
    print(class(txdb))
  } else {
    txdb <- txdbmaker::makeTxDbFromGFF(gtffile, chrominfo = GenomeInfoDb::Seqinfo(genome = genome))
    AnnotationDbi::saveDb(txdb, file = file.path(dir, paste0(genome, "_txdb.sql")))
  }

  if (file.exists(file.path(dir, paste0(genome, "_metaFeatures.rds"))) && !overwrite) {
    metaFeatures <- readRDS(file.path(dir, paste0(genome, "_metaFeatures.rds")))
    #print(class(metaFeatures))
  } else {
    metaFeatures <- prepare_5parts_genomic_features(txdb = txdb, longest = TRUE, meta = TRUE, nbins = 100, fiveP = -1000, threeP = 1000)
    saveRDS(metaFeatures, file = file.path(dir, paste0(genome, "_metaFeatures.rds")))
  }

  if (file.exists(file.path(dir, paste0(genome, "_geneFeatures.rds"))) && !overwrite) {
    geneFeatures <- readRDS(file.path(dir, paste0(genome, "_geneFeatures.rds")))
    #print(class(metaFeatures))
  } else {
    geneFeatures <- prepare_5parts_genomic_features(txdb = txdb, longest = TRUE, meta = FALSE, nbins = 100, fiveP = -3000, threeP = 3000)
    saveRDS(geneFeatures, file = file.path(dir, paste0(genome, "_geneFeatures.rds")))
  }

  if (file.exists(file.path(dir, paste0(genome, "_metaFeature3s.rds"))) && !overwrite) {
     metaFeatures3 <- readRDS(file.path(dir, paste0(genome, "_metaFeature3s.rds")))
     #print(class(metaFeatures3))
  } else {
     metaFeatures3 <- prepare_3parts_genomic_features(txdb = txdb, longest = TRUE, meta = TRUE, nbins = 100, fiveP = -1000, threeP = 1000)
     saveRDS(metaFeatures3, file = file.path(dir, paste0(genome, "_metaFeature3s.rds")))
  }

  if (file.exists(file.path(dir, paste0(genome, "_geneFeature3s.rds"))) && !overwrite) {
     geneFeatures3 <- readRDS(file.path(dir, paste0(genome, "_geneFeature3s.rds")))
     #print(class(metaFeatures3))
  } else {
     geneFeatures3 <- prepare_3parts_genomic_features(txdb = txdb, longest = TRUE, meta = FALSE, nbins = 100, fiveP = -3000, threeP = 3000)
     saveRDS(geneFeatures3, file = file.path(dir, paste0(genome, "_geneFeature3s.rds")))
  }

  return(list(gtffile = gtffile, txdb = txdb, metaFeatures = metaFeatures, geneFeatures = geneFeatures, metaFeatures3 = metaFeatures3, geneFeatures3 = geneFeatures3))
}

map_gene_id <- function(inputVector, inputIdType = "gene_id", outputIdType = "gene_name", gtfFile) {
  gff <- RCAS::importGtf(saveObjectAsRds = TRUE, filePath = gtfFile)
  gene_info_table <- gr2df(gff) %>%
    filter(type == "gene") %>%
    select(chr, start, end, gene_name, strand, gene_id, gene_status) %>%
    mutate(gene_name_new = case_when(
      !duplicated(gene_name) ~ gene_name,
      duplicated(gene_name) ~ gene_id
    )) %>%
    mutate(gene_name = gene_name_new)

  output <- gene_info_table[gene_info_table[[inputIdType]] %in% inputVector, outputIdType]
  names(output) <- gene_info_table[gene_info_table[[inputIdType]] %in% inputVector, inputIdType]

  output <- output[inputVector] ## restore the order to that of input vector

  invisible(output)
}

#' @title Find the longest transcript the 3' UTR belongs to
#' @description Given a bed file contains the genomic coordinates of 3' UTRs, export another bed file which contains
#' the genomic coordinates of longest transcripts the 3' UTRs belong to.
#'
#' @param utr3File path to the 3' UTR file
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param importParams a list of parameters for \code{handle_input}
#'
#' @return name of the transcripts file
#'
#' @author  Shuye Pu
#'
#' @examples
#'
#' txdb <- AnnotationDbi::loadDb(system.file("data", paste0(genome, "_txdb.sql"), package = "GenomicPlotData"))
#' utr3File <- system.file("data", "test_B.bed", package = "GenomicPlotData")
#'
#' importParams <- list(offset = 0, fix_width = 0, fix_point = "center", useScore = FALSE, outRle = FALSE, norm = FALSE, useSizeFactor = FALSE, genome = "hg19")
#' utr3_to_transcript(utr3File, txdb, importParams)
#'
#' @export utr3_to_transcript
#'
utr3_to_transcript <- function(utr3File, txdb, importParams = NULL) {
  transcriptFile <- gsub("\\.bed", "_TSS\\.bed", basename(utr3File))

  longest_tx <- extract_longest_tx(txdb)
  utrInput <- handle_input(utr3File, importParams)
  utrGR <- utrInput[[1]]$query
  print(paste("number of 3' UTRs", length(utrGR)))

  w <- width(utrGR)

  feature <- transcripts(txdb, use.name = TRUE)
  feature_longest <- feature[names(feature) %in% longest_tx$tx_name]

  feature_overlap <- filter_by_overlaps_stranded(feature_longest, utrGR, minoverlap = min(w))
  print(paste("number of transcripts", length(feature_overlap)))

  export.bed(feature_overlap, transcriptFile)

  return(transcriptFile)
}


#' @title Extract genomic coordinates of retained introns from VAST output
#' @description Get genomic coordinates of introns that are up-, down- or unregulated, and down sample the unregualted group
#' to match the size of up- and down-regulated groups.
#'
#' @param inputfile path to 'Intron_retention_regulated_events.tab', which is produced by the VAST tools
#' @param length_filter an integer detonating the minimum length of introns in the output
#'
#' @return a list of bed file names for the up-, down- or unregulated introns
#'
#' @author Shuye Pu
#' @export extract_event_VAST
#' @note used to called "extract_intron_VAST"
#'
#' @examples
#' vast_dir <- "C:/GREENBLATT/Gio/RNAseq/vast_tools_analysis"
#' extract_event_VAST(file.path(vast_dir, "Intron_retention_regulated_events.tab"), length_filter = 250L)
#'
#'
extract_event_VAST <- function(inputfile, length_filter = 0, event = "Intron") {
  library(ggplot2)

  in_df <- read.delim(inputfile, header = TRUE)
    #%>%mutate(CHANGE = gsub("Intron_retention_", "", GROUP))

  outfiles <- list()
  out_beds <- list()
  out_lens <- list()

  for (change in c("up", "down", "unregulated")) {
    # change <- "_unregulated"
    print(change)

    outfiles[[change]] <- gsub("\\.tab", paste0("_", change, "\\.bed"), inputfile)

    sub_df <- dplyr::filter(in_df, CHANGE == change)

    fullcoor <- lapply(as.list(sub_df$FullCO), function(x) {
      unlist(strsplit(x, split = ":", fixed = TRUE))
    })
    fullcoor_strand <- unlist(lapply(fullcoor, function(x) x[3]))

    introncoor <- lapply(as.list(sub_df$COORD), function(x) {
      unlist(strsplit(x, split = ":|-", fixed = FALSE))
    })
    introncoor_chr <- unlist(lapply(introncoor, function(x) x[1]))
    introncoor_start <- unlist(lapply(introncoor, function(x) x[2]))
    introncoor_end <- unlist(lapply(introncoor, function(x) x[3]))


    out_bed <- data.frame("chr" = introncoor_chr, "start" = as.integer(introncoor_start), "end" = as.integer(introncoor_end), "id" = sub_df$EVENT, "score" = sub_df$MV.dPsi._at_0.95, "strand" = fullcoor_strand)
    # write.table(out_bed, outputfile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    out_beds[[change]] <- out_bed
    out_lens[[change]] <- out_bed$end - out_bed$start + 1
  }

  ## test if up and down have different intron length distribution from the unregulated

  len <- c(out_lens[["up"]], out_lens[["down"]], out_lens[["unregulated"]])
  groups <- factor(rep(c("Up", "Down", "Unregulated"), c(length(out_lens[["up"]]), length(out_lens[["down"]]), length(out_lens[["unregulated"]]))))

  len_df <- data.frame(Length = len, Group = groups)
  p1 <- draw_mean_se_barplot(len_df, xc = "Group", yc = "Length", Ylab = paste(event, "length"), comp = list(c(1, 2), c(1, 3), c(2, 3)))
  p2 <- draw_boxplot_by_factor(len_df, xc = "Group", yc = "Length", Ylab = paste(event, "length"), comp = list(c(1, 2), c(1, 3), c(2, 3)), stats = "wilcox.test")

  print(p1)
  print(p2)


  ## filter final output using length_filter
  for (change in c("up", "down", "unregulated")) {
    if (length_filter >= 0) {
      out_beds[[change]] <- data.frame(out_beds[[change]]) %>% filter((end - start + 1) > length_filter)
    }else{
      stop("length_filter must be a whole number >= 0!")
    }
    write.table(out_beds[[change]], outfiles[[change]], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }

  unlist(outfiles)
}


test_first_intron <- function(txdb, intronFiles, importParams) {
  introns <- get_genomic_feature_coordinates(txdb,
    featureName = "intron", featureSource = NULL, export = FALSE,
    longest = FALSE, protein_coding = FALSE
  )$GRangesList

  cl <- start_parallel(40)

  clusterExport(cl, varlist = c("runValue", "strand"), envir = .GlobalEnv)
  first_introns <- parallel::parLapply(cl, introns, function(x) {
    if (length(x) > 0) {
      if (unique(runValue(strand(x))) == "+") {
        return(x[1])
      } else {
        return(x[length(x)])
      }
    } else {
      return(NA)
    }
  })
  stop_parallel(cl)

  first_introns <- as(first_introns[!is.na(first_introns)], "GRangesList")

  sink("stats_of_first_introns.txt")
  for (intronFile in intronFiles) {
    names(intronFile) <- names(intronFiles)[match(intronFile, intronFiles)]
    bed_gr <- handle_input(intronFile, importParams)[[1]]$query

    overlap <- GenomicRanges::findOverlaps(bed_gr, first_introns, ignore.strand = FALSE)

    nquery <- length(unique(queryHits(overlap)))
    ratio <- round(nquery / length(bed_gr), 2)

    print(paste("name of input file", names(intronFile)))
    print(paste("number of query introns", length(bed_gr)))
    print(paste("number of first introns in query", nquery))
    print(paste("ratio of first introns in query", ratio))
    cat("\n\n")
  }
  sink()
}

#' @title Extract genomic coordinates of cassette (CE) exons from VAST output
#' @description Extract genomic coordinates of CE from dPSI file and for each target, partition CE into up-,
#' down- and unregulated based on information in inputfile
#'
#' @param dpsifile path to dPSI file
#' @param inputfile path to inputfile
#' @param targets a character vector providing the names of treatments
#' @param length_filter an integer detonating the minimum length of exons in the output
#'
#' @author Shuye Pu
#'
#' @export extract_exon_dPSI


extract_exon_dPSI <- function(dpsifile, inputfile, targets = NULL, length_filter = 0) {
  for (target in targets) {
    diff_df <- read.table(inputfile, header = TRUE) %>%
      filter(if_any(contains(target), ~ . == TRUE))

    dpsi_df <- read.delim(dpsifile, header = TRUE) %>%
      filter(TYPE == "CE") %>%
      select(c(1:7, contains(target))) %>%
      # filter(EVENT %in% rownames(diff_df)) %>%
      filter(LENGTH > length_filter)

    bed <- apply(dpsi_df, 1, function(x) {
      y <- as.character(x[5])
      chr <- strsplit(y, split = ":", fixed = TRUE)[[1]][1]
      fullco <- strsplit(y, split = ":", fixed = TRUE)[[1]][2]
      C1 <- strsplit(fullco, split = ",", fixed = TRUE)[[1]][1]
      D <- strsplit(C1, split = "+", fixed = TRUE)[[1]][1]
      exon <- strsplit(x[3], split = ":", fixed = TRUE)[[1]][2]
      C2 <- strsplit(fullco, split = ",", fixed = TRUE)[[1]][3]
      A <- strsplit(C2, split = "+", fixed = TRUE)[[1]][1]

      start <- strsplit(exon, split = "-", fixed = TRUE)[[1]][1]
      end <- strsplit(exon, split = "-", fixed = TRUE)[[1]][2]

      # print(paste(C1, starts, ends, C2))
      strand <- ifelse(as.numeric(D) < as.numeric(A), "+", "-")
      return(c(chr, start, end, x[2], x[8], strand)) # score column is the dPSI value
    })

    bed <- data.frame(matrix(bed, ncol = 6, byrow = TRUE))
    colnames(bed) <- c("chr", "start", "end", "ID", "score", "strand")

    target_bed <- bed %>%
      filter(ID %in% rownames(diff_df))
    target_bed_up <- target_bed %>%
      filter(as.numeric(score) > 0)
    target_bed_down <- target_bed %>%
      filter(as.numeric(score) < 0)
    nontarget_bed <- bed %>%
      filter(!ID %in% rownames(diff_df))

    write.table(target_bed_up, gsub(".tab", paste0("_", target, "_up.bed"), inputfile), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(target_bed_down, gsub(".tab", paste0("_", target, "_down.bed"), inputfile), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(target_bed, gsub(".tab", paste0("_", target, ".bed"), inputfile), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    write.table(nontarget_bed, gsub(".tab", paste0("_non_", target, ".bed"), inputfile), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
}


#' @title Perform differential analysis using DESeq2
#' @description Perform differential analysis using the DESeq2 package
#'
#' @param data_table a numeric dataframe, if the type is double, the numbers will be rounded to integers
#' @param s2c a dataframe providing grouping information of the samples, sample names must be the column names of data_table
#' @param project a string denoting the prefix of output files
#' @param contrasts a dataframe denoting the comparisons, each row is in the form of c("groups", "numerator", "denominator")
#' @param plot logical, indicating whether to generate PCA plot and MA plot
#' @param useSizeFactor logical, indicating whether to estimate size factors, which are used to normalize the data_table
#'
#' @return a list of dataframes containing results of the analysis
#'
#' @author Shuye Pu
#' @export run_DESeq2

run_DESeq2 <- function(data_table, s2c, project, contrasts, plot = FALSE, useSizeFactor = TRUE) {
  ddsMatrix <- DESeqDataSetFromMatrix(round(data_table),
    colData = s2c,
    design = ~groups
  )
  class(ddsMatrix)

  if (!useSizeFactor) {
    sizeFactors(ddsMatrix) <- rep(1, ncol(data_table))
  }

  dds <- DESeq(ddsMatrix, fitType = "local")

  print(paste(project, "size factors"))
  print(sizeFactors(dds))

  res_list <- NULL

  for (i in 1:nrow(contrasts)) {
    # i <- 1
    contrast <- as.vector(as.character(contrasts[i, ]))
    print(paste("Processing ", paste(contrast, collapse = " ")))
    res <- results(dds, contrast = contrast, independentFiltering = TRUE, pAdjustMethod = "BH")
    treatSF <- paste(project, contrast[2], "vs", contrast[3], sep = "_")

    summary(res)
    res <- na.omit(res)

    mcols(res)

    res <- res[order(-res$stat), ]
    res_list[[treatSF]] <- res

    filename <- paste(treatSF, "DESeq2_results_table.tab", sep = "_")
    write.table(res, filename, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")

    if (plot) {
      pdf(paste(treatSF, "MA_plot.pdf", sep = "_"))
      plot(log2(res$baseMean), res$log2FoldChange, ylim = c(-5, 5))
      abline(h = 0)
      dev.off()
    }
  }
  return(res_list)
}

#' @title Convert genomic coordinates between genome versions
#' @description Take a bed file in one genome version and produce another bed file in another genome version
#'
#' @param chainFile path to chain file
#' @param queryfile path to a bed file to be converted
#'
#' @return path to output bed file
#' @author Shuye Pu
#' @export liftover_bed

liftover_bed <- function(chainFile, queryfile, fromGenome = "GRCh38", toGenome = "hg19") {
   if (grepl("\\.bed$", queryfile)) {
      outFile <- gsub("\\.bed$", paste0("_", toGenome, "\\.bed"), queryfile)
   } else {
      outFile <- gsub("\\.narrowPeak$", paste0("_", toGenome, "\\.narrowPeak"), queryfile)
   }

  chain <- rtracklayer::import.chain(chainFile)
  bedGr <- handle_bed(queryfile, setImportParams(useScore = TRUE,
                                                 outRle = FALSE,
                                                 genome = fromGenome
  ))$query

  lifted <- unlist(rtracklayer::liftOver(bedGr, chain))
  # mcols(lifted) <- mcols(lifted)[,2:3] ## remove the name column of numbers
  plyranges::write_bed(lifted, outFile)
  return(outFile)
}

#' @title Sub-sampling of vector elements
#' @description Given a large numeric vector and small numeric vector, produce another
#' vector that are similar in size and distribution to the small vector by random sampling of the large vector.
#'
#' @param v1 large vector to be sub-sampled
#' @param v2 small vector to provide density distribution
#' @param pv wilcox.test p-value to indicate similarity between the sub-sample and v2
#' @param maxiter maximum number of trials before the program stops
#'
#' @return a numeric vector
#' @author Shuye Pu
#'
#' @examples
#'
#' v1 <- rnorm(10000)
#' v2 <- rnorm(100, mean = 0.5, sd = 2)
#' pv <- 0.9
#'
#' out <- sub_sample(v1, v2, pv)
#'
#' @export sub_sample
#'

sub_sample <- function(v1, v2, pv = 0.5, maxiter = 10000) {
  ds <- density(v2)
  dv1 <- approx(x = ds$x, y = ds$y, xout = v1, rule = 1)
  probs <- dv1$y
  probs[is.na(probs)] <- 0
  summary(probs)

  success <- FALSE
  wtest <- NULL
  n <- 0
  vsub <- NULL

  ## find subset of vi that has the same distribution as v2
  while (!success) {
    n <- n + 1
    # vsub <- sample(v1, round(length(probs[probs>0])*0.99), replace=FALSE, prob=probs)
    vsub <- sample(v1, length(v2), replace = FALSE, prob = probs)
    wtest <- wilcox.test(vsub, v2)
    ktest <- ks.test(vsub, v2)

    print(paste(n, wtest$p.value, ktest$p.value))
    if (wtest$p.value > pv) {
      # if(lr < 1){
      success <- TRUE
      print("length of vsub")
      print(length(vsub))
      print("density distribution of v2")
      print(summary(v2))
      print("subsampling of v1 based on density distribution of v2")
      print(summary(vsub))
    }

    if (n > maxiter) {
      warning("Unable to subsample based on the given distribution and p-value! try with a lower p-value.")
      break
    }
  }

  wtest <- wilcox.test(v1, v2)
  print("p-value for v1 and v2")
  print(wtest)

  return(vsub)
}


#' @examples
#' v1 <- rnorm(50)
#' names(v1) <- paste0("N", seq(50))
#' v2 <- rnorm(10)
#' names(v2) <- paste0("P", seq(10))
#' oneTOmany <- FALSE
#' sub_sample_by_matching(v1, v2, oneTOmany)
#' oneTOmany <- TRUE
#' sub_sample_by_matching(v1, v2, oneTOmany)

sub_sample_by_matching <- function(v1, v2, oneTOmany = FALSE) {

  v1 <- v1[!is.na(v1)]
  v2 <- v2[!is.na(v2)]

  v_working <- v1

  onematched <- list()
  othermatched <- list()
  for(i in seq_along(v2)){
     v <- v2[i]
    diff <- abs(v_working - v)
    min_diff <- min(diff)
    candidates <- diff[diff == min_diff]

    onematched[[i]] <- candidates[1]
    if(length(candidates) > 1){
      othermatched[[i]] <- candidates[2:length(candidates)]
    }else{
      othermatched[[i]] <- NULL
    }

    candidate <- candidates[1]
    if(oneTOmany) candidate <- candidates
    v_working <- v_working[!names(v_working) %in% candidate]
    #print(length(v_working))
  }

  onematched <- names(unlist(onematched, use.names = TRUE))
  othermatched <- names(unlist(othermatched, use.names = TRUE))

  matched <- onematched
  if(oneTOmany) matched <- c(onematched, othermatched)

  matched_val <- v1[matched]

  #message("length of V1: ", length(v1), " length of v2: ", length(v2), " length of matched: ", length(matched))

  return(unique(matched))
}

#' @title Identify m6Am sites
#' @description Genomic m6Am sites are identified from miCLIP peaks based on the following rules: 1, the center of the peaks
#' must contains the nucleotide 'A', 2, the peaks must be located within the first 1/4 fraction of the 5' UTR, 3, the peaks do
#' not contain DRACH motif (m6A sites).
#'
#' @param bedfile path to the "A" containing miCLIP peaks bed file
#' @param m6Afile path to the m6A sites bed file
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param m6Amfile path to the output m6Am bed file
#'
#' @return NULL
#' @author Shuye Pu
#'
#' @export derive_m6Am_sites

derive_m6Am_sites <- function(bedfile, m6Afile, txdb, m6Amfile, importParams=NULL) {
  ## get "A" containing peaks in 5' UTR
  utr5_filtered <- filter_bed_by_genomicFeature(bedfile, featureName = "utr5", txdb = txdb, resizeFraction = 0.25, longest = TRUE, importParams = importParams)
  m6A <- rtracklayer::import.bed(m6Afile)
  ## remove peaks that are overlapping with m6A sites
  m6Am <- filter_by_nonoverlaps_stranded(utr5_filtered$query, m6A)
  ## convert to 6-column bed format
  utr5_filtered_bed <- gr2df(m6Am) %>%
    select(c(chr, start, end, name, score, strand)) %>%
    mutate(start = as.integer(start) - 1) %>%
    mutate(chr = as.character(chr))
  print(dim(utr5_filtered_bed))
  write.table(utr5_filtered_bed, m6Amfile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}



#' @title Filter peaks by differential analysis of read counts
#'
#' @description  Find peaks that have or lack read count support by comparing treatBam and refBam. In the case of filtering CLIP
#' peaks that lack Input support, CLIP peaks whose read counts in treatBam is not significantly higher than in refBam is excluded.
#' Also can be used to find the peaks that disappear based on read counts change after a treatment.
#' @param peakFile, a string denoting the peak file name, only .bed format is allowed. The peakfile can be associated with either
#' the treatBam or the refBam.
#' @param treatBam, a vector of strings denoting the path to the bam files used as numerator
#' @param refBam, a vector of strings denoting the path to the bam files used as denominator
#' @param treatName a string denoting the name of numerator
#' @param refName a string denoting the name of numerator
#' @param refDir path to the rds file storing refBam data, used to reduce I/O time when the number of refBams is large
#' @param pcutoff a numeric for pvalue cutoff during differential analysis
#' @param importParams a list of parameters for \code{handle_input}
#' @param verbose logical, whether to output differential analysis results and bed files derived from peakfile split into
#' 'above', 'below' and 'noChange' based on differential analysis
#'
#' @return a list with the elements c("pos", "neg", "res"), where "pos" and "neg" indicate the index of the peaks that have
#' significantly higer or lower read counts in treatBam than in refBam, respectively.
#' @author Shuye Pu
#'
#'
#' @export filter_peaks_by_differential
#'
filter_peaks_by_differential <- function(peakFile, treatBam, refBam, treatName = NULL, refName = NULL, refDir = NULL, importParams = NULL, pcutoff = 0.5, verbose = FALSE) {
  bamparam <- importParams
  bamparam$CLIP_reads <- FALSE
  bamparam$fix_width <- 0
  bamparam$useScore <- FALSE
  bamparam$outRle <- FALSE
  bamparam$norm <- FALSE
  bamparam$useSizeFactor <- FALSE

  if (!is.null(refDir)) {
    if (file.exists(file.path(refDir, "refBams.rds"))) {
      refBams <- readRDS(file.path(refDir, "refBams.rds"))
    } else {
      refBams <- handle_input(refBam, bamparam)
      saveRDS(refBams, file = file.path(refDir, "refBams.rds"))
    }
  } else {
    refBams <- handle_input(refBam, bamparam)
  }

  treatBams <- handle_input(treatBam, bamparam)

  bedparam <- importParams
  bedparam$CLIP_reads <- FALSE
  bedparam$useScore <- FALSE
  bedparam$outRle <- FALSE
  bedparam$norm <- FALSE
  bedparam$useSizeFactor <- FALSE
  peak <- handle_input(peakFile, bedarams)

  overlap_list <- lapply(c(treatBams, refBams), function(x) {
    ol <- countOverlaps(peak[[1]]$query, x$query, type = "any", maxgap = -1L)
  })

  overlap_mat <- do.call(cbind, overlap_list)

  head(overlap_mat)

  samples <- c(treatBam, refBam)
  groups <- factor(c(rep("Treat", length(treatBam)), rep("Reference", length(refBam))))
  s2c <- data.frame(samples, groups)
  rownames(s2c) <- samples
  ddsMatrix <- DESeqDataSetFromMatrix(overlap_mat,
    colData = s2c,
    design = ~groups
  )

  dds <- DESeq(ddsMatrix, fitType = "local")
  count_table <- DESeq2::counts(dds, normalized = TRUE)

  sf <- sizeFactors(dds)
  print(sf)

  ## Treat is the numerator, Reference is the denominator in the fold change
  res <- DESeq2::results(dds, contrast = c("groups", "Treat", "Reference"), independentFiltering = TRUE, pAdjustMethod = "BH")
  rownames(res) <- seq_along(res$pvalue)
  pos <- rownames(res)[!is.na(res$pvalue) & res$pvalue < pcutoff & sign(res$log2FoldChange) == 1]
  neg <- rownames(res)[!is.na(res$pvalue) & res$pvalue < pcutoff & sign(res$log2FoldChange) == -1]
  noChange <- rownames(res)[!rownames(res) %in% c(pos, neg)]
  nonpos <- rownames(res)[!rownames(res) %in% pos]

  out_table <- cbind(as.data.frame(res), as.data.frame(count_table))

  if (verbose) {
    res_list <- list("above" = pos, "below" = neg, "noChange" = noChange, "notAbove" = nonpos)

    peakbed <- read.delim(peakFile, header = FALSE, sep = "\t")

    print(dim(peakbed))
    lapply(names(res_list), function(x) {
      ind <- res_list[[x]]
      bed <- peakbed[ind, ]
      print(dim(bed))
      outfile <- file.path(dirname(peakFile), paste0(treatName, "_", x, "_", refName, ".bed"))
      write.table(bed, outfile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    })

    outfile <- file.path(dirname(peakFile), paste0(treatName, "_against_", refName, "_differential_results.tab"))
    write.table(out_table, outfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

  invisible(list("pos" = pos, "neg" = neg, "res" = out_table)) # return the order of significant peaks
}


#' @title Split GRangesList by chromosomes
#' @description Only work for GRangesList that have uniform chromosome in each list element, i.e. all intervals in each GRanges object must have identical sequence names.
#'
#' @param GRL a GRangesList object
#'
#' @return a list of GRangesLists
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' txdb <- AnnotationDbi::loadDb(system.file("data", paste0(genome, "_txdb.sql"), package = "GenomicPlotData"))
#' cds <- get_genomic_feature_coordinates(txdb, featureName = "cds", featureSource = "gtf", export = FALSE, longest = TRUE, protein_coding = TRUE)$GRangesList
#'
#' out <- split_grl_by_chr(cds)
#'
#' @export split_grl_by_chr

split_grl_by_chr <- function(GRL) {
  chrs <- seqlevels(GRL)
  chrList <- sapply(seqnames(GRL), runValue)
  out <- lapply(chrs, function(x) GRL[names(chrList)[chrList == x]])
  names(out) <- chrs
  invisible(out)
}

#' @title Parallel execution of scoreMatrixBin on multiple target windows
#'
#' @description Function for parallel execution of scoreMatrixBin. The 'windows' parameter of the scoreMatrixBin method is provided as a list of
#' GRangesList objects, allowing simultaneous processing of multiple windows objects.
#'
#' @param windowRs, a list of GRangesList objects.
#' @param queryRegions, a RleList object or Granges object providing input for the 'target' parameter of the scoreMatrixBin method
#' @param bin_num, number of bins the windows should be divided into
#' @param bin_op, operation on the signals in a bin, a string in c("mean", "max", "min", "median", "sum") is accepted.
#' @param weight_col, if the queryRegions is a GRanges object, a numeric column in meta data part can be used as weights.
#' @param stranded, logical, indicating if the strand of the windows should be considered to determine upstream and downstream
#'
#' @return a list of numeric matrices
#' @author Shuye Pu
#'
#'
#' @export parallel_apply_scoreMatrixBin
#'
#'
parallel_apply_scoreMatrixBin <- function(windowRs, queryRegions, bin_num, bin_op, weight_col, stranded) {
  nc <- length(windowRs)
  cl <- start_parallel(nc)
  clusterExport(cl, varlist = c("ScoreMatrixBin"), envir = environment())
  clusterExport(cl, varlist = c("queryRegions", "bin_op", "bin_num", "weight_col", "stranded"), envir = environment())
  call_scoreMatrixBin <- function(windowR) {
    smc <- ScoreMatrixBin(target = queryRegions, windows = windowR, bin.num = bin_num, bin.op = bin_op, weight.col = weight_col, strand.aware = stranded)
    smc@.Data
  }

  smcList <- parLapply(cl, windowRs, call_scoreMatrixBin)
  stop_parallel(cl)

  invisible(smcList)
}


#' @title Filter a bed file based on its overlap with a genomic feature
#' @description Takes a bed file and a genomic feature name, produce a list of
#' 2 GRanges objects
#'
#' @param bedfile path to a bed file
#' @param featureName a genomic feature name in c("utr3", "utr5", "cds",
#' "intron", "exon", "transcript", "gene")
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param importParams a list of parameters for \code{handle_input}
#' @param longest logical, indicating whether the output should be limited to
#' the longest transcript of each gene
#' @param protein_coding logical, indicating whether the output should be
#' limited to the protein_coding genes
#' @param resizeFraction a number between 0 and 1 denoting what fraction of the
#' feature should be used for overlap
#' @param featureGr a GRanges object providing a custom feature
#' @param featureBed a bed file providing a custom feature
#' @param maxgap minimum gap for defining overlap
#' @param verbose logical, indicate whether the results should be write to a file
#' @param stranded logical, indicate whether the strand should be considered
#'
#' @return a list of GRanges objects, the first is the overlapping bed, the
#' second is the overlapping feature
#'
#' @author Shuye Pu
#'
#' @export filter_bed_by_genomicFeature
#'
filter_bed_by_genomicFeature <- function(bedfile, featureName = "custom",
                                         txdb = NULL, importParams = NULL,
                                         resizeFraction = NULL, longest = TRUE,
                                         protein_coding = TRUE, featureGr = NULL,
                                         featureBed = NULL, maxgap = -1L,
                                         verbose = FALSE, stranded = TRUE) {
  if (featureName %in% c("utr3", "utr5", "cds", "intron", "exon", "transcript",
                         "gene")) {
    feature <- get_genomic_feature_coordinates(txdb, featureName,
                                               longest = longest)
    featureGr <- feature$GRanges
  } else if (is.null(featureGr)) {
    if (!is.null(featureBed)) {
      featureGr <- handle_bed(featureBed, importParams)$query
    } else {
      stop("Both featureGr and featureBed are null!")
    }
  } else {
    stop("Neither featureName nor featureGr nor featureBed is provided!")
  }

  if (!is.null(resizeFraction)) {
    new_w <- as.integer(width(featureGr) * resizeFraction)
    featureGr <- resize(featureGr, width = new_w, fix = "start",
                        use.names = TRUE, ignore.strand = FALSE)
  }

  bedGr <- handle_bed(bedfile, importParams)$query

  # filtered <- filter_by_overlaps_stranded(bedGr, featureGr, maxgap=maxgap)
  overlaps <- GenomicRanges::findOverlaps(bedGr, featureGr, maxgap = maxgap,
                                          ignore.strand = !stranded)

  print(paste(length(unique(overlaps@from)), "queries overlap",
              length(unique(overlaps@to)), "subjects"))
  if (verbose) {
    df <- gr2df(bedGr[unique(overlaps@from)])
    bedName <- names(bedfile)
    write.table(df, paste0(bedName, "_in_", featureName, ".bed"),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }

  return(list(query = bedGr[unique(overlaps@from)],
              feature = featureGr[unique(overlaps@to)]))
}


#' @title Get genomic coordinates of TSS that starts with a specific nucleotide
#' @description Get genomic coordinates of transcription start sites that may or may not start with a specific nucleotide. If 'longest' is true,
#' the TSS is limit to longest transcripts, otherwise all transcript. if protein-coding is true, limit to transcripts of protein-coding genes.
#'
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param protein_coding logical, indicating whether the output should be limited to protein-coiding genes only
#' @param bsg a `BSgenome` object
#' @param startWith a character vector in c("A", "C", "G", "T"), indicating whether the output should be limited to a specific nucleotides
#' @param export logical, indicating whether to write a bed file
#'
#' @return a GRanges object
#'
#' @author Shuye Pu
#' @export get_gr_at_TSS

get_gr_at_TSS <- function(txdb, longest = TRUE, protein_coding = TRUE,
                          bsg = BSgenome.Hsapiens.UCSC.hg19,
                          startWith = c("A"), export = FALSE) {
  tx <- get_genomic_feature_coordinates(txdb, featureName = "transcript",
                                        featureSource = NULL, export = FALSE,
                                        longest = longest,
                                        protein_coding = protein_coding)

  tss_tx <- resize(tx$GRanges, width = 1, fix = "start", use.names = FALSE,
                   ignore.strand = FALSE)
  TSS <- lapply(startWith, function(A) {
    seq_tss_tx <- BSgenome::getSeq(bsg, tss_tx, as.character = TRUE)
    A_tss_tx <- tss_tx[which(seq_tss_tx == A)]
    A_tss_tx
  })

  TSSgrl <- as(TSS, "GRangesList")
  TSSgr <- unlist(TSSgrl)

  if(export){
    outfile <- paste0("TSS_start_with_",
                      paste(startWith, collapse = ""), ".bed")
    mc <- mcols(TSSgr)
    mc <- subset(mc, select = -c(blocks))
    mcols(TSSgr) <- mc
    rtracklayer::export.bed(TSSgr, outfile)
    #TSSdf <- as.matrix(gr2df(TSSgr))
    ## write.table complains that the dataframe
    ## is a list, as.matrix is a work around
    #write.table(TSSdf, outfile, row.names = FALSE, col.names = FALSE,
     #           sep = "\t", quote = FALSE)
  }

  return(TSSgr)
}

#' @title Count reads or peaks overlapping TSS
#' @description Count the number of reads or peaks that overlapping with TSS
#'
#' @param queryfiles path to bam or bed files
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param protein_coding logical, indicating whether the output should be limited to protein-coiding genes only
#' @param startWith a character in c("A", "C", "G", "T"), indicating whether the output should be limited to transcripts starting with a specific nucleotide
#' @param TSSflank an integer denoting the number of nucleotide to include upstream and downstream of TSS
#' @param importParams a list of parameters for \code{handle_input} of queryfiles
#' @param norm logical, indicating if the intensities should be normalized by the sample library sizes
#'
#' @return a dataframe
#' @author Shuye Pu
#'
#' @examples
#'
#' txdb <- AnnotationDbi::loadDb(system.file("data", paste0(genome, "_txdb.sql"), package = "GenomicPlotData"))
#' queryfiles <- system.file("data", "test_clip.bam", package = "GenomicPlotData")
#' names(queryfiles) <- "query"
#' inputfiles <- system.file("data", "test_clip_input.bam", package = "GenomicPlotData")
#' names(inputfiles) <- "input"
#'
#' queryInputParams <- list(offset = 0, fix_width = 100, norm = TRUE, useScore = FALSE, outRle = FALSE, useSizeFactor = TRUE, genome = "hg19")
#'
#' res <- count_overlap_at_TSS(c(queryfiles, inputfiles), txdb, TSSflank = 10, importParams = queryInputParams, norm = TRUE, longest = TRUE, protein_coding = TRUE, startWith = "N")
#'
#' @export count_overlap_at_TSS

count_overlap_at_TSS <- function(queryfiles, txdb, TSSflank = 0, importParams = NULL, norm = FALSE, longest = TRUE, protein_coding = TRUE, startWith = "N") {
  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles

  input_list <- handle_input(inputFiles = queryfiles, importParams)

  print("Getting coordinates of TSS")
  A_tss <- get_gr_at_TSS(txdb, longest = longest, protein_coding = protein_coding, startWith = startWith)
  gr <- unique(flank(A_tss, width = TSSflank, both = TRUE))
  trim(gr)

  print("Counting overlaps")
  count_list <- lapply(queryfiles, function(queryfile) {
    queryRegions <- input_list[[queryfile]]$query
    overlap_count <- countOverlaps(gr, queryRegions, type = "any")
  })

  count_df <- bind_cols(count_list)
  if (norm) {
    print("Normalizing ...")
    sizes <- sapply(input_list, "[[", "size")
    count_df <- count_df * 1e6 / sizes
    print(sizes)
  }
  colnames(count_df) <- querylabels
  rownames(count_df) <- format_genomic_oordinates(gr)

  return(count_df)
}


#'
#' @title Plot count of reads in genomic features and/or peaks
#'
#' @description Plot the outputs of count_overlap_in_features. If sample grouping information is available, log fold change of (treat vs. control) evaluated by DESeq2 will be plotted as well.
#'
#' @param queryfiles paths of bam files for aligned reads
#' @param peakfile path of a bed file for the peaks
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param queryInputParams a list of parameters for \code{handle_input} of queryfiles
#' @param peakInputParams a list of parameters for \code{handle_input} of peakfile
#' @param TSSflank an integer denoting the number of nucleotides encompassing upstream and downstream of the TSS
#' @param sampleNames a character vector denoting sample names of queryfiles
#' @param subsets a list of vectors of sample labels associated with queryfiles, providing information for differential analysis
#' @param longest logical, indicating whether to use the longest transcript of genes
#' @param norm logical, indicating whether the queryfiles should be normalized by library size
#' @param input_type an arbitrary string indicating the type of queryfiles, default = "reads"
#'
#' @return NULL
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' txdb <- AnnotationDbi::loadDb(system.file("data", paste0(genome, "_txdb.sql"), package = "GenomicPlotData"))
#' peakfile <- system.file("data", "test_A.bed", package = "GenomicPlotData")
#' names(peakfile) <- "TestA"
#' queryfiles <- c(system.file("data", "test_clip.bam", package = "GenomicPlotData"), system.file("data", "test_clip_input.bam", package = "GenomicPlotData"))
#' names(queryfiles) <- c("query", "input")
#'
#' queryInputParams <- list(offset = -1, fix_width = 0, norm = TRUE, useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19")
#' peakInputParams <- list(offset = 0, fix_width = 21, fix_point = "center", norm = FALSE, useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19")
#'
#' plot_feature_overlap_count(queryfiles, peakfile = NULL, txdb, sampleNames = c("query", "input"), longest = TRUE, queryInputParams = queryInputParams, peakInputParams = peakInputParams, TSSflank = 5, norm = TRUE, subsets = NULL, input_type = "reads")
#'
#' plot_feature_overlap_count(queryfiles, peakfile, txdb, sampleNames = c("query", "input"), longest = TRUE, queryInputParams = queryInputParams, peakInputParams = peakInputParams, TSSflank = 5, norm = TRUE, subsets = NULL, input_type = "reads")
#'
#' @export plot_feature_overlap_count

plot_feature_overlap_count <- function(queryfiles, peakfile = NULL, txdb,
                                       sampleNames, longest = TRUE,
                                       queryInputParams = NULL,
                                       peakInputParams = NULL, TSSflank = 0,
                                       norm = TRUE, subsets = NULL,
                                       input_type = "reads", hw = c(8,8)) {

  functionName <- as.character(match.call()[[1]])
  params <- GenomicPlot:::plot_named_list(as.list(environment()))
  force(params)

  peaklabel <- ifelse(is.null(peakfile), "genomic_features", names(peakfile))
  input_type <- paste0(input_type, "_in_", peaklabel)

  feature_count_list <- count_overlap_in_features(queryfiles, peakfile, txdb,
                                                  longest = longest,
                                                  queryInputParams = queryInputParams,
                                                  peakInputParams = peakInputParams,
                                                  TSSflank = TSSflank)
  if (norm) {
    rpm <- feature_count_list[["inputsize"]] / 1000000
  } else {
    rpm <- rep(1, length(feature_count_list[["inputsize"]]))
  }
  names(rpm) <- names(queryfiles)
  print(rpm)

  print("plot counts in features")
  pdf(paste("overlap_count_of_", input_type, ".pdf", sep = ""), height = hw[1],
      width = hw[2])
  for (featureName in names(feature_count_list$count)) {
    count_table_feature <- feature_count_list[["count"]][[featureName]]

    count_table_feature <- count_table_feature[apply(
      count_table_feature, 1, sum) > 0, ]
    print(dim(count_table_feature))
    if (is.null(count_table_feature) || nrow(count_table_feature) == 0) next

    count_table_feature <- t(t(count_table_feature) / rpm)

    sum_table_feature <- sapply(sampleNames, function(sampleN) {
      sub_df <- as.data.frame(count_table_feature[, grepl(sampleN, colnames(count_table_feature))])
      asum <- apply(sub_df, 1, sum)
    })
    sum_table_feature <- as.data.frame(log2(sum_table_feature + 1))

    print(featureName)
    print(head(sum_table_feature))

    sum_df <- pivot_longer(sum_table_feature,
                           cols = colnames(sum_table_feature),
                           names_to = "Group", values_to = "Count") %>%
      mutate(Group = as.factor(Group)) %>%
      group_by(Group) %>%
      arrange(Count)

    if(length(levels(sum_df$Group)) < 2) next
    comps <- combn(seq_along(levels(sum_df$Group)), 2, simplify = FALSE)

    outp <- draw_combo_plot(sum_df, xc = "Group", yc = "Count",
                            Xlab = "", title = featureName,
                            Ylab = "log2(Count)", comp = comps)

    print(outp)

  }
  print(params)
  dev.off()


  if (!is.null(subsets)) {
    print("plot lfc for feautures")
    pdf(paste("lfc_at_TSNflank", TSSflank, "_and_transcript_for_", input_type,
              ".pdf", sep = ""), height = hw[1], width = hw[2])
    for (featureName in names(feature_count_list$count)) {
      lfc_list <- list()
      for (subject in names(subsets)) {
        subset <- subsets[[subject]]
        s2c <- data.frame(samples = subset,
                          groups = c(rep(c("GFP", subject), times = c(2, 2))))
        contrasts <- data.frame(groups = c(rep("groups", 1)), treat = subject,
                                control = c(rep("GFP", 1)))

        project <- paste(input_type, TSSflank, sep = "")

        data_table_feature <- feature_count_list[["count"]][[featureName]][, subset]

        if (is.null(count_table_feature) || nrow(data_table_feature) == 0) next
        data_table_feature <- data_table_feature[apply(
          data_table_feature, 1, sum) > 0, ]

        print(rpm[subset])
        data_table_feature <- round(t(t(data_table_feature) / rpm[subset]))
        print(head(data_table_feature))
        data_table_feature[is.na(data_table_feature)] <- 0

        res_feature <- NULL
        X <- try(
          res_feature <- run_DESeq2(data_table_feature, s2c, featureName,
                                    contrasts, plot = FALSE,
                                    useSizeFactor = FALSE)
        )
        if (!grepl("try-error", class(X))) {
          lfc_df <- data.frame(Group = rep(
            subject, length(res_feature[[1]]$log2FoldChange)),
            log2FoldChange = res_feature[[1]]$log2FoldChange)
          lfc_list[[subject]] <- lfc_df
        }
      }

      out_df <- as.data.frame(Reduce(rbind, lfc_list))
      if (nrow(out_df) == 0) next

      colnames(out_df) <- c("Group", "log2FoldChange")
      out_df <- mutate(out_df, Group = as.factor(Group)) %>%
        group_by(Group) %>%
        arrange(log2FoldChange)

      if (length(levels(out_df$Group)) < 2) next

      comps <- combn(seq_along(levels(out_df$Group)), 2, simplify = FALSE)

      outp <- draw_combo_plot(out_df, xc = "Group", yc = "log2FoldChange",
                              Xlab = "", title = featureName, comp = comps)

      print(outp)

    }
    print(params)
    dev.off()
  }
}

plot_feature_overlap_count_old <- function(queryfiles, peakfile = NULL, txdb, sampleNames, longest = TRUE, queryInputParams = NULL, peakInputParams = NULL, TSSflank = 0, norm = TRUE, subsets = NULL, input_type = "reads") {
  peaklabel <- ifelse(is.null(peakfile), "genomic_features", names(peakfile))
  input_type <- paste0(input_type, "_in_", peaklabel)

  feature_count_list <- count_overlap_in_features(queryfiles, peakfile, txdb, longest = longest, queryInputParams = queryInputParams, peakInputParams = peakInputParams, TSSflank = TSSflank)
  if (norm) {
    rpm <- feature_count_list[["inputsize"]] / 1000000
  } else {
    rpm <- rep(1, length(feature_count_list[["inputsize"]]))
  }
  names(rpm) <- names(queryfiles)
  print(rpm)

  print("plot counts in features")
  pdf(paste("overlap_count_of_", input_type, ".pdf", sep = ""), height = 5, width = 10)
  for (featureName in names(feature_count_list$count)) {
    count_table_feature <- feature_count_list[["count"]][[featureName]]

    count_table_feature <- count_table_feature[apply(count_table_feature, 1, sum) > 0, ]
    print(dim(count_table_feature))
    if (is.null(count_table_feature) || nrow(count_table_feature) == 0) next

    count_table_feature <- t(t(count_table_feature) / rpm)

    sum_table_feature <- sapply(sampleNames, function(sampleN) {
      sub_df <- as.data.frame(count_table_feature[, grepl(sampleN, colnames(count_table_feature))])
      asum <- apply(sub_df, 1, sum)
    })
    sum_table_feature <- as.data.frame(log2(sum_table_feature + 1))

    print(featureName)
    print(head(sum_table_feature))

    sum_df <- pivot_longer(sum_table_feature, cols = colnames(sum_table_feature), names_to = "Group", values_to = "Count") %>%
      mutate(Group = as.factor(Group)) %>%
      group_by(Group) %>%
      arrange(Count) %>%
      mutate(cumCount = cumsum(Count)) %>%
      mutate(Rank = order(cumCount)) %>%
      mutate(Fraction = cumCount / max(cumCount), Rank = Rank / max(Rank))

    p1 <- ggplot(data = sum_df, aes(x = Rank, y = Fraction, color = Group)) +
      geom_line() +
      ggtitle("Fingerprint plot") +
      labs(x = "Rank(Count)", y = paste("Fraction over total count"))


    median_df <- aggregate(sum_df$Count, list(sum_df$Group), median)
    colnames(median_df) <- c("Group", "median")

    if (length(levels(sum_df$Group)) < 2) next

    p <- ggplot(sum_df, aes(x = Group, y = Count, fill = Group)) +
      geom_boxplot(notch = FALSE) +
      theme_classic() +
      ggtitle(paste(featureName, "n =", nrow(count_table_feature))) +
      theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 14)) +
      theme(legend.position = "none") +
      labs(y = expression(paste(log[2], " (Count)"))) +
      theme(axis.text.x = element_text(face = "bold", size = 14, color = "black")) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "black") +
      geom_signif(comparisons = combn(levels(sum_df$Group), 2, simplify = FALSE), test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1)


    d <- ggplot(sum_df, aes(x = Count, color = Group, fill = Group)) +
      geom_density(alpha = 0.3) +
      ggtitle("Density plot") +
      labs(x = expression(paste(log[2], " (Count)")), y = "Density") +
      geom_vline(data = median_df, aes(xintercept = median, color = Group), linetype = "dotted", size = 1) +
      theme_classic()

    outp <- plot_grid(p, d, p1, ncol = 3, rel_widths = c(1, 1, 1))
    print(outp)
  }
  dev.off()


  if (!is.null(subsets)) {
    print("plot lfc for feautures")
    pdf(paste("lfc_at_TSNflank", TSSflank, "_and_transcript_for_", input_type, ".pdf", sep = ""), height = 8, width = 12)
    for (featureName in names(feature_count_list$count)) {
      lfc_list <- list()
      for (subject in names(subsets)) {
        subset <- subsets[[subject]]
        s2c <- data.frame(samples = subset, groups = c(rep(c("GFP", subject), times = c(2, 2))))
        contrasts <- data.frame(groups = c(rep("groups", 1)), treat = subject, control = c(rep("GFP", 1)))

        project <- paste(input_type, TSSflank, sep = "")

        data_table_feature <- feature_count_list[["count"]][[featureName]][, subset]

        if (is.null(count_table_feature) || nrow(data_table_feature) == 0) next
        data_table_feature <- data_table_feature[apply(data_table_feature, 1, sum) > 0, ]

        print(rpm[subset])
        data_table_feature <- round(t(t(data_table_feature) / rpm[subset]))
        print(head(data_table_feature))
        data_table_feature[is.na(data_table_feature)] <- 0

        res_feature <- NULL
        X <- try(
          res_feature <- run_DESeq2(data_table_feature, s2c, featureName, contrasts, plot = FALSE, useSizeFactor = FALSE)
        )
        if (!grepl("try-error", class(X))) {
          lfc_df <- data.frame(Group = rep(subject, length(res_feature[[1]]$log2FoldChange)), log2FoldChange = res_feature[[1]]$log2FoldChange)
          lfc_list[[subject]] <- lfc_df
        }
      }

      out_df <- as.data.frame(Reduce(rbind, lfc_list))
      if (nrow(out_df) == 0) next

      colnames(out_df) <- c("Group", "log2FoldChange")
      out_df <- mutate(out_df, Group = as.factor(Group)) %>%
        group_by(Group) %>%
        arrange(log2FoldChange) %>%
        mutate(cumCount = cumsum(log2FoldChange)) %>%
        mutate(Rank = order(cumCount)) %>%
        mutate(Fraction = cumCount / max(cumCount), Rank = Rank / max(Rank))

      p1 <- ggplot(data = sum_df, aes(x = Rank, y = Fraction, color = Group)) +
        geom_line() +
        ggtitle("Fingerprint plot") +
        labs(x = "Rank(log2FoldChange)", y = paste("Fraction over total log2FoldChange"))

      if (length(levels(out_df$Group)) < 2) next

      if (length(levels(out_df$Group)) > 1) {
        median_df <- aggregate(out_df$log2FoldChange, list(out_df$Group), median)
        colnames(median_df) <- c("Group", "median")

        count_df <- aggregate(out_df$log2FoldChange, list(out_df$Group), length)
        colnames(count_df) <- c("Group", "count")
        p <- ggplot(out_df, aes(x = Group, y = log2FoldChange, fill = Group)) +
          geom_boxplot(notch = FALSE) +
          theme_classic() +
          ggtitle(featureName) +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold", size = 14)) +
          theme(legend.position = "none") +
          labs(y = expression(paste(log[2], " (treat/control signal intensity)"))) +
          scale_x_discrete(limits = names(subsets)) + # labels=c("TSN" = expression(paste(m^6,"Am")), GroupName = expression(paste("5'-UTR ", m^6, "A")))) +
          theme(axis.text.x = element_text(face = "bold", size = 14, color = "black")) +
          geom_text(data = count_df, aes(x = Group, y = min(out_df$log2FoldChange) * 1.1, label = paste("n=", count, sep = ""))) +
          geom_signif(comparisons = combn(levels(out_df$Group), 2, simplify = FALSE), test = "t.test", map_signif_level = TRUE, step_increase = 0.1)

        # p1 <- add_pval(p, pairs = list(c(1, 2)), test='t.test')

        d <- ggplot(out_df, aes(x = log2FoldChange, color = Group, fill = Group)) +
          geom_density(alpha = 0.3) +
          ggtitle("Density plot") +
          ylab("Density") +
          geom_vline(data = median_df, aes(xintercept = median, color = Group), linetype = "dotted", size = 1) +
          theme_classic()

        outp <- plot_grid(p, d, p1, ncol = 3, rel_widths = c(1, 1, 1))
        print(outp)
      }
    }
    dev.off()
  }
}

#'
#'
#' @title Count reads in genomic features and/or peaks
#' @description Count the number of reads covering genomic features. If a peakfile is present, count the number of reads that cover all peaks (unrestricted)
#' and only peaks that are overlapping genomic features.
#'
#' @param queryfiles paths of bam files for aligned reads
#' @param peakfile path of a bed file for the peaks
#' @param txdb a TxDb object defined in the GenomicFeatures package
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param queryInputParams a list of parameters for \code{handle_input} of queryfiles
#' @param peakInputParams a list of parameters for \code{handle_input} of peakfile
#' @param TSSflank an integer denoting the number of nucleotides encompassing upstream and downstream of the TSS
#'
#' @return a list with two or three elements, 'count' is a list of dataframes, 'inputsize' is the total number of aligned reads
#' in the bam files, if peakfile is present, 'peaksize' is the number of peaks
#'
#' @author Shuye Pu
#'
#' @examples
#'
#' txdb <- AnnotationDbi::loadDb(system.file("data", paste0(genome, "_txdb.sql"), package = "GenomicPlotData"))
#' peakfile <- system.file("data", "test_B.bed", package = "GenomicPlotData")
#' names(peakfile) <- "TestB"
#' queryfiles <- system.file("data", "test_clip.bam", package = "GenomicPlotData")
#' names(queryfiles) <- "query"
#'
#' queryInputParams <- list(offset = -1, fix_width = 0, norm = TRUE, useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19")
#' peakInputParams <- list(offset = 0, fix_width = 21, fix_point = "center", norm = FALSE, useScore = FALSE, outRle = FALSE, useSizeFactor = FALSE, genome = "hg19")
#'
#' out <- count_overlap_in_features(queryfiles, peakfile, txdb, longest = TRUE, queryInputParams = queryInputParams, peakInputParams = peakInputParams, TSSflank = 2)
#' print(names(out$count))
#' print(names(out$count$cds))
#' print(sapply(out$count, colMeans))
#' @export count_overlap_in_features

count_overlap_in_features <- function(queryfiles, peakfile = NULL, txdb, longest = TRUE, queryInputParams = NULL, peakInputParams = NULL, TSSflank = 2) {
  querylabels <- names(queryfiles)
  names(querylabels) <- queryfiles
  input_list <- handle_input(inputFiles = queryfiles, importParams = queryInputParams)

  if (!is.null(peakfile)) {
    peaklabel <- names(peakfile)
    names(peaklabel) <- peakfile
    peak_list <- handle_input(inputFiles = peakfile, peakInputParams)
  }

  feature_count_list <- list()


  for (featureName in c("utr5", "cds", "utr3", "transcript", "TSA", "unrestricted")) {
    if (featureName == "TSA") {
      ## process transcription start sites with 'A'
      A_tss <- get_gr_at_TSS(txdb, longest = longest)
      feature_grl <- flank(A_tss, width = TSSflank, both = TRUE)
    } else if (featureName == "unrestricted") {
      if (!is.null(peakfile)) {
        feature_grl <- peak_list[[1]]$query
      } else {
        feature_grl <- NULL
      }
    } else {
      alist <- get_genomic_feature_coordinates(txdb, featureName, longest = longest)
      original_grl <- alist$GRangesList
      if (featureName %in% c("transcript")) {
        trimmed <- factR::trimTranscripts(original_grl, start = TSSflank, end = 0)
        feature_grl <- trimmed
      } else {
        feature_grl <- original_grl
      }
    }

    if (!is.null(peakfile)) {
      peakRegions <- peak_list[[names(peakfile)]]$query
      if (any(unique(runValue(strand(peakRegions))) %in% c("*", ".", ""))) {
        feature_grl <- filter_by_overlaps(peakRegions, feature_grl) ## feature_grl becomes peaks overlapping with feature
      } else {
        feature_grl <- filter_by_overlaps_stranded(peakRegions, feature_grl)
      }
    }

    if (!is.null(feature_grl)) {
      featureCount_df <- sapply(names(queryfiles), function(queryfile) {
        queryRegions <- input_list[[queryfile]]$query
        overlap_count <- countOverlaps(feature_grl, queryRegions, type = "any")
      })

      feature_count_list[[featureName]] <- as.data.frame(featureCount_df)
    }
  }

  bamsizes <- sapply(input_list, function(x) x$size)

  if (!is.null(peakfile)) {
    peaksizes <- sapply(peak_list, function(x) x$size)
    return(list("count" = feature_count_list, "inputsize" = bamsizes, "peaksize" = peaksizes))
  } else {
    return(list("count" = feature_count_list, "inputsize" = bamsizes))
  }
}


#' @title Identifies peak-targeted genes
#'
#' @description  Produce a table of genes targeted by peaks. A gene is targeted
#' if there is at least one peak overlaps with a user-specified feature, such as
#' Promoter, Exon, or Transcript
#' @param annotation a string denoting the peak annotation file name, which is a
#' output from \code{plot_peak_annotation}
#' @param feature a string denoting the feature in c("Promoter", "5'UTR", "CDS",
#' "3'UTR", "TTS", "Intron", "Exon", "Transcript"). Use "Exon" to find genes
#' that are targeted in mRNA, and use "Transcript" to find genes that are
#' targeted in both exon and intron.
#' @param verbose logical, whether to write the results into a file
#'
#' @return a dataframe with peak and gene information
#' @author Shuye Pu
#'
#'
#' @export peak_targeted_gene
#'
#'
peak_targeted_gene <- function(annotationFile, feature = "Exon",
                               verbose = FALSE, filter = FALSE) {
  stopifnot(file.exists(annotationFile))
  if(is.null(names(annotationFile))){
    names(annotationFile) <- gsub("_targeted_annotated_gene.tab", "",
                                  annotationFile, fixed = TRUE)
  }
  if(feature == "3'UTR") feature <- "X3.UTR"
  if(feature == "5'UTR") feature <- "X5.UTR"
  annot_table <- NULL

  annot_table <- read.delim2(annotationFile, header = TRUE) %>%
    dplyr::mutate(overlap_status = ifelse(.data[[feature]] > 0, "Targets",
                            "Nontargets"))

  if (verbose) {
    write.table(annot_table, annotationFile, sep = "\t", row.names = FALSE,
                quote = FALSE)
  }

  if (filter) {
    filtered_table <- annot_table %>%
      filter(overlap_status == "Targets") %>%
      mutate(score = 0) %>%
      select(chr, start, end, gene_name, score, strand)
    targetFile <- paste(names(annotationFile), feature, "target_gene.bed",
                        sep = "_")
    write.table(filtered_table, targetFile, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)

    filtered_table <- annot_table %>%
      filter(overlap_status != "Targets")%>%
      mutate(score = 0) %>%
      select(chr, start, end, gene_name, score, strand)
    targetFile <- paste(names(annotationFile), feature, "nontarget_gene.bed",
                        sep = "_")
    write.table(filtered_table, targetFile, sep = "\t", row.names = FALSE,
                col.names = FALSE, quote = FALSE)
  }

  invisible(annot_table)
}

#' @title Identifies peak_targeted enhancers
#'
#' @description  Produce a table of enhancers targeted by peaks, each enhancer is associated with a gene
#' @param peakfile, a string denoting the peak file name, only .bed format is allowed
#' @param enhancerfile, a string denoting the path to the enhancer file
#' @param enhancerGenefile, a string denoting the path to the enhancer-gene association file
#' @param maxg an integer denoting the distance that define overlap
#' @param importParams a list of parameters for \code{handle_input}
#'
#' @return a dataframe with peak and enhancer/gene information
#' @author Shuye Pu
#'
#'
#' @export peak_targeted_enhancer
#'
peak_targeted_enhancer <- function(peakfile, enhancerfile, enhancerGenefile, maxg = -1L, importParams = NULL) {
  EG <- read.delim(enhancerGenefile, header = FALSE)
  EG_list <- lapply(EG[, 1], function(x) {
    str1 <- unlist(strsplit(x, split = ":|_|\\$"))
    invisible(c(str1[1], unlist(strsplit(str1[2], split = "-")), str1[3:7]))
  })
  names(EG_list) <- paste0("EG", seq_along(EG_list))
  EG_df <- as.data.frame(do.call(rbind, EG_list))
  colnames(EG_df) <- c("chrEnhancer", "startEnhancer", "endEnhancer", "geneId", "geneName", "chrGene", "coordGene", "strandGene")

  return_table <- NULL
  peaklabel <- names(peakfile)

  importParams$useScore <- FALSE
  importParams$outRle <- FALSE
  importParams$fix_width <- 0

  peak <- handle_bed(inputFile = peakfile, importParams)$query
  enhancer <- handle_bed(inputFile = enhancerfile, importParams)$query

  print(length(peak))
  print(length(enhancer))
  Overlaps <- findOverlaps(peak, enhancer)
  peak_gr <- peak[Overlaps@from]
  peak_df <- gr2df(peak_gr) %>%
    mutate(
      chrPeak = as.character(chr),
      startPeak = as.integer(start) - 1,
      endPeak = as.integer(end),
      widthPeak = as.integer(score),
      strandPeak = "*",
      .keep = "unused"
    )
  enhancer_gr <- enhancer[Overlaps@to]
  names(enhancer_gr) <- seq_along(enhancer_gr)
  enhancer_df <- gr2df(enhancer_gr) %>%
    mutate(
      chrEnhancer = as.character(chr),
      startEnhancer = as.integer(start) - 1,
      endEnhancer = as.integer(end),
      widthEnhancer = as.integer(score),
      strandEnhancer = "*",
      .keep = "unused"
    )

  ot <- cbind(peak_df, enhancer_df) %>%
    select(chrPeak, startPeak, endPeak, widthPeak, strandPeak, chrEnhancer, startEnhancer, endEnhancer, widthEnhancer, strandEnhancer)

  return_table <- merge(ot, EG_df)
  return_table <- return_table[, colnames(return_table)[c(4:8, 1:3, 9:15)]]
  uni_peaks <- unique(return_table[, 1:5])
  uni_enhancers <- unique(return_table[, 6:10])
  uni_genes <- unique(return_table[11:15])

  fraction_peaks <- nrow(uni_peaks) / length(peak)
  fraction_enhancers <- nrow(uni_enhancers) / length(enhancer)

  sink(paste(peaklabel, "targeted_enhancer_gene_stats.tab", sep = "_"))
  cat(paste0("GenomicPlot.R::peak_targeted_enhancer()\t", date()))
  cat(paste0("\nSummary stats for\t", peaklabel, "_targeted_enhancer_gene.tab"))
  cat(paste0("\nTotal number of summit peaks\t", length(peak)))
  cat(paste0("\nTotal number of enhancers in HEK293\t", length(enhancer)))
  cat(paste0("\nNumber of peaks targeting enhancer\t", nrow(uni_peaks)))
  cat(paste0("\nFraction of peaks targeting enhancer\t", fraction_peaks))
  cat(paste0("\nNumber of enhancers targeted by peaks\t", nrow(uni_enhancers)))
  cat(paste0("\nFraction of enhancer targeted by peaks\t", fraction_enhancers))
  cat(paste0("\nNumber of genes associated with targeted enhancer\t", nrow(uni_genes)))
  sink()

  write.table(return_table, paste(peaklabel, "targeted_enhancer_gene.tab", sep = "_"), sep = "\t", row.names = FALSE, quote = FALSE)

  invisible(return_table)
}







# fiveP and threeP can be both positive and negative integers,
# positive means shift downstream, negative means shift upstream,
# so, a positive fiveP and a negative threeP will shrink the range.
expand_grange <- function(gr, fiveP, threeP, ignore.strand = FALSE) {
  if (ignor.strand) {
    start(gr) <- start() + fiveP
    end(gr) <- end(gr) + threeP
  } else {
    grPlus <- gr[as.vector(strand(gr)) %in% c("+", "*")]
    grMinus <- gr[as.vector(strand(gr)) == "-"]

    start(grPlus) <- start(grPlus) + fiveP
    end(grPlus) <- end(grPlus) + threeP
    start(grMinus) <- start(grMinus) - threeP
    end(grMinus) <- end(grMinus) - fiveP
    gr <- c(grPlus, grMinus)
  }
  return(gr)
}

## find genes whose TTS do not overlap with promoters of other genes
isolated_gene <- function(allGene, fiveP, threeP) {
  promoter <- GenomicRanges::promoters(allGene, upstream = fiveP, downstream = 0, use.names = TRUE)
  TTS <- GenomicRanges::flank(allGene, width = threeP, both = FALSE, start = FALSE, ignore.strand = FALSE)

  oc <- countOverlaps(TTS, promoter, type = "any", ignore.strand = TRUE)
  noc <- names(oc[oc == 0])

  iso <- allGene[noc]
  return(iso)
}


## partition genes along the length into three parts
partition_gene <- function(gr, fiveP, fiveIn, threeP, threeIn, export = FALSE) {
  # gr <- get_genomic_feature_coordinates(txdb, "gene", longest=TRUE, protein_coding = TRUE)[["GRanges"]]
  gr <- gr[width(gr) > (fiveIn + threeIn)]

  promoter <- GenomicRanges::promoters(gr, upstream = fiveP, downstream = fiveIn, use.names = TRUE)
  TTS <- GenomicRanges::resize(gr, width = 1, fix = "end", ignore.strand = FALSE)
  TTS <- GenomicRanges::promoters(TTS, upstream = threeIn, downstream = threeP, use.names = TRUE)

  grPlus <- gr[as.vector(strand(gr)) %in% c("+", "*")]
  grMinus <- gr[as.vector(strand(gr)) == "-"]
  grPlus <- narrow(grPlus, start = fiveIn, end = -threeIn)
  grMinus <- narrow(grMinus, start = threeIn, end = -fiveIn)

  centerGr <- c(grPlus, grMinus)

  if (export) {
    promoter_df <- gr2df(promoter)
    TTS_df <- gr2df(TTS)
    center_df <- gr2df(centerGr)
    write.table(TTS_df, paste0("TTS_", threeIn, "_", threeP, ".bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    write.table(promoter_df, paste0("Promoter_", fiveIn, "_", fiveP, ".bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    write.table(center_df, paste0("Body_", fiveIn, "_", threeIn, ".bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }

  return(list("TTS" = TTS, "Promoter" = promoter, "Body" = centerGr))
}

prepare_nonoverlapping_genes <- function(txdb, fiveP, fiveIn, threeP, threeIn, export = FALSE) {
  allGene <- get_genomic_feature_coordinates(txdb, "gene", longest = TRUE, protein_coding = TRUE)[["GRanges"]]
  gr <- isolated_gene(allGene, fiveP, threeP)
  res <- partition_gene(gr, fiveP, fiveIn, threeP, threeIn, export = export)
  invisible(res)
}

stalling_stopping_ratio <- function(experiment, op = "mean", dataType = "Signal intensity") {
  segment <- c("Promoter", "gene", "TTS")

  signals <- lapply(segment, function(x) {
    signal <- read.table(paste(dataType, experiment, x, "matrix.tab", sep = "_"), header = TRUE, sep = "\t")
    if (op == "sum") {
      signal_op <- apply(signal, 1, sum)
    } else if (op == "mean") {
      signal_op <- apply(signal, 1, mean)
    } else {
      stop("op is not supported!")
    }
  })
  names(signals) <- segment

  lapply(signals, length)
  common_gene <- Reduce(intersect, lapply(signals, names))
  total_signals <- signals$Promoter[common_gene] + signals$gene[common_gene] + signals$TTS[common_gene]
  stalling_ratio <- signals$Promoter[common_gene] / total_signals
  stopping_ratio <- signals$TTS[common_gene] / total_signals

  return(list("stalling" = stalling_ratio, "stopping" = stopping_ratio))
}


MACS_to_bed <- function(macs.out) {
  basef <- gsub("\\.macs\\.out\\.txt\\.tsv", "", basename(macs.out))
  macs <- read.delim(macs.out, header = TRUE, comment.char = "#")
  macs <- macs %>%
    mutate(name = paste0("chipPeak", seq(nrow(macs)))) %>%
    mutate(summitStart = as.integer(start + summit)) %>%
    mutate(summitEnd = as.integer(start + summit + 1)) %>%
    mutate(strand = "*")

  narrowPeak <- macs %>%
    select(chr, start, end, name, tags, strand)
  summitPeak <- macs %>%
    select(chr, summitStart, summitEnd, name, tags, strand)

  write.table(narrowPeak, paste0(basef, ".narrowPeak"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  write.table(summitPeak, paste0(basef, ".bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

MACS2_to_bed <- function(macs2.out, qCutoff = 3) {
   standard_chr <- paste0("chr", c(seq(1:22), "M", "X", "Y"))
   basef <- gsub("\\.xls|\\.macs\\.out\\.txt\\.tsv", "", basename(macs2.out))
   macs <- read.delim2(macs2.out, header = TRUE, comment.char = "#")
   macs <- macs %>%
      mutate(summitStart = as.integer(abs_summit - 1)) %>%
      mutate(summitEnd = as.integer(abs_summit)) %>%
      mutate(strand = "*") %>%
      filter(X.log10.qvalue. > qCutoff) %>%
      filter(chr %in% standard_chr)

   narrowPeak <- macs %>%
      select(chr, start, end, name, X.log10.qvalue., strand)
   summitPeak <- macs %>%
      select(chr, summitStart, summitEnd, name, X.log10.qvalue., strand)

   write.table(narrowPeak, paste0(basef, "_filtered.narrowPeak"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
   write.table(summitPeak, paste0(basef, "_summits_filtered.bed"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

add_chr_to_bed <- function(inbed) {
   if (grepl("\\.bed$", inbed)) {
      basef <- gsub("\\.bed$", "", basename(inbed))
      outbed <- paste0(basef, "_fixed.bed")
   } else {
      basef <- gsub("\\.narrowPeak$", "", basename(inbed))
      outbed <- paste0(basef, "_fixed.narrowPeak")
   }

   df <- read.delim2(inbed, header = FALSE, comment.char = "#")

   if(!grepl("chr", df[1,1])){
      df <- df %>%
         mutate(V1 = paste0("chr", V1))
   }

   write.table(df, outbed, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
   invisible(outbed)
}


#' @title Format genomic coordinates in GRanges or GRrangesList as strings used in igv
#' @description This function takes a GRanges or GRangesList object, and transform each range into a string
#' @param x a GRanges or GRangesList object
#' @return a vector of strings in the format of 'chr:start-end(strand)'
#' @author Shuye Pu
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr2", IRanges::IRanges(3, 6))
#' gr2 <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(7, 13), width = 3),
#'   strand = c("+", "-")
#' )
#' gr3 <- GenomicRanges::GRanges(c("chr1", "chr2"), IRanges::IRanges(c(1, 4), c(3, 9)),
#'   strand = "-"
#' )
#'
#' grl <- GenomicRanges::GRangesList(gr1 = gr1, gr2 = gr2, gr3 = gr3)
#' grl
#'
#' out <- format_genomic_coordinates(grl)
#' cat(out)
#'
#' @export format_genomic_coordinates
#'
format_genomic_coordinates <- function(x) {
  if (grepl("GRangesList", class(x))) x <- unlist(x) ## convert grl to gr

  chr <- as.vector(seqnames(x)) %>% as.character()
  start <- start(x) %>% as.numeric()
  end <- end(x) %>% as.numeric()
  strand <- strand(x) %>% as.character()

  out <- paste0(chr, ":", start, "-", end, "(", strand, ")")
  invisible(out)
}


#' @title Parallel execution of get_genomic_feature_coordinates
#' @description Get genomic coordinates for multiple features at once
#'
#' @param txdb a TxDb object
#' @param featureNames  a vector of gene features in c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")
#' @param longest logical, indicating whether the output should be limited to the longest transcript of each gene
#' @param protein_coding logical, indicating whether to limit to protein_coding genes
#' @param nc number of cores to use
#'
#' @return a list of lists(Granges, GRangesList)
#' @author Shuye Pu
#'
#' @examples
#' txdb <- AnnotationDbi::loadDb(system.file("extdata", "txdb_chr19.sql",
#'   package = "GenomicPlot"
#' ))
#'
#' coor <- GenomicPlot:::parallel_feature_coordinates(txdb, featureNames = c("utr5", "utr3"))
#'
#' @keywords internal
#'
#' @note not working because txdb cannot be exported to nodes, need to find a solution
#'

parallel_feature_coordinates <- function(txdb,
                                         featureNames,
                                         longest = TRUE,
                                         protein_coding = TRUE,
                                         nc = 2) {
  txdb <- force(txdb)

  cl <- start_parallel(min(length(featureNames), nc))

  clusterExport(cl, varlist = c("get_genomic_feature_coordinates"), envir = environment())
  clusterExport(cl, varlist = c("txdb", "longest", "protein_coding"), envir = environment())
  alist <- clusterApply(cl, featureNames, get_genomic_feature_coordinates, txdb = txdb, longest = longest, protein_coding = protein_coding)
  stop_parallel(cl)

  invisible(alist)
}

#' @title Parallel execution of binnedAverage
#'
#' @description Function for parallel computation of binnedAverage function in the GenomicRanges package
#'
#' @param Rle_list a list of RleList objects.
#' @param tileBins a GRanges object of tiled genome
#' @param nc integer, number of cores for parallel processing
#'
#' @return a list of numeric vectors
#' @author Shuye Pu
#'
#'
#' @export parallel_binnedAverage
#'
#'
parallel_binnedAverage <- function(Rle_list,
                                   tileBins,
                                   nc = 2) {
  # print(system.time({
  cl <- start_parallel(min(length(Rle_list), nc))

  clusterExport(cl, varlist = c("binnedAverage", "tileBins"), envir = environment())
  score_list <- parLapply(cl, Rle_list, function(x) {
    binAverage <- binnedAverage(tileBins, x, varname = "binned_score", na.rm = FALSE)
    binAverage$binned_score
  })

  stop_parallel(cl)
  # }))

  invisible(score_list)
}



#' @title Plot promoter, gene body and TTS
#
#' @description Plot reads or peak Coverage/base/gene of samples given in the query files around genes. The upstream and downstream windows flanking genes
#' can be given separately, the parameter 'meta' controls if gene or metagene plots are generated. If Input files are provided, ratio over Input
#' is computed and displayed as well.
#'
#' @param queryFiles a vector of sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param gFeatures genomic features as output of the function 'prepare_3parts_genomic_features'
#' @param inputFiles a vector of input sample file names. The file should be in .bam, .bed, .wig or .bw format, mixture of formats is allowed
#' @param importParams a list of parameters for \code{handle_input}
#' @param stranded logical, indicating whether the strand of the feature should be considered
#' @param scale logical, indicating whether the score matrix should be scaled to the range 0:1, so that samples with different baseline can be compared
#' @param smooth logical, indicating whether the line should smoothed with a spline smoothing algorithm
#' @param heatmap logical, indicating whether a heatmap of the score matrix should be generated
#' @param heatRange a numerical vector of two elements, defining range for heatmap color ramp generation
#' @param rmOutlier a numeric value serving as a multiplier of the MAD in Hampel filter for outliers identification, 0 indicating not removing outliers. For Gaussian distribution, use 3, adjust based on data distribution
#' @param outPrefix a string specifying output file prefix for plots (outPrefix.pdf)
#' @param transform a string in c("log", "log2", "log10"), default = NA indicating no transformation of data matrix
#' @param Ylab a string for y-axis label
#' @param verbose logical, whether to output additional information (data used for plotting or statistical test results)
#' @param hw a vector of two elements specifying the height and width of the output figures
#' @param nc integer, number of cores for parallel processing
#'
#' @return a dataframe containing the data used for plotting
#' @author Shuye Pu
#'
#' @export plot_3parts_metagene


plot_3parts_metagene <- function(queryFiles,
                                 gFeatures,
                                 inputFiles = NULL,
                                 scale = FALSE,
                                 verbose = FALSE,
                                 Ylab = "Coverage/base/gene",
                                 importParams = NULL,
                                 smooth = FALSE,
                                 stranded = TRUE,
                                 outPrefix = NULL,
                                 heatmap = FALSE,
                                 rmOutlier = 0,
                                 heatRange = NULL,
                                 transform = NA,
                                 hw = c(8, 8),
                                 nc = 2) {
   functionName <- as.character(match.call()[[1]])
   params <- plot_named_list(as.list(environment()))
   print(params)

   if (is.null(inputFiles)) {
      inputLabels <- NULL
      queryInputs <- handle_input(inputFiles = queryFiles, importParams, verbose = verbose, nc = nc)
   } else {
      inputLabels <- names(inputFiles)
      queryLabels <- names(queryFiles)
      if (length(queryFiles) == length(inputFiles)) {
         queryInputs <- handle_input(inputFiles = c(queryFiles, inputFiles), importParams, verbose = verbose, nc = nc)
      } else if (length(inputFiles) == 1) {
         queryInputs <- handle_input(inputFiles = c(queryFiles, inputFiles), importParams, verbose = verbose, nc = nc)
         queryInputs <- queryInputs[c(queryLabels, rep(inputLabels, length(queryLabels)))] # expand the list

         inputLabels <- paste0(names(inputFiles), seq_along(queryFiles)) # make each inputLabels unique
         names(queryInputs) <- c(queryLabels, inputLabels)
      } else {
         stop("the number of inputFiles must be 1 or equal to the number of queryFiles!")
      }
   }
   queryLabels <- names(queryInputs)

   if (!is.null(outPrefix)) {
      pdf(paste(outPrefix, "pdf", sep = "."), height = hw[1], width = hw[2])
   }

   windowRs <- gFeatures$windowRs
   featureNames <- names(windowRs)
   if (verbose) {
      message("Feature name: ", paste(featureNames, collapse = " "), "\n")
      message("Window sizes: ", paste(vapply(windowRs, length, numeric(1)), collapse = " "), "\n")
   }

   nbins <- gFeatures$nbins
   scaled_bins <- gFeatures$scaled_bins
   meta <- gFeatures$meta
   fiveP <- gFeatures$fiveP
   threeP <- gFeatures$threeP

   ## start overlapping

   if (verbose) message("Computing coverage for queryFiles...\n")
   scoreMatrix_list <- list()

   for (queryLabel in queryLabels) {
      myInput <- queryInputs[[queryLabel]]
      libsize <- myInput$size
      queryRegions <- myInput$query
      fileType <- myInput$type
      weight_col <- myInput$weight

      for (w in featureNames) {
         if (verbose) message("Feature name: ", w, "\n")
         windowR <- as(windowRs[[w]], "GRangesList")
         bin_num <- scaled_bins[w]

         bin_op <- "mean"
         fullMatrix <- parallel_scoreMatrixBin(queryRegions, windowR, bin_num, bin_op, weight_col, stranded, nc = nc)

         scoreMatrix_list[[queryLabel]][[w]] <- fullMatrix
      }
   }

   if (verbose) message("Plotting coverage for queryFiles...\n")
   mplot_df <- NULL
   vx <- c(1, cumsum(scaled_bins[seq_len(length(scaled_bins) - 1)]) + 1) ## x axis points for vlines that demarcate the genomic features
   names(vx) <- featureNames

   Ylab <- ifelse(!is.na(transform) && is.null(inputFiles), paste0(transform, " (", Ylab, ")"), Ylab)

   if (heatmap) heatmap_list <- list()

   processed_matrix <- list()
   for (queryLabel in queryLabels) {
      plot_df <- NULL

      dims <- vapply(scoreMatrix_list[[queryLabel]], dim, numeric(2))

      if (any(dims[1, ] != dims[1, 1])) {
         message(paste(dims[1, ], collapse = " "), "\n")
         stop("Number of genes are not equal among features, make sure all feature windows are within chromosome lengths of query regions,
            as genomation will remvove all feature windows outside chromosome boundaries")
      } else {
         featureMatrix <- as.matrix(bind_cols(scoreMatrix_list[[queryLabel]]))
         if (is.null(inputFiles)) {
            fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier, transform = transform, verbose = verbose)
         } else {
            fullMatrix <- process_scoreMatrix(fullMatrix, scale = FALSE, rmOutlier = rmOutlier, transform = NA, verbose = verbose)
         }
         processed_matrix[[queryLabel]] <- featureMatrix

         colm <- apply(featureMatrix, 2, mean)
         colsd <- apply(featureMatrix, 2, sd)
         colse <- colsd / sqrt(nrow(featureMatrix))
         querybed <- rep(queryLabel, ncol(featureMatrix))
         collabel <- list()
         featuretype <- list()
         for (w in featureNames) {
            if (verbose) message("Feature name: ", w, "\n")
            if (scaled_bins[w] > 0) {
               bin_num <- scaled_bins[w]
               collabel[[w]] <- seq(vx[w], vx[w] + bin_num - 1)
               featuretype[[w]] <- rep(w, bin_num)
            }
         }
         collabel <- unlist(collabel)
         featuretype <- unlist(featuretype)
         names(collabel) <- featuretype

         if (heatmap) {
            dataname <- paste(Ylab, queryLabel, "gene", sep = ":")
            heatmap_list[dataname] <- draw_matrix_heatmap(featureMatrix, dataName = dataname, labels_col = collabel, levels_col = featureNames, ranges = heatRange, verbose = verbose)
         }

         plot_df <- data.frame("Intensity" = colm, "sd" = colsd, "se" = colse, "Position" = collabel, "Query" = querybed, "Feature" = featuretype)
      }

      if (smooth) {
         plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df = as.integer(nbins / 5))$y)
         plot_df$se <- as.vector(smooth.spline(plot_df$se, df = as.integer(nbins / 5))$y)
      }

      mplot_df <- rbind(mplot_df, plot_df)
   }

   mplot_df <- mutate(mplot_df, lower = Intensity - se, upper = Intensity + se)

   xmax <- max(mplot_df$Position)
   pp <- draw_region_landmark(featureNames, vx, xmax)
   ppp <- draw_region_name(featureNames, scaled_bins, xmax)
   marker <- plot_grid(pp, ppp, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 2))
   ## plot individual sample lines with error band
   plot_list <- list()
   for (queryLabel in queryLabels) {
      aplot_df <- mplot_df %>%
         filter(Query == queryLabel)
      p <- draw_region_profile(plot_df = aplot_df, cn = "Query", vx = vx, Ylab = Ylab)
      outp <- plot_grid(p, marker, ncol = 1, align = "v", axis = "lr", rel_heights = c(10, 1))
      plot_list[[queryLabel]] <- outp
   }
   rowp <- plot_grid(plotlist = plot_list, nrow = 1, align = "h")
   # print(rowp)

   if (heatmap) {
      groblist <- lapply(heatmap_list, function(x) grid.grabExpr(draw(x, heatmap_legend_side = "left")))
      heatp <- plot_grid(plotlist = groblist, nrow = 1, align = "h")

      composite <- plot_grid(rowp, heatp, ncol = 1, align = "v")
      print(composite)
   } else {
      print(rowp)
   }

   ## plot multi-sample lines with error band
   p <- draw_region_profile(plot_df = mplot_df, cn = "Query", vx = vx, Ylab = Ylab)
   outp <- plot_grid(p, marker, ncol = 1, align = "v", axis = "lr", rel_heights = c(10, 1))
   print(outp)

   if (!is.null(inputFiles)) {
      if (verbose) message("Computing coverage for ratio over input...\n")
      Ylabr <- ifelse(is.na(transform), "Ratio-over-Input", paste0(transform, " (Ratio-over-Input)"))
      if (heatmap) heatmap_list_ratio <- list()

      inputMatrix_list <- processed_matrix[inputLabels]
      ratiolabels <- queryLabels[!queryLabels %in% inputLabels]
      ratioMatrix_list <- processed_matrix[ratiolabels]

      for (i in seq_along(ratiolabels)) {
         fullMatrix <- ratio_over_input(ratioMatrix_list[[ratiolabels[i]]], inputMatrix_list[[inputLabels[i]]], verbose)
         fullMatrix <- process_scoreMatrix(fullMatrix, scale, rmOutlier, transform, verbose)

         ratioMatrix_list[[ratiolabels[i]]] <- fullMatrix
      }


      if (verbose) message("Plotting coverage for ratio over input...\n")
      mplot_df_ratio <- NULL
      for (ratiolabel in ratiolabels) {
         plot_df <- NULL

         featureMatrix <- ratioMatrix_list[[ratiolabel]]

         colm <- apply(featureMatrix, 2, mean)
         colsd <- apply(featureMatrix, 2, sd)
         colse <- colsd / sqrt(nrow(featureMatrix))
         querybed <- rep(ratiolabel, ncol(featureMatrix))
         collabel <- list()
         featuretype <- list()
         for (w in featureNames) {
            if (verbose) message("Feature name: ", w, "\n")
            if (scaled_bins[w] > 0) {
               bin_num <- scaled_bins[w]
               collabel[[w]] <- seq(vx[w], vx[w] + bin_num - 1)
               featuretype[[w]] <- rep(w, bin_num)
            }
         }
         collabel <- unlist(collabel)
         featuretype <- unlist(featuretype)
         names(collabel) <- featuretype

         if (heatmap) {
            dataname <- paste(Ylabr, ratiolabel, "gene", sep = ":")
            heatmap_list_ratio[dataname] <- draw_matrix_heatmap(featureMatrix, dataName = dataname, labels_col = collabel, levels_col = featureNames, ranges = heatRange, verbose = verbose)
         }
         plot_df <- data.frame("Intensity" = colm, "sd" = colsd, "se" = colse, "Position" = collabel, "Query" = querybed, "Feature" = featuretype)

         if (smooth) {
            plot_df$Intensity <- as.vector(smooth.spline(plot_df$Intensity, df = as.integer(nbins / 5))$y)
            plot_df$se <- as.vector(smooth.spline(plot_df$se, df = as.integer(nbins / 5))$y)
         }

         mplot_df_ratio <- rbind(mplot_df_ratio, plot_df)
      }

      mplot_df_ratio <- mutate(mplot_df_ratio, lower = Intensity - se, upper = Intensity + se)

      ## plot individual sample lines with error band
      plot_list <- list()
      for (ratiolabel in ratiolabels) {
         aplot_df <- mplot_df_ratio %>%
            filter(Query == ratiolabel)
         p <- draw_region_profile(plot_df = aplot_df, cn = "Query", vx = vx, Ylab = Ylabr)
         outp <- plot_grid(p, marker, ncol = 1, align = "v", axis = "lr", rel_heights = c(10, 1))
         plot_list[[ratiolabel]] <- outp
      }
      rowp <- plot_grid(plotlist = plot_list, nrow = 1, align = "h")
      # print(rowp)

      if (heatmap) {
         groblist <- lapply(heatmap_list_ratio, function(x) grid.grabExpr(draw(x, heatmap_legend_side = "left")))
         heatp <- plot_grid(plotlist = groblist, nrow = 1, align = "h")
         composite <- plot_grid(rowp, heatp, ncol = 1, align = "v")
         print(composite)
      } else {
         print(rowp)
      }

      ## plot multi-sample lines with error band
      p <- draw_region_profile(plot_df = mplot_df_ratio, cn = "Query", vx = vx, Ylab = Ylabr)
      outp <- plot_grid(p, marker, ncol = 1, align = "v", axis = "lr", rel_heights = c(10, 1))
      print(outp)
   }

   if (!is.null(outPrefix)) {
      print(params)
      on.exit(dev.off(), add = TRUE)
   }

   invisible(mplot_df)
}

gene_type_distribution <- function(gff, gene_list){

  gene_info_table <- gr2df(gff) %>%
    filter(type == "gene")
  ensembl <- vapply(strsplit(gene_info_table$gene_id,
                               split=".", fixed = TRUE), function(x) x[1],
                      FUN.VALUE = character(1))
  gene_info_table <- gene_info_table %>%
    mutate(gene_ID = ensembl)

  selected <- gene_info_table %>%
    filter(gene_ID %in% gene_list) %>%
    select(gene_ID, gene_name, gene_biotype)

  df <- selected %>%
    group_by(gene_biotype) %>%
    summarize(count = dplyr::n_distinct(gene_ID)) %>%
    rename_with(~ c("Feature", "Count")) %>%
    arrange(desc(Count)) %>%
    top_n(n = 3, wt = Count)

  p <- ggplot(df, aes(x = reorder(Feature, -Count), y = Count, fill = Feature)) +
    geom_bar(stat = "identity") +
    theme(legend.position = "none") +
    xlab(label = "Gene type")
  print(p)

  return(selected)
}

#' @examples
#' bedA <- system.file("extdata", "test_chip_peak_chr19.bed",
#'                      package = "GenomicPlot"
#'                      )
#'
#' bedB <- system.file("extdata", "test_chip_peak_chr19.narrowPeak",
#'                      package = "GenomicPlot"
#'                      )
#'
#' res <- bed_overlap(bedA, bedB, maxgap = -1, stranded = FALSE)
#'
bed_overlap <- function(bedA, bedB, maxgap = -1, minoverlap = 0, ext = 20,
                         genome = "hg19", stranded = TRUE, export = FALSE){
  if(is.null(names(bedA))) names(bedA) <- basename(bedA)
  if(is.null(names(bedB))) names(bedB) <- basename(bedB)
  importParams <- GenomicPlot::setImportParams(outRle = FALSE, fix_point = "center",
                                               fix_width = ext, genome = genome)
  grA <- GenomicPlot::handle_bed(bedA, importParams = importParams)$query
  grB <- GenomicPlot::handle_bed(bedB, importParams = importParams)$query

  if(stranded){
    AoverlapB <- GenomicPlot::filter_by_overlaps_stranded(grA, grB, maxgap = -1,
                                                          minoverlap = 0,
                                                          ignore.order = FALSE)
    AnotoverlapB <- GenomicPlot::filter_by_nonoverlaps_stranded(grA, grB, maxgap = -1,
                                                             minoverlap = 0,
                                                             ignore.order = FALSE)
  }else{
    AoverlapB <- plyranges::filter_by_overlaps(grA, grB, maxgap = -1,
                                                 minoverlap = 0)
    AnotoverlapB <- plyranges::filter_by_non_overlaps(grA, grB, maxgap = -1,
                                                 minoverlap = 0)
  }

  AoverlapB_df <- GenomicPlot::gr2df(AoverlapB)
  AnotoverlapB_df <- GenomicPlot::gr2df(AnotoverlapB)
  if(export){
    write.table(AoverlapB_df,
                paste0(names(bedA), "_overlap_", names(bedB), ".bed"),
                row.names = FALSE, col.names = FALSE, sep = "\t",
                quote = FALSE)
    write.table(AnotoverlapB_df,
                paste0(names(bedA), "_notOverlap_", names(bedB), ".bed"),
                row.names = FALSE, col.names = FALSE, sep = "\t",
                quote = FALSE)
  }

  invisible(list(overlap = AoverlapB, notoverlap = AnotoverlapB))
}


compare_methylation_levels <- function(bedfile, GLORI_file, methylation_df,
                                       genome, slop=0, title = "title"){
  genome_df <- as.data.frame(GenomeInfoDb::Seqinfo(genome=genome)) %>%
    tibble::rownames_to_column(var = "chrom") %>%
    mutate(size = seqlengths, .keep = "unused")
  valr::check_genome(genome_df)

  CITS <- valr::read_bed(bedfile)
  if(slop > 0) CITS <- bed_slop(CITS, genome = genome_df, both = slop)
  bed <- names(bedfile)
  GLORI_m6A <- valr::read_bed(GLORI_file)
  GLORI_in_CITS <- bed_intersect(CITS, GLORI_m6A)

  GLORI_in_CITS_df <- data.frame(score = as.numeric(GLORI_in_CITS$score.y),
                                 Cluster_info = bed)

  combined_df <- rbind(methylation_df, GLORI_in_CITS_df) %>%
    mutate(Methylation_level = score) %>%
    mutate(Cluster_info = factor(Cluster_info,
                                 levels=c("Clustered", "Non-cluster", bed)))

  p1 <- draw_combo_plot(combined_df, xc = "Cluster_info",
                        yc = "Methylation_level",
                        comp = combn(c(1,2,3), m=2, simplify = FALSE))

  p2 <- ggplot(combined_df, aes(x=Methylation_level, color=Cluster_info)) +
    geom_density(linewidth=2) +
    theme_classic()
  print(cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(2, 1),
                           labels = c(title, " ")))
}
