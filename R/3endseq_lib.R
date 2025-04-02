#' Find files by regex pattern
#' 
#' This function searches for files in a specified directory that match a regex pattern
#' and filters them based on a list of sample names.
#' 
#' @param path Character string specifying the directory path to search in
#' @param pattern Character string specifying the regex pattern for file filtering
#' @param samples Vector of character strings with sample names to match
#' @param start Logical, if TRUE prepends "^" to each sample name for exact start matching (default: FALSE)
#' @return Named vector of file paths that match the criteria
match_samples <- function(path, pattern, samples, start=FALSE){
    
    queryfiles <- list.files(path=path, pattern=pattern)
    queryfiles.idx <- unlist(sapply(samples, function(x){
        if(start) x <- paste0("^", x)
        grep(x, queryfiles, fixed = FALSE)
    }))
    queryfiles <- queryfiles[queryfiles.idx]
    queryfiles <- file.path(path, queryfiles)
    names(queryfiles) <- names(queryfiles.idx)
    
    return(queryfiles)
}


#' Define clusters of cleavage sites (CS)
#' 
#' This function clusters cleavage sites based on their genomic proximity.
#' It processes a GRanges object containing cleavage sites and groups them
#' into clusters based on specified distance parameters.
#' 
#' @param CSgr GRanges object containing cleavage sites
#' @param separatingDist Numeric, minimum distance in bp to separate clusters (default: 75)
#' @param addonDist Numeric, distance threshold in bp to merge sites into existing clusters (default: 5)
#' @return GRanges object with clustered cleavage sites
CS2peak <- function(CSgr, separatingDist = 75, addonDist = 5) {
    CSgr <- sortSeqlevels(CSgr)
    CSgrl <- split(CSgr, paste(seqnames(CSgr), strand(CSgr)))
    CSgrl <- CSgrl[sapply(CSgrl, function(x) length(x) > 0)]
    
    cl <- start_parallel(nc = length(CSgrl), verbose = TRUE)
    clusterExport(cl, varlist = c("CS2peakChr", "separatingDist", "addonDist"), envir = environment())
    overall <- parLapply(cl, CSgrl, CS2peakChr, separatingDist, addonDist)
    # out <- lapply(CSgrl, CS2peakChr, separatingDist, addonDist)
    stop_parallel(cl)
    
    overall <- overall[!sapply(overall, is.null)]
    overall_gr <- unlist(as(overall, "GRangesList"), use.names = FALSE)
    overall_gr <- sort(overall_gr, by = ~name)
    
    return(overall_gr)
}


#' Cluster cleavage sites on a single chromosome
#' 
#' This function processes cleavage sites on a single chromosome and clusters them
#' based on proximity. It is based on the method from Rot et al., 2017, Cell Reports 19, 1056–1067.
#' 
#' @param CHRgr GRanges object containing cleavage sites from a single chromosome
#' @param separatingDist Numeric, minimum distance in bp to separate clusters (default: 75)
#' @param addonDist Numeric, distance threshold in bp to merge sites into existing clusters (default: 5)
#' @return GRanges object with clustered cleavage sites for the chromosome, or NULL if input is empty
CS2peakChr <- function(CHRgr, separatingDist = 75, addonDist = 5) {
    if (length(CHRgr) == 0) {
        return(NULL)
    }
    library(GenomicRanges)
    CHRgr <- sort(CHRgr, by = ~name)
    mcols(CHRgr)$status <- "eliminated"
    kept <- CHRgr[1]
    CHRgr[1]$status <- "kept"
    
    for (i in 2:length(CHRgr)) {
        current <- CHRgr[i]
        nearest_namex <- nearest(current, kept)
        
        if (is.na(nearest_namex)) {
            kept <- c(kept, current)
            print(paste(i, "kept"))
            CHRgr[i]$status <- "kept"
        } else {
            name <- kept[nearest_namex]$name
            min_dist <- abs(GenomicRanges::distance(current, kept[nearest_namex]))
            if (min_dist < separatingDist) {
                if (min_dist < addonDist) {
                    print(paste(i, "added"))
                    CHRgr[i]$status <- "added"
                    score(CHRgr[CHRgr$name == name]) <- score(CHRgr[CHRgr$name == name]) + score(current)
                    start(CHRgr[CHRgr$name == name]) <- min(start(CHRgr[CHRgr$name == name]), start(current))
                    end(CHRgr[CHRgr$name == name]) <- max(end(CHRgr[CHRgr$name == name]), end(current))
                } else {
                    print(paste(i, "eliminated"))
                    CHRgr[i]$status <- "eliminated"
                }
            } else {
                kept <- c(kept, current)
                print(paste(i, "kept"))
                CHRgr[i]$status <- "kept"
            }
        }
    }
    
    return(CHRgr)
}

#' Plot cleavage site width distribution as a histogram
#' 
#' This function reads a BED file containing cleavage sites and plots
#' the distribution of their widths as a histogram.
#' 
#' @param CSbed Character string, path to the BED file containing cleavage sites
#' @return Histogram object showing the distribution of cleavage site widths
explore_CS <- function(CSbed) {
    cs_df <- read.delim(CSbed, header = FALSE)
    colnames(cs_df) <- c("chr", "start", "end", "id", "score", "strand")
    
    cs_df <- cs_df %>%
        mutate(width = end - start)
    brks <- c(1, 5, 10, 15, 20, 25, 1004)
    h <- hist(cs_df$width, breaks = brks)
}

#' Analyze nucleotide occurrence in FASTA sequences
#' 
#' This function analyzes the occurrence of a specific nucleotide in sequences from a FASTA file,
#' typically used to examine the sequence context around cleavage sites.
#' 
#' @param fa Character string, path to the FASTA file
#' @param ch Character, the nucleotide to count (e.g., "T")
#' @return Histogram showing the distribution of the specified nucleotide's occurrences
character_occurance_fasta <- function(fa, ch) {
    library(stringr)
    library(seqinr)
    library(parallel)
    
    fseq <- read.fasta(fa, seqtype = "DNA", as.string = TRUE, seqonly = TRUE)
    
    cl <- start_parallel(nc = 20)
    counts_char <- parLapply(cl, fseq, str_count, pattern = ch)
    stop_parallel(cl)
    
    counts_char <- unlist(counts_char)
    breaks <- seq(min(counts_char) - 1, max(counts_char), 1)
    
    print(hist(counts_char, breaks = breaks, xlab = paste("Occurance of", ch), ylab = "Number of reads", main = names(fa)))
}

#' Calculate Relative Cleavage Efficiency (RCE)
#' 
#' This function calculates the relative cleavage efficiency between control and experimental
#' conditions for cleavage sites, using Fisher's exact test for statistical significance.
#' 
#' @param count_df Data frame containing count data for cleavage sites
#' @param ctl_col Column name in count_df containing control counts
#' @param exp_col Column name in count_df containing experimental counts
#' @param RCE_cutoff Numeric, minimum mean RCE to include in results (default: 0.02)
#' @return Data frame with RCE values and statistical test results
cal_RCE <- function(count_df, ctl_col, exp_col, RCE_cutoff = 0.02) {
    cl <- start_parallel(nc = 20)
    clusterEvalQ(cl, expr = library(dplyr))
    clusterExport(cl, varlist = c("count_df", "ctl_col", "exp_col"), envir = environment())
    result_list <- parLapply(cl, seq(1, nrow(count_df)), function(i) {
        arow <- count_df[i, ]
        tx <- as.character(arow["tx_name"])
        sub_df <- count_df %>%
            filter(tx_name == tx)
        ctl_count_gene <- sum(sub_df[, ctl_col])
        ctl_count <- sum(arow[ctl_col])
        exp_count_gene <- sum(sub_df[, exp_col])
        exp_count <- sum(arow[exp_col])
        
        rce_ctl <- ctl_count / ctl_count_gene
        rce_exp <- exp_count / exp_count_gene
        rce_foldChange <- rce_exp / rce_ctl
        rce_mean <- (rce_ctl + rce_exp) / 2
        
        count_mat <- matrix(c(exp_count, exp_count_gene - exp_count, ctl_count, ctl_count_gene - ctl_count), ncol = 2, byrow = FALSE)
        f <- fisher.test(count_mat)
        p <- f$p.value
        e <- f$estimate
        
        new_value <- c(rce_ctl, rce_exp, rce_mean, rce_foldChange, p)
        names(new_value) <- c("RCE_ctl", "RCE_exp", "RCE_mean", "RCE_foldChange", "fisherP")
        return(c(arow, new_value))
    })
    stop_parallel(cl)
    
    result <- bind_rows(result_list)
    
    res <- result %>%
        filter(RCE_mean > RCE_cutoff) %>%
        mutate(fisherPadj = p.adjust(fisherP)) %>%
        mutate(name = as.factor(name)) %>%
        mutate(gene_name = as.factor(gene_name))
    
    invisible(res)
}

#' Determine UTR length changes between conditions
#' 
#' This function analyzes alternative polyadenylation (APA) events to determine
#' if genes show 3' UTR shortening or lengthening between conditions.
#' 
#' @param df Data frame containing RCE and statistical data for cleavage sites
#' @param gene_col Column name in df containing gene identifiers (default: "gene_name")
#' @param start_col Column name in df containing genomic start positions (default: "start")
#' @param strand_col Column name in df containing strand information (default: "strand")
#' @param sorting_col Column name in df for sorting sites (default: "RCE_mean")
#' @param padj_col Column name in df containing adjusted p-values (default: "fisherPadj")
#' @param direction_col Column name in df containing fold change direction (default: "RCE_foldChange")
#' @param padj_cutoff Numeric, significance threshold for adjusted p-values (default: 0.1)
#' @param norm Numeric, normalization factor for fold change (default: 1)
#' @return Data frame with UTR change annotations for each cleavage site
UTR_change <- function(df, gene_col = "gene_name", start_col = "start", strand_col = "strand", sorting_col = "RCE_mean", padj_col = "fisherPadj", direction_col = "RCE_foldChange", padj_cutoff = 0.1, norm = 1) {
    genes <- unique(df[[gene_col]])
    
    cl <- start_parallel(nc = 20)
    clusterEvalQ(cl, expr = library(dplyr))
    clusterExport(cl, varlist = c("df", "gene_col", "start_col", "strand_col", "sorting_col", "padj_col", "direction_col", "padj_cutoff", "norm"), envir = environment())
    change_df <- parLapply(cl, genes, function(gene) {
        gene_df <- df %>%
            filter(.data[[gene_col]] == gene) %>%
            arrange(desc(.data[[sorting_col]])) %>%
            mutate(position = NA)
        
        strand <- gene_df[1, strand_col]
        start <- gene_df[1, start_col]
        change <- "noChange"
        
        if (nrow(gene_df) > 1) {
            if (is.na(gene_df[1, direction_col])) {
                print(paste(gene, "has no change!"))
            } else {
                found <- FALSE
                for (i in 1:nrow(gene_df)) {
                    if (gene_df[i, padj_col] < padj_cutoff && gene_df[i, direction_col] < norm) {
                        for (j in 1:nrow(gene_df)) {
                            if (gene_df[j, direction_col] > norm) {
                                if (strand == "-") {
                                    change <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "lengthening", "shortening")
                                    gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                    gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                } else {
                                    change <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "shortening", "lengthening")
                                    gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                    gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                }
                                found <- TRUE
                                break
                            }
                        }
                    } else if (gene_df[i, padj_col] < padj_cutoff && gene_df[i, direction_col] > norm) {
                        for (j in 1:nrow(gene_df)) {
                            if (gene_df[j, direction_col] < norm) {
                                if (strand == "-") {
                                    change <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "shortening", "lengthening")
                                    gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                    gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                } else {
                                    change <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "lengthening", "shortening")
                                    gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                    gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                }
                                found <- TRUE
                                break
                            }
                        }
                    }
                    if (found) break
                }
                
                ## deal with no change
                if(!found){
                    for (i in 1:nrow(gene_df)) {
                        if (gene_df[i, direction_col] < norm) {
                            for (j in 1:nrow(gene_df)) {
                                if (gene_df[j, direction_col] > norm) {
                                    if (strand == "-") {
                                        gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                        gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                    } else {
                                        gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                        gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                    }
                                    found <- TRUE
                                    break
                                }
                            }
                        } else if (gene_df[i, direction_col] > norm) {
                            for (j in 1:nrow(gene_df)) {
                                if (gene_df[j, direction_col] < norm) {
                                    if (strand == "-") {
                                        gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                        gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                    } else {
                                        gene_df[j, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "proximal", "distal")
                                        gene_df[i, "position"] <- ifelse(gene_df[j, start_col] < gene_df[i, start_col], "distal", "proximal")
                                    }
                                    found <- TRUE
                                    break
                                }
                            }
                        }
                        if (found) break
                    }
                }
            }
        }
        # print(paste(gene, change))
        gene_df <- gene_df %>%
            mutate(change = change)
    })
    
    stop_parallel(cl)
    
    output <- bind_rows(change_df)
    
    return(output)
}


#' Assign Alternative Polyadenylation (APA) site types
#' 
#' This function assigns APA site types to cleavage sites:
#' S = constitutive (single site per gene)
#' P = proximal (most upstream site)
#' D = distal (all other sites, numbered from upstream to downstream)
#' 
#' @param df Data frame containing cleavage site information
#' @param gene_col Column name in df containing gene identifiers (default: "gene_name")
#' @param start_col Column name in df containing genomic start positions (default: "start")
#' @param strand_col Column name in df containing strand information (default: "strand")
#' @return Data frame with APA site type annotations for each cleavage site
APA_assign <- function(df, gene_col = "gene_name", start_col = "start", strand_col = "strand") {
    genes <- unique(df[[gene_col]])
    
    cl <- start_parallel(nc = 20)
    clusterEvalQ(cl, expr = library(dplyr))
    clusterExport(cl, varlist = c("df", "gene_col", "start_col", "strand_col"), envir = environment())
    
    apa <- parLapply(cl, genes, function(gene) {
        gene_df <- df %>%
            filter(.data[[gene_col]] == gene) %>%
            arrange(.data[[start_col]])
        strand <- gene_df[1, strand_col]
        start <- gene_df[1, start_col]
        prox <- vector()
        if (nrow(gene_df) == 1) {
            prox[1] <- "S"
        } else if (strand == "+") {
            prox[1] <- "P"
            for (i in 2:nrow(gene_df)) {
                prox[i] <- paste0("D", i - 1)
            }
        } else if (strand == "-") {
            prox[nrow(gene_df)] <- "P"
            for (i in 1:(nrow(gene_df) - 1)) {
                prox[i] <- paste0("D", nrow(gene_df) - i)
            }
        }
        
        gene_df$APA <- prox
        return(gene_df)
    })
    stop_parallel(cl)
    
    names(apa) <- genes
    apa_df <- bind_rows(apa)
    
    return(apa_df)
}


#' Detect polyadenylation signals (PAS) upstream of cleavage sites
#' 
#' This function identifies canonical and non-canonical polyadenylation signals
#' in the upstream region of cleavage sites.
#' 
#' @param BSgenome BSgenome object containing the reference genome
#' @param gr GRanges object containing cleavage sites
#' @param upstream Numeric, length of upstream region to search for PAS (default: 100)
#' @return GRanges object with PAS status annotations (PAS_strong, PAS_weak, or PAS_less)
detect_PAS_signal <- function(BSgenome, gr, upstream = 100) {
    USR <- flank(gr, width = upstream, start = TRUE, both = FALSE, use.names = TRUE, ignore.strand = FALSE)
    seqUSR <- toupper(getSeq(BSgenome, USR, as.character = TRUE))
    seqUSR <- gsub("T", "U", seqUSR, fixed = TRUE)
    
    PAS_strong <- "AAUAAA|AUUAAA"
    PAS_weak <- "AGUAAA|UAUAAA|CAUAAA|GAUAAA|AAUAUA|AAUACA|AAUAGA|ACUAAA|AAGAAA|AAUGAA"
    
    strong <- grepl(PAS_strong, seqUSR, perl = TRUE)
    weak <- grepl(PAS_weak, seqUSR, perl = TRUE)
    
    status <- rep("PAS_less", length(gr))
    for (i in seq_along(status)) {
        if (strong[i]) {
            status[i] <- "PAS_strong"
        } else if (weak[i]) {
            status[i] <- "PAS_weak"
        }
    }
    gr$status <- status
    return(gr)
}


#' Process DEXSeq results and create BED files for different APA categories
#' 
#' This function processes DEXSeq differential expression results for APA analysis
#' and creates BED files for different categories of cleavage sites (proximal/distal,
#' shortening/lengthening, etc.).
#' 
#' @param factor Character string, experimental factor name (default: "siZNF281")
#' @param dxrFile Character string, path to DEXSeq results file
#' @param outdir Character string, output directory for BED files (default: NULL)
#' @return List of gene IDs categorized by APA change (shortening, lengthening, stable)
results_from_dexseq_relative <- function(factor = "siZNF281", dxrFile, outdir = NULL) {
    message("[results_from_dexseq_relative] output bed files ...")
    dxr <- read.delim2(dxrFile, header = TRUE) %>%
        mutate(name = paste(groupID, APA, sep = ":")) %>%
        mutate(exonBaseMean=as.numeric(exonBaseMean))
    
    group_abundance <- tapply(dplyr::select(dxr, dplyr::starts_with("countData.gAAVS1")), dxr$groupID, FUN=sum)
    
    if (!is.null(outdir)) {
        abed <- dxr %>%
            select(genomicData.seqnames, genomicData.start, genomicData.end, name, exonBaseMean, genomicData.strand)
        bedfile <- file.path(outdir, paste0(factor, "_all_cleavage_site.bed"))
        write.table(abed, bedfile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
        
        noChange_forShortening_geneName <- NULL ## for noChange, make sure proximal and distal are from the same gene after random sampling
        noChange_forLengthening_geneName <- NULL
        for (apa in c("proximal", "distal")) {
            bed <- dxr %>%
                filter(grepl(apa, position)) %>%
                select(genomicData.seqnames, genomicData.start, genomicData.end, name, exonBaseMean, genomicData.strand)
            bedfile <- file.path(outdir, paste0(factor, "_", apa, "_cleavage_site.bed"))
            write.table(bed, bedfile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
            
            shortening_abundance <- NULL
            lengthening_abundance <- NULL
            for (ch in c("shortening", "lengthening", "noChange")) {
                abed <- dxr %>%
                    filter(grepl(apa, position) & change == ch) %>%
                    select(genomicData.seqnames, genomicData.start, genomicData.end, groupID, exonBaseMean, genomicData.strand)
                
                if(ch == "shortening" && apa == "proximal"){
                    shortening_abundance <- group_abundance[abed$groupID]
                }else if(ch == "lengthening" && apa == "proximal"){
                    lengthening_abundance <- group_abundance[abed$groupID]
                }
                
                if(ch == "noChange"){
                    if(apa == "proximal"){
                        noChange_abundance <- group_abundance[abed$groupID]
                        noChange_forShortening_geneName <- sub_sample_by_matching(noChange_abundance, shortening_abundance, oneTOmany = TRUE)
                        noChange_forLengthening_geneName <- sub_sample_by_matching(noChange_abundance, lengthening_abundance, oneTOmany = TRUE)
                        noChange_forShortening_geneName <- noChange_forShortening_geneName[sample(length(noChange_forShortening_geneName), length(shortening_abundance))]
                        noChange_forLengthening_geneName <- noChange_forLengthening_geneName[sample(length(noChange_lengthening_geneName), length(lengthening_abundance))]
                    }
                    
                    
                    abed_forShortening <- abed[abed$groupID %in% noChange_forShortening_geneName,]
                    abed_forLengthening <- abed[abed$groupID %in% noChange_forLengthening_geneName,]
                    
                    bedfileS <- file.path(outdir, paste0(factor, "_noChangeS_", apa, "_cleavage_site.bed"))
                    bedfileL <- file.path(outdir, paste0(factor, "_noChangeL_", apa, "_cleavage_site.bed"))
                    bedfile <- file.path(outdir, paste0(factor, "_noChange_", apa, "_cleavage_site.bed"))
                    write.table(abed_forShortening, bedfileS, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
                    write.table(abed_forLengthening, bedfileL, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
                    write.table(abed, bedfile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
                
                }else{
                    bedfile <- file.path(outdir, paste0(factor, "_", ch, "_", apa, "_cleavage_site.bed"))
                    write.table(abed, bedfile, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
                }
            }
        }
    }
    
    df <- dxr
    
    apa_gene_list <- list(
        shortening = df[df$change == "shortening", "groupID"],
        lengthening = df[df$change == "lengthening", "groupID"],
        stable = df[df$change == "noChange", "groupID"]
    )
    return(apa_gene_list)
}

#' Calculate distance between proximal and distal sites
#' 
#' This function calculates the distance between proximal and distal cleavage sites.
#' 
#' @param p Numeric vector, positions of proximal sites
#' @param d Numeric vector, positions of distal sites
#' @return Numeric, average distance between proximal and distal sites
#' Calculate average distance between proximal and distal sites
#'
#' This function calculates the average distance between proximal and distal cleavage sites.
#'
#' @param p Numeric vector, positions of proximal sites
#' @param d Numeric vector, positions of distal sites
#' @return Numeric, average distance between proximal and distal sites
pd_distance <- function(p=c(1,2), d=c(3,4)){
    d <- round(mean(d)) - round(mean(p))
}
