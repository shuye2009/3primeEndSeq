# 3' End Sequencing Analysis Guide

This document provides instructions for running the 3' end sequencing analysis pipeline using the `run_3endseq.R` script.

## Overview

The 3' end sequencing analysis pipeline processes raw cleavage sites and internal priming sites from 3' end sequencing data, clusters the sites, and performs differential usage analysis to identify genes with alternative polyadenylation (APA) leading to 3' UTR shortening or lengthening.

## Prerequisites

- R (version 3.6 or higher)
- Required R packages:
  - GenomicRanges
  - DEXSeq
  - dplyr
  - ggplot2
  - pheatmap
  - GenomicPlot
  - and other dependencies

## Directory Structure

Before running the analysis, ensure you have the following directory structure:

```
working_directory/
├── merged_CSs.bed        # Merged cleavage sites
├── merged_IPs.bed        # Merged internal priming sites
├── sample_info_*.txt     # Sample information files (one per factor)
└── other input files     # Other required input files
```

## Running the Analysis

The analysis is executed using the `run_3endseq.R` script with the following command:

```bash
Rscript run_3endseq.R <working_directory> <factor> <option>
```

### Parameters

- `working_directory`: Path to the directory containing input files
- `factor`: Experimental factor (e.g., "gU2AF1", "gU170K", "gSF1")
- `option`: Analysis option (one of "global", "peak", "read", "other")

### Analysis Options

1. **global**: Processes single nucleotide resolution cleavage sites and internal priming sites, clusters them into consolidated peaks, and runs DEXSeq analysis to define shortening and lengthening genes.

2. **peak**: Integrates 3' end sequencing data with iCLIP peaks to analyze the relationship between RNA-binding proteins and cleavage sites.

3. **read**: Integrates 3' end sequencing data with iCLIP reads for more detailed analysis.

4. **other**: Performs additional analyses such as plotting iCLIP peaks with respect to gene features.

## Analysis Workflow

### Global Analysis

The global analysis performs the following steps:

1. Processes raw cleavage sites and internal priming sites
2. Clusters cleavage sites into consolidated peaks
3. Associates cleavage site clusters with genes
4. Performs DEXSeq analysis for different genomic features (3'UTR, Intron, Transcript)
5. Identifies genes with 3' UTR shortening or lengthening

### Integration with iCLIP Data

The pipeline can integrate 3' end sequencing data with iCLIP data to analyze the relationship between RNA-binding proteins and cleavage sites:

1. Plots iCLIP peaks around proximal and distal CPA sites
2. Analyzes the distribution of RNA-binding proteins around cleavage sites
3. Quantifies distances between proximal and distal CPA sites

## Output Files

The analysis generates various output files, including:

- Metagene plots for cleavage sites and internal priming sites
- Plots showing the distribution of sites around 3' UTR regions
- DEXSeq results tables
- BED files for different categories of cleavage sites
- Plots showing the relationship between RNA-binding proteins and cleavage sites

## Example Usage

```bash
# Run global analysis for factor gU2AF1
Rscript run_3endseq.R /path/to/working/directory gU2AF1 global

# Integrate 3' end sequencing data with iCLIP peaks
Rscript run_3endseq.R /path/to/working/directory gU2AF1 peak

# Integrate 3' end sequencing data with iCLIP reads
Rscript run_3endseq.R /path/to/working/directory gU2AF1 read
```

## Sample Information File Format

The sample information file (`sample_info_*.txt`) should be a tab-delimited file with the following columns:

- `sample`: Sample name
- `condition`: Condition (e.g., "control" or "treat")
- Other metadata columns as needed

Example:
```
sample  condition
sample1 control
sample2 control
sample3 treat
sample4 treat
```

## Additional Resources

For more information on the functions used in this pipeline, refer to the documentation in the following files:

- `3endseq.R`: Contains the main processing functions
- `3endseq_lib.R`: Contains utility functions for the pipeline
