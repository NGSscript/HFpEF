# HFpEF
This repository contains the scripts used for upstream raw-data processing and downstream analysis in this study. The workflow is organized as follows:

### 1. Upstream Raw-Data Processing
**RNA-seq Pipeline:**
RNAseq_rawdata_pre_processing.sh: Processes RNA-seq raw FASTQ files to generate the final gene expression matrix. It orchestrates the pipeline by calling the following two scripts:
1.Quantification_featurecounts.R: An R wrapper around featureCounts (Rsubread) that quantifies gene-level counts for a single sorted, deduplicated BAM. Called once per sample within the main script.
2.Merge_sample_quantification.py: Merges the per-sample count files into a single gene × sample expression matrix (with gene symbol annotation from the GTF) using the sample metadata. Executed once after all samples finish processing.
      
**ATAC-seq Pipeline:**
ATACseq_rawdata_pre_processing.sh: Processes ATAC-seq raw FASTQ files to generate final bigWig tracks and peak files.
**Chip-seq Pipeline:**
Chipseq_peak_processing.sh: Performs peak counting and bivalent peak definition for H3K4me3 and H3K27me3.
*Note: These upstream steps utilize tools including fastp, HISAT2, samblaster, featureCounts, bowtie2, sambamba, MACS2, and bamCoverage, as appropriate for each assay.*

### 2. Downstream Bioinformatics Analysis (R)
### 2. Downstream Bioinformatics Analysis (R)
Bioinformatics_analysis.R implements the comprehensive downstream bioinformatics workflow across multiple omics layers, including Bulk RNA-seq, ChIP-seq, ATAC-seq, and Metabolomics. The script is organized into the following major analytical modules:
1. Bulk RNA-seq Analysis (LNAFv & LNAFg datasets)
Differential Expression (DESeq2): Identifies differentially expressed genes (DEGs) across multiple pairwise comparisons (e.g., StageA/B/C vs. Control, Foxo6kd vs. Control).
Gene Clustering & Heatmaps: Clusters DEGs into distinct expression patterns (C1–C5) based on Z-scores and visualizes them using ComplexHeatmap. 
Batch Correction: Integrates independent RNA-seq datasets (LNAFv and LNAFg) using ComBat to remove batch effects and generates a combined PCA plot.
Functional Enrichment: Performs GO (Gene Ontology) biological process enrichment for the identified clusters using clusterProfiler and visualizes selected pathways with dot plots and line charts.
2. ChIP-seq Analysis (H3K4me3 & H3K27me3)
Peak Quantification: Calculates RPKM values for H3K4me3 peaks based on raw counts and peak widths. Classifies peaks into "High" and "Low" H3K4me3 occupancy based on the 85th percentile threshold.
Peak Annotation: Annotates ChIP-seq peaks (H3K4me3, H3K27me3, and bivalent peaks) using ChIPseeker to determine associated genes.
Permutation Testing: Evaluates the statistical significance of overlaps between specific epigenetic states (e.g., H3K4me3 High/Low, H3K27me3, bivalent peaks) and DEGs (Up/Down) using 100,000 random simulations.
3. ATAC-seq Analysis
Differential Accessibility: Utilizes DiffBind to identify differentially accessible regions (DARs) between different groups (StageC vs. Control, Foxo6 vs. StageC, Foxo6 vs. Control).
Peak Annotation: Annotates differential ATAC-seq peaks using ChIPseeker (mm10 genome) and extracts promoter-associated genes for downstream integration.
4. Multi-omics Data Integration
Integrates RNA-seq, ChIP-seq, and ATAC-seq data to identify key regulatory genes. 
5. Metabolomics Analysis
Uses OPLS-DA/limma to identify significantly altered metabolites between Control and StageC groups.
Generates heatmaps for specific differential metabolites and a large-scale heatmap across five different experimental groups.
