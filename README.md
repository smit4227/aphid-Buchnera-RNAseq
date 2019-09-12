# aphid-Buchnera-RNAseq
Scripts for processing and analysis of aphid and Buchnera transcripts

Thomas E Smith, University of Texas at Austin, Department of Intregrative Biology, Moran lab

These scripts are associated with the publication, "...", doi:... . In this study, we looked at differences in gene expression among seven different A. pisum genotypes displaying a wide range of Buchnera titer, or abundance. This trait is dependent on the host-genotype; the Buchnera are genetically identical across all host genotypes. We found that Buchnera gene expression differs significantly among host genotypes, as the bacteriocytes that house Buchnera exhibit distinct cellular environments in low titer versus high titer aphid genotypes.



What you will need:

1) RNA-Seq reads are available for download from NCBI (PRJNA564856)
2) The "Reference Files" folder contains the pea aphid and Buchnera genome annotation files, but you will need to download the genome sequence files from NCBI. We used "GCF_000142985.2_Acyr_2.0_genomic.fna" for A. pisum and "GCF_000009605.1_ASM960v1_genomic.fna" for Buchnera.
3) bash programs:
  fastqc v0.11.5
  FASTX-Toolkit v0.0.13
  cutadapt v1.16
  hisat2 v 2.1.0
  bwa v0.7.17
  htseq-count v0.9.1
4) R programs:
  R v3.6.1
  RStudio v1.1.463
  DESeq2 v1.24.0
  AnnotationHub v2.16.0
  WGCNA v1.68
  flashClust v1.01-2
  stats v3.6.1
  dendextend v1.12.0
  pheatmap v1.0.12
  RColorBrewer v1.1-2
  ggplot2 v3.2.1
  reshape v0.8.8
  clusterProfiler v3.12.0
  DOSE v3.10.2
  ggraph v1.0.2
  igraph v1.2.4.1
  enrichplot v1.4.0
  pathview v1.24.0
  png v0.1-7
  KEGGREST v1.24.0
  ggfortify v0.4.7
  clustsig v1.1
  ape v5.3
  ggtree v1.16.5
  
  
  
What to do:

1) Raw reads are processed using bash scripts in the "Processing" folder. PRJNA564856 reads have been processed up to script "5_removeAdapters.sh".
2) Follow the remaining bash scripts in the "Processing" folder in order from "6_readQuality.sh" to "9_count.sh". This will map and count the reads for both the A. pisum and Buchnera genomes.
3) Run the "deseq2.R" R script in the "DESeq2" folder to calculate differential expression
4) To correlate host and symbiont gene expression, follow the R scripts in the "HostSymbiontGeneCorrelation" folder from "1_CorrelationAnalysis.R" to "4_EnrichmentTables".
5) Many figures were constructed in R using custom scripts located in the "Plots" folder. Follow R scripts from "1_PCA.R" to "3_SIMPROF.R" to generate them. 
6) For WGCNA, follow the R scripts located in foler "WGCNA", from "1_DataInputAndCleaning.R" to "5_FunctionalEnrichmentBa.R".

Note: Functional enrichment analyses using clusterProfiler and pathview R packages utilize customized "cnetplot" and "pathview" functions. These have "TES" appended to the name and are located in the "CustomRFunctions" folder.


Updated Sep 2019
