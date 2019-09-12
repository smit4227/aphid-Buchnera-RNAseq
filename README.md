# aphid-Buchnera-RNAseq
Scripts for processing and analysis of aphid and Buchnera transcripts

Thomas E Smith, University of Texas at Austin, Department of Intregrative Biology, Moran lab

These scripts are associated with the publication, "...", doi:... . In this study, we looked at differences in gene expression among seven different A. pisum genotypes displaying a wide range of Buchnera titer, or abundance. This trait is dependent on the host-genotype; the Buchnera are genetically identical across all host genotypes. We found that Buchnera gene expression differs significantly among host genotypes, as the bacteriocytes that house Buchnera exhibit distinct cellular environments in low titer versus high titer aphid genotypes.

What you will need:
1) RNA-Seq reads are available for download from NCBI (PRJNA564856)
2) Download the "Reference Files" folder containing pea aphid and Buchnera genomes and annotation files
3) bash programs:
  fastqc
  fastx
  cutadapt
  hisat2
  bwa
  htseq-count
4) R programs:
  DESeq2
  rapportools
  AnnotationHub
  WGCNA
  flashClust
  stats
  dendextend
  pheatmap
  RColorBrewer
  ggplot2
  reshape
  clusterProfiler
  DOSE
  ggraph
  igraph
  enrichplot
  pathview
  png
  KEGGREST
  ggfortify
  clustsig
  ape
  ggtree
  
  
  
What to do:
1) Raw reads are processed using bash scripts in the "Processing" folder. PRJNA564856 reads have been processed up to script "5_removeAdapters.sh".
2) Follow the remaining bash scripts in the "Processing" folder in order from "6_readQuality.sh" to "9_count.sh". This will map and count the reads for both the A. pisum and Buchnera genomes.
3) Run the "deseq2.R" R script in the "DESeq2" folder to calculate differential expression
4) To correlate host and symbiont gene expression, follow the R scripts in the "HostSymbiontGeneCorrelation" folder from "1_CorrelationAnalysis.R" to "4_EnrichmentTables".
5) Many figures were constructed in R using custom scripts located in the "Plots" folder. Follow R scripts from "1_PCA.R" to "3_SIMPROF.R" to generate them. 
6) For WGCNA, follow the R scripts located in foler "WGCNA", from "1_DataInputAndCleaning.R" to "5_FunctionalEnrichmentBa.R".


Updated Sep 2019
