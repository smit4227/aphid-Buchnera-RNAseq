#!/bin/bash
# script to count mapped reads

# requires gene annotation files for A. pisum and Buchnera
# aphid gene annotation file = GCF_000142985.2_Acyr_2.0_genomic.gff
# Buchnera genome is ASM960v1 from NCBI, file = GCF_000009605.1_ASM960v1_genomic.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-1_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-1_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-1_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-1_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-2_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-2_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-2_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-2_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-3_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-3_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-3_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-3_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-4_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-4_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-4_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-4_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-5_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-5_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-5_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-5_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-6_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-6_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-6_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-6_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-7_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-7_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-7_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-7_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-8_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-8_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-8_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-8_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-9_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-9_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-9_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-9_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-10_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-10_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-10_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-10_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-11_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-11_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-11_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-11_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-12_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-12_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-12_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-12_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-13_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-13_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-13_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-13_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-14_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-14_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-14_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-14_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-15_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-15_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-15_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-15_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-16_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-16_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-16_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-16_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-17_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-17_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-17_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-17_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-18_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-18_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-18_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-18_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-19_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-19_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-19_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-19_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-20_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-20_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-20_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-20_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-21_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-21_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-21_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-21_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-22_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-22_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-22_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-22_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-23_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-23_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-23_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-23_cat_trim_filt_adapt_bwa_Ba_gene.gff


htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-24_cat_trim_filt_adapt_hisat2_Ap.sam GCF_000142985.2_Acyr_2.0_genomic.gff > TES-24_cat_trim_filt_adapt_hisat2_Ap_gene.gff

htseq-count -m intersection-nonempty --stranded reverse -i ID -t gene TES-24_cat_trim_filt_adapt_bwa_Ba.mem.sam GCF_000009605.1_ASM960v1_genomic.gff > TES-24_cat_trim_filt_adapt_bwa_Ba_gene.gff


