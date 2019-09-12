#!/bin/bash
# script to map reads to aphid genome using hisat2
# and to Buchnera genome with bwa

hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-1_cat_trim_filt_adapt.fastq -S TES-1_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile 
TES-1_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-1_cat_trim_filt_adapt.fastq > TES-1_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-2_cat_trim_filt_adapt.fastq -S TES-2_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-2_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-2_cat_trim_filt_adapt.fastq > TES-2_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-3_cat_trim_filt_adapt.fastq -S TES-3_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-3_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-3_cat_trim_filt_adapt.fastq > TES-3_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-4_cat_trim_filt_adapt.fastq -S TES-4_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-4_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-4_cat_trim_filt_adapt.fastq > TES-4_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-5_cat_trim_filt_adapt.fastq -S TES-5_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-5_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-5_cat_trim_filt_adapt.fastq > TES-5_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-6_cat_trim_filt_adapt.fastq -S TES-6_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-6_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-6_cat_trim_filt_adapt.fastq > TES-6_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-7_cat_trim_filt_adapt.fastq -S TES-7_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-7_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-7_cat_trim_filt_adapt.fastq > TES-7_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-8_cat_trim_filt_adapt.fastq -S TES-8_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-8_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-8_cat_trim_filt_adapt.fastq > TES-8_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-9_cat_trim_filt_adapt.fastq -S TES-9_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-9_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-9_cat_trim_filt_adapt.fastq > TES-9_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-10_cat_trim_filt_adapt.fastq -S TES-10_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-10_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-10_cat_trim_filt_adapt.fastq > TES-10_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-11_cat_trim_filt_adapt.fastq -S TES-11_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-11_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-11_cat_trim_filt_adapt.fastq > TES-11_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-12_cat_trim_filt_adapt.fastq -S TES-12_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-12_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-12_cat_trim_filt_adapt.fastq > TES-12_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-13_cat_trim_filt_adapt.fastq -S TES-13_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-13_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-13_cat_trim_filt_adapt.fastq > TES-13_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-14_cat_trim_filt_adapt.fastq -S TES-14_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-14_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-14_cat_trim_filt_adapt.fastq > TES-14_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-15_cat_trim_filt_adapt.fastq -S TES-15_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-15_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-15_cat_trim_filt_adapt.fastq > TES-15_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-16_cat_trim_filt_adapt.fastq -S TES-16_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-16_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-16_cat_trim_filt_adapt.fastq > TES-16_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-17_cat_trim_filt_adapt.fastq -S TES-17_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-17_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-17_cat_trim_filt_adapt.fastq > TES-17_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-18_cat_trim_filt_adapt.fastq -S TES-18_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-18_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-18_cat_trim_filt_adapt.fastq > TES-18_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-19_cat_trim_filt_adapt.fastq -S TES-19_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-19_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-19_cat_trim_filt_adapt.fastq > TES-19_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-20_cat_trim_filt_adapt.fastq -S TES-20_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-20_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-20_cat_trim_filt_adapt.fastq > TES-20_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-21_cat_trim_filt_adapt.fastq -S TES-21_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-21_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-21_cat_trim_filt_adapt.fastq > TES-21_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-22_cat_trim_filt_adapt.fastq -S TES-22_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-22_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-22_cat_trim_filt_adapt.fastq > TES-22_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-23_cat_trim_filt_adapt.fastq -S TES-23_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-23_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-23_cat_trim_filt_adapt.fastq > TES-23_cat_trim_filt_adapt_bwa_Ba.mem.sam


hisat2 -x GCF_000142985.2_Acyr_2.0_genomic_index.fna -U TES-24_cat_trim_filt_adapt.fastq -S TES-24_cat_trim_filt_adapt_hisat2_Ap.sam --phred33 --novel-splicesite-outfile TES-24_cat_trim_filt_adapt_hisat2_Ap.junctions --rna-strandness R --dta -t

bwa mem GCF_000009605.1_ASM960v1_genomic_index.fna TES-24_cat_trim_filt_adapt.fastq > TES-24_cat_trim_filt_adapt_bwa_Ba.mem.sam

