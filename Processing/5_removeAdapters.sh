#!/bin/bash
# script to removes reads with adaptor sequences in them from filtered reads

# necessary to perform cutadapt twice to remove all adaptor sequences

# removal of GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCCAATCTCGTATGCCGTCTTCTGCTTG
# is for

# removal of GATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
# is for

cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-1_cat_trim_filt_adapt1.fastq TES-1_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-1_cat_trim_filt_adapt.fastq TES-1_cat_trim_filt_adapt1.fastq
rm TES-1_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCACAATCTCGTATGCCGTCTTCTGCTTG -o TES-2_cat_trim_filt_adapt1.fastq TES-2_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-2_cat_trim_filt_adapt.fastq TES-2_cat_trim_filt_adapt1.fastq
rm TES-2_cat_trim_filt_adapt1.fastq

cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-3_cat_trim_filt_adapt1.fastq TES-3_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-3_cat_trim_filt_adapt.fastq TES-3_cat_trim_filt_adapt1.fastq
rm TES-3_cat_trim_filt_adapt1.fastq

cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCACAATCTCGTATGCCGTCTTCTGCTTG -o TES-4_cat_trim_filt_adapt1.fastq TES-4_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-4_cat_trim_filt_adapt.fastq TES-4_cat_trim_filt_adapt1.fastq
rm TES-4_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-5_cat_trim_filt_adapt1.fastq TES-5_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-5_cat_trim_filt_adapt.fastq TES-5_cat_trim_filt_adapt1.fastq
rm TES-5_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATCAATCTCGTATGCCGTCTTCTGCTTG -o TES-6_cat_trim_filt_adapt1.fastq TES-6_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-6_cat_trim_filt_adapt.fastq TES-6_cat_trim_filt_adapt1.fastq
rm TES-6_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-7_cat_trim_filt_adapt1.fastq TES-7_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-7_cat_trim_filt_adapt.fastq TES-7_cat_trim_filt_adapt1.fastq
rm TES-7_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGACAATCTCGTATGCCGTCTTCTGCTTG -o TES-8_cat_trim_filt_adapt1.fastq TES-8_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-8_cat_trim_filt_adapt.fastq TES-8_cat_trim_filt_adapt1.fastq
rm TES-8_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-9_cat_trim_filt_adapt1.fastq TES-9_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-9_cat_trim_filt_adapt.fastq TES-9_cat_trim_filt_adapt1.fastq
rm TES-9_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTCAATCTCGTATGCCGTCTTCTGCTTG -o TES-10_cat_trim_filt_adapt1.fastq TES-10_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-10_cat_trim_filt_adapt.fastq TES-10_cat_trim_filt_adapt1.fastq
rm TES-10_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACCAATCTCGTATGCCGTCTTCTGCTTG -o TES-11_cat_trim_filt_adapt1.fastq TES-11_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-11_cat_trim_filt_adapt.fastq TES-11_cat_trim_filt_adapt1.fastq
rm TES-11_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTACAATCTCGTATGCCGTCTTCTGCTTG -o TES-12_cat_trim_filt_adapt1.fastq TES-12_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-12_cat_trim_filt_adapt.fastq TES-12_cat_trim_filt_adapt1.fastq
rm TES-12_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG -o TES-13_cat_trim_filt_adapt1.fastq TES-13_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-13_cat_trim_filt_adapt.fastq TES-13_cat_trim_filt_adapt1.fastq
rm TES-13_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-14_cat_trim_filt_adapt1.fastq TES-14_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-14_cat_trim_filt_adapt.fastq TES-14_cat_trim_filt_adapt1.fastq
rm TES-14_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCACAATCTCGTATGCCGTCTTCTGCTTG -o TES-15_cat_trim_filt_adapt1.fastq TES-15_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-15_cat_trim_filt_adapt.fastq TES-15_cat_trim_filt_adapt1.fastq
rm TES-15_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-16_cat_trim_filt_adapt1.fastq TES-16_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-16_cat_trim_filt_adapt.fastq TES-16_cat_trim_filt_adapt1.fastq
rm TES-16_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTAGAGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-17_cat_trim_filt_adapt1.fastq TES-17_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-17_cat_trim_filt_adapt.fastq TES-17_cat_trim_filt_adapt1.fastq
rm TES-17_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-18_cat_trim_filt_adapt1.fastq TES-18_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-18_cat_trim_filt_adapt.fastq TES-18_cat_trim_filt_adapt1.fastq
rm TES-18_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACAATCTCGTATGCCGTCTTCTGCTTG -o TES-19_cat_trim_filt_adapt1.fastq TES-19_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-19_cat_trim_filt_adapt.fastq TES-19_cat_trim_filt_adapt1.fastq
rm TES-19_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-20_cat_trim_filt_adapt1.fastq TES-20_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-20_cat_trim_filt_adapt.fastq TES-20_cat_trim_filt_adapt1.fastq
rm TES-20_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-21_cat_trim_filt_adapt1.fastq TES-21_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-21_cat_trim_filt_adapt.fastq TES-21_cat_trim_filt_adapt1.fastq
rm TES-21_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-22_cat_trim_filt_adapt1.fastq TES-22_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-22_cat_trim_filt_adapt.fastq TES-22_cat_trim_filt_adapt1.fastq
rm TES-22_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGCAATCTCGTATGCCGTCTTCTGCTTG -o TES-23_cat_trim_filt_adapt1.fastq TES-23_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-23_cat_trim_filt_adapt.fastq TES-23_cat_trim_filt_adapt1.fastq
rm TES-23_cat_trim_filt_adapt1.fastq


cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGTAGCCAATCTCGTATGCCGTCTTCTGCTTG -o TES-24_cat_trim_filt_adapt1.fastq TES-24_cat_trim_filt.fastq
cutadapt -m 22 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o TES-24_cat_trim_filt_adapt.fastq TES-24_cat_trim_filt_adapt1.fastq
rm TES-24_cat_trim_filt_adapt1.fastq

