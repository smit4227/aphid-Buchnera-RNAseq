#!/bin/bash
# script to build hisat2 index files from aphid genome
# and to build bwa index files from Buchnera genome

# aphid genome is Acyr2.0 from NCBI, file = GCF_000142985.2_Acyr_2.0_genomic.fna
# Buchnera genome is ASM960v1 from NCBI, file = GCF_000009605.1_ASM960v1_genomic.fna

hisat2 -build GCF_000142985.2_Acyr_2.0_genomic.fna
GCF_000142985.2_Acyr_2.0_genomic_index.fna

bwa -index -a bwtsw GCF_000009605.1_ASM960v1_genomic.fna GCF_000009605.1_ASM960v1_genomic_index.fna