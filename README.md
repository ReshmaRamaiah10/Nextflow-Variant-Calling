# Nextflow-Variant-Calling
Simple but functional genomics pipeline for Variant Calling using GATK tutorial

## 1. Input File Format: sample_list.txt
Your sample_list.txt will look like this:

| Sample1    /path/to/data/Sample1_R1_001.fastq.gz    /path/to/data/Sample1_R2_001.fastq.gz
Sample2    /path/to/data/Sample2_R1_001.fastq.gz    /path/to/data/Sample2_R2_001.fastq.gz
Where each line contains:

Sample name (e.g., Sample1)
Path to R1 (e.g., /path/to/data/Sample1_R1_001.fastq.gz)
Path to R2 (e.g., /path/to/data/Sample1_R2_001.fastq.gz)
