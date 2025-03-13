# Nextflow-Variant-Calling
Simple but functional genomics pipeline for Variant Calling using GATK.

This Nextflow pipeline performs genotype variant calling from paired-end FASTQ files, processing multiple samples in parallel. The workflow trims adapters, aligns reads, marks duplicates, downsamples BAM files, calls variants, filters them, and merges VCFs for final analysis.

The pipeline is designed for high-throughput sequencing data and leverages bwa-mem2, samtools, picard, and bcftools inside Docker containers for reproducibility.

## 📂 Input Files

* **`sample_list.txt`** → A tab-delimited file containing sample names and paths to paired-end FASTQ files:
```text
  Sample1    /path/to/Sample1_R1.fastq.gz    /path/to/Sample1_R2.fastq.gz
  Sample2    /path/to/Sample2_R1.fastq.gz    /path/to/Sample2_R2.fastq.gz
```
* Reference Genome (`params.ref_fasta`) → FASTA file for alignment and variant calling.

## 📤 Output Files

For each sample (`Sample1`, `Sample2`, etc.), the following files will be generated:

Intermediate Files
* `Sample1_val_1.fq.gz`, `Sample1_val_2.fq.gz` → Trimmed FASTQ files
* `Sample1_sort.bam` → Sorted BAM file
* `Sample1_rmDup.bam` → BAM file after duplicate removal
* `Sample1_75pcReads.bam`, Sample1_50pcReads.bam, Sample1_25pcReads.bam → Downsampled BAMs
* `Sample1_bcftools.vcf` → Raw variant calls

##  Final Outputs
✅ Filtered Sample VCFs:
```
Sample1_rmDup_variants.vcf
Sample2_rmDup_variants.vcf
```
...
✅ Final Merged VCF:
```
merged_variants.vcf
```


