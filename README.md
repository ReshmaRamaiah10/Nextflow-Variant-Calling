# Nextflow-Variant-Calling
Simple but functional genomics pipeline for Variant Calling.

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
**Filtered Sample VCFs:**
```
Sample1_rmDup_variants.vcf
Sample2_rmDup_variants.vcf
```
**Final Merged VCF:**
```
merged_variants.vcf
```
## 🛠️ Workflow Steps

1. **Reference Indexing** (`bwa-mem2 index`) → Prepares the genome for alignment.
2. **Adapter Trimming** (`Trim Galore`) → Removes sequencing adapters from raw FASTQ files.
3. **Read Alignment** (`bwa-mem2 mem`) → Maps trimmed reads to the reference genome.
4. **Sorting & Duplicate Removal** (`samtools sort`, `picard MarkDuplicates`) → Sorts and removes PCR duplicates.
5. **Downsampling** (`picard DownsampleSam`) → Creates 75%, 50%, and 25% read subsets.
6. **Variant Calling** (`bcftools mpileup` & `call`) → Calls variants from full and downsampled BAM files.
7. **Filtering Variants** (`bcftools norm` & `view`) → Normalizes and filters variants.
8. **Merging VCFs** (`bcftools merge`) → Combines all final filtered VCFs into a single dataset.


## 🚀 How to Run the Workflow
1. Install Nextflow
If you haven’t installed Nextflow, do so using:
```
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/
```
2. Prepare Input Files
Ensure:
* `sample_list.txt` contains the correct sample names and paths to FASTQ files.
* `params.ref_fasta` points to the FASTA reference genome.
3. Run the Pipeline
Use **Docker** for reproducibility:
```
nextflow run main.nf -profile docker --input sample_list.txt --ref_fasta /path/to/genome.fa
```

## Notes
* Modify `params.ref_fasta` → Ensure the correct reference genome is provided.
* Check `sample_list.txt` → Ensure FASTQ paths are absolute or relative.
* Storage Considerations → BAM files and VCFs can be large; ensure enough disk space.
* **Parallel Execution** → The pipeline automatically processes multiple samples in parallel.
