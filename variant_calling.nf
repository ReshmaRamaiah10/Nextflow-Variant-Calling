nextflow.enable.dsl=2

params.input = "sample_list.txt"  // Path to the sample list file
params.ref_fasta = "/path/to/reference/genome.fa"  // Update with the actual path to reference genome

process BWA_INDEX {
    tag "Indexing Reference Genome"
    container 'community.wave.seqera.io/library/bwa-mem2:2.2.1--1842774b9b0b4729'

    input:
    path ref_genome

    output:
    path "hg38_bw2.*", emit: index_files

    script:
    """
    bwa-mem2 index -p hg38_bw2 ${ref_genome}
    """
}

process TRIM_GALORE {
    tag "$sample - Trimming"
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_val_1.fq.gz"), path("${sample}_val_2.fq.gz"), emit: trimmed_reads

    script:
    """
    trim_galore --paired --basename ${sample} -j 8 ${read1} ${read2}
    """
}

process BWA_ALIGN {
    tag "$sample - Aligning"
    container 'community.wave.seqera.io/library/bwa-mem2:2.2.1--1842774b9b0b4729'

    input:
    tuple val(sample), path(trimmed_read1), path(trimmed_read2), path(index_files)

    output:
    tuple val(sample), path("${sample}_sort.bam"), emit: sorted_bam

    script:
    """
    bwa-mem2 mem -t 8 hg38_bw2 ${trimmed_read1} ${trimmed_read2} | \
    samtools sort -@ 8 -o ${sample}_sort.bam -
    """
}

process MARK_DUPLICATES {
    tag "$sample - Removing Duplicates"
    container 'broadinstitute/picard:2.27.4'

    input:
    tuple val(sample), path(sorted_bam)

    output:
    tuple val(sample), path("${sample}_rmDup.bam"), emit: dedup_bam

    script:
    """
    picard MarkDuplicates \
        VERBOSITY=WARNING \
        INPUT=${sorted_bam} \
        OUTPUT=${sample}_rmDup.bam \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        METRICS_FILE=tmp_${sample}_rmDup.log
    """
}

process DOWNSAMPLE {
    tag "$sample - Downsampling $percent%"
    container 'broadinstitute/picard:2.27.4'

    input:
    tuple val(sample), path(dedup_bam)
    val percent

    output:
    tuple val(sample), path("${sample}_${percent}pcReads.bam"), emit: downsampled_bam

    script:
    """
    picard DownsampleSam \
        INPUT=${dedup_bam} \
        OUTPUT=${sample}_${percent}pcReads.bam \
        RANDOM_SEED=50 \
        PROBABILITY=0.${percent} \
        VALIDATION_STRINGENCY=SILENT
    samtools index ${sample}_${percent}pcReads.bam
    """
}

process VARIANT_CALLING {
    tag "$sample - Calling Variants"
    container 'biocontainers/bcftools:v1.17-1-deb_cv1'

    input:
    tuple val(sample), path(dedup_bam), path(downsampled_bams), path(ref_genome)

    output:
    tuple val(sample), path("${sample}_bcftools.vcf"), emit: raw_vcf

    script:
    """
    bcftools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --annotate DP,AD \
        -f ${ref_genome} ${dedup_bam} ${downsampled_bams} | \
    bcftools call --multiallelic-caller --variants-only -Ov -o ${sample}_bcftools.vcf
    """
}

process VCF_FILTER {
    tag "$sample - Filtering Variants"
    container 'biocontainers/bcftools:v1.17-1-deb_cv1'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(raw_vcf), path(ref_genome)

    output:
    tuple val(sample), path("${sample}_bcftools-final.vcf"), path("${sample}_rmDup_variants.vcf"), emit: filtered_vcf
    path "${sample}_rmDup_variants.vcf", emit: rmdup_vcf

    script:
    """
    bcftools norm -Ou -m-any ${raw_vcf} | \
    bcftools norm -Ov -f ${ref_genome} -o ${sample}_bcftools-final.vcf

    bcftools view ${raw_vcf} -c 1 > ${sample}_rmDup_variants.vcf
    """
}

process MERGE_VCF {
    tag "Merging All VCFs"
    container 'biocontainers/bcftools:v1.17-1-deb_cv1'
    publishDir params.outdir, mode: 'copy'

    input:
    path vcf_files

    output:
    path "merged_variants.vcf", emit: merged_vcf

    script:
    """
    bcftools merge -m all -Ov ${vcf_files.join(' ')} > merged_variants.vcf
    """
}

workflow VariantCalling {
    samples = Channel
        .fromPath(params.input) // Read input list from sample_list.txt
        .splitCsv(header: false, sep: "\t")  // Split based on tab separator
        .map { row -> tuple(row[0], file(row[1]), file(row[2])) } // Create a tuple with sample name, R1, and R2

    index_files = BWA_INDEX(file(params.ref_fasta))

    trimmed_reads = TRIM_GALORE(samples)

    sorted_bams = BWA_ALIGN(trimmed_reads, index_files)

    dedup_bams = MARK_DUPLICATES(sorted_bams)

    downsampled_bams = Channel.of(75, 50, 25)
        .cross(dedup_bams)
        .map { percent, bam -> DOWNSAMPLE(bam, percent) }

    raw_vcfs = VARIANT_CALLING(dedup_bams, downsampled_bams, file(params.ref_fasta))

    final_vcfs = VCF_FILTER(raw_vcfs, file(params.ref_fasta))

    rmdup_vcf_files = final_vcfs.map { it[1] } // Extract only the "${sample}_rmDup_variants.vcf"

    MERGE_VCF(rmdup_vcf_files.collect())
}
