/*
 * Indexing reference genome
 */

process BWA_INDEX {
    tag "Indexing reference genome"
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
