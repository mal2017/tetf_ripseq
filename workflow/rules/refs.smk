
rule copy_unmasked_genome:
    input:
        refs_wf("results/flybase-anno/genome.fasta.gz")
    output:
        "results/refs/unmasked-genome.fasta"
    shell:
        """
        cat {input} | gunzip -c > {output}
        """

rule copy_masked_genome_plus_tes:
    """
    star needs features to be exon
    """
    input:
        refs_wf("results/combined-anno/genome-plus-tes.fasta.gz")
    output:
        "results/refs/masked-genome-plus-tes.fasta"
    shell:
        """
        cat {input} | gunzip -c | awk '{{ if ($2 == "TIDAL") $3 = "exon"; print }}' > {output}
        """

rule copy_masked_genome_plus_tes_gtf:
    """
    star needs features to be exon
    """
    input:
        refs_wf("results/combined-anno/transcripts-plus-tes.gtf")
    output:
        "results/refs/masked-genome-plus-tes.gtf"
    shell:
        """
        cat {input} | awk '{{ if ($2 == "TIDAL") $3 = "exon"; print }}' > {output}
        """

rule copy_unmasked_genome_gtf:
    input:
        refs_wf("results/flybase-anno/transcriptome.gtf.gz")
    output:
        "results/refs/unmasked-genome.gtf"
    shell:
        """
        cat {input} | gunzip -c > {output}
        """


rule star_index:
    """
    fasta must be unzipped
    """
    input:
        fasta="results/refs/{ref}.fasta",
        gtf = "results/refs/{ref}.gtf"
    output:
        directory("results/star-indices/{ref}/"),
    params:
        extra="--genomeSAindexNbases 12",
    threads:
        24
    resources:
        time=60,
        mem=48000,
        cpus=24
    log:
        "results/logs/star_index_{ref}.log",
    wrapper:
        "v3.3.3/bio/star/index"