rule star_se:
    input:
        fq1="results/trimmed/{sample}.fastp.fastq.gz",
        # path to STAR reference genome index
        idx="results/star-indices/{ref}/",
    output:
        # see STAR manual for additional output files
        aln="results/star-se/{ref}/{sample}/se_aligned.bam",
        log="results/star-se/{ref}/{sample}/Log.out",
        log_final="results/star-se/{ref}/{sample}/Log.final.out",
        reads_per_gene="results/star-se/{ref}/{sample}/ReadsPerGene.out.tab", # need quantMode GeneCounts for this to work, but we want to quantify after removing multimappers
        unmapped="results/star-se/{ref}/{sample}/unmapped.fastq",
    log:
        "results/logs/star-se/{sample}_{ref}.log",
    threads:
        24
    resources:
        time=60,
        mem=48000,
        cpus=24
    params:
        # optional parameters
        extra="--outSAMtype BAM Unsorted --quantMode GeneCounts --outFilterMultimapNmax 1 --outSAMmultNmax 1",
    wrapper:
        "v3.3.3/bio/star/align"


rule samtools_prep:
    """
    initial filtering steps
    as with bwa rule, meant to mimic the most important steps of nf.core/chipseq
    samtools view excludes non-primary mappings and unmapped reads

    exclude multimappers (q < 1)
    """
    input:
        "results/star-se/{ref}/{sample}/se_aligned.bam",
    output:
        temp("results/bams/{ref}/{sample}.initial_filt.bam")
    params:
        args = "-F 0x004 -F 0x0100 -q 1"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.17--hd87286a_1"
    shell:
        """
        samtools view {params.args} -O BAM {input} > {output}
        """
        

# probably don't want to remove duplicates for expression data
# rule picard_sort_by_name:
#     input:
#         "results/bams/{ref}/{sample}.initial_filt.bam"
#     output:
#         bam = temp("results/bams/{ref}/{sample}.nsort.bam"),
#     singularity:
#         "docker://quay.io/biocontainers/picard:3.0.0--hdfd78af_1"
#     resources:
#         time=40,
#         mem=24000,
#         cpus=1
#     shell:
#         """
#         picard SortSam -I {input} -O {output.bam} -SO queryname 
#         """

# rule picard_markdup:
#     """
#     namesorted input is important so that multiple alignments of the same read are adjacent and therefore removable
#     see picard markduplicates docs for more details

#     https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
#     """
#     input:
#         rules.picard_sort_by_name.output.bam
#     output:
#         bam = temp("results/bams/{ref}/{sample}.markdup.bam"),
#         stats = "results/bams/{ref}/{sample}.markdup.stats"
#     singularity:
#         "docker://quay.io/biocontainers/picard:3.0.0--hdfd78af_1"
#     resources:
#         time=40,
#         mem=24000,
#         cpus=1
#     shell:
#         """
#         picard MarkDuplicates --REMOVE_DUPLICATES true -I {input} -M {output.stats} -O {output.bam} 
#         """


rule samtools_sort:
    input:
        "results/bams/{ref}/{sample}.initial_filt.bam",
    output:
        bam = "results/bams/{ref}/{sample}.srt.bam",
        bai = "results/bams/{ref}/{sample}.srt.bam.bai"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.17--hd87286a_1"
    resources:
        time=30,
        mem=24000,
        cpus=18
    threads:
        18
    priority: 2
    shell:
        """
        samtools sort -@ {threads} -m 1G {input} -o {output.bam} &&
        samtools index -@ {threads} {output.bam}
        """