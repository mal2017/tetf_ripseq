
rule samtools_idxstats:
    input:
        "results/bams/{ref}/{sample}.srt.bam"
    output:
        "results/stats/{ref}/{sample}.idxstats.txt"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.17--hd87286a_1"
    threads: 16
    resources:
        time=20,
        mem=24000,
        cpus=16
    shell:
        """
        samtools idxstats {input} > {output}
        """

rule samtools_stats:
    input:
        "results/bams/{ref}/{sample}.srt.bam"
    output:
        "results/stats/{ref}/{sample}.stats.txt"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.17--hd87286a_1"
    threads: 16
    resources:
        time=20,
        mem=24000,
        cpus=16
    shell:
        """
        samtools stats -@ {threads} {input} > {output}
        """


localrules: collect_samtools_stats, collect_samtools_idxstats, 

rule collect_samtools_stats:
    input:
        expand("results/stats/{{ref}}/{sample}.stats.txt", sample=SAMPLES)
    output:
        "results/stats/{ref}/all_stats.txt"
    shell:
        """
        echo "" > {output};
        for file in {input}; 
            do
                paste <(echo $file) <(grep 'MQ0' $file) | cut -f 1,3,4,5 >> {output};
                paste <(echo $file) <(grep 'reads mapped:' $file) | cut -f 1,3,4,5 >> {output};
            done
        """

rule collect_samtools_idxstats:
    """
    https://unix.stackexchange.com/questions/117568/adding-a-column-of-values-in-a-tab-delimited-file
    """
    input:
        expand("results/stats/{{ref}}/{sample}.idxstats.txt", sample=SAMPLES)
    output:
        "results/stats/{ref}/all_idxstats.txt"
    shell:
        """
        echo "" > {output};
        for file in {input}; 
            do
                paste $file <(yes $file | head -n $(cat $file | wc -l)) >> {output}; 
            done
        """