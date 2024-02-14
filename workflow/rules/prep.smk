def get_subsample_uri(wc):
    """
    Returns the uri for a given subsample's fastq.
    """
    df = pep.subsample_table[pep.subsample_table.sample_name == wc.sample]
    return list(set(df[df.subsample_name == wc.subsample].accession))

rule get_fq:
    """
    assumes provided accession code is actually a uri...

    alternative when an actually accession code: curl $(ffq --ftp {params.acc} | jq '.[0] | .url' | head -n 1 | tr -d '"') > {output}
    """
    output:
        temp("results/fastq/{sample}_subsample{subsample}.fastq.gz")
    params:
        acc = lambda wc: get_subsample_uri(wc)
        #acc = lambda wc: list(set(pep.sample_table[pep.sample_table.sample_name == wc.sample].accession))
    retries: 10
    resources:
        time=20,
        mem=24000,
        cpus=1
    priority: 100
    shell:
        """
        curl {params.acc} > {output}
        """

flatten = lambda t: [item for sublist in t for item in sublist]

def get_subsamples(wc):
    """
    Returns the experiment accessions for a given sample (experiment).
    """
    return list(set(pep.subsample_table[pep.subsample_table.sample_name == wc.sample].subsample_name))


rule concat_runs:
    """
    We concatenated fastqs from the same experiment (during sample table creation we enforce only R1, to keep a level playing field)
    """
    input:
        fqs = lambda wc: expand("results/fastq/{{sample}}_subsample{subsample}.fastq.gz",subsample = get_subsamples(wc))
    output:
        temp("results/concat/{sample}.fastq.gz")
    threads:
        1
    resources:
        time=120,
        mem=24000,
        cpus=1
    priority: 2
    shell:
        "cat {input.fqs} > {output}"


rule fastp:
    """
    in general this is true for anything pulled from SRA
    """
    input:
        rules.concat_runs.output
    output:
        fq = "results/trimmed/{sample}.fastp.fastq.gz",
        h = "results/trimmed/{sample}.fastp.html",
        j = "results/trimmed/{sample}.fastp.json",
    singularity:
        "docker://quay.io/biocontainers/fastp:0.23.2--h5f740d0_3"
    threads:
        4
    resources:
        time=20,
        mem=24000,
        cpus=4
    shell:
        """
        fastp --in1 {input} --out1 {output} -w {threads} -h {output.h} -j {output.j} 
        """


