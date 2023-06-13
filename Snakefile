s="FAW82928_pass_barcode01_2687e1a7_f06207b4_1"

# glob fastq files from folder
# align
#file_pattern="Fastq-passed/{sample}.fastq.gz"
file_pattern="alignment/barcode01/{sample}.bam"
SAMPLES,=glob_wildcards(file_pattern)



rule sort:
    input:
        "aligned.sam"
    output:
        "aligned.sorted.bam"
    container: "library://vi.ya/default/biotools"
    shell:
        "samtools sort -o {output} -O bam {input}"

rule align:
    input:
        "all.fastq"
    output:
        "aligned.sam"
    params:
        prefix="Fasta/barcode01.final"
    container: "docker://dceoy/bwa-mem2"
    threads: 10
    shell:
        "bwa-mem2 mem -t {threads} {params.prefix} {input} > {output}"

rule sort_bam:
    input:
        "tmp/fakemerged.bam"
    output:
        "tmp/fakemergedsorted.bam"
    container: "library://vi.ya/default/biotools"
    shell:
        "samtools sort -o {output} -O bam {input}"

rule merge_bam:
    input:
        expand(file_pattern, sample=SAMPLES)
    output:
        "tmp/fakemerged.bam"
    container: "library://vi.ya/default/biotools"
    shell:
        "python scripts/merge_bams.py {input} {output}"

rule merge_fastq:
    input:
        expand(file_pattern, sample=SAMPLES)
    output:
        "all.fastq"
    shell:
        "zcat {input} > {output}"


rule index_fasta:
    input:
        fasta="Fasta/barcode01.final.fasta"
    params:
        prefix="Fasta/barcode01.final"
    output:
        "Fasta/barcode01.final.bwt.2bit.64"
    container: "docker://dceoy/bwa-mem2"
    shell:
        "bwa-mem2 index -p {params.prefix} {input}" 

rule mpileup:
    input:
        fasta="Fasta/barcode01.final.fasta",
        bam="tmp/fakemergedsorted.bam"
    output:
        "tmp/mpileup"
    container: "library://vi.ya/default/biotools:latest"
    shell:
        "samtools mpileup --redo-BAQ -6 --output-MQ -f {input.fasta} {input.bam} > {output}"

rule makeFileList:
    input:
    output: "tmp/filelist"
    shell: "ls -d -1 alignment/barcode01/*bam > {output}"

#rule mpileup:
#    input:
#        fasta="Fasta/barcode01.final.fasta",
#        filelist="tmp/filelist"
#    output:
#        "tmp/mpileup"
#    container: "library://vi.ya/default/biotools:latest"
#    shell:
#        "samtools mpileup --redo-BAQ -6 -f {input.fasta} --output-MQ -b {input.filelist} > {output}"

rule plot_box:
        input:
            pileup="tmp/mpileup"
        output:
            png="plots/dots.png"
        params:
            start=1,
            end=3000
        container: "library://vi.ya/rnaseq-dbs/analysis_environment:latest"
        script:
            "scripts/plot.R"
            
rule fastqc:
        input:
            bam="alignment/barcode01/{sample}.bam"
        output:
            html="fastqc/{sample}_fastqc.html"
        params:
            outdir='fastqc'
        container: "library://vi.ya/rnaseq-dbs/qc:latest"
        script:
            "fastqc {input} -o {params.outdir}"
