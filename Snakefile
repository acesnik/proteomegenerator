shell.executable("/bin/bash")

singularity: "docker://continuumio/miniconda3:4.5.4"

import os, multiprocessing
SMDir = os.getcwd()

# PAR initializes to template parameter file unless one is provided
MQ = SMDir + "/MaxQuant/bin/MaxQuantCmd.exe"
PAR = SMDir + "/MaxQuant/mqpar_template.xml"

### Begin User Variables

# Directories
WD = "/mnt/f/ProjectsActive/Spritz/ProteomeGenerator"
workdir: WD
TMP = os.path.join(WD, 'tmp')
UNIPROT="/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/reference/UP000005640.fasta"

# References,
GTF='/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/reference/gencode.v21.annotation.gtf'
INDEX = '/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/reference/STAR'
FASTA='/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/reference/GRCh38.genome.fa'
# # Comment out or delete the "GTF" variable to force de novo mode

# Samples,  K052
# SAMPLES = "SRR6425178"
# R1="/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/RNASeq/K052/{sample}_1.fastq"
# R2="/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/RNASeq/K052/{sample}_2.fastq"

# References, Jurkat
# Comment out or delete the "GTF" variable to force de novo mode
SAMPLES = ["SRR791578", "SRR791579", "SRR791580", "SRR791581", "SRR791582", "SRR791583", "SRR791584", "SRR791585", "SRR791586"]
R1="/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/RNASeq/Jurkat/{sample}_1.fastq"
R2="/mnt/f/ProjectsActive/Spritz/ProteomeGenerator/RNASeq/Jurkat/{sample}_2.fastq"

# MaxQuant
RAW = os.path.join(WD, "160116_K052_OffLRP_RP_f01.raw")
THREADS = str(multiprocessing.cpu_count())

### End User Variables

if os.path.exists(SMDir + "/MaxQuant/mqpar.xml") :
    PAR = SMDir + "/MaxQuant/mqpar.xml"

if 'GTF' in locals():
    ruleorder: STAR_GTF > STAR_denovo
    MODELS = 'merged reference'.split()
else:
    ruleorder: STAR_denovo > STAR_GTF
    MODELS = 'merged'

snakemake.utils.makedirs(TMP)
snakemake.utils.makedirs('out/benchmarks')

rule all:
    input: expand("out/all-merge/{model}/proteome.unique.fasta", model=MODELS)

include: "pgm.smk"

rule STAR_GTF:
    input: r1=R1, r2=R2, index=INDEX
    # output: "out/{sample}.Aligned.sortedByCoord.out.bam"
    output: "out/{sample}.Aligned.sortedByCoord.out.bam"
    benchmark: "out/benchmarks/{sample}.align.json"
    log: "out/logs/{sample}.align.txt"
    conda: "envs/myenv.yaml"
    threads: 12
    params: R="'span[hosts=1] rusage[mem=20]'", J="align", o="out/logs/align.out", eo="out/logs/align.err"
    shell: "STAR \
        --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix out/{wildcards.sample}. \
        --outSAMattributes NH HI XS \
        --outSAMattrRGline ID:{wildcards.sample} LB:1 PL:illumina PU:1 SM:{wildcards.sample} \
        --runThreadN {threads} \
        --outSAMtype BAM SortedByCoordinate \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --outReadsUnmapped None \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --alignMatesGapMax 1000000 \
        --alignIntronMax 1000000 \
        --outFilterType Normal \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 1 \
        --outSJfilterReads Unique \
        --outFilterMultimapNmax 10 \
        --outBAMcompression 10 \
        --sjdbOverhang 100 \
        --sjdbGTFfile {GTF} 2> {log}"

rule STAR_denovo:
    input: r1=R1, r2=R2, index=INDEX
    output: "out/{sample}.Aligned.sortedByCoord.out.bam"
    benchmark: "out/benchmarks/{sample}.align.json"
    log: "out/logs/{sample}.align.txt"
    conda: "envs/myenv.yaml"
    threads: 12
    params: R="'span[hosts=1] rusage[mem=20]'", J="align", o="out/logs/align.out", eo="out/logs/align.err"
    shell: "STAR \
        --genomeDir {input.index} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix out/{wildcards.sample}. \
        --outSAMattributes NH HI XS \
        --outSAMattrRGline ID:{wildcards.sample} LB:1 PL:illumina PU:1 SM:{wildcards.sample} \
        --runThreadN {threads} \
        --outSAMtype BAM SortedByCoordinate \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --outReadsUnmapped None \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --alignMatesGapMax 1000000 \
        --alignIntronMax 1000000 \
        --outFilterType Normal \
        --outBAMcompression 10 \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 1 \
        --outSJfilterReads Unique \
        --outFilterMultimapNmax 10 \
        --sjdbOverhang 100 2> {log}"


#snakemake --snakefile Snakefile-K0562 --cluster "bsub -J {params.J} -n {params.n} -R {params.R} -W 4:00 -o {params.o} -eo {params.eo}" --jn {rulename}.{jobid}.sj -j 50 -k --latency-wait 60 --use-conda --use-singularity --singularity-args "--bind /data:/data,/lila:/lila" --ri -n
# --rulegraph | dot -Tpdf > dag.pdf
