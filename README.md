# ProteomeGenerator

## Introduction

ProteomeGenerator is an open, modular, and scalable framework for reference guided and de novo proteogenomic database generation written in the Snakemake workflow management system. The workflow consists of a base file that defines the project details and input samples, and pgm, a modular set of rules sourced as an include.


## Installation

### Dependencies

ProteomeGenerator depends on multiple free and open source tools. The following software and all dependencies must be installed:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [STAR](https://github.com/alexdobin/STAR)
* [Picard](http://broadinstitute.github.io/picard/)
* [samtools](http://samtools.sourceforge.net)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)
* [gffread](https://github.com/gpertea/gffread)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [R](https://conda.io/docs/user-guide/tasks/use-r-with-conda.html)

The Bioconductor suite as well as the Biostrings R package are also required.

* [Bioconductor](https://bioconductor.org/)
* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

Optionally, many of the outputs from ProteomeGenerator may be viewed using the [Integrative Genomics Viewer](http://software.broadinstitute.org/software/igv/).

### Installing ProteomeGenerator

ProteomeGenerator is installed to your preferred computing environment using git clone:

```bash
git clone https://github.com/jtpoirier/proteomegenerator
```


This command downloads Snakefile-K0562, the configuraton file used in Cifani et al., as well as the required includes and accessory scripts to run ProteomeGenerator on your data set.

## Configuration

ProteomeGenerator is configured by setting a series of variables to the explicit location of each of the required tools above in order to avoid conflicting versions that may be available in the path. To run ProteomeGenerator on our own data, create a new Snakefile using Snakefile-K0562 as a template. Next, you will need to edit the "User Variables" section to reflect your computing environment. This section is broken into four parts: "Directories", "References", "Dependencies", and "Samples".

### Directories

The "WD" variable sets the path that will serve as the relative path for all analyses. In a cluster environment, this path should be accessible to all nodes on the network. The "TMP" variable should point to a local scratch storage directory that should be mounted locally to each node on a cluster, as opposed to a high performance network storage location. This folder is not required to be accessible to all nodes.

### References

ProteomeGenerator requires a FASTA formated genome reference and accepts a GTF transcript model with coordinates based on the same reference. STAR and BLAST indexes are automatically generated by ProteomeGenerator based on the provided references.

### Dependencies

Variables pointing to the explicit paths of the above listed dependies are located in this section and must be edited by the user to reflect the intended computing environment.

### Samples

The "Samples" section should be edited to reflect the type, number, and naming structure of the input FASTQ files for a given experiment. In addition, the user must edit the input for the rule STAR to reflect the specific naming convention used.

## Running ProteomeGenerator

ProteomeGenerator can be run locally or in a variety of cluster environments. The below example demonstrates how to run ProteomeGenerator on an LSF cluster head node from within screen in case the connection to the head node is lost.

```bash
screen -S pg
snakemake --snakefile Snakefile-K0562 --cluster \
"bsub -J {params.J} -n {params.n} -R {params.R} -W 4:00 -o {params.o} -eo {params.eo}" \
--jn {rulename}.{jobid}.sj -j 50 -k --latency-wait 60 --ri
```

## Expected Output

ProteomeGenerator will generate an indexed bam filed of mapped and filtered reads of the format {sample}.Aligned.trimmed.out.bam, a sample-specific GTF of the format {sample}-stringtie.gtf, and a proteogenomic database called proteome.unique.fasta. A GFF3 corresponding to each entry in the fasta database is also generated with the predicted spliced peptide sequences mapped onto genome space for easy viewing in the [Integrative Genomics Viewer](http://software.broadinstitute.org/software/igv/).