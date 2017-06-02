# IPSA-nf

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-blue.svg)](http://nextflow.io)
[![CircleCI](https://circleci.com/gh/guigolab/ipsa-nf/tree/master.svg?style=shield)](https://circleci.com/gh/guigolab/ipsa-nf/tree/master)

An Integrative Pipeline for Splicing Analyses (IPSA) written in the Nextflow DSL.

The pipeline allows to perform the following splicing analyses:

* Quantification of splice junctions and splice sites
* Calculation of exon-centric and intron-centric splicing metrics
* Identification of micro-exons

## Quickstart

Install nextflow with the following command:
```
curl -fsSL get.nextflow.io | bash
```

Pull the docker image:
```
docker pull guigolab/ipsa-nf@sha256:29750072f2b42ee8ea094a331f53bc5906183591a71ed26619ea76b12b6be3ed
```

Launch the test pipeline with the following command:
```
./nextflow run guigolab/ipsa-nf
```

## Pipeline usage

Launching the pipeline with the `--help ` parameter shows the help message:

```
nextflow run ipsa-nf --help
```

```
N E X T F L O W  ~  version 0.24.4
Launching `guigolab/ipsa-nf` [dreamy_aryabhata] - revision: v4.0

I P S A ~ Integrative Pipeline for Splicing Analyses
----------------------------------------------------
Run IPSA on a set of data.

Usage: 
    ipsa-nf [options]

Options:
--index INDEX_FILE        the index file in TSV format
--genome GENOME_FILE      the genome file in FASTA format
--annot ANNOTATION_FILE   the annotation file in gtf format
--merge MERGE             prefix for merged output files (default: all)
--dir DIRECTORY           the output directory
--sjcount-params PARAMS   additional `sjcount` paramters
--margin MARGIN           margin for aggregate (default: 5)
--mincount MIN_COUNT      minimum number of counts for denominators when calculationg fractions (default: 10)
--deltaSS DELTA           distance threshold for splice sites (default: 10)
--entropy ENTROPY         entropy lower threshold (default: 1.5)
--status STATUS           annotation status lower threshold (default: 0)
--microexons              include microexons, default=false
```

## Input format

IPSA-nf takes as input a tab separated file containing information about the input data. The file must be specified using the `--index` parameter. The format of the index file is as follows:

1. sample identifier
2. path to the bam file to be processed
3. library type:
	* `Single-End`
	* `Paired-End`
4. strandnesss of the data:
	* `NONE` for unstranded
	* `SENSE`/`ANTISENSE` for `Single-End` stranded
	* `MATE1_SENSE`/`MATE2_SENSE` for `Paired-End` stranded

Here is an index file example:

```
E14_rep1	E14AlnRep1.sub.bam	Paired-End	MATE2_SENSE
E14_rep2	E14AlnRep2.sub.bam	Paired-End	MATE2_SENSE
E18_rep1	E18AlnRep1.sub.bam	Paired-End	MATE2_SENSE
E18_rep2	E18AlnRep2.sub.bam	Paired-End	MATE2_SENSE
```

## Pipeline results

Analyses results are saved into the folder specified with the `--dir` parameter. By default it is the `data` directory within the current working folder.

Output files are organinzed into folders corresponding to the different pipeline endpoints:

* `A01` - splice junctions and splice sites counts
* `A02` - aggregated splice junctions and splice sites counts
* `A03` - aggregated junctions with annotation status and splice site nucleotides
* `A04` - aggregated junctions with selected strand and constrained splice sites
* `A06` - aggregated counts filtered by annotation status and entropy
* `A07` - splicing indices
* `E06` - BED files with splicing indices from A06

And if `--microexons` is used:

* `D01` - aggregated 2-splits
* `D02` - aggregated constrained 2-splits
* `D06` - extracted microexons from constrained 2-splits

## Requirements

IPSA-nf is configured to run using the [Docker](https://www.docker.com/) container engine by default. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details. 

In order to run the pipeline with Doecker the following dependencies have to be met:

* Java 7/8
* [Nextflow](https://www.nextflow.io) 0.24.x (or higher)
* [Docker](https://www.docker.com/) 1.10 (or higher)

The pipeline can also be used without Docker by installing the following software components on your system (natively or by using [Environemnt Modules](http://modules.sourceforge.net/)):

* [maptools](https://github.com/pervouchine/maptools)
* [sjcount](https://github.com/pervouchine/sjcount)
* [samtools](https://github.com/samtools/samtools)
* [Perl](http://www.perl.org/) standard modules
