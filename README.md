# IPSA-nf

[![CircleCI status](https://circleci.com/gh/guigolab/ipsa-nf.svg?style=shield)](https://circleci.com/gh/guigolab/ipsa-nf/tree/master)

An Integrative Pipeline for Splicing Analyses (IPSA) written in the Nextflow DSL.

The pipeline performs the following splicing analyses:

* Quantification of splice junctions and splice boundaries
* Calculation of splicing indices, exon- and intron-centric

## Quickstart

Install nextflow with the following command:
```
curl -fsSL get.nextflow.io | bash
```

Pull the docker image:
```
docker pull guigolab/ipsa-nf@sha256:88e680da318023d2a577893d5c4f0324ad720f83b13830b4e29f2d03f77490bb
```

Launch the test pipeline with the following command:
```
./nextflow run guigolab/ipsa-nf
```