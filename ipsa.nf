/*
 * Copyright (c) 2014-2017, Centre for Genomic Regulation (CRG)
 *
 * Copyright (c) 2014-2017, Dmitri Pervouchine
 *
 * Copyright (c) 2016-2017, Emilio Palumbo
 * 
 * This file is part of 'IPSA-nf': 
 * Integrative Pipeline for Splicing Analyses (IPSA) in Nextflow
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// library imports
import IPSA

// parameters for test run
params.index = "${baseDir}/data/index.tsv"
params.genome = "${baseDir}/data/genome.fa"
params.annot = "${baseDir}/data/annotation.gtf"


// parameters
params.merge = "all"
params.dir = 'data'
params.sjcountParams = ''
params.margin = 5
params.mincount = 10
params.deltaSS = 10
params.entropy = 1.5
params.status = 0
params.microexons = false
params.help = false

//print usage
if (params.help) {
  log.info ''
  log.info 'I P S A ~ Integrative Pipeline for Splicing Analyses'
  log.info '----------------------------------------------------'
  log.info 'Run IPSA on a set of data.'
  log.info ''
  log.info 'Usage: '
  log.info '    ipsa-nf [options]'
  log.info ''
  log.info 'Options:'
  log.info '--index INDEX_FILE        the index file in TSV format'
  log.info '--genome GENOME_FILE      the genome file in FASTA format'
  log.info '--annot ANNOTATION_FILE   the annotation file in gtf format'
  log.info '--merge MERGE             prefix for merged output files (default: all)'
  log.info '--dir DIRECTORY           the output directory'
  log.info '--sjcount-params PARAMS   additional `sjcount` paramters'
  log.info '--margin MARGIN           margin for aggregate (default: 5)'
  log.info '--mincount MIN_COUNT      minimum number of counts for denominators when calculationg fractions (default: 10)'
  log.info '--deltaSS DELTA           distance threshold for splice sites (default: 10)'
  log.info '--entropy ENTROPY         entropy lower threshold (default: 1.5)'
  log.info '--status STATUS           annotation status lower threshold (default: 0)'
  log.info '--microexons              include microexons, default=false'
  exit 1
}

// check mandatory options
if (!params.genome) {
    exit 1, "Reference genome not specified"
}

if (!params.annot) {
    exit 1, "Annotation not specified"
}

log.info ""
log.info "I P S A ~ Integrative Pipeline for Splicing Analyses"
log.info ""
log.info "General parameters"
log.info "------------------"
log.info "Index file                         : ${params.index}"
log.info "Genome                             : ${params.genome}"
log.info "Annotation                         : ${params.annot}"
log.info "Merge output name                  : ${params.merge}"
log.info "Output dir                         : ${params.dir}"
log.info "Sjcount parameters                 : ${params.sjcountParams ?: '-'}"
log.info "Include microexons                 : ${params.microexons}"
log.info ""
log.info "Thresholds and limits"
log.info "---------------------"
log.info "Margin for aggregate               : ${params.margin}"
log.info "Minimum counts for denominator     : ${params.mincount}"
log.info "Splice sites distance threshold    : ${params.deltaSS}"
log.info "Entropy lowewr threshold           : ${params.entropy}"
log.info "Annotation status lower threshold  : ${params.status}"
log.info ""

if (params.genome =~ /.fa$/) {
  process genomeIndex {
    input:
    file genome from file(params.genome)

    output:
    set file("${prefix}.dbx"), file("${prefix}.idx") into genomeIdx
    
    script:
    prefix = genome.name.replace(/.fa/, '')
    """
    transf -dir ./${genome} -dbx ${prefix}.dbx -idx ${prefix}.idx -exactdir
    """
  }
} else {
  Channel.value([file("${params.genome}.dbx"), file("${params.genome}.idx")])
    .set{ genomeIdx }
}

if (params.annot =~ /.g[tf]f$/) {
  process txElements {
    input:
    file annotation from file(params.annot)

    output:
    file "${prefix}.gfx" into txIdxAnnotate, txIdxZeta, txIdxZetaMex, txIdxPsicas

    script:
    prefix = annotation.name.replace(/.gtf/,'')
    """
    transcript_elements.pl - < ${annotation} | sort -k1,1 -k4,5n > ${prefix}.gfx
    """
  }
} else {
  Channel.value(file("${params.annot}"))
    .into { txIdxAnnotate; txIdxZeta; txIdxZetaMex; txIdxPsicas }
}

Channel
  .from(IPSA.parseIndexFile(file(params.index)))
  .set { bams } 

process preprocBams {
  input:
  set id, file(bam), readType, readStrand from bams

  output:
  set id, file(bam), readType, readStrand, stdout into bamsWreadLength

  script:
  prefix = bam.name.replace(/.bam/,'')
  """
  samtools view -F4 ${bam} | head -1 | awk '\$0=length(\$10)' | tr -d '\\n'
  """
}

process sjcount {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(bam), readType, readStrand, readLength from bamsWreadLength

  output:
  set id, file("${prefix}.ssc.tsv"), readLength into A01ssc
  set id, file("${prefix}.ssj.tsv"), readLength into A01ssj
  set id, file("${prefix}.ssj.tsv") into A01mex

  script:
  endpoint = 'A01'
  prefix = "${id}.${endpoint}"
  def strandParams
  switch (readStrand) {
    case'SENSE':
      strandParams = '-read1 0'
      break
    case'MATE1_SENSE':
      strandParams = '-read1 0 -read2 1'
      break
    case 'NONE':
      strandParams = '-unstranded'
      break
    default:
      strandParams = ''
      break
  }
  """
  sjcount -bam ${bam} \
          -ssc ${prefix}.ssc.tsv \
          -ssj ${prefix}.ssj.tsv \
          -nbins ${readLength} \
          ${strandParams} \
          ${params.sjcountParams ?: ''} \
          -quiet
  """
}

process aggregateSSC {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(tsv), readLength from A01ssc

  output:
  set id, file("${prefix}.tsv") into sscA02

  script:
  degree = 0
  endpoint = 'A02'
  prefix = "${id}.${endpoint}.ssc"
  """
  aggregate.awk -v degree=${degree} -v readLength=${readLength} -v margin=${params.margin} -v prefix= -v logfile=${prefix}.log ${tsv} > ${prefix}.tsv
  """
}

process aggregateSSJ {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(tsv), readLength from A01ssj

  output:
  set id, file("${prefix}.tsv") into ssjA02

  script:
  degree = 1
  endpoint = 'A02'
  prefix = "${id}.${endpoint}.ssj"
  """
  aggregate.awk -v degree=${degree} -v readLength=${readLength} -v margin=${params.margin} -v prefix= -v logfile=${prefix}.log ${tsv} > ${prefix}.tsv
  """
}

process aggregateMex {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(tsv) from A01mex

  when:
  params.microexons

  output:
  set id, file("${prefix}.tsv") into D01

  script:
  degree = 2
  endpoint = 'D01'
  prefix = tsv.name.replace(/.tsv/,'').replace(/A01.ssj/,'D01')
  """
  aggregate.awk -v degree=${degree} -v logfile=${prefix}.log ${tsv} > ${prefix}.tsv
  """
}

process annotate {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set file(genomeDBX), file(genomeIDX) from genomeIdx
  file annotation from txIdxAnnotate
  set id, file(ssj) from ssjA02

  output:
  set id, file("${prefix}.tsv") into A03

  script:
  endpoint = 'A03'
  prefix = "${id}.${endpoint}.ssj"
  """
  annotate.pl -annot ${annotation} -dbx ${genomeDBX} -idx ${genomeIDX} -deltaSS ${params.deltaSS} -in ${ssj} > ${prefix}.tsv
  """
}

process chooseStrand {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssj) from A03

  output:
  set id, file("${prefix}.tsv") into ssjA04, ssj4constrain, ssj4constrainMult

  script:
  endpoint = 'A04'
  prefix = "${id}.${endpoint}.ssj"
  """
  choose_strand.awk ${ssj} > ${prefix}.tsv
  """
}

ssj4constrain.combine(sscA02, by: 0)
  .map { 
    [it[0]] + it[1..-1].sort { it.baseName }
  }.set{ constrain }

if ( params.microexons ) {
  ssj4constrainMult.combine(D01, by: 0)
    .map { 
      [it[0]] + it[1..-1].sort { it.baseName }
    }.set { constrainMult }
} else {
  Channel.empty()
    .set { constrainMult }
}

process constrainSSC {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssc), file(ssj) from constrain

  output:
  set id, file("${prefix}.tsv") into sscA04

  script:
  endpoint = 'A04'
  prefix = "${id}.${endpoint}.ssc"
  """
  constrain_ssc.awk -v jncfile=${ssj} ${ssc} > ${prefix}.tsv
  """
}

process constrainMex {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssj), file(ssjMex) from constrainMult

  output:
  set id, file("${prefix}.tsv") into D02

  script:
  endpoint = 'D02'
  prefix = "${id}.${endpoint}"
  """
  constrain_mex.awk -v jncfile=${ssj} ${ssjMex} > ${prefix}.tsv
  """
}

process extractMex {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssjMex) from D02

  output:
  set id, file("${prefix}.tsv") into D06

  script:
  endpoint = 'D06'
  prefix = "${id}.${endpoint}"
  """
  extract_mex.awk ${ssjMex} > ${prefix}.tsv 
  """
}

process sscA06 {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssc) from sscA04

  output:
  file "${prefix}.tsv" into sscA06
  set id, file("${prefix}.tsv") into ssc4merge, ssc4allA06

  script:
  endpoint = 'A06'
  prefix = "${id}.${endpoint}.ssc"
  """
  awk '\$4>=${params.entropy}' ${ssc} > ${prefix}.tsv
  """
}

process ssjA06 {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssj) from ssjA04

  output:
  file "${prefix}.tsv" into ssjA06, ssj4gffA06
  set id, file("${prefix}.tsv") into ssj4merge, ssj4allA06, ssj4psicasA06

  script:
  endpoint = 'A06'
  prefix = "${id}.${endpoint}.ssj"
  """
  awk '\$4>=${params.entropy} && \$5>=${params.status}' ${ssj} > ${prefix}.tsv
  """
}

ssj4merge.toSortedList { a,b -> a[0] <=> b[0] }
  .map { list ->
    ids = []
    ssjs = []
    list.each { ids << it[0]; ssjs << it[1] }
    [ids, ssjs]
  }.set { ssj4mergeSplit }

ssc4merge.toSortedList { a,b -> a[0] <=> b[0] }
  .map { list ->
    ids = []
    sscs = []
    list.each { ids << it[0]; sscs << it[1] }
    [ids, sscs]
  }.set { ssc4mergeSplit }

if ( params.microexons ) {
  ssj4allA06.combine(ssc4allA06, by: 0).combine(D06, by: 0)
    .map {
      [it[0]] + it[1..-1].sort { it.baseName }
    }.set { allMex }
  Channel.empty().set { allA06 }
} else {
  ssj4allA06.combine(ssc4allA06, by: 0)
    .map {
      [it[0]] + it[1..-1].sort { it.baseName }
    }.set { allA06 }
  Channel.empty().set { allMex }
}

process mergeTsvSSJ {
  publishDir "${params.dir}"
  
  input:
  set ids, file(ssjs) from ssj4mergeSplit

  output:
  file "${prefix}.tsv"

  shell:
  by = 1
  val = 2
  sep = '_'
  input = [ssjs.toList(), ids].transpose().flatten().collect { "'$it'" }.join(',')
  prefix = "${params.merge}.counts.ssj"
  template 'merge_tsv.pl'
}

process mergeTsvSSC {
  publishDir "${params.dir}"
  
  input:
  set ids, file(sscs) from ssc4mergeSplit

  output:
  file "${prefix}.tsv"

  shell:
  by = 1
  val = 2
  sep = '_'
  input = [sscs.toList(), ids].transpose().flatten().collect { "'$it'" }.join(',')
  prefix = "${params.merge}.counts.ssc"
  template 'merge_tsv.pl'
}

process zeta {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  file annotation from txIdxZeta
  set id, file(ssc), file(ssj) from allA06

  output:
  set id, file("${prefix}.gff") into A07

  script:
  endpoint = 'A07'
  prefix = ssj.name.replace(/.tsv/,'').replace(/A06.ssj/, endpoint)
  """
  zeta.pl -annot ${annotation} -ssc ${ssc} -ssj ${ssj} -mincount ${params.mincount} > ${prefix}.gff 
  """
}

process psicas {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  file annotation from txIdxPsicas
  set id, file(ssj) from ssj4psicasA06

  output:
  set id, file("${prefix}.gff") into B07

  script:
  endpoint = 'B07'
  prefix = ssj.name.replace(/.tsv/,'').replace(/A06.ssj/, endpoint)
  """
  psicas.pl -ssj ${ssj} -annot ${annotation} -mincount ${params.mincount} > ${prefix}.gff 
  """
}

process zetaMex {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  file annotation from txIdxZetaMex
  set id, file(ssc), file(ssj), file(exons) from allMex

  output:
  set id, file("${prefix}.gff") into A07mex

  script:
  endpoint = 'A07'
  prefix = ssj.name.replace(/.tsv/,'').replace(/A06.ssj/, endpoint)
  """
  zeta.pl -annot ${annotation} -ssc ${ssc} -ssj ${ssj} -exons ${exons} -mincount ${params.mincount} > ${prefix}.gff 
  """
}

if ( params.microexons ) {
  A07mex.toSortedList { a,b -> a[0] <=> b[0] }
    .map { list ->
      ids = []
      gffs = []
      list.each { ids << it[0]; gffs << it[1] }
      [ids, gffs]
    }.set { A074merge }
} else {
  A07.toSortedList { a,b -> a[0] <=> b[0] }
    .map { list ->
      ids = []
      gffs = []
      list.each { ids << it[0]; gffs << it[1] }
      [ids, gffs]
    }.set { A074merge }
}

B07.toSortedList { a,b -> a[0] <=> b[0] }
  .map { list ->
    ids = []
    gffs = []
    list.each { ids << it[0]; gffs << it[1] }
    [ids, gffs]
  }.set { B074merge }

process mergeGFFzeta {
  publishDir "${params.dir}"
  
  input:
  set ids, file(sscs) from A074merge

  output:
  file "${prefix}.psi.tsv"
  file "${prefix}.psi3.tsv"
  file "${prefix}.psi5.tsv"
  file "${prefix}.psit.tsv"
  file "${prefix}.cosi.tsv"
  file "${prefix}.cosi3.tsv"
  file "${prefix}.cosi5.tsv"
  file "${prefix}.cosit.tsv"
  file "${prefix}.exc.tsv"
  file "${prefix}.inc.tsv"
  file "${prefix}.ret.tsv"

  shell:
  prefix = "${params.merge}.A"
  input = [sscs.toList(), ids].transpose().flatten().collect { "'$it'" }.join(',')
  features = ['cosi', 'cosi3', 'cosi5', 'cosit', 'psi', 'psi3', 'psi5', 'psit', 'exc', 'inc', 'ret']
  output = features.collect { "'${it}', '${prefix}.${it}.tsv'" }.join(',')
  percent = 0.25
  transf = 'log'
  template 'merge_gff.pl'
}

process mergeGFFpsicas {
  publishDir "${params.dir}"
  
  input:
  set ids, file(sscs) from B074merge

  output:
  file "${prefix}.psicas.tsv"
  file "${prefix}.psiloc.tsv"

  shell:
  prefix = "${params.merge}.B"
  input = [sscs.toList(), ids].transpose().flatten().collect { "'$it'" }.join(',')
  features = ['psicas', 'psiloc']
  output = features.collect { "'${it}', '${prefix}.${it}.tsv'" }.join(',')
  percent = 0.25
  transf = 'log'
  template 'merge_gff.pl'
}

process ssjTsv2bed {

  publishDir "${params.dir}/${endpoint}"

  input:
  file ssj from ssjA06

  output:
  file "${prefix}.bed" into E06

  script:
  endpoint = 'E06'
  prefix = ssj.name.replace(/.tsv/,'').replace(/A06/, endpoint)
  """
  tsv2bed.pl < ${ssj} -extra 2,3,4,5,6,7 > ${prefix}.bed
  """
}

process sscTsv2bed {

  publishDir "${params.dir}/${endpoint}"

  input:
  file ssc from sscA06

  output:
  file "${prefix}.bed" into E06ssc

  script:
  endpoint = 'E06'
  prefix = ssc.name.replace(/.tsv/,'').replace(/A06/, endpoint)
  """
  tsv2bed.pl < ${ssc} -extra 2 -ssc > ${prefix}.bed
  """
}

process tsv2gff {

  publishDir "${params.dir}/${endpoint}"

  input:
  file ssj from ssj4gffA06

  output:
  file "${prefix}.gff" into E06ssj

  script:
  endpoint = 'E06'
  prefix = ssj.name.replace(/.tsv/,'').replace(/A06/, endpoint)
  """
  tsv2gff.pl < ${ssj} -o count 2 -o stagg 3 -o entr 4 -o annot 5 -o nucl 6 -o IDR 7 > ${prefix}.gff
  """
}

workflow.onComplete {
    log.info """
Pipeline execution summary
---------------------------
Completed at: ${workflow.complete}
Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
exit status : ${workflow.exitStatus}
Error report: ${workflow.errorReport ?: '-'}
"""
}
