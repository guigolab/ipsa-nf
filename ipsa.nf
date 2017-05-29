// library imports
import TsvIndexFile

// parameters for test run
params.index = "${baseDir}/data/index.tsv"
params.genome = "${baseDir}/data/genome.fa"
params.annot = "${baseDir}/data/annotation.gtf"


// parameters
params.deltaSS = 10
params.dir = 'data'
params.entropy = 1.5
params.group = 'labExpId'
params.idr = 0.1
params.margin = 5
params.merge = "all"
params.mincount = 10
params.smpid = 'labExpId'
params.status = 0
params.microexons = false
params.help = false
params.param = ''
params.repository = ''

//print usage
if (params.help) {
  log.info ''
  log.info 'I P S A ~ Integrative Pipeline for Splicing Analyses'
  log.info '----------------------------------------------------'
  log.info 'Run IPSA on a set of data.'
  log.info ''
  log.info 'Usage: '
  log.info '    ipsa.nf [options]'
  log.info ''
  log.info 'Options:'
  log.info '--index INDEX_FILE        the index file'
  log.info '--genome GENOME_FILE      the genome file (FASTA)'
  log.info '--annot ANNOTATION_FILE   the annotation file (gtf)'
  log.info '--deltaSS DELTA           distance threshold for splice sites, default=10'
  log.info '--dir DIRECTORY           the output directory, obligatory'
  log.info '--entropy ENTROPY         entropy lower threshold, default=1.5'
  log.info '--group GROUP             the grouping field for IDR, default=labExpId'
  log.info '--idr IDR                 IDR upper threshold, default=0.1'
  log.info '--margin MARGIN           margin for aggregate, default=5'
  log.info '--merge MERGE             the name of the output to merge in case if blocks are missing, default=all'
  log.info '--mincount MIN_COUNT      min number of counts for the denominator, default=10'
  log.info '--param PARAMS            parameters passed to sjcount'
  log.info '--repository REPOSITORY   the repository subdirectory for bam files'
  log.info '--smpid SAMPLE_ID_FIELD   sample id field, default=labExpId'
  log.info '--status STATUS           annotation status lower threshold, default=0'
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
log.info "Splice sites distance threshold    : ${params.deltaSS}"
log.info "Output dir                         : ${params.dir}"
log.info "Entropy lowewr threshold           : ${params.entropy}"
log.info "Grouping field                     : ${params.group}"
log.info "IDR upper threshold                : ${params.idr}"
log.info "Margin for aggregate               : ${params.margin}"
log.info "Merge output name                  : ${params.merge}"
log.info "Minimum counts for denominator     : ${params.mincount}"
log.info "Sjcount parameters                 : ${params.param ?: '-'}"
log.info "BAM files repository               : ${params.repository ?: '-'}"
log.info "Sample id field                    : ${params.smpid}"
log.info "Annotation status lower threshold  : ${params.status}"
log.info "Include microexons                 : ${params.microexons}"
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
  genomeIdx = Channel.value [file("${params.genome}.dbx"), file("${params.genome}.idx")]
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
  Channel.value file("${params.annot}")
    .into { txIdxAnnotate; txIdxZeta; txIdxZetaMex; txIdxPsicas }
}

Channel
  .from(TsvIndexFile.parse(file(params.index)))
  .into { bams } 

process preprocBams {
  input:
  set sample, id, file(bam), type, view, readType, readStrand from bams

  output:
  set sample, id, file(bam), type, view, readType, readStrand, stdout into bamsWreadLength

  script:
  prefix = bam.name.replace(/.bam/,'')
  """
  samtools view -F4 ${bam} | head -1 | awk '\$0=length(\$10)' | tr -d '\\n'
  """
}

process sjcount {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set sample, id, file(bam), type, view, readType, readStrand, readLength from bamsWreadLength

  output:
  set id, val(1), file("${prefix}.A01.ssj.tsv"), readLength into A01ssc
  set id, val(0), file("${prefix}.A01.ssc.tsv"), readLength into A01ssj
  set id, val(2), file("${prefix}.A01.ssj.tsv") into A01mex

  script:
  endpoint = 'A01'
  prefix = bam.name.replace(/.bam/,'')
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
          -ssc ${prefix}.A01.ssc.tsv \
          -ssj ${prefix}.A01.ssj.tsv \
          -nbins ${readLength} \
          ${strandParams} \
          ${params.param ?: ''} \
          -quiet
  """
}

process aggregate {
  
  publishDir "${params.dir}/${endpoint}"

  input:
  set id, splits, file(tsv), readLength from A01ssc.mix(A01ssj)

  output:
  set id, file("${prefix}.tsv") into A02

  script:
  endpoint = 'A02'
  prefix = tsv.name.replace(/.tsv/,'').replace(/A01/,'A02')
  """
  awk '\$2==${splits}' ${tsv} | agg.pl -readLength ${readLength} -margin ${params.margin} -logfile ${prefix}.log > ${prefix}.tsv 
  """
}

process aggregateMex {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, splits, file(tsv) from A01mex

  when:
  params.microexons

  output:
  set id, file("${prefix}.tsv") into D01

  script:
  endpoint = 'D01'
  prefix = tsv.name.replace(/.tsv/,'').replace(/A01.ssj/,'D01')
  """
  awk '\$2==${splits}' ${tsv} | agg.pl -logfile ${prefix}.log > ${prefix}.tsv 
  """
}

sscA02 = Channel.create()
ssjA02 = Channel.create()

A02.choice(sscA02, ssjA02) {
    f = it[1]
    f.name =~ /ssc/ ? 0 : 1
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
  prefix = ssj.name.replace(/.tsv/,'').replace(/A02/,'A03')
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
  prefix = ssj.name.replace(/.tsv/,'').replace(/A03/,'A04')
  """
  choose_strand.pl - < ${ssj} -logfile ${prefix}.log > ${prefix}.tsv
  """
}

constrain = ssj4constrain.combine(sscA02, by: 0)
.map { 
  [it[0]] + it[1..-1].sort { it.baseName }
}

if ( params.microexons ) {
  constrainMult = ssj4constrainMult.combine(D01, by: 0)
  .map { 
    [it[0]] + it[1..-1].sort { it.baseName }
  }
} else {
  constrainMult = Channel.empty()
}

process constrainSSC {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssc), file(ssj) from constrain

  output:
  set id, file("${prefix}.tsv") into sscA04

  script:
  endpoint = 'A04'
  prefix = ssc.name.replace(/.tsv/,'').replace(/A02/,'A04')
  """
  constrain_ssc.pl -ssj ${ssj} < ${ssc} > ${prefix}.tsv  
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
  prefix = ssjMex.name.replace(/.tsv/,'').replace(/D01/,'D02')
  """
  constrain_mult.pl -ssj ${ssj} < ${ssjMex} > ${prefix}.tsv
  """
}

process extractMex {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssjMex) from D02

  output:
  set id, file("${prefix}.tsv") into D03

  script:
  endpoint = 'D03'
  prefix = ssjMex.name.replace(/.tsv/,'').replace(/D02/,'D03')
  """
  extract_mex.pl < ${ssjMex} > ${prefix}.tsv
  """
}

process ssjIDR {

  publishDir "${params.dir}/${endpoint}"
  
  input:
  set id, file(tsv) from ssjA04

  output:
  set id, file("${prefix}.tsv") into ssjA05

  script:
  endpoint = 'A05'
  prefix = "${id}.${endpoint}.ssj"
  """
  idr4sj.pl ${tsv} > ${prefix}.tsv
  """
}

process sscIDR {

  publishDir "${params.dir}/${endpoint}"
  
  input:
  set id, file(tsv) from sscA04

  output:
  set id, file("${prefix}.tsv") into sscA05

  script:
  endpoint = 'A05'
  prefix = "${id}.${endpoint}.ssc"
  """
  idr4sj.pl ${tsv} > ${prefix}.tsv
  """
}

process idrMex {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(tsv) from D03

  output:
  set id, file("${prefix}.tsv") into D06

  script:
  endpoint = 'D06'
  prefix = "${id}.${endpoint}.mex"
  """
  idr4sj.pl ${tsv} > ${prefix}.tsv
  """
}

process ssjA06 {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssj) from ssjA05

  output:
  file "${prefix}.tsv" into ssjA06, ssj4gffA06
  set id, file("${prefix}.tsv") into ssj4merge, ssj4allA06, ssj4psicasA06

  script:
  endpoint = 'A06'
  prefix = ssj.name.replace(/.tsv/,'').replace(/A05/,'A06')
  """
  awk '\$4>=${params.entropy} && \$5>=${params.status} && \$7<${params.idr}' ${ssj} > ${prefix}.tsv
  """
}

process sscA06 {

  publishDir "${params.dir}/${endpoint}"

  input:
  set id, file(ssc) from sscA05

  output:
  file "${prefix}.tsv" into sscA06
  set id, file("${prefix}.tsv") into ssc4merge, ssc4allA06

  script:
  endpoint = 'A06'
  prefix = ssc.name.replace(/.tsv/,'').replace(/A05/,'A06')
  """
  awk '\$4>=${params.entropy} && \$7<${params.idr}' ${ssc} > ${prefix}.tsv
  """
}

ssj4merge.toSortedList { a,b -> a[0] <=> b[0] }
.map { list ->
  ids = []
  ssjs = []
  list.each { ids << it[0]; ssjs << it[1] }
  [ids, ssjs]
}.into {ssj4mergeSplit}

ssc4merge.toSortedList { a,b -> a[0] <=> b[0] }
.map { list ->
  ids = []
  sscs = []
  list.each { ids << it[0]; sscs << it[1] }
  [ids, sscs]
}.into {ssc4mergeSplit}

if ( params.microexons ) {
  ssj4allA06.combine(ssc4allA06, by: 0).combine(D06, by: 0).map {
    [it[0]] + it[1..-1].sort { it.baseName }
  }.into { allMex }
  Channel.empty().into { allA06 }
} else {
  ssj4allA06.combine(ssc4allA06, by: 0).map {
    [it[0]] + it[1..-1].sort { it.baseName }
  }.into { allA06 }
  Channel.empty().into { allMex }
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
  }.into {A074merge}
} else {
  A07.toSortedList { a,b -> a[0] <=> b[0] }
  .map { list ->
    ids = []
    gffs = []
    list.each { ids << it[0]; gffs << it[1] }
    [ids, gffs]
  }.into {A074merge}
}

B07.toSortedList { a,b -> a[0] <=> b[0] }
.map { list ->
  ids = []
  gffs = []
  list.each { ids << it[0]; gffs << it[1] }
  [ids, gffs]
}.into {B074merge}

process mergeGFFzeta {
  publishDir "${params.dir}"
  
  input:
  set ids, file(sscs) from A074merge

  output:
  file "${prefix}.psi.tsv"
  file "${prefix}.psi3.tsv"
  file "${prefix}.psi5.tsv"
  file "${prefix}.cosi.tsv"
  file "${prefix}.cosi3.tsv"
  file "${prefix}.cosi5.tsv"

  shell:
  prefix = "${params.merge}.A"
  input = [sscs.toList(), ids].transpose().flatten().collect { "'$it'" }.join(',')
  features = ['cosi', 'cosi3', 'cosi5', 'psi', 'psi3', 'psi5']
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

  shell:
  prefix = "${params.merge}.B"
  input = [sscs.toList(), ids].transpose().flatten().collect { "'$it'" }.join(',')
  features = ['psicas']
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
