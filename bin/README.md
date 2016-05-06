agg.pl
This utility aggregates the output of sjcount (STDIN) by the 5th column (offset) and outputs a TSV (STDOUT) with three extra columns being (5) total count, (6) staggered read count, (7) entropy
	-logfile ..., name of the log file
	-margin ..., the margin for offset, default=0
	-maxintron ..., max intron length; only applied to splits of degree 1, default=0
	-minintron ..., min intron length; only applied to splits of degree 1, default=4
	-prefix ..., prefix in the chromosome name
	-readLength ..., the read length, default=0
=====================================================
annotate.pl
This utility takes an aggregated TSV file, the genomic annotation, and the genome, and outputs a TSV with two more columns: (8) annotation status and (9) splice sites
If the input was strandless ('.' in column 4) then each line will be reported twice, one for '+' and and for '-' strand
	-annot ..., the annotation (gtf), obligatory
	-dbx ..., the genome (dbx), obligatory
	-deltaSS ..., distance threshold for splice sites, default=0
	-idx ..., the genome (idx), obligatory
	-in ..., the input tsv file, obligatory
	-logfile ..., name of the log file
=====================================================
choose_strand.pl
This utility takes a TSV4+3+2 file (STDIN) where column (8) is the annotation status and column (9) is splice sites, and selects for each pair of beg/end the strand based on these two columns (STDOUT)
	-annot ..., annotation column, default=5
	-logfile ..., name of the log file
	-sites ..., splice site column, default=6
=====================================================
constrain_ssc.pl
This utility takes a BED6 ssc file (STDIN) and constraints its content to splice sites which are present in BED6 ssj file (STDOUT)
	-ssj ..., input ssj with constraints, obligatory
=====================================================
idr4sj.pl
Generic npIDR routine for more than two bioreplicas
Died at Perl/idr4sj.pl line 10.
=====================================================
make.pl
This utility creates a makefile for the sj pipeline, taking the index file from STDIN and printing the makefile to STDOUT
	-annot ..., the annotation (gtf), obligatory
	-deltaSS ..., distance threshold for splice sites, default=10
	-dir ..., the output directory, obligatory
	-entropy ..., entropy lower threshold, default=1.5
	-genome ..., the genome (without .dbx or .idx), obligatory
	-group ..., the grouping field for IDR, default=labExpId
	-idr ..., IDR upper threshold, default=0.1
	-in ..., the index file, obligatory
	-margin ..., margin for aggregate, default=5
	-merge ..., the name of the output to merge in case if blocks are missing, default=all
	-mincount ..., min number of counts for the denominator, default=10
	-param ..., parameters passed to sjcount
	-repository ..., the repository subdirectory for bam files
	-smpid ..., sample id field, default=labExpId
	-status ..., annotation status lower threshold, default=0
=====================================================
make_large.pl
This utility creates a makefile for the sj pipeline, taking the index file from STDIN and printing the makefile to STDOUT
	-annot ..., the annotation (gtf), obligatory
	-deltaSS ..., distance threshold for splice sites, default=10
	-dir ..., the output directory, obligatory
	-entropy ..., entropy lower threshold, default=1.5
	-genome ..., the genome (without .dbx or .idx), obligatory
	-group ..., the grouping field for IDR, default=labExpId
	-idr ..., IDR upper threshold, default=0.1
	-in ..., the index file, obligatory
	-margin ..., margin for aggregate, default=5
	-merge ..., the name of the output to merge in case if blocks are missing, default=all
	-mincount ..., min number of counts for the denominator, default=10
	-other ..., subset of exons and introns
	-param ..., parameters passed to sjcount
	-repository ..., the repository subdirectory for bam files
	-smpid ..., sample id field, default=labExpId
	-status ..., annotation status lower threshold, default=0
=====================================================
merge_gff.pl
This utility merges gtf/gff files into rectangular tables with colnames
Example: Perl/merge_gff.pl -i <gtf_1> <label_1> -i <gtf_2> <label_2> ... -o feature_1 <tsv1> -o feature_2 <tsv2>
	-i ..., input gtf file name and label, array=hash
	-o ..., feature and output tsv file name, array=hash
	-percent ..., do not report rows that have more than <value> NAs, default=0.25
	-transf ..., transformation (log, logit)
=====================================================
merge_large_gff.pl
	-dir ..., directory name, obligatory
	-ext ..., file extention, obligatory
	-f ..., feature
	-group ..., sample id field, default=labExpId
	-out ..., feature and output tsv file name, array=hash
	-percent ..., do not report rows that have more than <value> NAs, default=0.25
	-subset ..., file to subset id on
=====================================================
merge_large_tsv.pl
	-by ..., comma separated list of key columns, default=1
	-dir ..., directory name, obligatory
	-ext ..., file extention, obligatory
	-group ..., sample id field, default=labExpId
	-sep ..., separator, default=_
	-subset ..., file to subset id on
	-val ..., column with values, default=2
=====================================================
merge_tsv.pl
Merger script for tsv files with SJ counts
	-by ..., comma separated list of key columns, default=1
	-i ..., input gtf file name and label, array=hash
	-sep ..., separator, default=_
	-val ..., column with values, default=2
=====================================================
mk_stat.pl
	-i ..., input gtf file name and label, array=hash
=====================================================
mk_stat_large.pl
	-dir ..., directory name, obligatory
	-ext ..., file extention, obligatory
	-group ..., sample id field, default=labExpId
=====================================================
psicas.pl
This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs
	-annot ..., the annotation (GTF) file, obligatory
	-mincount ..., the min value of the denominator of the fraction, default=10
	-ssj ..., the input ssj  (tsv) file, obligatory
=====================================================
transcript_elements.pl
This utility takes a gtf annotation (STDIN) and reformats it into a more compact, quickly readable file (STDOUT). Only exons are taken into account.
	-features ..., list of features to output, default=transcript_id
	-lim ..., line limit for debugging, default=0
	-source ..., the content of the source field, default=SJPIPE
=====================================================
tsv2bed.pl
This utility transforms tsv to bed (tsv is supplied at STDIN)
	-extra ..., comma-separated list of column numbers to make BED6+n
	-minus ..., base to subtract from the coordinates, default=1
	-name, use junction id as name
	-ssc, if exon-intron junction
=====================================================
tsv2gff.pl
This utility transforms tsv to bed (tsv is supplied at STDIN)
	-minus ..., base to subtract from the coordinates, default=0
	-o ..., hash, attr name >> column number}, array=hash
	-source ..., the source field, default=IPSA
=====================================================
zeta.pl
This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs
	-annot ..., the annotation (GTF) file
	-exons ..., exons file with exons to quantify
	-mincount ..., the min value of the denominator of the fraction, default=10
	-ssc ..., the input ssc (BED) file
	-ssj ..., the input ssj (BED) file, obligatory
	-stranded ..., 1(yes) or 0(no), default=1
=====================================================

