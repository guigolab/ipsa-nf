#!/usr/bin/perl
use Perl::utils;

@suffixes = ('ssj','ssc');

if(@ARGV==0) {
    print STDERR "This utility creates a makefile for the sj pipeline, taking the index file from STDIN and printing the makefile to STDOUT\n";
}

parse_command_line(     in      => {description=>'the index file', ifabsent=>'index file not specified'},
                        dir     => {description=>'the output directory', ifabsent=>'output directory not specified'},
                        repository => {description=>'the repository subdirectory for bam files'},
                        param   => {description=>'parameters passed to sjcount'},
                        group   => {description=>'the grouping field for IDR', default=>'labExpId'},
                        smpid   => {description=>'sample id field', default=>'labExpId'},
                        margin  => {description=>'margin for aggregate', default=>5},
                        entropy => {description=>'entropy lower threshold', default=>1.5},
                        deltaSS => {description=>'distance threshold for splice sites', default=>10},
                        status  => {description=>'annotation status lower threshold', default=>0},
                        mincount=> {description=>'min number of counts for the denominator', default=>10},
                        idr     => {description=>'IDR upper threshold', default=>0.1},
                        annot   => {description=>'the annotation (gtf)', ifabsent=>'annotation not specified'},
                        genome  => {description=>'the genome (without .dbx or .idx)', ifabsent=>'genome not specified'},
                        merge   => {description=>'the name of the output to merge in case if blocks are missing', default=>"all"},
			other	=> {description=>'subset of exons and introns'},
                        SJCOUNTDIR =>{variable=>T,ifabsent=>'SJCOUNTDIR not specified'});

@group = split /\,/, $group;
@smpid = split /\,/, $smpid;

open FILE, $in || die("Can't read $in exiting\n");
while($line=<FILE>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);

    $grp = join("_", @attr{@group});	# grouping id for IDR
    $smp = join("_", @attr{@smpid});    # sample id

    next unless($grp && $smp);

    if($attr{'type'} eq "bam" && $attr{'view'}=~/^Alignments/) {

	if($attr{'readType'}=~/(\d+)(D*)$/) {
	    $readLength = $1;
	    $stranded   = ($2 eq "D" ? "" : "-unstranded");
	}
	else {
	    die("Read length not specified");
	}
	die("Incorrect read length") unless($readLength>0);

        @a = split /\//, $file;
        $target = pop(@a);

        if($file=~/^http/ || $file=~/^ftp/) {
            make(script=>"wget", before=>$file, output=>{-O=>"$repository$target"}, mkdir=>T, endpoint=>'download');
            $file = "$repository$target";
        }

        $target  =~ s/\.bam$//;
        $name = "$target";

	$PARAMS{'touch'}=>T;

	make(script=>$SJCOUNTDIR."sjcount", input=>{-bam=>$file}, output=>{-ssj=>fn($name,A01,ssj,tsv), -ssc=>fn($name,A01,ssc,tsv), -log=>fn($name,A01,ssj,'log')},
             after=>"-nbins $readLength $param $stranded -quiet", mkdir=>T, endpoint=>A01);

	$prm = "-readLength $readLength -margin $margin";
	make(script=>"awk '\$\$2==1'",input=>{''=>fn($name,A01,ssj,tsv)}, output=>{'>'=>fn($name,A02,ssj,tsv), -logfile=>fn($name,A02,ssj,'log')}, between=>"|perl Perl/agg.pl $prm", endpoint=>A02);
        make(script=>"awk '\$\$2==0'",input=>{''=>fn($name,A01,ssc,tsv)}, output=>{'>'=>fn($name,A02,ssc,tsv), -logfile=>fn($name,A02,ssc,'log')}, between=>"|perl Perl/agg.pl $prm", endpoint=>A02);

	make(script=>"offset.r", input=>{-t=>fn($name,A02,ssj,'log')}, output=>{-p=>fn($name,A02,ssj,'pdf')}, endpoint=>QC1);

	make(script=>"annotate.pl", input=>{-in=>fn($name,A02,ssj,tsv), -annot=>$annot, -dbx=>"$genome.dbx", -idx=>"$genome.idx"}, output=>{'>'=>fn($name,A03,ssj,tsv)}, after=>"-deltaSS $deltaSS", endpoint=>A03);
	make(script=>"choose_strand.pl", input=>{'<'=>fn($name,A03,ssj,tsv)}, output=>{'>'=>fn($name,A04,ssj,tsv), -logfile=>fn($name,A04,ssj,'log')}, before=>"-", endpoint=>A04);
	make(script=>"constrain_ssc.pl", input=>{'<'=>fn($name,A02,ssc,tsv),-ssj=>fn($name,A04,ssj,tsv)}, output=>{'>'=>fn($name,A04,ssc,tsv)}, endpoint=>A04);	

	make(script=>"awk", before=>"'\$\$4>=$entropy && \$\$5>=$status'", input=>{''=>fn($name,A04,ssj,tsv)}, output=>{'>'=>fn($name,F04,ssj,tsv)}, endpoint=>F04);

        push @{$IDR{$grp}{ssj}}, fn($name,A04,ssj,tsv);
	push @{$IDR{$grp}{ssc}}, fn($name,A04,ssc,tsv);
    }
}
close FILE;

foreach $grp(keys(%IDR)) {
    make(script=>"idr4sj.pl", input=>{''=>join(" ", @{$IDR{$grp}{ssj}})}, output=>{'>'=>fn($grp,A05,ssj,tsv)}, endpoint=>A05);
    make(script=>"idr4sj.pl", input=>{''=>join(" ", @{$IDR{$grp}{ssc}})}, output=>{'>'=>fn($grp,A05,ssc,tsv)}, endpoint=>A05);

    make(script=>"awk", before=>"'\$\$4>=$entropy && \$\$5>=$status && \$\$7<$idr'", input=>{''=>fn($grp,A05,ssj,tsv)}, output=>{'>'=>fn($grp,A06,ssj,tsv)}, endpoint=>A06);
    make(script=>"awk", before=>"'\$\$4>=$entropy && \$\$7<$idr'", input=>{''=>fn($grp,A05,ssc,tsv)}, output=>{'>'=>fn($grp,A06,ssc,tsv)}, endpoint=>A06);

    make(script=>'tsv2bed.pl', input=>{'<'=>fn($grp,A06,ssj,tsv)}, output=>{'>'=>fn($grp,E06,ssj,bed)}, between=>"-extra 2,3,4,5,6,7", endpoint=>E06);
    make(script=>'tsv2gff.pl', input=>{'<'=>fn($grp,A06,ssj,tsv)}, output=>{'>'=>fn($grp,E06,ssj,gff)}, between=>"-o count 2 -o stagg 3 -o entr 4 -o annot 5 -o nucl 6 -o IDR 7", endpoint=>E06);

    $prm = "-mincount $mincount";
    make(script=>"zeta.pl", input=>{-ssj=>fn($grp,A06,ssj,tsv), -ssc=>fn($grp,A06,ssc,tsv), -annot=>$annot}, output=>{'>'=>fn($grp,A07,gff)}, between=>$prm, endpoint=>A07);
    make(script=>"psicas.pl", input=>{-ssj=>fn($grp,A06,ssj,tsv), -annot=>$annot}, output=>{'>'=>fn($grp,B07,gff)}, endpoint=>B07);
}

#######################################################################################################################################################################

print "all :: A\n";

foreach $x([psi, cosi, A07, exon], [psi5, psi3, A07, intron], [cosi5, cosi3, A07, intron]) {
    ($one, $two, $l, $f) = @{$x};
    make2(script=>'merge_large_gff.pl', inputs=>{'<'=>{''=>$in}}, before=>"-ext .$l.gff -group $group -dir $dir$l/ -f $f $other", 
	 outputs=>{-out=>{$one=>fn($merge,$one,$l,tsv),$two=>fn($merge,$two,$l,tsv)}}, endpoint=>master);
}


make(script=>'merge_large_tsv.pl', input=>{'<'=>$in}, before=>"-ext .A06.ssj.tsv -group $group -dir $dir"."A06/ $other", output=>{'>'=>fn($merge,counts,ssj,tsv)}, endpoint=>master);
make(script=>'merge_large_tsv.pl', input=>{'<'=>$in}, before=>"-ext .A06.ssc.tsv -group $group -dir $dir"."A06/ $other", output=>{'>'=>fn($merge,counts,ssc,tsv)}, endpoint=>master);

make2(script=>"mk_stat_large.pl",   inputs=>{'<'=>{''=>$in}}, before=>"-ext .A06.ssj.tsv -group $group -dir $dir"."A06/", outputs=>{-out=>{'>'=>fn($merge,stats,tsv)}}, endpoint=>stats);
make(script=>"mk_stat.r", input=>{-i=>fn($merge,stats,tsv)}, output=>{-o=>fn($merge,stats,pdf)}, endpoint=>stats);


sub fn {
    return(@_[1]=~/^[A-Z]\d+$/ ? join(undef, $dir, @_[1], "/", join('.', @_)) : join(undef, $dir, join('.', @_)));
    
}


