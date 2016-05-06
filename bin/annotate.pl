#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility takes an aggregated TSV file, the genomic annotation, and the genome, and outputs a TSV with two more columns: ";
    print STDERR "(8) annotation status and (9) splice sites\n";
    print STDERR "If the input was strandless ('.' in column 4) then each line will be reported twice, one for '+' and and for '-' strand\n";
}

parse_command_line(	in	=> {description=>'the input tsv file', ifunreadable=>'input not specified'},
			annot	=> {description=>'the annotation (gtf)', ifunreadable=>'annotation not specified'}, 
			dbx	=> {description=>'the genome (dbx)', ifunreadable=>'dbx not specified'},
                        idx     => {description=>'the genome (idx)', ifunreadable=>'idx not specified'},
			deltaSS => {description=>'distance threshold for splice sites', default=>0},
			logfile => {description=>'name of the log file'},
			MAPTOOLSDIR  =>{variable=>T, ifabsent=>'MAPTOOLSDIR not specified'});


read_junctions($in);
read_annotation($annot);
index_ss();

$program = $MAPTOOLSDIR."bin/getsegm2 -limit 4 -margins -1 0 -spacer 0 -inp_type 2 -out_type 1";
%seq = split /[\t\n]/, `cat $in | $program -dbx $dbx -idx $idx`;
#print STDERR "[cat $in | $program -dbx $dbx -idx $idx]\n";

open FILE, $in || die;
while($line=<FILE>) {
    chomp $line;
    ($id, $count, $rest) = split /\t/, $line, 3;
    ($chr, $beg, $end, $strand) = split /\_/, $id;
    foreach $str("+", "-") {
	if($strand eq $str || $strand eq '.') {	
	    $status = annot_status($chr, $beg, $end, $str);
	    $id = join("_", $chr, $beg, $end, $str);
	    $nucl   = $seq{$id};
	    $nucl   = "NA" unless($nucl);
	    $nucl =~ tr/[a-z]/[A-Z]/;
	    print join("\t", $id, $count, $rest, $status, $nucl), "\n";
	    $stat1{$status}{$nucl}++;
	    $stat2{$status}{$nucl}+=$count;
	}
    }
}
close FILE;

if($logfile) {
    open FILE, ">$logfile";
    foreach $status(sort keys(%stat1)) {
	foreach $nucl(sort keys(%{$stat1{$status}})) {
	    print FILE join("\t", $logfile, $status, $nucl, $stat1{$status}{$nucl}, $stat2{$status}{$nucl});
	}
    }
    close FILE;
}

#################################################################################################################################################################

sub read_junctions {
    my $N=0;
    print STDERR "[<@_[0]";
    open FILE, @_[0] || die ('Cannot read junctions');
        while($line=<FILE>) {
        chomp $line;
        ($chr, $beg, $end, $strand) = split /[\_\t]/, $line;
        foreach $str("+", "-") {
            if($strand eq $str || $strand eq '.') {
                $SS{"b;$chr;$str"}{$beg} = f;
                $SS{"e;$chr;$str"}{$end} = f;
            }
        }
        $N++;
    }
    close FILE;
    print STDERR ", $N junctions]\n";
}

sub read_annotation {
    my $N=0;
    print STDERR "[<@_[0]";
    open FILE, @_[0] || die ('Cannot read annotation');
        while($line=<FILE>) {
        chomp $line;
        ($chr, $source, $feature, $beg, $end, $score, $str, $frame, $group) = split /\t/, $line;
        if($feature eq "intron") {
            $SJ{$chr}{$beg}{$end}{$str}++;
            $SS{"b;$chr;$str"}{$beg} = t;
            $SS{"e;$chr;$str"}{$end} = t;
        }
        $N++;
    }
    close FILE;
    print STDERR ", $N annotated introns]\n";
}

sub index_ss {
    foreach $key(keys(%SS)) {
        @arr = sort {$a<=>$b} keys(%{$SS{$key}});
        $SS_next{$key}{0} = $arr[1];
        for($i=1;$i<@arr-1;$i++) {
            $SS_next{$key}{$arr[$i]} = $arr[$i+1];
            $SS_prev{$key}{$arr[$i]} = $arr[$i-1];
        }
        $SS_prev{$key}{$arr[-1]} = $arr[-2];
    }
}

#################################################################################################################################################################

sub call_site {
    my ($type, $chr, $pos, $str) = @_;
    my $key = "$type;$chr;$str";
    my %res = ();
    for(my $x = $pos; abs($x-$pos)<=$deltaSS; $x = $SS_next{$key}{$x}) {
        last unless($x);
        $res{$x} = abs($x-$pos) if($SS{$key}{$x} eq t);
    }
    for(my $x = $pos; abs($x-$pos)<=$deltaSS; $x = $SS_prev{$key}{$x}) {
        last unless($x);
        $res{$x} = abs($x-$pos) if($SS{$key}{$x} eq t);
    }
    my @out = sort {$res{$a} <=> $res{$b}} keys(%res);
    return($SS{$key}{$pos} eq t ? 2 : (@out > 0 ? 1 : 0));
}

sub annot_status {
    my ($chr, $beg, $end, $str) = @_;
    return($SJ{$chr}{$beg}{$end}{$str} ? 3 : $SITE_MATRIX[call_site(b, $chr, $beg, $str)][call_site(e, $chr, $end, $str)]);
}

