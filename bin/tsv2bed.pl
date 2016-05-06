#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility transforms tsv to bed (tsv is supplied at STDIN)\n";
}

parse_command_line(extra=>{description => 'comma-separated list of column numbers to make BED6+n'},
		   minus=>{description => 'base to subtract from the coordinates', default=>1},
		   name =>{description => 'use junction id as name', store=>T},
		   ssc  =>{description => 'if exon-intron junction', store=>T});

@extra = split /\,/, $extra;
while($line=<STDIN>) {
    chomp $line;
    @array = split /\t/, $line;
    unshift @array, 0;
    if($ssc) {
    	($chr, $beg, $str) = split /\_/, $array[1];
	$end = $beg;
    }
    else {
	($chr, $beg, $end, $str) = split /\_/, $array[1];
    }
    $score = int(100*log($array[2])/log(2));
    $score = 1000 if($score>1000);
    print join("\t", $chr, $beg - $minus, $end - $minus, $name eq T ? $id : '.', $score, $str, @array[@extra]), "\n";
}
