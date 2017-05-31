#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility transforms tsv to bed (tsv is supplied at STDIN)\n";
}

parse_command_line(minus =>{description => 'base to subtract from the coordinates', default=>0},
		   source=>{description => 'the source field',default=>'IPSA'},
		   o     =>{description => 'hash, attr name >> column number}', array=>hash});

%fields = @o;

while($line=<STDIN>) {
    chomp $line;
    @array = split /\t/, $line;
    unshift @array, undef;
    ($chr, $beg, $end, $str) = split /\_/, $array[1];
    $score = int(100*log($array[2])/log(2));
    $score = 1000 if($score>1000);
    %output = ();
    foreach $key(sort keys(%fields)) {
	$output{$key} = $array[$fields{$key}];
    }
    print join("\t", $chr, $source, 'SJ', $beg - $minus, $end - $minus, $score, $str, '.', set_attributes(%output)), "\n";
}
