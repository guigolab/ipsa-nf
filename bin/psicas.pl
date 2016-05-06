#!/usr/bin/perl
use Perl::utils;


if(@ARGV==0) {
    print STDERR "This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs\n";
}

parse_command_line(annot    => {description=>'the annotation (GTF) file', ifunreadable=>'annotation not specified'},
                   ssj      => {description=>'the input ssj  (tsv) file', ifunreadable=>'input not specified'},
                   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'},
                  );

open FILE, $annot || die();
print STDERR "[<$annot";
while($line=<FILE>) {
    chomp $line;
    ($chr, $trash, $element, $beg, $end, $trash, $str, $frame, $attr) = split /\t/, $line;
    $elem{$element}{$chr}{$str}{$beg}{$end}++;
    $mele{$element}{$chr}{$str}{$end}{$beg}++;
}
close FILE;
print STDERR "]\n";

print STDERR "[<$ssj";
open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($chr, $beg, $end, $str, $count) = split /[\t\_]/, $line;
    $data{$chr}{$str}{$beg}{$end} = $count;
}
close FILE;
print STDERR "]\n";

foreach $chr(sort keys(%{$elem{exon}})) {
    foreach $str(sort keys(%{$elem{exon}{$chr}})) {
	foreach $x(sort{$a<=>$b} keys(%{$elem{exon}{$chr}{$str}})) {
	    foreach $y(sort{$a<=>$b} keys(%{$elem{exon}{$chr}{$str}{$x}})) {
		$a = (sort {$a<=>$b} keys(%{$mele{intron}{$chr}{$str}{$x}}))[-1];
		$b = (sort {$a<=>$b} keys(%{$elem{intron}{$chr}{$str}{$y}}))[0];
		next unless($elem{intron}{$chr}{$str}{$a}{$b});
		$inc = $data{$chr}{$str}{$a}{$x} + $data{$chr}{$str}{$y}{$b} + 0;
		$exc = $data{$chr}{$str}{$a}{$b} + 0;
		$psi  = frac($inc, 2*$exc);
		next unless($psi =~ /\d/);
        	print join("\t", $chr, 'SJPIPE', 'exon', $x, $y, int(1000*$psi), $str, '.', set_attributes(psicas=>$psi,inc=>$inc,exc=>$exc)), "\n";
	    }
	}
    }
}



