#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
}

parse_command_line( ssj => {description=>'the splice junction file', ifunreadable=>'bed not specified'},
                    minstaggered=>{description=>'he minimum umber of staggered reads', default=>2},
                    nucleotides =>{description=>'the splice site nucleotides', default=>GTAG});

open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($id, $total, $staggered, $entropy, $annot, $nuc) = split /\t/, $line;
    $jnc{$id}++ if($nuc eq $nucleotides && $staggered>=$minstaggered);
}
close FILE;

while($line=<STDIN>) {
    chomp $line;
    ($id, $rest) = split /\t/, $line, 2;
    @arr = split /\_/, $id;
    $chr = shift(@arr);
    $strand = pop(@arr);
    foreach $str("+", "-") {
        if($strand eq $str || $strand eq '.') {
	    $flag = 1;
	    for($i = 0; $i < @arr; $i += 2) {
		$flag = undef unless($jnc{join("_", $chr, $arr[$i], $arr[$i+1], $str)});
	    }
    	    print join("\t", join("_", $chr, @arr, $str), $rest), "\n" if($flag);
	}
    }
}
