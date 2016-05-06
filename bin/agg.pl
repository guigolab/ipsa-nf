#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility aggregates the output of sjcount (STDIN) by the 5th column (offset) and outputs a TSV (STDOUT) with three extra columns being ";
    print STDERR "(5) total count, (6) staggered read count, (7) entropy\n";
}

parse_command_line(margin	=>{default=>0, 	description=>'the margin for offset'}, 
		   readLength	=>{default=>0, 	description=>'the read length'},
		   logfile	=>{description=>'name of the log file'},
		   minintron    =>{default=>4,  description=>'min intron length; only applied to splits of degree 1'},
                   maxintron    =>{default=>0,  description=>'max intron length; only applied to splits of degree 1'},
		   prefix	=>{description=>'prefix in the chromosome name', default=>""}
		  );

while($line=<STDIN>) {
    chomp $line;
    ($id, $degree, $offset, $count) = split /\t/, $line;
    if($degree==1) {
        ($chr, $beg, $end, $str) = split /\_/, $id;
         next if($minintron  && ($end - $beg - 2 < $minintron) || $maxintron && ($end - $beg - 2 > $maxintron));
    }
    $id ="$prefix$id";
    $stat[$offset]+=$count;
    next if($readLength && $margin && ($offset < $margin || $offset >= $readLength - $margin));
    push @{$data{$id}}, $count;
}

foreach $id(sort keys(%data)) {
    @stats = aggstat(@{$data{$id}});
    print join("\t", $id, @stats), "\n";
}

if($logfile) {
    open FILE, ">$logfile";
    for($offset=0;$offset<$readLength;$offset++) {
	print FILE join("\t", $logfile, $offset, 0+$stat[$offset]), "\n";
    }
    close FILE;
}
