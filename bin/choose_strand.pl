#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility takes a TSV4+3+2 file (STDIN) where column (8) is the annotation status and column (9) is splice sites, and selects for each pair of beg/end the strand based on these two columns (STDOUT)\n";
}

parse_command_line(annot   => {default=>5, description=>'annotation column'},
		   sites   => {default=>6, description=>'splice site column'},
		   logfile => {description=>'name of the log file'});

while($line=<STDIN>) {
    chomp $line;
    @array = split /\t/, $line;
    ($chr, $beg, $end, $str) = split /\_/, $array[0];
    $id = join("_", $chr, $beg, $end);
    push @{$data{$id}}, [@array];
}

foreach $id(sort keys(%data)) {
    @array = sort my_sort @{$data{$id}};
    print join("\t", @{$array[0]}),"\n";
    $stat1{$array[0]->[4]}{$array[0]->[5]}++;
    $stat2{$array[0]->[4]}{$array[0]->[5]}+=$count;
}

if($logfile) {
    open FILE, ">$logfile";
    foreach $status(sort keys(%stat1)) {
        foreach $nucl(sort keys(%{$stat1{$status}})) {
            print FILE join("\t", $logfile, $status, $nucl, $stat1{$status}{$nucl}, $stat2{$status}{$nucl}),"\n";
        }
    }
    close FILE;
}


sub  my_sort {
    return(-1) if($b->[$annot-1] <  $a->[$annot-1]);
    return(1)  if($b->[$annot-1] >  $a->[$annot-1]);
    return(-1) if($b->[$sites-1] lt $a->[$sites-1]);
    return(1)  if($b->[$sites-1] gt $a->[$sites-1]);
    return(0)
}


