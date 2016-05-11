#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
}


while($line=<STDIN>) {
    chomp $line;
    ($id, $c, $s, $e) = split /\t/, $line;
    @arr = split /\_/, $id;
    next unless(@arr==6);
    ($chr, $x, $y, $z, $t, $str) = split /\_/, $id;
    $id = join("_", $chr, $y, $z, $str);
    $count{$id}+= $c;
    $stag{$id} = $s if($s > $stag{$id}+0);
    $entr{$id} = $e if($e > $entr{$id}+0);
}

foreach $id(sort keys(%count)) {
    print join("\t", $id, $count{$id}, $stag{$id}, $entr{$id}), "\n";
}
