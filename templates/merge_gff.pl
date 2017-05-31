#!/usr/bin/env perl
use Perl::utils;

%input = (!{input});
%output = (!{output});
$percent = !{percent};
$transf = "!{transf}";

die unless(keys(%output)>0 && keys(%input)>0);

foreach $file(keys(%input)) {
    $name = $input{$file};
    print STDERR "[",++$N,":$file $name";
    open FILE, $file;
    while($line = <FILE>) {
	($chr, $source, $feature, $beg, $end, $score, $str, $frame, $attr) = split /\t/, $line;
        $id = "$chr\_$beg\_$end\_$str";
	%attr = get_attributes($attr);
	foreach $key(keys(%output)) {
	    next unless($attr{$key}=~/\w/);
	    $data{$key}{$name}{$id} = $attr{$key};
	    $rows{$key}{$id}++;
	    $cols{$key}{$name}++;
	}
    }
    print STDERR "]\n"; 
    close FILE;
}

foreach $key(sort keys(%output)) {
    @c = sort keys(%{$cols{$key}});
    @r = sort keys(%{$rows{$key}});
    print STDERR "[$key=>$output{$key}";
    open FILE, ">$output{$key}" || die("can't write to this file");
    print FILE join("\t", @c),"\n";
    foreach $id(@r) {
	@arr = ($id);
	$num = 0;
    	foreach $name(@c) {
	    $value = $data{$key}{$name}{$id};
            push @arr, $value =~/\d/ ? $value : "NA";
            $num++ if($value =~/\d/);
        }
        print FILE join("\t", @arr), "\n" if($num>=$percent*@c);
    }
    close FILE;
    print STDERR "]\n";
}