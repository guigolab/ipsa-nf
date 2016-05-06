#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "Merger script for tsv files with SJ counts\n";
}

parse_command_line(i   => {description=>'input gtf file name and label', array=>hash},
		   by  => {description=>'comma separated list of key columns', default=>"1"},
		   val => {description=>'column with values', default=>'2'},
		   sep => {description=>'separator',default=>'_'});

%input  = @i;
@key = split /\,/, $by;

die unless(keys(%input)>0);

foreach $file(keys(%input)) {
    $name = $input{$file};
    print STDERR "[",++$N," $file $name";
    open FILE, $file;
    while($line = <FILE>) {
	@array = split /\t/, $line;
	unshift(@array, undef);
        $id = join($sep, @array[@key]);
	$data{$name}{$id} += $array[$val];
	$rows{$id}++;
	$cols{$name}++;
    }
    print STDERR "]\n";
    close FILE;
}

@c = sort keys(%cols);
@r = sort keys(%rows);
print join("\t", @c),"\n";
foreach $id(@r) {
    @arr = ($id);
    foreach $name(@c) {
        push @arr, 0 + $data{$name}{$id};
    }
    print join("\t", @arr), "\n";
}


