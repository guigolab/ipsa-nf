#!/usr/bin/env perl

%input = (!{input});
@key = (!{by});
$out = "!{prefix}.tsv";

die unless(keys(%input)>0);

foreach $file(keys(%input)) {
    $name = $input{$file};
    print STDERR "[",++$N," $file $name";
    open FILE, $file;
    while($line = <FILE>) {
      @array = split /\t/, $line;
      unshift(@array, undef);
      $id = join(!{sep}, @array[@key]);
      $data{$name}{$id} += $array[!{val}];
      $rows{$id}++;
      $cols{$name}++;
    }
    print STDERR "]\n";
    close FILE;
}

@c = sort keys(%cols);
@r = sort keys(%rows);
open OUTFILE, '>', $out;
print OUTFILE join("\t", @c), "\n";
foreach $id(@r) {
    @arr = ($id);
    foreach $name(@c) {
        push @arr, 0 + $data{$name}{$id};
    }
    print OUTFILE join("\t", @arr), "\n";
}
close OUTFILE;