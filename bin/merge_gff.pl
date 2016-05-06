#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility merges gtf/gff files into rectangular tables with colnames\n";
    print STDERR "Example: $0 -i <gtf_1> <label_1> -i <gtf_2> <label_2> ... -o feature_1 <tsv1> -o feature_2 <tsv2>\n";
}

parse_command_line(i => {description=>'input gtf file name and label', array=>hash},
		   o => {description=>'feature and output tsv file name', array=>hash},
		   percent=>{description=>'do not report rows that have more than <value> NAs', default=>0.25},
		   transf=>{description=>'transformation (log, logit)'});

################################################################################################################################

%input  = @i;
%output = @o;

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
	    next unless($attr{$key});
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
        print FILE join("\t", @arr), "\n" if($num>$percent*@c);
    }
    close FILE;
    print STDERR "]\n";
}


