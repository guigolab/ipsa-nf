#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
}


parse_command_line(dir    => {description=>'directory name', ifabsent=>'need directory name'},
                   ext    => {description=>'file extention', ifabsent=>'need file extention'},
                   group  => {description=>'sample id field', default=>'labExpId'},
		   subset => {description=>'file to subset id on'},
                   by     => {description=>'comma separated list of key columns', default=>"1"},
                   val    => {description=>'column with values', default=>'2'},
                   sep    => {description=>'separator',default=>'_'});

@group = split /\,/, $group;
die unless(@group>0);

@key = split /\,/, $by;

if($subset) {
    open FILE, $subset || die("Can't read $subset\n");
    while($line=<FILE>) {
    	chomp $line;
    	$subset{$line}++;
    }
    close FILE;
}

while($line=<STDIN>) {
    chomp $line;
    ($file, $attr) = split /\t/, $line;
    %attr = get_features($attr);
    $name = join("_", @attr{@group});    # sample id
    next unless($name);
    if($attr{'type'} eq "bam" && $attr{'view'}=~/^Alignments/) {
        $file = "$dir$name$ext";
	$input{$file} = $name;
    }
}

foreach $file(keys(%input)) {
    $name = $input{$file};
    print STDERR "[",++$N," $file $name";
    open FILE, $file;
    while($line = <FILE>) {
	@array = split /\t/, $line;
        unshift(@array, undef);
        $id = join($sep, @array[@key]);
	next if($subset && $subset{$id} eq undef);
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
