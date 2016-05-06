#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility computes (a) the global exon inclusion and processing rates for a given set of annotated exons and (b) inclusion and processing rates of SJs\n";
}

parse_command_line(annot    => {description=>'the annotation (GTF) file'},
		   exons    => {description=>'exons file with exons to quantify'},
                   ssj      => {description=>'the input ssj (BED) file', ifunreadable=>'input not specified'},
                   ssc      => {description=>'the input ssc (BED) file'},
		   mincount => {default=>10,  description=>'the min value of the denominator of the fraction'},
                   stranded => {default=>1,   description=>'1(yes) or 0(no)'}
                  );

#######################################################################################################################

# Read SJs from ssj bed file and memorize splice sites
print STDERR "[<$ssj";
open FILE, $ssj || die();
while($line=<FILE>) {
    chomp $line;
    ($id, $count) = split /\t/, $line;
    ($chr, $beg, $end, $s) = split /\_/, $id;
    $str = strand_c2i($s) * $stranded;
    $site{$chr}{$str}{$beg} = 1;
    $site{$chr}{$str}{$end} = 1;
}
close FILE;
if(-r $annot) {
    print STDERR ",$annot";
    open FILE, $annot;
    while($line=<FILE>) {
        chomp $line;
        ($chr, $source, $feature, $beg, $end, $score, $s) = split /\t/, $line;
	next unless($feature eq "exon");
	$str = strand_c2i($s) * $stranded;
	$site{$chr}{$str}{$beg} = 1;
    	$site{$chr}{$str}{$end} = 1;
    }
    close FILE;
}
print STDERR ", indexing";
foreach $chr(sort keys(%site)) {
    foreach $str(sort keys(%{$site{$chr}})) {
	foreach $pos(sort {$a<=>$b} keys(%{$site{$chr}{$str}})) {
	    $index{$chr}{$str}{$pos} = ++$n;
	}
    }
}
print STDERR ", n=$n]\n";

#######################################################################################################################

# Read exon table and create two hash tables, lb and rb
# lb is an array which appoints to each exon boundary l |-> list of exons whose left boundary is l
# rb ... r |-> list of exons whose right boundary is r

if(-r $annot) {
    print STDERR "[<$annot";
    open FILE, $annot;
    while($line=<FILE>) {
    	chomp $line;
    	($chr, $source, $feature, $beg, $end, $score, $s) = split /\t/, $line;
	next unless($feature eq "exon");
    	$e = join("_", $chr, $beg, $end, $s);
	next if($been{$e}++);
    	$str = strand_c2i($s) * $stranded;
    	$l = $index{$chr}{$str}{$beg};
    	$r = $index{$chr}{$str}{$end};
    	push @{$lb[$l]}, $e;
    	push @{$rb[$r]}, $e;
    }
    close FILE;
    print STDERR "] ", 0+keys(%been), "\n";
}

if(-r $exons) {
    print STDERR "[<$exons";
    open FILE, $exons || die();
    while($line=<FILE>) {
    	chomp $line;
    	($e) = split /\t/, $line;
	($chr, $beg, $end, $s) = split /\_/, $e;
	next if($been{$e}++);
	$str = strand_c2i($s) * $stranded;
	$l = $index{$chr}{$str}{$beg};
	$r = $index{$chr}{$str}{$end};
	push @{$lb[$l]}, $e;
	push @{$rb[$r]}, $e;
    }
    close FILE;
    print STDERR "] ", 0+keys(%been), "\n";
}

#######################################################################################################################

# Read ssj table again and gather inclusion and exclusion counts
print STDERR "[<$ssj";
open FILE, $ssj || die("Can't open $ssj\n");
while($line=<FILE>) {
    chomp $line;
    ($id, $count) = split /\t/, $line;
    ($chr, $beg, $end, $s) = split /\_/, $id;
    $str = strand_c2i($s) * $stranded;
    $start = $index{$chr}{$str}{$beg}; 		# index of the left side of SJ
    $stop  = $index{$chr}{$str}{$end}; 		# index of the right side of SJ
    next unless($start && $stop && $start < $stop);
    foreach $e(@{$rb[$start]}) {			# for each exon whose right boundary is the left side of the SJ
	$inclusion_r{$e} += $count;			# increment inclusion counter
    }
    foreach $e(@{$lb[$stop]}) {				# for each exon whose left boundary is the right side of the SJ
	$inclusion_l{$e} += $count;			# increment inclusion counter
    }
    for($k = $start + 1; $k<$stop; $k++) {		# for each splice site in between left and right sides of the SJ
	foreach $e(@{$lb[$k]}) {			#     for each exon with that boundary
	    ($chr1, $beg1, $end1, $s1) = split /\_/, $e;
	    $str1 = strand_c2i($s1) * $stranded;		# if the exon is fully contain in the interval defined by SJ
	    if($chr1 eq $chr && $str1 eq $str && $beg1>$beg && $end1<$end) {
		$exclusion{$e} += $count;		# increment the exclusion counter
	    }
    	}
    }
    ($beg, $end) = reverse ($beg, $end) if($str<0);
    $count53{$chr}{$str}{$beg}{$end}+=$count;
    $count5X{$chr}{$str}{$beg}+=$count;
    $countX3{$chr}{$str}{$end}+=$count;
}
close FILE;
print STDERR "]\n";

if(-e $ssc) {
    print STDERR "[<$ssc";
    open FILE, $ssc || die();
    while($line=<FILE>) {
    	chomp $line;
	($id, $count) = split /\t/, $line;
    	($chr, $pos, $s) = split /\_/, $id;
	$str = strand_c2i($s) * $stranded;
        $loc = $index{$chr}{$str}{$pos};
	next unless($loc);
        foreach $e(@{$rb[$loc]}) {
            $retention_r{$e} += $count;
        }
        foreach $e(@{$lb[$loc]}) {
            $retention_l{$e} += $count;
    	}
        $count00{$chr}{$str}{$pos}+=$count;
    }
    print STDERR "]\n";
    close FILE;
}

#######################################################################################################################

print STDERR "[>stdout";
for($j=1;$j<=$n;$j++) {
    foreach $e(@{$lb[$j]}) {
	$inc = $inclusion_l{$e} + $inclusion_r{$e} + 0;
	$exc = $exclusion{$e} + 0;
	$ret = $retention_l{$e} + $retention_r{$e} + 0;
	$psi  = frac($inc, 2*$exc);
	$cosi = frac($inc + 2*$exc, $ret);
	($chr, $beg, $end, $s) = split /\_/, $e;
	next unless(join(undef, $psi, $cosi) =~ /\d/);
	print join("\t", $chr, 'SJPIPE', 'exon', $beg, $end, int(1000*$psi), $s, '.', set_attributes(psi=>$psi, cosi=>$cosi, inc=>$inc, exc=>$exc, ret=>$ret)), "\n";
    }
}
foreach $chr(sort keys(%count53)) {
    foreach $str(sort {$a<=>$b} keys(%{$count53{$chr}})) {
        foreach $beg(sort {$a<=>$b} keys(%{$count53{$chr}{$str}})) {
            foreach $end(sort {$a<=>$b} keys(%{$count53{$chr}{$str}{$beg}})) {
                $nDA = $count53{$chr}{$str}{$beg}{$end};
                $nDX = $count5X{$chr}{$str}{$beg} - $count53{$chr}{$str}{$beg}{$end};
                $nXA = $countX3{$chr}{$str}{$end} - $count53{$chr}{$str}{$beg}{$end};
                $nD = $count00{$chr}{$str}{$beg} + 0;
                $nA = $count00{$chr}{$str}{$end} + 0;
                $psi5  = frac($nDA, $nDX);
                $psi3  = frac($nDA, $nXA);
                $cosi5 = frac($nDA + $nDX, $nD);
                $cosi3 = frac($nDA + $nXA, $nA);
		next unless(join(undef, $psi5,$psi3,$cosi5,$cosi3) =~ /\d/);
		($x, $y) = sort {$a<=>$b} ($beg, $end);
                print join("\t", $chr, 'SJPIPE', 'intron', $x, $y, $psi5=~/\d/ ? int(500*($psi5 + $psi3)) : '.', strand_i2c($str), '.',
                        set_attributes(psi5=>$psi5, psi3=>$psi3, cosi5=>$cosi5, cosi3=>$cosi3, nDA=>$nDA, nDX=>$nDX, nXA=>$nXA, nD=>$nD, nA=>$nA)), "\n";
            }
        }
    }
}
print STDERR "]\n";



