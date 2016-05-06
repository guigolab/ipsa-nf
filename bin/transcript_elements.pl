#!/usr/bin/perl
use Perl::utils;

if(@ARGV==0) {
    print STDERR "This utility takes a gtf annotation (STDIN) and reformats it into a more compact, quickly readable file (STDOUT). Only exons are taken into account.\n";
}

parse_command_line(source  =>{description=>'the content of the source field', default=>'SJPIPE'},
		   features=>{description=>'list of features to output', default=>'transcript_id'},
		   lim     =>{description=>'line limit for debugging', default=>0});

foreach $item (split /\,/, $features) {
    ($f, $g) = split /\:/, $item;
    $g = 'uniq' unless($g);
    $feature_list{$f} = \&$g;
    print STDERR "[function $g($f)]\n"; 
}

print STDERR "[<stdin";
while($line=<STDIN>) {
    chomp $line;
    ($chr, $src, $element, $beg, $end, $trash, $str, $trash, $attr) = split /\t/, $line;
    %attr =  get_attributes($attr);
    $tid = $attr{'transcript_id'};
    $CHR{$tid}{$chr} = $STR{$tid}{$str} = 1;
    push @{$exons{$tid}}, [$beg, $end] if($element eq "exon");
    push @{$cds{$tid}}, [$beg, $end]   if($element eq "CDS");
    foreach $f(keys(%feature_list)) {
            $tr_feature{$f}{$tid} = $attr{$f};
    }
    $tr_feature{source}{$tid} = $src;
    $n++;
    last if($n>$lim && $lim>0);
}
print STDERR "]\n";

foreach $tid(keys(%exons)) {
    $N++;
    progressbar($N, 0+keys(%exons), "Processing ");
    if(keys(%{$CHR{$tid}})==1 && keys(%{$STR{$tid}})==1) {
	($chr, $str) = (keys(%{$CHR{$tid}}), keys(%{$STR{$tid}}));
	@exon_list = sort {$a->[0]<=>$b->[0]} @{$exons{$tid}};
	@cds_list  = sort {$a->[0]<=>$b->[0]} @{$cds{$tid}};
	for($i=0,$j=0; $i<@exon_list; $i++) {
	    for(; $j<@cds_list && $cds_list[$j]->[0]<$exon_list[$i]->[0]; $j++) {}
	    $coding_part = min($exon_list[$i]->[1], $cds_list[$j]->[1]) - max($exon_list[$i]->[0],$cds_list[$j]->[0]) + 1;
	    $key = join("\t", $chr, $source, 'exon',   $exon_list[$i]->[0],   $exon_list[$i]->[1], '.', $str, '.');
	    push @{$feature{$key}{transcript_id}}, $tid;
	    push @{$feature{$key}{position}}, relpos($i, @exon_list-1);
	    push @{$feature{$key}{coding}}, max($coding_part/($exon_list[$i]->[1] - $exon_list[$i]->[0] + 1), 0);
	    foreach $f(keys(%feature_list)) {
                push @{$feature{$key}{$f}}, $tr_feature{$f}{$tid} if($tr_feature{$f}{$tid});
            }
	    if($i>0) {
		$key = join("\t", $chr, $source, 'intron', $exon_list[$i-1]->[1], $exon_list[$i]->[0], '.', $str, '.');
		push @{$feature{$key}{transcript_id}}, $tid;   
		push @{$feature{$key}{position}}, relpos($i-1, @exon_list-2);
		foreach $f(keys(%feature_list)) {
                    push @{$feature{$key}{$f}}, $tr_feature{$f}{$tid} if($tr_feature{$f}{$tid});
            	}
	    }
	}
    }
    else {
	$trans_spliced++;
    }
}

print STDERR "[WARNING: $trans_spliced trans spliced transcripts excluded]" if($trans_spliced);

print STDERR "[>stdout";
foreach $key(sort keys(%feature)) {
    my %res=();
    foreach $f(sort keys(%feature_list)) {
        $res{$f} = $feature_list{$f}->(@{$feature{$key}{$f}}) if(@{$feature{$key}{$f}}>0);
    }
    print $key, "\t", set_attributes(%res), "\n";
}
print STDERR "]\n";

sub relpos {
    return('NA') unless(@_[0]>=0 && @_[1]>0);
    return(strand_c2i($str)<0 ? 1 - @_[0]/@_[1] : @_[0]/@_[1]);
}

sub abspos {
    return('I') if(@_[0]>0 && @_[0]<@_[1]);
    return(strand_c2i($str)>0 ? '5' : '3') if(@_[0]==0);
    return(strand_c2i($str)>0 ? '3' : '5');
}



