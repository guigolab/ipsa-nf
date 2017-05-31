#!/usr/bin/awk -f
# This script is used to select 2-splits that have canonical SJs and at least 2 staggered reads
# Perhaps this is to be changed to annotated status
BEGIN {
    if(length(minstaggered) == 0) minstaggered = 2;
    if(length(nucleotides)  == 0) nucleotides  = "GTAG";
    while (getline < jncfile >0) {
        if($3>=minstaggered && $6==nucleotides) info[$1]++;
    }
    strand["+"] = 1;
    strand["-"] = -1;
    OFS = "\t";
};
{
    n = split($1, a, "_");
    for(s in strand) {
        if(a[n] == "." || a[n] == s) {
	    flag = 1;
	    id = a[1];
	    for(i = 2; i <= n - 1; i += 2) {
		if(! info[a[1] "_" a[i] "_" a[i+1] "_" s]) flag = 0;
		id = id "_" a[i] "_" a[i+1]
	    }
	    $1 = id "_" s;
            if(flag) print $0;
        }
    }
}
