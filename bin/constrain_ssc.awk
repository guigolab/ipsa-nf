#!/usr/bin/awk -f
# This script chooses the strand for 0-splits based on filtered SJs
BEGIN {
    if(length(jncfile) == 0) {
	print "junction file missing" > "/dev/stderr"
        exit 1 
    }
    while (getline < jncfile >0) {
	split($1, a, "_");
	info[a[1]" "a[2]" "a[4]] = info[a[1]" "a[3]" "a[4]] = 1;
    }
    strand["+"] = strand["-"] = 1;
    OFS = "\t";
};
{
    split($1, a, "_");
    for(s in strand) {
	if((a[3] == "." || a[3] == s) && info[a[1]" "a[2]" "s]) {
	    $1 = a[1] "_" a[2] "_" s;
	    print $0;
	}
    }
}


