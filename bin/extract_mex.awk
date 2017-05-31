#!/usr/bin/awk -f
# This script extracts mini-exon information from 2-splits and reports them
function max(a, b) {
    return(a>b ? a : b);
};
{
    n = split($1,a,"_");
    if(n != 6) next;
    id = a[1] "_" a[3] "_" a[4] "_" a[6];
    count[id] += $2;
    stag[id] = max(stag[id], $3);
    entr[id] = max(entr[id], $4);
    OFS = "\t";
};
END{
    for(id in count) {
	print id, 0 + count[id], 0 + stag[id], 0 + entr[id];
    }
}

