#!/usr/bin/awk -f
# This script is used to choose strand for SJs on opposite chains of DNA
function comp_func(v1, v2) {
    split(v1, a1, "\t");
    split(v2, a2, "\t");
    if(a1[annot] < a2[annot]) return(-1);
    if(a1[annot] == a2[annot] && a1[sites] < a2[sites]) return(-1); #// needs to be adjusted to GT/AG or other tabulated values
    return(1);
}
BEGIN {
    if(length(annot) == 0) annot = 5;
    if(length(sites) == 0) sites = 6;
};
{   split($1, array, "_");
    key = array[1] array[2] array[3];
    data[key, ++n[key]] = $0;
};

END {
    for(key in n) {
	print (n[key]==1) ? data[key,1] : (comp_func(data[key,1], data[key,2]) > 0) ? data[key,1] : data[key,2];
    }
};
