#!/usr/bin/awk -f
# This script aggregates the output of sjcount by offsets
BEGIN {
    if(length(degree) == 0)    degree = 1;
    if(length(minintron) == 0) minintron = 4;
    if(length(maxintron) == 0) maxintron = 0;
    OFS = "\t";
};
{
    if($2 != degree) next;
    id = prefix $1;
    split(id, array, "_");
    if($2 == 1 && (array[3] - array[2] - 2 < minintron || (maxintron > 0 && array[3] - array[2] - 2 > maxintron)))  next;
    stat[$3] += $4;

    if(readLength && margin && ($3 < margin || $3 >= readLength - margin)) next;
    data[id, ++n[id]] = $4;
};
END {
    for(id in n) {
	s = 0;
    	c = 0;
    	l = 0;
	for(i=1; i<=n[id]; i++) {
	    value = data[id, i];
	    s += value;
            c += 1;
            l += value*log(value);
        }
    	h = log(s) - l/s;
	printf "%s\t%i\t%i\t%2.2lf\n", id, s, c, h/log(2);	
    };
    if(length(logfile) > 0) {
	for(i in stat) {
	    print i, stat[i] > logfile;
	}
    };
};
