#!/usr/bin/awk

BEGIN {
	FS = "\t";
}


BEGIN {
	FS = "\t";
}


NR == FNR {
	strand[$12] = 1;
	start[$10] = 1;
	end[$11] = 1;

	if($2 == "0") {
		keep_strand[$12] = 1;
		keep_start[$10] = 1;
		keep_end[$11] = 1;
	}

}

NR > FNR && FNR == 1 {
	print $0;
}
NR > FNR && FNR > 1 {
	if(!strand[$3] || !start[$1] || !end[$2]) {
		printf("%s\t\n",$0);
	}

	if(keep_strand[$3] && keep_start[$1] && keep_end[$2]) {
		printf("%s\tno_overlap join\n",$0);
	}
}
