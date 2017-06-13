#!/usr/bin/awk

BEGIN {
	FS = "\t";

	rev["+"] = "-";
	rev["-"] = "+";
}

NR == 1 {
	printf("%s", $0);
}

NR > 1 {
	printf("\n%s", $0);

	# Fake complement coordinates to IT with no id_grp
	id_grp = $15;
	if(! id_grp) {
		start = $1;
		end = $2;
		strand = $3;
		printf("\n%s\t%s\t%s\tFake\t", start, end, rev[strand]);
	}
}