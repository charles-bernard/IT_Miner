#!/usr/bin/awk

BEGIN {
	FS = "\t";
	# The variable INTRA_FILE is sent to this script.
}

NR == 1 {
	printf("%s", $0);
	printf("%s", $0) > INTRA_FILE;
}

NR > 1 {
	upDist = $10;
	dwDist = $14;

	if(upDist > 0 && dwDist > 0) {
		printf("\n%s", $0);
	} else {
		printf("\n%s", $0) > INTRA_FILE;
	}
}
