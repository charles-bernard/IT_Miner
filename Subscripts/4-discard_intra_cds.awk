#!/usr/bin/awk

BEGIN {
	FS = "\t";
}

NR == 1 {
	printf("%s", $0);
}

NR > 1 {
	upDist = $10;
	dwDist = $14;

	if(upDist > 0 && dwDist > 0) {
		printf("\n%s", $0);
	}
}
