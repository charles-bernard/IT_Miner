#!/usr/bin/awk

BEGIN {
	FS = "\t";

	# The variable 'INTRA_FILE' is sent to this script.
	# It is needed to print the ITs that are intra CDS
	# into another file. This information can be
	# indeed valuable for the user.
}

NR == 1 {
	printf("%s", $0);
	printf("%s", $0) > INTRA_FILE;
}

NR > 1 {
	upDist = $10;
	dwDist = $14;

	if((upDist == "" || upDist > 0) && (dwDist == "" || dwDist > 0)) {
		printf("\n%s", $0);
	} else {
		printf("\n%s", $0) > INTRA_FILE;
	}
}
