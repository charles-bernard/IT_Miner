#!/usr/bin/awk

BEGIN {
	FS = "\t";
	# The variable 'cutoff' is sent to this script.
}

NR == 1 {
	printf("%s", $0);
}

NR > 1 && $10 <= cutoff {
	printf("\n%s", $0);
}
