#!/usr/bin/awk

function abs(x) {
	return v < 0 ? -x : x
}

BEGIN {
	FS = "\t";
	t = -1; # terminator counter
	# The variable 'dev' is sent to this script

	# These variables will be needed for further comparisons
	old_strand = "";
	old_gene = old_start = 0;
	best_bit = 0;
}

NR == 1 {
	header = $0;
}

NR > 1 {
	start = $1;
	strand = $3;
	bit = $5;
	gene = $7;

	# If the current terminator is the same as the previous one(s)
	if(strand == old_strand && gene == old_gene \
		&& abs(start-old_start) <= dev) {

		# If the current terminator has the best bit score,
		# its line is stored
		if(bit > best_bit) {
			best_bit = bit;
			best_line = $0;
		}

	} else {

		# Since the current terminator is different from the previous one(s),
		# the line of the previous best terminator will be further printed.  
		line[t] = best_line;
		t++;

		best_bit = bit;
		best_line = $0;
	}

	old_start = start;
	old_strand = strand;
	old_gene = gene;
}

END {
	printf("%s", header);
	# We don't retrieve line[-1] since it stores the first terminator,
	# which might not be the best duplicate
	for(i = 0; i < t; i ++) {
		printf("\n%s", line[i]);
	}
}