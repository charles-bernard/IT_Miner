#!/usr/bin/awk

BEGIN {
	FS = "\t";

	# file counter
	file_idx = 0;

	# terminator counter
	t = 0; 
}

FNR == 1 {
	file_idx++;

	if(file_idx == 1) {
		# First Input File corresponds to RNIE Genome output
		mode = "Genome";
	} else {
		# Second to RNIE Gene output
		mode = "Gene";
	}
}

{
	# Important Fields:
	start[t] = $4;
	end = $5;
	bit_score = $6;
	strand = $7;

	line[t] = start[t] "\t" end "\t" strand "\t" mode "\t" bit_score;

	t++;
}

END {
	header = "Start\tEnd\tStrand\tRNIE Mode\tBit Score";
	printf("%s", header);

	# Sort Terminators based on 'Start'
	# 'sorted_order' stores the new order
	n = asorti(start, sorted_order, "@val_num_asc");

	# Print each line in the sorted order
	for (i = 1; i <= n; i++) {
		printf("\n%s", line[sorted_order[i]]);
	}
}