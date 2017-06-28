#!/usr/bin/awk

BEGIN {
	FS = "\t";
}

NR == 1 {
	header = $0;
	t = 0;
}

NR > 1 {
	line[t] = $0;
	start[t] = $10;
	overlap[t] = $2;
	t++;
}

END {
	printf("%s", header);

	asorti(start, new_order, "@val_num_asc");

	best_line = line[new_order[1]];
	best_overlap = overlap[new_order[1]];
	old_start = start[new_order[1]];

	for(i = 2; i <= t; i++) {

		cur_line = line[new_order[i]];
		cur_overlap = overlap[new_order[i]];
		cur_start = start[new_order[i]];

		if(cur_start && cur_start == old_start) {
			if(cur_overlap > best_overlap) {
				best_overlap = cur_overlap;
				best_line = cur_line;
			}
		} else {
			printf("\n%s", best_line);
			best_overlap = cur_overlap;
			best_line = cur_line;
		}

		old_start = cur_start;
	}
}