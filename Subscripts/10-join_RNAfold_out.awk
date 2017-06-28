#!/usr/bin/awk

BEGIN {
	file_idx = 0;
}

FNR == 1 {
	file_idx++;
	t = 0;

	if(file_idx == 1) {
		FS = " ";
	} else {
		FS = "\t";
	}
}

file_idx == 1 && FNR % 2 == 0 {
	ss[t] = $1;
	free_energy[t] = $2;
	gsub(/[()]/, "", free_energy[t]);
	t++;
}

file_idx == 2 && FNR == 1 {
	header = $0 "\tFree Energy (kcal/mol)\tDot-Bracket Secondary Structure";
	printf("%s", header);
}

file_idx == 2 && FNR > 1 {
	printf("\n%s\t%s\t%s", $0, free_energy[t], ss[t]);
	t++;
}
