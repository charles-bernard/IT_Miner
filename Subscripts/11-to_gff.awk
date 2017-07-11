#/usr/bin/awk

BEGIN {
	FS = "\t";
	source = "IT_Miner";
	# the variables 'organism' is sent to this script;
	if(!organism) {
		organism = "whole genome";
	};
	feature = "terminator";
	strength = ".";

	header = "##<organism>" \
		"\t<source>" \
		"\t<feature>" \
		"\t<start>" \
		"\t<end>" \
		"\t<RNIE_bit_score>" \
		"\t<strand>" \
		"\t<.>" \
		"\t<comment>"
	print header; 
}

NR > 1 {
	start = $1;
	end = $2;
	strand = $3;
	score = $5;
	sequence = $6;
	energy = $20;
	structure = $21;

	if(strand == "+") {
		id = "IT-" end;
	} else {
		id = "IT-" start;
	}

	printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t", 
		organism,
		source,
		feature,
		start,
		end,
		score,
		strand)
	printf("ID=%s;sequence=%s;pred_structure=%s;pred_free_energy=%s;pred_strength=%s\n",
		id,
		sequence,
		structure,
		energy,
		strength)
}
