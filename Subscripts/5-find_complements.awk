#!/usr/bin/awk

function free_dictionary(id,  s,  k,  strand) {
	for(s = 0; s < 2; s++) {
		k = 0;
		if(s == 0) { strand = "+"; } else { strand = "-"; }
		while(complement[id][strand][k] != "") {
			delete complement[id][strand][k];
			k++;
		}
	}
}

function remove_grp(end_idx, n,  i) {
	for(i = end_idx - n + 1; i <= end_idx; i++) {
		delete id_grp[i];
	}
}

function get_overlap_perc(ref_start, ref_end, compl_start, compl_end) {
	if(ref_start >= compl_start) {
		overlap_start = ref_start;
	} else {
		overlap_start = compl_start;
	}

	if(ref_end > compl_end) {
		overlap_end = compl_end;
	} else {
		overlap_end = ref_end;
	}

	ref_len = ref_end - ref_start + 1;
	overlap_len = overlap_end - overlap_start + 1;
	overlap_perc = overlap_len / ref_len * 100;

	return overlap_perc;
}

function get_best_complement(grp_id, idx) {
	ref_prefix = prefix[idx];
	split(ref_prefix, ref_fields, "\t");

	# Get info on the current terminator
	ref_start = ref_fields[1];
	ref_end = ref_fields[2];
	ref_strand = ref_fields[3];

	best_compl_idx = "";
	best_overlap = -1;

	# Visit dictionary of reverse complements, searching for best overlap
	k = 0;
	while(complement[grp_id][rev[ref_strand]][k] != "") {
		
		# Get info on the current reverse complement
		compl_idx = complement[grp_id][rev[ref_strand]][k];
		compl_prefix = prefix[compl_idx];
		split(compl_prefix, compl_fields, "\t");
		compl_start = compl_fields[1];
		compl_end = compl_fields[2];

		# best complement is the one which has the best overlap
		overlap = get_overlap_perc(ref_start, ref_end, compl_start, compl_end);
		if(overlap > best_overlap) {
			best_compl_idx = compl_idx;
			best_compl_start = compl_start;
			best_compl_end = compl_end;
			best_overlap = overlap;
		}
		k++;
	}

	# Write down some basic info on the reverse complement
	compl_info = "Start=" best_compl_start ";End=" best_compl_end \
		";Strand=" rev[ref_strand] ";% Overlapping=" best_overlap;

	# N.B: Keep in mind that the notion of best reverse compl is 
	# not necessarily reciprocal btw two complements

	return best_compl_idx "\t" compl_info;
}

function get_suffix(ref_i, compl_i) {
	split(prefix[ref_i], ref_fields, "\t");
	split(prefix[compl_i], compl_fields, "\t");

	# Tell which one in the couple is the most distant
	# from its upstream gene
	ref_upDist = ref_fields[10];
	compl_upDist = compl_fields[10];
	if(ref_upDist > compl_upDist) {
		dist_class = "Most Distant";
	} else {
		dist_class = "Least Distant";
	}

	# Tell if couple is convergent or not
	# Legend:
	# ----->: Gene
	# ||    : Intrinsic Terminator
	#
	# COVERGENCE is defined by the following genomic landscape
	# Scheme: 
	# ----->_||_<-----
	# Details:
	# + ___----->__||_____----->_
	# - _<-----____||__<-----____
	# 
	# CO-DIRECTIONAL
	# Scheme:
	# ----->_||_----->
	# Details:
	# + ___----->__||__----->____
	# - _<-----____||____<-----__
	#
	# DIVERGENCE
	# Scheme: 
	# <-----_||_----->
	# Details:
	# + _----->____||__----->____
	# - ___<-----__||_____<-----_

	ref_strand = ref_fields[3];

	ref_upStart = ref_fields[8];
	ref_upEnd = ref_fields[9];
	ref_dwStart = ref_fields[12];
	ref_dwEnd = ref_fields[13];

	compl_upStart = compl_fields[8];
	compl_upEnd = compl_fields[9];
	compl_dwStart = compl_fields[12];
	compl_dwEnd = compl_fields[13];

	if(ref_strand == "+") {
		if(ref_upEnd >= compl_dwEnd && ref_dwStart >= compl_upStart) {
			geno_class = "Convergence";
		} else if(ref_upEnd < compl_dwEnd && ref_dwStart < compl_upStart) {
			geno_class = "Divergence";
		} else {
			geno_class = "Co-directionality";
		}
	} else {
		if(ref_upStart <= compl_dwStart && ref_dwEnd <= compl_upEnd) {
			geno_class = "Convergence";
		} else if(ref_upStart > compl_dwStart && ref_dwEnd > compl_upEnd) {
			geno_class = "Divergence";
		} else {
			geno_class = "Co-directionality";
		}
	}

	ref_suffix = dist_class "\t" geno_class;
	return ref_suffix;
}

BEGIN {
	FS = "\t";

	# Counter for terminators
	t = 0;

	# These variables will be used for furthered comparisons
	old_start = old_end = -1;
	old_strand = "";

	# Required counters for building complements dictionary
	cur_c = -1;
	grp_count = for_c_count = rev_c_count = 0;

	rev["+"] = "-";
	rev["-"] = "+";

	# 'id_prefix' is sent to this script. 
	# If not sent, then the 'id_grp' will be basically the current 'grp_count'.
}

NR == 1 {
	# If input file has no line that has been searched before, the new header is stored 
	if(NF == 14) {
		header = $0 "\tID complementary couple\tCompl. Dist. Class\tCompl. Genomic Class\tInfo on the rev. compl.";
	} else {
		header = $0;
	}
}

# Only records with no complements fields (NF < 18 or no id) will be searched.
NR > 1 && (NF == 14 || $15 == "") {

	start = $1;
	end = $2;
	strand = $3;

	# prefix is built from fields ranging from 1 to 14
	prefix[t] = $1; for(i = 2; i <= 14; i++) { prefix[t] = prefix[t] "\t" $i; }


	# If the current terminator overlaps the previous
	if((start >= old_start && start < old_end) || (end > old_start && end <= old_end)) {
		# cur_c == -1 means that it is the first time the condition
		# is satisfied in the current group of reverse complements.
		if(cur_c == -1) {
			# Then the group counter needs to be incremented.
			grp_count++;

			if(old_start > -1) {
				# Previous terminator is also part of this group.
				id_grp[t-1] = grp_count;

				# The deduplication step discarded only exact duplicates.
				# Some duplicates with low overlap % might still remain.
				# A counter for such "partial" duplicates is required 
				# for further comparisons. One counter per strand.
				if(old_strand == "+") {
					for_c_count++;
					cur_c = for_c_count - 1; #index should start with O
				} else {
					rev_c_count++;
					cur_c = rev_c_count - 1;
				}

				# A dictionary of all terminators that are either
				# partial duplicates or complements is built.
				# Element is accessed by strand and strand_counter
				complement[grp_count][old_strand][cur_c] = t-1;
			}
		}

		id_grp[t] = grp_count;

		if(strand == "+") {
			for_c_count++;
			cur_c = for_c_count - 1; 
		} else {
			rev_c_count++;
			cur_c = rev_c_count - 1;
		}

		complement[grp_count][strand][cur_c] = t;

		# N.B: A dictionary might contain partial duplicates
		# and no complements. Such dictionary will be furthered removed.

	# Here, the current terminator closes the previous dictionary
	} else {
		# Dictionary containing only partial duplicates is removed
		if(cur_c > -1 && (for_c_count == 0 || rev_c_count == 0)) {
			free_dictionary(grp_count);
			remove_grp(k-1, cur_c+1)
			grp_count--;
		}

		cur_c = -1;
		for_c_count = rev_c_count = 0;
		id_grp[t] = "";
	}

	old_start = start;
	old_end = end;
	old_strand = strand;

	t++;
}

# If the current IT is a real complement:
NR > 1 && NF >= 18 && $15 != "" {
	# if previous line was also a real complement
	# concatenate the current line to the previous
	if(already_line[t]) {
		already_line[t] = already_line[t] "\n" $0;
	} else {
		already_line[t] = $0;
	}
} 

END {
	printf("%s", header);

	for(i = 0; i < t; i++) {

		if(already_line[i]) {
			printf("\n%s", already_line[i]);
		}

		# Each time a complement is encountered
		if(id_grp[i]) {

			# Find idx of the best reverse complement
			best_compl_output = get_best_complement(id_grp[i], i);
			split(best_compl_output, fields, "\t");
			compl_i = fields[1];
			compl_info = fields[2];


			suffix = get_suffix(i, compl_i);
			line = prefix[i] "\t" id_prefix id_grp[i] "\t" suffix "\t" compl_info;
		} else {
			line = prefix[i] "\t\t\t\t";
		}

		printf("\n%s", line);
	}
	# last line might be an already line
	if(already_line[i]) {
		printf("\n%s", already_line[i]);
	}
}