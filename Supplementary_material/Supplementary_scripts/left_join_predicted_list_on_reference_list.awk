#!/usr/bin/awk

# This script generate a list containging all records from the reference list,
# and the match records from the predicted list

# Usage example:
# awk -f left_join_predicted_list_on_reference_list.awk ref.csv pred.csv > left_join.csv


function abs(v) {return v < 0 ? -v : v}


function compare_locations(ref_start, ref_end, rnie_start, rnie_end, strand,\
  overlap_perc,  distance,  ref_len) {

	# Make the "-" strand position negative, in order to make the function
	# treat the two strands as a single case 
	if(strand == "-") {
		old_ref_start = ref_start;
		old_rnie_start = rnie_start;
		ref_start = ref_end * (-1);
		ref_end = old_ref_start * (-1);
		rnie_start = rnie_end * (-1);
		rnie_end = old_rnie_start * (-1);
	}

	diff_start = rnie_start - ref_start;
	diff_end = rnie_end - ref_end;
	

	if(diff_start >= 0) {
		overlap_start = rnie_start;
	} else {
		overlap_start = ref_start;
	}

	if(diff_end >= 0) {
		overlap_end = ref_end;
	} else {
		overlap_end = rnie_end;
	}

	ref_len = ref_end - ref_start + 1;
	overlap_len = overlap_end - overlap_start + 1;

	if(rnie_start <= ref_end && rnie_end >= ref_start) {
		overlap_perc = overlap_len / ref_len * 100;
	} else if ( rnie_start <= ref_end && rnie_end < ref_start ) {
		overlap_len = 0;
		overlap_perc = 0;
	} else {
		overlap_len = 0;
		overlap_perc = 0;
	}

	location_fields = overlap_perc "\t" diff_start "\t" diff_end;
	return location_fields;
}

function get_best_predicted_term(n, overlap_perc, diff_start, rnie_mode, g,\
  best_idx,  i,  j) {

	if(n == 1) {
		return 1;
	}	

	max_overlap = 0;
	min_diff = 1*10^9;
	is_o = 0;

	for(i = 1; i <= n; i++) {

		if (overlap_perc[i] > max_overlap) {
			max_overlap = overlap_perc[i];
			o_idx[max_overlap] = i;
			is_o = 1;
		} else if ( overlap_perc[i] == max_overlap && is_o == 1 ) {
			o_idx[max_overlap] = o_idx[max_overlap] "\t" i;
		}

		cur_diff = abs(diff_start[i])
		if(cur_diff <= min_diff) {
			min_diff = cur_diff;
			d_idx[min_diff] = i;
		} else if ( cur_diff == min_diff && is_o == 0 ) {
			d_idx[min_diff] = d_idx[min_diff] "\t" i;
		}

	}

	if(is_o == 1) {
		n_best_idx = split(o_idx[max_overlap], best_idx, "\t");
	} else {
		n_best_idx = split(d_idx[min_diff], best_idx, "\t");
	}

	if(n_best_idx == 1) {
		return best_idx[1];
	} else {
		for (j = 1; j <= n_best_idx; j++) {
			if(rnie_mode[g][best_idx[j]] == "Genome") {
				return best_idx[j];
			}
		}
		return best_idx[1];
	}

}

BEGIN {

	FS = "\t";
	k = 0;
	file_idx = 0;

}

FNR == 1 {
	file_idx++;
}

file_idx == 1 && FNR == 1 {
	ref_NF = NF;
	ref_header = $0;
}

# GET THE REFERENCE LIST
file_idx == 1 && FNR > 1 {

	# g stands for "gene_name" ; o for "gene_current_occurence"
	g = $4;
			
	if(n_occ[g] == "") {
		n_occ[g] = 1;
	} else {
		n_occ[g]++;
	}
	o = n_occ[g];

	gene_name[k] = g;
	occ[g][k] = o;

	strand[g][o] = $1;
	ref_start[g][o] = $2;
	ref_end[g][o] = $3;

	if(ref_NF >= 4) {
		ref_suffix[g][o] = $4;
		for(i = 5; i <= ref_NF; i++) {
			ref_suffix[g][o] = ref_suffix[g][o] "\t" $i;
		}
	}

	k++;

}

file_idx == 2 && FNR == 1 {
	pref_NF = NF;
	pred_header = $0;
}

file_idx == 1 && FNR > 1 {

	g = $7;

	if(ref_start[g][1] != "") {

		if(pred_n_occ[g] == "") {				
			pred_n_occ[g] = 1;
		} else {
			pred_n_occ[g]++;
		}
		o = pred_n_occ[g];

		rnie_mode[g][o] = $4;
		pred_start[g][o] = $1;
		pred_end[g][o] = $2;
		pred_suffix[g][o] = $5;
		for(i = 6; i <= pred_NF; i++) {
			pred_suffix[g][o] = pred_suffix[g][o] "\t" $i;
		}
	} 
}

	}
}

END {

	header = ref_header "\t" pred_header;
	printf("%s", header);

	for(i = 0; i < k; i++) {


		g = gene_name[i];
		o = occ[g][i];

		if(rnie_mode[gene_name[i]][1] != "") {

			tot_occ = pred_n_occ[g]
			for(j = 1; j <= tot_occ; j++) {

				list_location_fields[j] = compare_locations(ref_start[g][o], ref_end[g][o], \
					pred_start[g][j], pred_end[g][j], strand[g][j]);
				split(list_location_fields[j], field, "\t", seps);
				list_overlap_perc[j] = field[1];
				list_diff_start[j] = field[2];

			}
			best_pred_idx = get_best_predicted_term(tot_occ, list_overlap_perc, list_diff_start, rnie_mode, g)
			location_fields = list_location_fields[best_pred_idx];
			b = best_pred_idx;

		} else {
			b = 0;
			pred_start[g][b] = pred_end[g][b] = rnie_mode[g][b] = "\t";
			location_fields = "\t\t";
			pred_suffix[g][b] = "";
			for(nf = 1; nf <= pred_NF - 7; nf++) {
				pred_suffix[g][b] = pred_suffix[g][b] "\t";
			}
		}
		delete list_location_fields;
		delete list_overlap_perc;
		delete list_distance;

		
		printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
			g, o, strand[g][o], \
			ref_start[g][o], ref_end[g][o], ref_suffix[g][o], \
			pred_start[g][b], pred_end[g][b], \
			rnie_mode[g][b], location_fields, \
			pred_suffix[g][b])
	}

}