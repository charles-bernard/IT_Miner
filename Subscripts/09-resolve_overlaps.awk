#!/usr/bin/awk

function abs(x) {
	return x < 0 ? -x : x
}


function compute_median(vector, n) {
	m = int(n/2) + 1;
	asorti(vector, sorted_order, "@val_num_asc");
	median = vector[sorted_order[m+1]];
	delete sorted_order;

	return median;
}


function compare(g, a, b) {
	if(nb_overlap[a] > nb_overlap[b]) {
		return a;
	} else if(nb_overlap[a] < nb_overlap[b]) {
		return b;
	}

	if(abs(start[g][a] - start[g][b] <= dev)) {
		if(bit[g][a] > bit[g][b]) {
			return a;
		} else if(bit[g][a] < bit[g][b]) {
			return b;
		}
	}

	if(abs(dist[g][a] - median) < abs(dist[g][b] - median)) {
		return a;
	} else if(abs(dist[g][a] - median) > abs(dist[g][b] - median)) {
		return b;
	}

	if(dist[g][a] <= dist[g][b]) {
		return a;
	} else {
		return b;
	}
}


function resolve_overlap(median, gene, n,\
	  i,  j,  k) {
	nb_overlap[1] = sum_overlap_len[1] = 0;

	for(i = 1; i <= n; i++) {
		ref_end = end[g][i];
		for(j = i + 1; j <= n; j++) {
			# N.B: comp_start always > ref_start since the input file is sorted
			comp_start = start[g][j];

			# Break the nested loop if the current IT doesn't overlap with 
			# the ref IT. Get the length of the overlap otherwise
			if(abs(comp_start > ref_end)) {
				break; 
			} else {
				cur_overlap_len = ref_end - comp_start + 1;

				sum_overlap_len[i] = overlap_len[i] + cur_overlap_len;
				nb_overlap[i]++;

				sum_overlap_len[j] = cur_overlap_len;
				nb_overlap[j] = 1;

				list_neighbor[i][nb_overlap[i]] = j;
				list_neighbor[j][1] = i;

				are_neighbors[i][j] = are_neighbors[j][i] = 1;
			}
			j++;
		}
	}

	asorti(sum_overlap_len, new_order, "@val_num_desc");

	for(i = 1; i <= n; i++) {
		# access to the IT by sum_overlap_len in desc order;
		ref_idx = new_order[i];

		# if IT has not been compared
		if(!is_discarded[ref_idx] && !is_best[ref_idx]) {

			cur_best = ref_idx;
			cur_best_len = sum_overlap_len[ref_idx];

			k = i + 1;
			comp_idx = new_order[k];
			comp_len = sum_overlap_len[comp_idx];

			while(cur_best_len && cur_best_len == comp_len) {

				# if IT to compare has not been discarded and is neighbor of the cur best
				if(!is_discarded[comp_idx] && are_neighbors[cur_best][comp_idx]) {
					cur_best = compare(cur_best, comp_idx);
				}

				k++;
				comp_idx = new_order[k];
				comp_len = sum_overlap_len[b];

			}

			# Update the best IT among all compared
			is_best[cur_best] = 1;

			# And discard all its neighbors
			if(nb_overlap[cur_best] > 0) {
				ngb = 1;
				ngb_idx = list_neighbor[cur_best][ngb];
				while(ngb_idx) {
					is_discarded[ngb_idx];
					ngb++;
					ngb_idx = list_neighbor[cur_best][ngb];
				}
			}

			retained_line[r] = line[g][cur_best];
			r++;
		}
	}
	delete is_best; 
	delete is_discarded;
	delete new_order;
	delete sum_overlap_len;
	delete nb_overlap;
	delete are_neighbors;
	delete list_neighbor;
}


BEGIN {
	FS = "\t";

	# counter for terminator
	t = 0;

	# counter for unique gene
	ng = 0;

	# array to count gene occurence
	og[""] = 0; 

	# counter for unnamed genes		
	unknown_ng = 0;	

	# counter for relevant distance from upgene
	d = 0;
	relevant_upgene_dist[0] = 0;

	# The variable 'dev' is sent to this script

	# counter for line that will be retained
	l = 0;
}

NR == 1 {
	header = $0;
}

NR > 1 {
	g = gene[t] = $7;

	# If gene has a name and current IT terminates the same gene 
	# as an anterior one, then the occurence of gene is incremented. 
	# Else initialized to 1.
	if(g && og[g] > 0) {
		og[g]++;
	} else {
		if(!g) {
			# Unnamed genes are assigned a generic id
			gene[t] = "Unknown n°" unknown_ng;
			unknown_ng++;
		} 
		og[g] = 1;
		ng++;
	}
	o = og[g];

	# Save key variables (accessed by gene name, and gene occurence)
	line[g][o] = $0;
	start[g][o] = $1;
	end[g][o] = $2;
	strand[g][o] = $3;
	bit[g][o] = $5;
	upgene_dist[g][o] = $10;

	# 'relevant_upgene_dist' will be furthered used to compute the median
	# distance from upstream gene
	if(upgene_dist[g][o] <= 300) {
		relevant_upgene_dist[d] = upgene_dist[g][o];
		d++;
	}
}

END {
	printf("%s", header);

	# compute median distance from upstream gene
	median = compute_median(relevant_upgene_dist, d);

	for(i = 0; i < ng; i++) {

		# Retrieve the gene and nb of time it occured in total. 
		g = gene[i];
		tot_o = og[g];

		# if total occurence is 1, then print the line
		if(tot_o == 1) {
			retained_line[r] = line[g][tot_o];
			r++;
		} else if(!already_examined[g]) {
			# retained_line and r are updated within the function.
			resolve_overlap(median, g, tot_o);
		}
	}
}
