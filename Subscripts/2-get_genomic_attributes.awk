#!/usr/bin/awk

function get_gene_name(attribute) {
	# Extract gene tag matching Regex within attribute string
	match(attribute, /;gene=[a-zA-Z0-9.+:_-]+;/, gene_tag);

	# Get char indexes corresponding to gene name in the attribute string
	string_start = gene_tag[0, "start"] + length(";gene=");
	string_len = gene_tag[0, "length"] - length(";gene=" ";");

	# Get gene name
	gene_name = substr(attribute, string_start, string_len);
	return gene_name;
}

function write_seq(strand, start, end) {
	# The terminator sequence is just a substring of the genome sequence
	seq = substr(genome_seq, start, end-start+1);

	if(strand == "+") {
		return seq;
	} else {
		# If the terminator is on the "-" strand
		# the sequence needs to be reversed
		rev_seq = "";

		# The sequence is read from the end to the start
		for(i = length(seq); i > 0; i--) {
			# Then each nt will be furthered reversed
			# and concatenated to the previous
			nt = substr(seq, i, 1);

			# Simple Check to see if a nt can be reversed
			if(nt ~ /[ATCG]/) {
				rev_seq = rev_seq rev[nt];
			# A gap character "-" will be use if not
			} else {
				rev_seq = rev_seq "-";
			}

		}
		return rev_seq;
	}
}

function get_genomic_landscape(strand, start, end,  g) {
	if(strand == "+") {
		g = old_for_idx;
		# As long as the terminator is downstream of a CDS
		# Get the next CDS
		while(start > for_gene_start[g]) {
			g++; 
		}

		# The Current CDS is the first downstream from the terminator
		dw_gene_start = for_gene_start[g];
		dw_gene_end = for_gene_end[g];
		dw_gene_name = for_gene_name[g];
		if(dw_gene_start != "") { 
			dw_dist = dw_gene_start - end;
		} else { 
			dw_dist = ""; 
		}

		up_gene_start = for_gene_start[g-1];
		up_gene_end = for_gene_end[g-1];
		up_gene_name = for_gene_name[g-1];
		if(up_gene_start != "") { 
			up_dist = start - up_gene_end;
		} else {
			up_dist = "";
		}

		old_for_idx = g;
	} else {
		g = old_rev_idx;
		# For reverse terminators
		# functional start corresponds to annotation end
		while(end >= rev_gene_end[g]) {
			g++;
		}

		# The Current CDS is the first upstream from the terminator
		up_gene_start = rev_gene_start[g];
		up_gene_end = rev_gene_end[g];
		up_gene_name = rev_gene_name[g];
		if(up_gene_start != "") { 
			up_dist = up_gene_start - end;
		} else {
			up_dist = "";
		}

		dw_gene_start = rev_gene_start[g-1];
		dw_gene_end = rev_gene_end[g-1];
		dw_gene_name = rev_gene_name[g-1];
		if(dw_gene_start != "") {
			dw_dist = start - dw_gene_start;
		} else {
			dw_dist = "";
		}

		old_rev_idx = g;
	}
	genomic_landscape = up_gene_name "\t" up_gene_start "\t" up_gene_end "\t" up_dist \
		"\t" dw_gene_name "\t" dw_gene_start "\t" dw_gene_end "\t" dw_dist;

	return genomic_landscape;
}

BEGIN {
	file_idx = 0;
	FS = "\t";
	
	genome_seq = "";

	for_count = 0;
	rev_count = 0;

	rev["A"] = "T";
	rev["T"] = "A";
	rev["C"] = "G";
	rev["G"] = "C";

	t = 0; # terminator counter
}

FNR == 1 {
	file_idx++;
}

# The input file 1 corresponds to the genome in fasta
file_idx == 1 && FNR > 1 {
	# Concatenate each line in the fasta to get
	# the whole genome sequence as a single string
	genome_seq = genome_seq $0;
}

# The input file 2 is the annotation of the genome
# Lines starting with a comment are ignored
file_idx == 2 && !/^#.*$/ {
	region = $3;
	if(region == "CDS") {
		# Only the Coding Sequences matter for us

		# Important Fields
		start = $4;
		end = $5;
		strand = $7;
		attribute = $9;

		if(strand == "+") {
			# First Gene Dictionary refers to the "+" strand
			for_gene_start[for_count] = start;
			for_gene_end[for_count] = end;
			for_gene_name[for_count] = get_gene_name(attribute);
			for_count++;
		} else {
			# Second to "-"
			rev_gene_start[rev_count] = start;
			rev_gene_end[rev_count] = end;
			rev_gene_name[rev_count] = get_gene_name(attribute);
			rev_count++;
		}
	}
}

# The input file 3 is the list of terminators
file_idx == 3 && FNR == 1 {
	header_prefix = $0;
	header_suffix = "IT Seq\t" \
		"UpGene Name\tUpGene Start\tUpGene Stop\tUpGeneEnd-termStart\t" \
		"DwGene Name\tDwGene Start\tDwGene Stop\tDwGeneStart-termEnd";
}

file_idx == 3 && FNR > 1 {
	# Genomic landscape of each terminator will be searched
	# In the Gene Dictionaries.
	# As the terminator list is already sorted, the two
	# following counters are there to avoid reading parts
	# of the dictionaries which have been already read 
	old_for_idx = 0;
	old_rev_idx = 0;

	line_prefix = $0;

	start = $1;
	end = $2;
	strand = $3;

	# Get the sequence of the terminator
	term_seq = write_seq(strand, start, end);

	# The fct get_genomic_landscape will increment the dictionary index
	line_suffix = get_genomic_landscape(strand, start, end);

	line[t] = line_prefix "\t" term_seq "\t" line_suffix;
	t++;
}

END {
	printf("%s\t%s", header_prefix, header_suffix);

	for(i = 0; i < t; i++) {
		printf("\n%s", line[i]);
	}
}