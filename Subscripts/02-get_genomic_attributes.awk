#!/usr/bin/awk

function get_gene_name(attribute) {
	# The Gene tag is the match of the Regex 
	# against the 'attribute' string
	match(attribute, /; *(gene_)?[nN]ame[= ]+"*[a-zA-Z0-9.+:_-]+[ "]*;/, gene_tag);
	if(!gene_tag[0]) {
		# If gene has no name, get its ID
		match(attribute, /; *(gene_)?[iI]d[= ]+[a-zA-Z0-9.+:_-]+[ "]*;/, gene_tag);
	}
	gsub(/^; *(gene_)?[nNiI](ame|d)[= ]+"*/, "", gene_tag[0]);
	gsub(/[ "]*;$/, "", gene_tag[0]);
	gene_name = gene_tag[0];
	if(!gene_name) {
		gene_name = "unknown n°" ug;
		ug++;
	}
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
			# Then each nt will be reversed
			# and concatenated to the previous
			nt = substr(seq, i, 1);

			# Check if nt is A,C,T or G 
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
	# The env variables 'old_for_idx' and 'old_rev_idx',
	# tell the function where to begin to read the gene dictionary

	if(strand == "+") {
		g = old_for_idx;

		# As long as the terminator is downstream of a CDS
		# Get the next CDS
		while(start > for_gene_start[g] && for_gene_start[g]) {
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
		# For reverse IT, functional start corresponds to annotation end
		while(end >= rev_gene_end[g] && rev_gene_end[g]) {
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

	# The genomic_landscape is a summary of all the information 
	# on the direct upstream and downstream genes of the IT
	genomic_landscape = up_gene_name "\t" up_gene_start "\t" up_gene_end "\t" up_dist \
		"\t" dw_gene_name "\t" dw_gene_start "\t" dw_gene_end "\t" dw_dist;

	return genomic_landscape;
}

BEGIN {
	file_idx = 0;
	FS = "\t";
	
	genome_seq = "";

	# Counter for gene/CDS (one per strand)
	for_count = 0;
	rev_count = 0;

	rev["A"] = "T";
	rev["T"] = "A";
	rev["C"] = "G";
	rev["G"] = "C";

	 # terminator counter
	t = 0;
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
	# Only the Genes/Coding Sequences matter for us
	if(region == "gene") {
		

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
	# N.B: If some lines of the current list have already been annotated
	# in the past, then the number of fields of the header should be >= 14.

	if(NF == 5) {
		header = $0 "\tIT Seq\t" \
			"UpGene Name\tUpGene Start\tUpGene End\tDist. upSTOPcodon-ITfunctionalStart\t" \
			"DwGene Name\tDwGene Start\tDwGene End\tDist. ITfunctionalEnd-dwSTARTcodon";
	} else {
		header = $0;
	}
	
	# Genomic landscape of each terminator will be searched
	# In the Gene Dictionaries.
	# As the terminator list is already sorted, the two
	# following counters are there to avoid reading parts
	# of the dictionaries which have been already read 
	old_for_idx = 0;
	old_rev_idx = 0;
}

# Only term with no genomic attribute (nb fields < 14) will be annotated
file_idx == 3 && FNR > 1 && NF == 5 {

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

file_idx == 3 && FNR > 1 && NF >= 14 {
	line[t] = $0; t++;
}

END {
	printf("%s", header);

	for(i = 0; i < t; i++) {
		printf("\n%s", line[i]);
	}
}