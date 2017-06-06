#!/usr/bin/awk

# Purpose:
# This script is meant to join the output of the two modes of RNIE
# The Genome mode will serve as a reference list while the gene mode will allow to
# add to theis list the terminators that are not already present in the reference list 

# USAGE EXAMPLE (in bash)
# RNIE_GENOME="/home/charles/clones/Investigation_of_readthrough_signatures_in_bacteria/Build_reference_list_of_classified_terminators_in_E.coli/Dataset/1-List_of_Intrinsic_Terminators_predicted_by_RNIE/E_coli_K12_r2_(U00096.2)/01_Genome_(specific)_mode/2-list_predicted_intrinsic_terminators_with_seq_&_upstream_gene_redundant_mode.csv"
# RNIE_GENE="/home/charles/clones/Investigation_of_readthrough_signatures_in_bacteria/Build_reference_list_of_classified_terminators_in_E.coli/Dataset/1-List_of_Intrinsic_Terminators_predicted_by_RNIE/E_coli_K12_r2_(U00096.2)/02_Gene_(sensitive)_mode/2-list_predicted_intrinsic_terminators_with_seq_&_upstream_gene_redundant_mode.csv"
# RNIE_JOIN="/home/charles/clones/Investigation_of_readthrough_signatures_in_bacteria/Build_reference_list_of_classified_terminators_in_E.coli/Dataset/1-List_of_Intrinsic_Terminators_predicted_by_RNIE/E_coli_K12_r2_(U00096.2)/Concatenated_list.csv"
# SCRIPT="/home/charles/clones/Investigation_of_readthrough_signatures_in_bacteria/Build_reference_list_of_classified_terminators_in_E.coli/Scripts/2-treat_list_infered_by_RNIE/2_0-concatenate_the_two_modes_outputs_and_count_gene_occurence.awk"
# awk -f "$SCRIPT" "$RNIE_GENOME" "$RNIE_GENE" > "$RNIE_JOIN"

function get_line(mode, gene_occ) {
	line = $1 "\t" $2 "\t" mode "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" gene_occ "\t" $9 "\t" $10 "\t" $11 "\t" $12;
	return line;
}

BEGIN {
	FS = "\t";
	k = 0;
}

{
	if ( NR == FNR ) {

		# GET THE REFERENCE LIST

		if ( NR == 1 ) {

			header = get_line("Detection Mode", "Gene Occurence");
			printf("%s\n", header);

		} else {

			gene_name[k] = $8;

			if (counter[gene_name[k]] == "") {
				counter[gene_name[k]] = gene_occ[k] = 1;
			} else {
				counter[gene_name[k]]++;
				gene_occ[k] = counter[gene_name[k]];
			}

			gene_line[k] = get_line("Genome", gene_occ[k]);
			gene_start[k] = $4;

			k++;

		}

	} else if ( FNR != 1 ) {

		# EXTEND THE LIST WITH NEW DETECTED TERMINATORS IN GENE MODE

		gene_name[k] = $8;

		if (counter[gene_name[k]] == "") {
			counter[gene_name[k]] = gene_occ[k] = 1;
		} else {
			counter[gene_name[k]]++;
			gene_occ[k] = counter[gene_name[k]];
		}

		gene_line[k] = get_line("Gene", gene_occ[k]);
		gene_start[k] = $4;

		k++;

	}
}

END {
	N = asorti(gene_start, sorted_genes, "@val_num_asc");

	for (i = 1; i <= N; i++) {
		printf("%s\n", gene_line[sorted_genes[i]]);
	}
}