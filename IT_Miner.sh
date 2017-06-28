#!/bin/bash
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
TMP_TOOL_STDERR=$(mktemp);

###########################################################
# WELCOME MESSAGE
###########################################################
cat <<WELCOME_MESSAGE
###########################################################
You are using IT Miner!

Let us find the intrinsic terminators in your bacterial 
genome!

WELCOME_MESSAGE

###########################################################
# I. FUNCTIONS
###########################################################
###### I.1 Handle Errors ##################################
###########################################################
function error_exit {
	# exit with error message and error code sent as args
	local ERR_MSG="$1"; local ERR_CODE="$2"; local LOG="${3:-$(mktemp)}";
	printf "ERROR:\n""$ERR_MSG""\n" | tee -a "$LOG" >&2; 
	exit "${ERR_CODE:-1}";				
}

function check_tool_stderr {
	# check if a tool has printed in stderr, if yes exit.
	local TOOL_EXEC_STDERR_FILE="$1";
	local TOOL_NAME="$2"; local LOG="$3";

	if [ -s "$TOOL_EXEC_STDERR_FILE" ]; then
		local TOOL_STDERR="The tool "\""$TOOL_NAME"\"" exited with the following error:\n" 
		TOOL_STDERR="$TOOL_STDERR"`cat "$TOOL_EXEC_STDERR_FILE"`
		rm "$TOOL_EXEC_STDERR_FILE";
		error_exit "$TOOL_STDERR" 4 "$LOG";
	fi
}

###########################################################
###### I.2 Get nb of terminators ##########################
###########################################################
function get_nb_term {
	local FILE="$1"; local DECREMENT="${2:-1}";
	N_IT=`awk -v d=$DECREMENT 'END { printf("%s", NR-d); }' "$FILE"`;
	printf $N_IT;
}

###########################################################
###### I.3 Test Pipeline Parameters #######################
###########################################################
function check_param {
	local PARAM="$1"; local NAME="$2"; local OPTION="$3";
	if [ ! "$PARAM" ]; then
		ERROR="Parameter "\""$NAME"\"" is required!\n";
		ERROR="$ERROR""Use option ""$OPTION"
		error_exit "$ERROR" 1;
	fi
}

function check_outdir {
	# Check if output dir exists; create it if needed 
	local DIR="$1";
	if [ ! -d "$DIR" ]; then
		printf "The ouput dir \"$DIR\" does not exist.\n"; 
		printf "Do you want to create it? (y|n)\n";
		while true; do
			read -p "> " YN
			case $YN in
				[Yy]* )	mkdir -p "$DIR"; break;;
				[Nn]* ) error_exit \
					"The output dir \"$DIR\" was not valid!" 3; break;;
				* )		printf "Please answer yes or no.\n";;
			esac
		done
		printf "\n";
	fi
}

function check_cutoff {
	CUTOFF="$1";
	if [ ! "$CUTOFF" ]; then
		CUTOFF="inclusive";
		printf "You did not mention a cutoff method.\n";
		printf "\tAn inclusive method will be used by default.\n";
	fi
}

function check_log {
	# Create log file if necessary
	local LOCAL_LOG="$1";
	if [ ! $LOCAL_LOG ]; then 
		LOG="$OUTPUT_DIR"/IT_Miner.log;
		printf "You did not mention a path for the log file.\n";
		printf "The log file will be automatically generated at:\n";
		printf "\t\"$LOG\"\n\n";
	fi
}

function check_file {
	# Check if file exists
	local FILE="$1";
	if [ ! -f "$FILE" ]; then
		echo "$FILE"
		error_exit \
		"The file \"$FILE\" does not exist!";
	fi
}

function check_rnie {
	# Check if rnie path is correct
	local RNIE_PATH="$1";
	if [[ ! "$RNIE_PATH" =~ ^.*rnie.pl$ ]]; then
		error_exit \
		"The RNIE path does not end with \"rnie.pl\"!" 2;
	fi
}

function check_fasta {
	# Check if input genome is in fasta
	local GENOME="$1";
	if [[ ! "$GENOME" =~ ^.*\.fa(sta)?$ ]]; then
		error_exit \
		"The genome file has no \".fa\" nor \".fasta\" extension!" 2;
	fi
	read -r -n 1 FIRST_CHAR < "$GENOME";
	if [[ "$FIRST_CHAR" != ">" ]]; then
		error_exit \
		"The fasta file does not begin with the required \">\" character!" 2;
	fi
}

function check_gff {
	# Check if input annotation is a gff file
	local ANNOTATION="$1";
	if [[ ! "$ANNOTATION" =~ ^.*\.gff(3)?$ ]]; then
		error_exit \
		"The annotation file has no \".gff\" extension!" 2;
	fi
}

###########################################################
###### I.4 Run RNIE #######################################
###########################################################
function run_rnie {
	# Run RNIE

	# Args
	local RNIE_PATH="$1"; local MODE="$2";
	local GENOME="$3"; local BIT_SCORE_THRESH="$4";
	local PREFIX="$5"; local COMMAND="$6";

	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "RNIE COMMAND:\n" >> "$COMMAND";
	(set -x;
		# perl "$RNIE_PATH" \
		# --fastafile "$GENOME" \
		# --"$MODE" \
		# --thresh "$BIT_SCORE_THRESH" \
		# --prefix "$PREFIX" \
		# 2>"$TMP_TOOL_STDERR"
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "RNIE/rnie.pl" "$LOG";
}

###########################################################
###### I.5 Concatenate ####################################
###########################################################
function concatenate {
	# Concatenate the two output lists of RNIE

	# Args
	local SCRIPT_PATH="$1";	
	local GENOME_GFF="$2"; local GENE_GFF="$3";	
	local OUT_LIST="$4"; local COMMAND="$5";

	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "CONCATENATE COMMAND:\n" >> "$COMMAND";
	(set -x; 
		awk -f "$SCRIPT_PATH" "$GENOME_GFF" "$GENE_GFF" \
		> "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";


	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "01-concatenate.awk" "$LOG";
}

###########################################################
###### I.6 Get Genomic Attributes #########################
###########################################################
function get_genomic_attributes {
	# Get Sequence of each IT, 
	# & up/down stream genes landscape

	# Args
	local SCRIPT_PATH="$1"; local GENOME="$2";
	local ANNOTATION="$3"; local INPUT_LIST="$4";
	local OUT_LIST="$5"; local COMMAND="$6";

	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "GET GENOMIC ATTRIBUTES COMMAND:\n" >> "$COMMAND";
	(set -x; 
		awk -f "$SCRIPT_PATH" "$GENOME" "$ANNOTATION" \
		"$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND"
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "02-get_genomic_attributes.awk" "$LOG";
}

###########################################################
###### I.7 Deduplicate ####################################
###########################################################
function deduplicate {
	# Deduplicate list

	# Args
	local SCRIPT_PATH="$1";
	local NT_DEV="$2"; local INPUT_LIST="$3";
	local OUT_LIST="$4"; local COMMAND="$5";


	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "DEDUPLICATE COMMAND:\n" >> "$COMMAND";
	(set -x;
		awk -v dev=$NT_DEV -f "$SCRIPT_PATH" \
		"$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "03-deduplicate.awk" "$LOG";
}

###########################################################
###### I.8 Discard ITs within CDS #########################
###########################################################
function discard_intra_cds {
	# Discard intra cds ITs

	# Args
	local SCRIPT_PATH="$1"; local INPUT_LIST="$2"; 
	local OUT_LIST="$3"; local INTRA_LIST="$4"; local COMMAND="$5";

	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "DISCARD ITs INTRA CDS COMMAND:\n" >> "$COMMAND";
	(set -x;
		awk -v INTRA_FILE="$INTRA_LIST" -f "$SCRIPT_PATH" \
		"$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "04-discard_intra_cds.awk" "$LOG";
}

###########################################################
###### I.9 Find Reverse Complementary ITs #################
###########################################################
function find_complements {
	# Find reverse complementary ITs 

	# Args
	local SCRIPT_PATH="$1"; local INPUT_LIST="$2"; 
	local OUT_LIST="$3"; local COMMAND="$4"; local ID_PREFIX="${5:-""}";

	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "FIND REVERSE COMPLEMENTARY ITS:\n" >> "$COMMAND";
	(set -x;
		awk -v id_prefix="$ID_PREFIX" -f "$SCRIPT_PATH" \
		"$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "05-identify_complements.awk" "$LOG";
}

function stats_complements {
	# Get Nb of complements per class..

	# Args
	local LIST="$1"; LOG="$2";
	AWK_TMP_OUT=$(mktemp);

	awk 'BEGIN { 
		FS="\t"; old_id = n = con = div = cod = 0; 
	}
	NR > 1 {
		if($15 != "") {
			id = $15; n++;
			if(id > old_id) {
				if($17 == "Convergence") {
					con++;
				} else if($17 == "Divergence") {
					div++;
				} else {
					cod++;
				}
			}
			old_id = id;
		}
	}
	END {
		printf("%d %d %d %d %d", n, id, con, div, cod);

	}' "$LIST" > "$AWK_TMP_OUT";

	read -a STATS < "$AWK_TMP_OUT";
	rm "$AWK_TMP_OUT";

cat $(mktemp) | tee -a "$LOG" <<COMPLEMENTS_MESSAGE

   ** ${STATS[0]} reverse complementary ITs detected!\n
           (belonging to ${STATS[1]} groups/couples):

      *** ${STATS[2]} groups are in gene convergence context

                     --->_||_<---

      *** ${STATS[3]} in divergence

                     <---_||_--->

      *** ${STATS[4]} in co-directionality

                     --->_||_--->

COMPLEMENTS_MESSAGE

}

###########################################################
###### I.10 Fake complements ##############################
###########################################################
function fake_complements {
	# Fake complement to ITs with no predicted complement

	# Args
	local FAKE_SCRIPT_PATH="$1"; local GENO_SCRIPT_PATH="$2";
	local FIND_COMPL_SCRIPT_PATH="$3"; local GENOME="$4";
	local ANNOTATION="$5"; local INPUT_LIST="$6";
	local OUT_LIST="$7"; local COMMAND="$8";
	
	#######################################################
	# Step 1: Fake complement coordinates
	#######################################################
	TMP_FAKE_OUT_LIST=$(mktemp);
	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "FAKE COMPLEMENTS:\n" >> "$COMMAND";
	(set -x;
		awk -f "$FAKE_SCRIPT_PATH" "$INPUT_LIST" \
		> "$TMP_FAKE_OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "06-fake_complements.awk" "$LOG";

	#######################################################
	# Step 2: Get genomic attributes of the faked complements
	#######################################################

	TMP_GENO_OUT_LIST=$(mktemp);
	get_genomic_attributes "$GENO_SCRIPT_PATH" \
		"$GENOME" "$ANNOTATION" \
		"$TMP_FAKE_OUT_LIST" "$TMP_GENO_OUT_LIST" "$COMMAND";

	#######################################################
	# Step 3: Get info on the couples formed by the fake complements
	#######################################################
	find_complements "$FIND_COMPL_SCRIPT_PATH" \
	"$TMP_GENO_OUT_LIST" "$OUT_LIST" "$COMMAND" "FakeCouple-";

	rm "$TMP_FAKE_OUT_LIST" "$TMP_GENO_OUT_LIST";
}

###########################################################
###### I.11 Compute cutoff ################################
###########################################################
function compute_cutoff {
	# Compute cutoff distance from upstream stop codon

	# Args
	local SCRIPT_PATH="$1"; local INPUT_LIST="$2";
	local CUTOFF="$3"; local OUTDIR="$4"; 
	local OUTFILE="$5"; local LOG="$6"; local COMMAND="$7";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "COMPUTE CUTOFF:\n" >> "$COMMAND";

	(set -x;
		Rscript "$SCRIPT_PATH" --table "$INPUT_LIST" \
		--cutoff "$CUTOFF" --fig_dir "$OUTDIR" --out_file "$OUTFILE" \
		2> "$TMP_TOOL_STDERR" | tee -a "$LOG" ;
	) 2>> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "07-compute_cutoff.R" "$LOG";
}

###########################################################
###### I.12 Filter ########################################
###########################################################
function filter {
	# Compute cutoff distance from upstream stop codon

	# Args
	local SCRIPT_PATH="$1"; local INPUT_LIST="$2";
	local OUT_LIST="$3"; local CUTOFF="$4"; local COMMAND="$5";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "FILTER:\n" >> "$COMMAND";

	(set -x;
		awk -v cutoff=$CUTOFF -f "$SCRIPT_PATH" \
		"$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "08-filter.awk" "$LOG";
}

###########################################################
###### I.13 Resolve overlaps ##############################
###########################################################
function resolve_overlaps {
	# Discard intra cds ITs

	# Args
	local SCRIPT_PATH="$1"; local INPUT_LIST="$2"; 
	local OUT_LIST="$3"; local COMMAND="$4";

	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "RESOLVE OVERLAPS COMMAND:\n" >> "$COMMAND";
	(set -x;
		awk -v dev=5 -f "$SCRIPT_PATH" \
		"$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "09-resolve_overlaps.awk" "$LOG";
}

###########################################################
###### I.14 Compute Free Energy ###########################
###########################################################
function compute_free_energy {
	# Compute free energy of ITs based on primary seq

	# Args
	local SCRIPT_PATH="$1"; local INPUT_LIST="$2"; 
	local OUT_LIST="$3"; local COMMAND="$4";

	RNAFOLD_IN="$(mktemp)";
	RNAFOLD_OUT="$(mktemp)";

	awk 'BEGIN { FS = "\t"; } NR > 1 { print $6; }' "$INPUT_LIST" > "$RNAFOLD_IN";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "COMPUTE FREE ENERGY COMMAND:\n" >> "$COMMAND";
	(set -x;
		RNAfold -i "$RNAFOLD_IN" --noPS > "$RNAFOLD_OUT" \
	2>"$TMP_TOOL_STDERR";) 2>> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "RNAfold" "$LOG";
	
	# Execute command while printing it into log
	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "JOIN RNIEfold OUT TO LIST OF ITs:\n" >> "$COMMAND";
	(set -x;
		awk -f "$SCRIPT_PATH" \
		"$RNAFOLD_OUT" "$INPUT_LIST" > "$OUT_LIST" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "10-join_RNAfold_out.awk" "$LOG";

	rm "$RNAFOLD_IN" "$RNAFOLD_OUT";
}

###########################################################
# II. PARAMETERS FOR THE PIPELINE
###########################################################
###### II.1 Read Options ##################################
###########################################################
TEMP=`getopt \
	-o \
		o:g:a:c:l:: \
	--long \
		output-dir:,rnie-path:,genome:,annotation:,cutoff:,log:: \
	-n 'report' -- "$@"`;
eval set -- "$TEMP";

###########################################################
###### II.2 Extract Options ###############################
###########################################################
while true ; do
	case "$1" in
		-o | --output-dir )
			case "$2" in
				"") error_exit "No output dir has not been provided!" 1; shift 2 ;;
				*) OUTPUT_DIR="$2"; shift 2 ;;
			esac ;;
		-g | --genome )
			case "$2" in
				"") error_exit "No genome (fasta) has been provided!" 1; shift 2 ;;
				*) GENOME="$2"; shift 2 ;;
			esac ;;
		-a | --annotation )
			case "$2" in
				"") error_exit "No annotation file (gff) has been provided!" 1; shift 2 ;;
				*) ANNOTATION="$2"; shift 2 ;;
			esac ;;
		-c | --cutoff )
			case "$2" in
				"") CUTOFF="conservative"; shift 2 ;;
				"conservative") CUTOFF="$2"; shift 2 ;;
				"inclusive") CUTOFF="$2"; shift 2;;
				"average") CUTOFF="$2"; shift 2;;
				*) error_exit "$2 is not an appropriate cutoff method!" 1: shift 2 ;;
			esac ;;
		-l | --log )
			case "$2" in
				*) LOG="$2"; shift 2 ;;
			esac ;;
		-- ) shift; break ;;
		*) error_exit "Internal error!" 1; 
	esac
done

###########################################################
###### II.3 Check pipeline parameters #####################
###########################################################
printf "###########################################################\n"
printf "PRE-PROCESSING)\n";
printf " * Checking the validity of the pipeline parameters\n\n";

RECQ_PARAMS=("$OUTPUT_DIR" "$GENOME" "$ANNOTATION");
PARAM_NAMES=("Output Directory" "Genome" "Annotation");
OPTIONS=("-o/--output-dir" "-g/--genome" "-a/--annotation");
for (( i=0; i<${#RECQ_PARAMS[@]}; i++ )); do 
	check_param "${RECQ_PARAMS[$i]}" "${PARAM_NAMES[$i]}" "${OPTIONS[$i]}";
done

check_outdir "$OUTPUT_DIR"; 
STEPS_DIR="$OUTPUT_DIR"/"Output_of_each_step"; mkdir -p "$STEPS_DIR";
FIG_DIR="$OUTPUT_DIR"/"Figures"; mkdir -p "$FIG_DIR";

check_cutoff "$CUTOFF";
check_log "$LOG"; > "$LOG";

INPUT_FILES=("$RNIE_PATH" "$GENOME" "$ANNOTATION");
for FILE in ${INPUT_FILES[@]}; do check_file "$FILE"; done

check_fasta "$GENOME";
check_gff "$ANNOTATION";

printf "Parameters Check was successful!\n\n";

COMMAND="$OUTPUT_DIR"/"commands.log"; > "$COMMAND";

###########################################################
###### II.4 Print pipeline parameters #####################
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "PARAMETERS LIST)\n" | tee -a "$LOG";
printf " * OUTPUT DIRECTORY:"\ \""$OUTPUT_DIR"\""\n" | tee -a "$LOG";
printf " * LOG FILE:        "\ \""$LOG"\""\n" | tee -a "$LOG";
printf " * GENOME:          "\ \""$GENOME"\""\n" | tee -a "$LOG";
printf " * ANNOTATION:      "\ \""$ANNOTATION"\""\n" | tee -a "$LOG";
printf " * CUTOFF:          "\ "$CUTOFF\n\n" | tee -a "$LOG";

###########################################################
# III. RUN RNIE
###########################################################
###### III.1 Define output parameters #####################
###########################################################
RNIE_GENOME_DIR="$STEPS_DIR"/"00-RNIE_Genome_Mode_Specific";
RNIE_GENE_DIR="$STEPS_DIR"/"00-RNIE_Gene_Mode_Sensitive";
mkdir -p "$RNIE_GENOME_DIR" "$RNIE_GENE_DIR";
PREFIX="IT_Miner";
BIT_SCORE_THRESH=14;

###########################################################
###### III.2 Run RNIE #####################################
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 0) Run RNIE\n" | tee -a "$LOG";

RNIE_PATH="$SCRIPT_PATH"/"Subscripts"/"RNIE"/"rnie.pl";
export RNIE=`dirname "$RNIE_PATH"`;

# Run genome mode
printf " * Genome (Specific) Mode\n" | tee -a "$LOG";
run_rnie "$RNIE_PATH" "genome" \
	"$GENOME" "$BIT_SCORE_THRESH" \
	"$RNIE_GENOME_DIR"/"$PREFIX" "$COMMAND";
RNIE_GENOME_OUT_GFF="$RNIE_GENOME_DIR"/"$PREFIX""-genomeMode-rnie.gff";

N_TERM=$(get_nb_term "$RNIE_GENOME_OUT_GFF" 0);
printf "   ** ""$N_TERM"" ITs predicted!\n" | tee -a "$LOG";

# Run gene mode
printf " * Gene (Sensitive) Mode\n" | tee -a "$LOG";
run_rnie "$RNIE_PATH" "gene" \
	"$GENOME" "$BIT_SCORE_THRESH" \
	"$RNIE_GENE_DIR"/"$PREFIX" "$COMMAND";
RNIE_GENE_OUT_GFF="$RNIE_GENE_DIR"/"$PREFIX""-geneMode-rnie.gff";

N_TERM=$(get_nb_term "$RNIE_GENE_OUT_GFF" 0);
printf "   ** ""$N_TERM"" ITs predicted!\n" | tee -a "$LOG";

###########################################################
# IV. Concatenate two RNIE mode outputs
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 1) Concatenate RNIE outputs\n" | tee -a "$LOG";

CONC_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"01-concatenate.awk";
CONC_OUT="$STEPS_DIR"/"Step01-Concatenated_list.csv";

concatenate "$CONC_SCRIPT_PATH" \
	"$RNIE_GENOME_OUT_GFF" "$RNIE_GENE_OUT_GFF" \
	"$CONC_OUT" "$COMMAND";
printf "   ** "$(get_nb_term "$CONC_OUT")" ITs retained!\n" | tee -a "$LOG";

###########################################################
# V. Get Genomic Attributes
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 2) Get Genomic Attributes of ITs\n" | tee -a "$LOG";
GENO_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"02-get_genomic_attributes.awk";
GENOMIC_OUT="$STEPS_DIR"/"Step02-Genomic_attributes_list.csv";

get_genomic_attributes "$GENO_SCRIPT_PATH" \
	"$GENOME" "$ANNOTATION" \
	"$CONC_OUT" "$GENOMIC_OUT" "$COMMAND";

###########################################################
# VI. Deduplicate
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 3) Deduplicate identical ITs\n" | tee -a "$LOG";

NT_DEV=3;
DEDUP_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"03-deduplicate.awk";
DEDUP_OUT="$STEPS_DIR"/"Step03-Deduplicated_list.csv";

deduplicate "$DEDUP_SCRIPT_PATH" $NT_DEV \
	"$GENOMIC_OUT" "$DEDUP_OUT" "$COMMAND";

printf "   ** "$(get_nb_term "$DEDUP_OUT")" ITs retained!\n" | tee -a "$LOG";

###########################################################
# VII. Discard ITs within CDS
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 4) Discard ITs intra CDS\n" | tee -a "$LOG";

DISCARD_INTRA_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"04-discard_intra_cds.awk";
DISCARD_INTRA_OUT="$STEPS_DIR"/"Step04-Without_intra_cds_list.csv";
LIST_INTRA="$STEPS_DIR"/"Step04-Only_intra_cds_list.csv";

discard_intra_cds "$DISCARD_INTRA_SCRIPT_PATH" \
	"$DEDUP_OUT" "$DISCARD_INTRA_OUT" "$LIST_INTRA" "$COMMAND";

printf "   ** "$(get_nb_term "$DISCARD_INTRA_OUT")" ITs retained!\n" | tee -a "$LOG";

###########################################################
# VIII. Find reverse complements
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 5) Find reverse complementary ITs\n" | tee -a "$LOG";

FIND_COMPL_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"05-identify_complements.awk";
FIND_COMPL_OUT="$STEPS_DIR"/"Step05-With_complements_list.csv";

find_complements "$FIND_COMPL_SCRIPT_PATH" \
	"$DISCARD_INTRA_OUT" "$FIND_COMPL_OUT" "$COMMAND";

stats_complements "$FIND_COMPL_OUT" "$LOG";

###########################################################
# IX. Fake complement for ITs with no predicted complement
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 6) Fake complement for ITs with no predicted complement\n" | tee -a "$LOG";

FAKE_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"06-fake_complements.awk";
FAKE_OUT="$STEPS_DIR"/"Step06-With_fake_complements_list.csv";

fake_complements "$FAKE_SCRIPT_PATH" \
	"$GENO_SCRIPT_PATH" "$FIND_COMPL_SCRIPT_PATH" \
	"$GENOME" "$ANNOTATION" \
	"$FIND_COMPL_OUT" "$FAKE_OUT" "$COMMAND";

###########################################################
# X. Find cut off distance to discard unlikely complements
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 7) Compute cut-off distance from upstream gene\n" | tee -a "$LOG"; 

CUTOFF_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"07-compute_cutoff.R";
CUTOFF_OUT="$STEPS_DIR"/"Step07-Cutoff.txt"

compute_cutoff "$CUTOFF_SCRIPT_PATH" "$FAKE_OUT" \
	"$CUTOFF" "$FIG_DIR" "$CUTOFF_OUT" "$LOG" "$COMMAND";

read -r CUTOFF_DISTANCE < "$CUTOFF_OUT";

###########################################################
# XI. Filter
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 8) Filter ITs above cutoff\n" | tee -a "$LOG"; 

FILTER_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"08-filter.awk";
TMP_FILTER_OUT=$(mktemp);
FILTER_OUT="$STEPS_DIR"/"Step08-Filtered_list.csv";

filter "$FILTER_SCRIPT_PATH" "$DISCARD_INTRA_OUT" \
	"$TMP_FILTER_OUT" "$CUTOFF_DISTANCE" "$COMMAND";

printf "   ** "$(get_nb_term "$TMP_FILTER_OUT")" ITs retained!\n" | tee -a "$LOG";

printf "\n .... ) Reinitializing stats on reverse complementary ITs\n"

find_complements "$FIND_COMPL_SCRIPT_PATH" \
	"$TMP_FILTER_OUT" "$FILTER_OUT" "$COMMAND";

stats_complements "$FILTER_OUT" "$LOG";
rm "$TMP_FILTER_OUT";

###########################################################
# XII. Resolve Overlapping ITs (partial duplicates)
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 9) Resolve overlapping ITs\n" | tee -a "$LOG"; 

RESOLVE_OVERLAPS_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"09-resolve_overlaps.awk";
RESOLVE_OVERLAPS_OUT="$STEPS_DIR"/"Step09-Resolved_overlaps_list.csv"

resolve_overlaps "$RESOLVE_OVERLAPS_SCRIPT_PATH" "$FILTER_OUT" \
	"$RESOLVE_OVERLAPS_OUT" "$COMMAND";

printf "   ** "$(get_nb_term "$RESOLVE_OVERLAPS_OUT")" ITs retained!\n" | tee -a "$LOG";

###########################################################
# XIII. Compute Free Energy and Secondary Structure of RNA
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "STEP 10) Compute free energy of ITs\n" | tee -a "$LOG"; 

JOIN_RNAFOLD_SCRIPT_PATH="$SCRIPT_PATH"/"Subscripts"/"10-join_RNAfold_out.awk";
FREE_ENERGY_OUT="$STEPS_DIR"/"Step10-With_free_energy_list.csv";

compute_free_energy "$JOIN_RNAFOLD_SCRIPT_PATH" \
	"$RESOLVE_OVERLAPS_OUT" "$FREE_ENERGY_OUT" "$COMMAND";