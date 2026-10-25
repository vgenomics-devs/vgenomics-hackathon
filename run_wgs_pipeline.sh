#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# ============================================================
# PIPELINE RUNTIME TIMER
# ============================================================
PIPELINE_START_EPOCH=$(date +%s)
PIPELINE_START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
CURRENT_STEP=""
CURRENT_STEP_START_EPOCH=0

# ============================================================
# CONTAINER-INTERNAL FIXED PATHS
# ============================================================
PIPELINE_DIR="/pipeline"
SCRIPT_DIR="${PIPELINE_DIR}/scripts"
EXOMISER_DIR="/opt/exomiser"
EXOMISER_JAR="${EXOMISER_DIR}/exomiser-cli-15.1.0.jar"

# ============================================================
# PYTHON SCRIPTS
# ============================================================
HPO_YAML_SCRIPT="${SCRIPT_DIR}/yml_hpo_parser.py"
VCF_EXTRACT_SCRIPT="${SCRIPT_DIR}/wgs_vcf_extract.py"
JSON_EXTRACT_SCRIPT="${SCRIPT_DIR}/exomiser_json_extract.py"
JSON_TSV_MERGE_SCRIPT="${SCRIPT_DIR}/exomiser_json_tsv_merge.py"
VCF_EXOMISER_MERGE_SCRIPT="${SCRIPT_DIR}/wgs_vcf_exomiser_merge.py"
PRIORITIZER_SCRIPT="${SCRIPT_DIR}/wgs_variants_prioritizer2.py"

# ============================================================
# FUNCTIONS
# ============================================================
usage() {
    cat <<EOF
Usage:
run_wgs_pipeline.sh \
    --config /config/wgs_pipeline.env \
    --input-vcf /input/sample.vcf.gz \
    --case-id CASE123 \
    --hpo-ids "HP:0001250,HP:0001263" \
    --output-dir /output

Required arguments:
    --config
        External pipeline configuration file.
    --input-vcf
        bgzip-compressed input VCF.
    --case-id
        Case/sample identifier.
    --hpo-ids
        Comma-separated HPO IDs.
    --output-dir
        Output directory.

Optional:
    --top N
        Number of final prioritized variants.
        Overrides TOP_N from configuration.
EOF
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

format_duration() {
    local total_seconds="$1"
    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))

    printf "%02dh:%02dm:%02ds" \
        "$hours" \
        "$minutes" \
        "$seconds"
}

start_step() {
    CURRENT_STEP="$1"
    CURRENT_STEP_START_EPOCH=$(date +%s)

    log "${CURRENT_STEP}"
}

end_step() {
    local step_end_epoch
    local step_elapsed

    step_end_epoch=$(date +%s)
    step_elapsed=$((step_end_epoch - CURRENT_STEP_START_EPOCH))

    log "${CURRENT_STEP} completed in $(format_duration "$step_elapsed")"

    CURRENT_STEP=""
    CURRENT_STEP_START_EPOCH=0
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

check_file() {
    [[ -f "$1" ]] || \
        die "Required file not found: $1"
}

check_dir() {
    [[ -d "$1" ]] || \
        die "Required directory not found: $1"
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || \
        die "Required command not found: $1"
}

check_variable() {
    local var_name="$1"
    [[ -n "${!var_name:-}" ]] || \
        die "Missing required configuration variable: ${var_name}"
}

on_error() {
    local exit_code=$?
    local line_number=$1
    local failure_epoch
    local total_elapsed

    failure_epoch=$(date +%s)
    total_elapsed=$((failure_epoch - PIPELINE_START_EPOCH))

    echo "ERROR: Pipeline failed at line ${line_number}" >&2
    echo "Exit code: ${exit_code}" >&2

    if [[ -n "${CURRENT_STEP}" && "${CURRENT_STEP_START_EPOCH}" -gt 0 ]]; then
        local step_elapsed
        step_elapsed=$((failure_epoch - CURRENT_STEP_START_EPOCH))

        echo "Failed step: ${CURRENT_STEP}" >&2
        echo "Step elapsed time: $(format_duration "$step_elapsed")" >&2
    fi

    echo "Pipeline elapsed time before failure: $(format_duration "$total_elapsed")" >&2

    exit "${exit_code}"
}

trap 'on_error $LINENO' ERR

# ============================================================
# ARGUMENTS
# ============================================================
CONFIG=""
INPUT_VCF=""
CASE_ID=""
HPO_IDS=""
OUTPUT_DIR=""

# Command-line --top overrides config TOP_N if supplied
CLI_TOP_N=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            [[ $# -ge 2 ]] || \
                die "Missing value for --config"
            CONFIG="$2"
            shift 2
            ;;
        --input-vcf)
            [[ $# -ge 2 ]] || \
                die "Missing value for --input-vcf"
            INPUT_VCF="$2"
            shift 2
            ;;
        --case-id)
            [[ $# -ge 2 ]] || \
                die "Missing value for --case-id"
            CASE_ID="$2"
            shift 2
            ;;
        --hpo-ids)
            [[ $# -ge 2 ]] || \
                die "Missing value for --hpo-ids"
            HPO_IDS="$2"
            shift 2
            ;;
        --output-dir)
            [[ $# -ge 2 ]] || \
                die "Missing value for --output-dir"
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --top)
            [[ $# -ge 2 ]] || \
                die "Missing value for --top"
            CLI_TOP_N="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

# ============================================================
# VALIDATE REQUIRED CLI ARGUMENTS
# ============================================================
[[ -n "$CONFIG" ]] || \
    die "--config is required"

[[ -n "$INPUT_VCF" ]] || \
    die "--input-vcf is required"

[[ -n "$CASE_ID" ]] || \
    die "--case-id is required"

[[ -n "$HPO_IDS" ]] || \
    die "--hpo-ids is required"

[[ -n "$OUTPUT_DIR" ]] || \
    die "--output-dir is required"
check_file "$CONFIG"

# ============================================================
# LOAD EXTERNAL CONFIGURATION
# ============================================================
# shellcheck source=/dev/null
source "$CONFIG"

# ============================================================
# REQUIRED CONFIG VARIABLES
# ============================================================
REQUIRED_CONFIG_VARS=(
    REFERENCE_FASTA
    EXOMISER_PROPERTIES
    EXOMISER_TEMPLATE
    SECONDARY_GENE_LIST
    EXOMISER_DATA_DIR
)

for variable in "${REQUIRED_CONFIG_VARS[@]}"; do
    check_variable "$variable"
done

# ============================================================
# CONFIG DEFAULTS
# ============================================================
THREADS="${THREADS:-8}"
MIN_DP="${MIN_DP:-10}"
MIN_GQ="${MIN_GQ:-20}"
HET_AB_MIN="${HET_AB_MIN:-0.20}"
HET_AB_MAX="${HET_AB_MAX:-0.80}"
TOP_N="${TOP_N:-25}"
# JAVA_OPTS="${JAVA_OPTS:--Xms4g -Xmx16g}"
JAVA_XMS="${JAVA_XMS:-4g}"
JAVA_XMX="${JAVA_XMX:-16g}"
EXOMISER_ASSEMBLY="${EXOMISER_ASSEMBLY:-hg38}"

# ------------------------------------------------------------
# CLI --top overrides config TOP_N
# ------------------------------------------------------------
if [[ -n "$CLI_TOP_N" ]]; then
    TOP_N="$CLI_TOP_N"
fi

# ============================================================
# BASIC VALUE VALIDATION
# ============================================================
[[ "$THREADS" =~ ^[0-9]+$ ]] || \
    die "THREADS must be an integer"

[[ "$MIN_DP" =~ ^[0-9]+$ ]] || \
    die "MIN_DP must be an integer"

[[ "$MIN_GQ" =~ ^[0-9]+$ ]] || \
    die "MIN_GQ must be an integer"

[[ "$TOP_N" =~ ^[0-9]+$ ]] || \
    die "TOP_N must be an integer"

# ============================================================
# VALIDATE SOFTWARE
# ============================================================
check_command bcftools
check_command tabix
check_command python3
check_command java

# ============================================================
# VALIDATE INPUT / CONFIG FILES
# ============================================================
check_file "$INPUT_VCF"
check_file "$REFERENCE_FASTA"
check_file "$EXOMISER_JAR"
check_file "$EXOMISER_PROPERTIES"
check_file "$EXOMISER_TEMPLATE"
check_file "$SECONDARY_GENE_LIST"
check_dir "$EXOMISER_DATA_DIR"

# ============================================================
# VALIDATE PIPELINE SCRIPTS
# ============================================================
check_file "$HPO_YAML_SCRIPT"
check_file "$VCF_EXTRACT_SCRIPT"
check_file "$JSON_EXTRACT_SCRIPT"
check_file "$JSON_TSV_MERGE_SCRIPT"
check_file "$VCF_EXOMISER_MERGE_SCRIPT"
check_file "$PRIORITIZER_SCRIPT"
[[ "$INPUT_VCF" == *.vcf.gz ]] || \
    die "Input VCF must end with .vcf.gz"


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"
WORK_DIR="${OUTPUT_DIR}/work"
RESULT_DIR="${OUTPUT_DIR}/results"
LOG_DIR="${OUTPUT_DIR}/logs"
PRIORITIZED_DIR="${OUTPUT_DIR}/prioritized"

mkdir -p \
    "$WORK_DIR" \
    "$RESULT_DIR" \
    "$LOG_DIR" \
    "$PRIORITIZED_DIR"

LOG_FILE="${LOG_DIR}/${CASE_ID}_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ============================================================
# START LOG
# ============================================================
log "============================================================"
log "Starting WGS pipeline"
log "============================================================"
log "Pipeline start       : ${PIPELINE_START_TIME}"
log "Case ID              : ${CASE_ID}"
log "Input VCF            : ${INPUT_VCF}"
log "Output directory     : ${OUTPUT_DIR}"
log "Config file          : ${CONFIG}"
log "Reference FASTA      : ${REFERENCE_FASTA}"
log "Exomiser properties  : ${EXOMISER_PROPERTIES}"
log "Exomiser template    : ${EXOMISER_TEMPLATE}"
log "Exomiser data        : ${EXOMISER_DATA_DIR}"
log "Secondary SF list    : ${SECONDARY_GENE_LIST}"
log "Threads              : ${THREADS}"
log "Top N                : ${TOP_N}"
log "============================================================"

# ============================================================
# OUTPUT FILE NAMES
# ============================================================
PRIMARY_VCF="${WORK_DIR}/${CASE_ID}_primary_chr.vcf.gz"
NORMALIZED_VCF="${WORK_DIR}/${CASE_ID}_normalized_split.vcf.gz"
PASS_VCF="${WORK_DIR}/${CASE_ID}_pass.vcf.gz"
NONREF_VCF="${WORK_DIR}/${CASE_ID}_nonref.vcf.gz"
QC_VCF="${RESULT_DIR}/${CASE_ID}_qc_filtered.vcf.gz"
CASE_YAML="${WORK_DIR}/${CASE_ID}_exomiser.yml"
VCF_TSV="${RESULT_DIR}/${CASE_ID}_vcf.tsv"
EXOMISER_PREFIX="${CASE_ID}_exomiser"
EXOMISER_JSON="${RESULT_DIR}/${EXOMISER_PREFIX}.jsonl"
EXOMISER_VARIANTS="${RESULT_DIR}/${EXOMISER_PREFIX}.variants.tsv"
EXOMISER_DETAILS="${RESULT_DIR}/${EXOMISER_PREFIX}_details.tsv"
EXOMISER_ENRICHED="${RESULT_DIR}/${EXOMISER_PREFIX}_enriched.tsv"
FINAL_MERGED="${RESULT_DIR}/${CASE_ID}_wgs_merged.tsv"

# ============================================================
# STEP 1
# PRIMARY CHROMOSOME FILTER
# ============================================================
# Check first 5 variant chromosome names
set +o pipefail
CHROMS=$(bcftools query -f '%CHROM\n' "$INPUT_VCF" | head -n 5)
set -o pipefail

# echo "Detected chromosome names:"
# echo "$CHROMS"

if echo "$CHROMS" | grep -q '^chr'; then
    echo "VCF already uses chr prefix."
    VCF_FOR_FILTER="$INPUT_VCF"
else
    echo "Adding chr prefix to chromosome names."

    CHR_RENAME_MAP="${OUTPUT_DIR}/${CASE_ID}.chr_rename_map.txt"

    {
        for chr in {1..22}; do
            printf "%s\tchr%s\n" "$chr" "$chr"
        done

        printf "X\tchrX\n"
        printf "Y\tchrY\n"
        printf "MT\tchrM\n"
    } > "$CHR_RENAME_MAP"

    NORMALIZED_CHROM_VCF="${OUTPUT_DIR}/${CASE_ID}.chrom_normalized.vcf.gz"

    bcftools annotate \
        --threads "$THREADS" \
        --rename-chrs "$CHR_RENAME_MAP" \
        -Oz \
        -o "$PRIMARY_VCF" \
        "$INPUT_VCF"

    tabix -f -p vcf "$NORMALIZED_CHROM_VCF"

    VCF_FOR_FILTER="$NORMALIZED_CHROM_VCF"
fi
end_step

# ============================================================
# STEP 2
# NORMALIZATION + SPLIT MULTIALLELIC
# ============================================================
start_step "STEP 2: Normalize and split multiallelic variants"
bcftools norm \
    --threads "$THREADS" \
    -f "$REFERENCE_FASTA" \
    -m -any \
    -c w \
    -Oz \
    -o "$NORMALIZED_VCF" \
    "$PRIMARY_VCF"

tabix -f -p vcf "$NORMALIZED_VCF"
end_step

# ============================================================
# STEP 3
# PASS FILTER
# ============================================================
start_step "STEP 3: Retain FILTER=PASS variants"
bcftools view \
    --threads "$THREADS" \
    -f PASS \
    -Oz \
    -o "$PASS_VCF" \
    "$NORMALIZED_VCF"

tabix -f -p vcf "$PASS_VCF"
end_step

# ============================================================
# STEP 4
# REMOVE HOMOZYGOUS REFERENCE
# ============================================================
start_step "STEP 4: Remove homozygous-reference variants"
bcftools view \
    --threads "$THREADS" \
    -e 'GT[0]="RR"' \
    -Oz \
    -o "$NONREF_VCF" \
    "$PASS_VCF"

tabix -f -p vcf "$NONREF_VCF"
end_step

# ============================================================
# STEP 5
# GENOTYPE QC
#
# DP >= MIN_DP
# GQ >= MIN_GQ
#
# HET:
# HET_AB_MIN <= AB <= HET_AB_MAX
#
# HOM_ALT:
# DP/GQ only
#
# HAPLOID_ALT:
# DP/GQ only
# ============================================================
start_step "STEP 5: DP/GQ/allele-balance filtering"

QC_EXPRESSION="\
FORMAT/DP[0] >= ${MIN_DP} \
&& FORMAT/GQ[0] >= ${MIN_GQ} \
&& ( \
    ( \
        GT[0]=\"het\" \
        && FORMAT/AD[0:0] + FORMAT/AD[0:1] > 0 \
        && FORMAT/AD[0:1] / \
           (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= ${HET_AB_MIN} \
        && FORMAT/AD[0:1] / \
           (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= ${HET_AB_MAX} \
    ) \
    || GT[0]=\"AA\" \
    || GT[0]=\"A\" \
)"

bcftools view \
    --threads "$THREADS" \
    -i "$QC_EXPRESSION" \
    -Oz \
    -o "$QC_VCF" \
    "$NONREF_VCF"

tabix -f -p vcf "$QC_VCF"
end_step
# ============================================================
# STEP 6
# GENERATE CASE-SPECIFIC EXOMISER YAML
# ============================================================
start_step "STEP 6: Generate Exomiser case YAML"
python3 "$HPO_YAML_SCRIPT" \
    --input "$EXOMISER_TEMPLATE" \
    --hpo-ids "$HPO_IDS" \
    --output "$CASE_YAML"

check_file "$CASE_YAML"
end_step

# ============================================================
# STEP 7
# RUN EXOMISER
# ============================================================
start_step "STEP 7: Run Exomiser 15.1.0"
export EXOMISER_DATA_DIRECTORY="$EXOMISER_DATA_DIR"

# JAVA_OPTS intentionally expanded as words.
# shellcheck disable=SC2086
java \
    "-Xms${JAVA_XMS}" \
    "-Xmx${JAVA_XMX}" \
    -Dspring.config.location="$EXOMISER_PROPERTIES" \
    -jar "$EXOMISER_JAR" \
    analyse \
    --analysis="$CASE_YAML" \
    --assembly="$EXOMISER_ASSEMBLY" \
    --vcf="$QC_VCF" \
    --output-directory="$RESULT_DIR" \
    --output-filename="$EXOMISER_PREFIX"

check_file "$EXOMISER_JSON"
check_file "$EXOMISER_VARIANTS"
end_step

# ============================================================
# STEP 8
# QC VCF → TSV
# ============================================================
start_step "STEP 8: Extract QC VCF data"
python3 "$VCF_EXTRACT_SCRIPT" \
    --input "$QC_VCF" \
    --output "$VCF_TSV"

check_file "$VCF_TSV"
end_step

# ============================================================
# STEP 9
# EXOMISER JSONL → TSV
# ============================================================
start_step "STEP 9: Extract Exomiser JSONL data"
python3 "$JSON_EXTRACT_SCRIPT" \
    --input "$EXOMISER_JSON" \
    --output "$EXOMISER_DETAILS"

check_file "$EXOMISER_DETAILS"
end_step

# ============================================================
# STEP 10
# EXOMISER variants.tsv + JSON details
# ============================================================
start_step "STEP 10: Merge Exomiser variants + JSON details"
python3 "$JSON_TSV_MERGE_SCRIPT" \
    --variants-tsv "$EXOMISER_VARIANTS" \
    --json-tsv "$EXOMISER_DETAILS" \
    --output "$EXOMISER_ENRICHED"

check_file "$EXOMISER_ENRICHED"
end_step

# ============================================================
# STEP 11
# VCF TSV + ENRICHED EXOMISER TSV
# ============================================================
start_step "STEP 11: Merge VCF + enriched Exomiser data"
python3 "$VCF_EXOMISER_MERGE_SCRIPT" \
    --vcf-tsv "$VCF_TSV" \
    --exomiser-tsv "$EXOMISER_ENRICHED" \
    --output "$FINAL_MERGED"

check_file "$FINAL_MERGED"
end_step

# ============================================================
# STEP 12
# FINAL VARIANT PRIORITIZATION
# ============================================================
start_step "STEP 12: Final variant prioritization"
python3 "$PRIORITIZER_SCRIPT" \
    --input "$FINAL_MERGED" \
    --secondary-gene-list "$SECONDARY_GENE_LIST" \
    --output-dir "$PRIORITIZED_DIR" \
    --prefix "$CASE_ID" \
    --top "$TOP_N"

end_step

# ============================================================
# PIPELINE END TIME / TOTAL RUNTIME
# ============================================================
PIPELINE_END_EPOCH=$(date +%s)
PIPELINE_END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
PIPELINE_TOTAL_SECONDS=$((PIPELINE_END_EPOCH - PIPELINE_START_EPOCH))
PIPELINE_TOTAL_RUNTIME=$(format_duration "$PIPELINE_TOTAL_SECONDS")

# ============================================================
# FINAL SUMMARY
# ============================================================
log "============================================================"
log "WGS PIPELINE COMPLETED SUCCESSFULLY"
log "Pipeline end         : ${PIPELINE_END_TIME}"
log "Total runtime        : ${PIPELINE_TOTAL_RUNTIME}"
log "Total seconds        : ${PIPELINE_TOTAL_SECONDS}"
log "============================================================"

echo
echo "Runtime:"
echo "  Start time     : ${PIPELINE_START_TIME}"
echo "  End time       : ${PIPELINE_END_TIME}"
echo "  Total seconds  : ${PIPELINE_TOTAL_SECONDS}"
echo "  Total runtime  : ${PIPELINE_TOTAL_RUNTIME}"

echo
echo "Case:"
echo "  ${CASE_ID}"
echo
echo "Configuration:"
echo "  ${CONFIG}"
echo
echo "Reference FASTA:"
echo "  ${REFERENCE_FASTA}"
echo
echo "QC VCF:"
echo "  ${QC_VCF}"
echo
echo "VCF TSV:"
echo "  ${VCF_TSV}"
echo
echo "Exomiser variants:"
echo "  ${EXOMISER_VARIANTS}"
echo
echo "Exomiser JSONL:"
echo "  ${EXOMISER_JSON}"
echo
echo "Exomiser details:"
echo "  ${EXOMISER_DETAILS}"
echo
echo "Enriched Exomiser:"
echo "  ${EXOMISER_ENRICHED}"
echo
echo "Final merged table:"
echo "  ${FINAL_MERGED}"
echo
echo "Final prioritized outputs:"
echo "  ${PRIORITIZED_DIR}/${CASE_ID}_top${TOP_N}.tsv"
echo "  ${PRIORITIZED_DIR}/${CASE_ID}_possible_compound_het.tsv"
echo "  ${PRIORITIZED_DIR}/${CASE_ID}_secondary_findings.tsv"
echo
echo "Pipeline log:"
echo "  ${LOG_FILE}"
echo
echo "============================================================"